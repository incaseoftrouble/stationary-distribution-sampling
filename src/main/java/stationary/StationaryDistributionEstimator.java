package stationary;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.naturals.set.NatBitSets;
import de.tum.in.probmodels.SelfLoopHandling;
import de.tum.in.probmodels.explorer.DefaultExplorer;
import de.tum.in.probmodels.generator.DtmcGenerator;
import de.tum.in.probmodels.graph.SccDecomposition;
import de.tum.in.probmodels.model.Distribution;
import de.tum.in.probmodels.model.MarkovChain;
import de.tum.in.probmodels.util.Sample;
import de.tum.in.probmodels.util.Util;
import de.tum.in.probmodels.util.Util.KahanSum;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntStack;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.function.IntConsumer;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import parser.State;
import prism.PrismException;
import stationary.component.AbsorbingComponent;
import stationary.component.BottomComponent;
import stationary.component.NontrivialApproximatingComponent;
import stationary.util.Bound;
import stationary.util.Explore;
import stationary.util.FrequencyRecord;

public final class StationaryDistributionEstimator {
  private record SamplingTarget(double weight, BottomComponent component, Distribution.WeightFunction weights) {}

  private static final Logger logger = Logger.getLogger(StationaryDistributionEstimator.class.getName());
  private static final int MAX_EXPLORES_PER_SAMPLE = 10;

  private final double precision;
  private final MarkovChain model;
  private final Mode mode;

  private final DefaultExplorer<State, MarkovChain> explorer;
  private final List<BottomComponent> components = new ArrayList<>();
  private final IntSet newStatesSinceComponentSearch = new IntOpenHashSet();
  private final Int2ObjectMap<BottomComponent> statesInBottomComponents = new Int2ObjectOpenHashMap<>();
  private final Int2ObjectMap<Int2ObjectMap<Bound>> componentReachability = new Int2ObjectOpenHashMap<>();
  private final Int2DoubleMap exitProbability = new Int2DoubleOpenHashMap();

  private int loopCount = 0;
  private int loopStopsUntilCollapse;
  private long mainLoopCount = 0L;
  private long sampledStatesCount = 0L;
  private long computedComponentUpdates = 0L;
  private long computedTransientUpdates = 0L;
  private long lastProgressUpdate = 0L;
  private long lastProgressUpdateLoopCount = 0L;
  private long abortedSamples = 0L;
  private final boolean preExplore;

  public StationaryDistributionEstimator(DtmcGenerator generator, double precision, Mode mode, boolean preExplore) {
    this.precision = precision;
    this.mode = mode;
    this.model = new MarkovChain();
    this.preExplore = preExplore;

    this.explorer = DefaultExplorer.of(model, generator, SelfLoopHandling.KEEP);
    loopStopsUntilCollapse = 10;

    // Fail-fast if these are accessed for a non-transient state
    exitProbability.defaultReturnValue(Double.NaN);
  }

  public State getState(int s) {
    return explorer.getState(s);
  }

  private void preExplore() throws PrismException {
    int initialState = model.getFirstInitialState();
    Explore.explore(explorer, initialState);
    model.findDeadlocks(true);

    List<NatBitSet> components = SccDecomposition.computeSccs(this.model::getSuccessors, IntSet.of(initialState), s -> true, false)
        .stream().filter(component -> SccDecomposition.isBscc(this.model::getSuccessors, component)).toList();
    if (logger.isLoggable(Level.FINE)) {
      logger.log(Level.FINE, String.format("Found %d BSCCs with %d states",
          components.size(), components.stream().mapToInt(NatBitSet::size).sum()));
    }
    components.forEach(this::createComponent);
  }

  public Int2ObjectMap<FrequencyRecord> solve() throws PrismException {
    int initialState = model.getFirstInitialState();
    lastProgressUpdate = System.currentTimeMillis();
    newStatesSinceComponentSearch.add(initialState);

    if (preExplore) {
      preExplore();
    }
    while (computeTotalError(initialState) > precision) {
      mainLoopCount += 1L;
      BottomComponent initialComponent = statesInBottomComponents.get(initialState);
      if (initialComponent == null) {
        sample(initialState);
      } else {
        computedComponentUpdates += 1;
        initialComponent.countVisit();
        initialComponent.update(initialState);
      }
      logProgress(false);
    }
    logProgress(true);

    Int2ObjectMap<FrequencyRecord> frequency = new Int2ObjectOpenHashMap<>(statesInBottomComponents.size());
    for (BottomComponent component : components) {
      Bound reachability = getComponentReachabilityBounds(component, initialState);
      double error = reachability.upper() * component.error();
      component.states().intStream().forEach((int s) ->
          frequency.put(s, new FrequencyRecord(component.frequency(s).lower() * reachability.lower(), error)));
    }

    assert Util.lessOrEqual((1.0 - getAbsorptionLowerBound(initialState)) + components.stream().mapToDouble(component ->
        getComponentReachabilityUpperBound(component, initialState) * component.error()).max().orElseThrow(), precision);

    return frequency;
  }

  private void logProgress(boolean force) {
    if (logger.isLoggable(Level.INFO)) {
      long now = System.currentTimeMillis();
      long timeDelta = now - this.lastProgressUpdate;
      if (!force && timeDelta < TimeUnit.SECONDS.toMillis(5L)) {
        return;
      }
      long loopDelta = mainLoopCount - this.lastProgressUpdateLoopCount;
      this.lastProgressUpdate = now;
      this.lastProgressUpdateLoopCount = mainLoopCount;

      int initialState = model.getFirstInitialState();

      String boundsString;
      if (components.isEmpty()) {
        boundsString = "None";
      } else {
        boundsString = components.stream()
            .sorted(Comparator.comparingDouble(component ->
                -getComponentReachabilityLowerBound(component, initialState) * component.error()))
            .limit(10L)
            .map(component -> "%s@[%.5f]".formatted(getComponentReachabilityBounds(component, initialState), component.error()))
            .collect(Collectors.joining(", "));
        if (components.size() > 10) {
          boundsString += " ... (%d)".formatted(components.size());
        }
      }
      logger.log(Level.INFO, String.format("Progress Report:%n"
              + "    %d loops (%d/sec), %d sampling steps, %d aborted%n"
              + "    loop condition %.5g, exit probability %.5g in initial state%n"
              + "    %d explored states, %d states in BSSCs (%d absorbing)%n"
              + "    %d component updates, %d transient state updates%n"
              + "    Bounds: %s",
          mainLoopCount, loopDelta * TimeUnit.SECONDS.toMillis(1L) / timeDelta, sampledStatesCount, abortedSamples,
          computeTotalError(initialState), getExitProbabilityUpperBound(initialState),
          explorer.exploredStateCount(), statesInBottomComponents.size(),
          statesInBottomComponents.values().stream().map(BottomComponent::states).mapToInt(IntSet::size).filter(i -> i == 1).count(),
          computedComponentUpdates, computedTransientUpdates,
          boundsString
      ));
    }
  }

  private double getExitProbabilityUpperBound(int state) {
    if (preExplore) {
      return 0.0;
    }

    assert !(exitProbability.containsKey(state) && statesInBottomComponents.containsKey(state));
    double bound = statesInBottomComponents.containsKey(state) ? 0.0 : exitProbability.getOrDefault(state, 1.0);
    assert Util.lessOrEqual(0.0, bound) && Util.lessOrEqual(bound, 1.0);
    //noinspection FloatingPointEquality
    assert explorer.isExploredState(state) || bound == 1.0;
    return bound;
  }

  private double getAbsorptionLowerBound(int state) {
    if (statesInBottomComponents.containsKey(state)) {
      return 1.0;
    }
    Int2ObjectMap<Bound> reachabilityMap = componentReachability.get(state);
    return reachabilityMap == null ? 0.0 : reachabilityMap.values().stream().mapToDouble(Bound::lower).sum();
  }

  private double computeTotalError(int state) {
    double error;
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      Int2ObjectMap<Bound> reachabilityMap = componentReachability.get(state);
      if (reachabilityMap == null) {
        error = 1.0;
      } else {
        KahanSum exploredWeightedError = new KahanSum();
        KahanSum absorptionProbability = new KahanSum();
        for (BottomComponent component : components) {
          double componentError = component.error();
          Bound bound = reachabilityMap.get(component.index);
          if (bound != null) {
            exploredWeightedError.add(componentError * bound.lower());
            absorptionProbability.add(bound.lower());
          }
        }
        assert Util.lessOrEqual(exploredWeightedError.get(), 1.0) : "Got unexpected error: %s".formatted(exploredWeightedError);
        assert Util.lessOrEqual(absorptionProbability.get(), 1.0) : "Got unexpected absorption: %s".formatted(absorptionProbability);
        error = exploredWeightedError.get() + (1.0 - absorptionProbability.get());
      }
    } else {
      error = stateComponent.error();
    }
    assert Util.lessOrEqual(0.0, error) && Util.lessOrEqual(error, 1.0) : "Got unexpected frequency error estimate %.5f".formatted(error);
    return error;
  }

  private Bound getComponentReachabilityBounds(BottomComponent component, int state) {
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      return computeComponentReachabilityBoundsTransient(component, state);
    }
    return stateComponent.equals(component) ? Bound.ONE : Bound.ZERO;
  }

  private Bound computeComponentReachabilityBoundsTransient(BottomComponent component, int state) {
    assert !statesInBottomComponents.containsKey(state);
    Int2ObjectMap<Bound> reachability = componentReachability.get(state);
    return reachability == null ? Bound.UNKNOWN : reachability.getOrDefault(component.index, Bound.UNKNOWN);
    /* if (reachability == null) {
      return Bound.UNKNOWN;
    }
    Bound bound = reachability.get(component.index);
    if (bound == null) {
      return Bound.UNKNOWN;
    }
    double sum = reachability.values().stream().mapToDouble(Bound::lower).sum();
    if (bound.upper() > 1 - sum + bound.lower()) {
      Bound newBound = Bound.of(bound.lower(), 1 - sum + bound.lower());
      reachability.put(component.index, newBound);
      return newBound;
    }
    return bound; */
  }

  private double getComponentReachabilityLowerBound(BottomComponent component, int state) {
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      return getComponentReachabilityLowerBoundTransient(component, state);
    }
    return stateComponent.equals(component) ? 1.0 : 0.0;
  }

  private double getComponentReachabilityLowerBoundTransient(BottomComponent component, int state) {
    assert !statesInBottomComponents.containsKey(state);
    return getComponentReachabilityBounds(component, state).lower();
  }

  private double getComponentReachabilityUpperBound(BottomComponent component, int state) {
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      return getComponentReachabilityUpperBoundTransient(component, state);
    }
    return stateComponent.equals(component) ? 1.0 : 0.0;
  }

  private double getComponentReachabilityUpperBoundTransient(BottomComponent component, int state) {
    assert !statesInBottomComponents.containsKey(state);
    return getComponentReachabilityBounds(component, state).upper();
  }

  private void updateComponentReachabilityBound(BottomComponent component, int state, Bound bound) {
    assert !statesInBottomComponents.containsKey(state);
    if (Util.isZero(bound.lower()) && Util.isOne(bound.upper())) {
      return;
    }
    Int2ObjectMap<Bound> reachability = componentReachability.computeIfAbsent(state, k -> new Int2ObjectOpenHashMap<>());
    assert reachability.getOrDefault(component.index, Bound.UNKNOWN).contains(bound)
        : "Updating reachability of %d from %s to %s".formatted(state, reachability.getOrDefault(component.index, Bound.UNKNOWN), bound);
    reachability.put(component.index, bound);
    // Bound oldBound = reachability.getOrDefault(component.index, Bound.UNKNOWN);
    // reachability.put(component.index, Bound.of(bound.lower(), Math.min(bound.upper(), oldBound.upper())));
  }

  private BottomComponent createComponent(NatBitSet states) {
    assert states.intStream().noneMatch(statesInBottomComponents::containsKey) : "States %s already in component".formatted(states);
    assert SccDecomposition.isBscc(this.model::getSuccessors, states);
    int index = components.size();
    BottomComponent component = states.size() == 1
        ? new AbsorbingComponent(index, states.firstInt())
        : new NontrivialApproximatingComponent(index, states, model::getTransitions);
    states.forEach((IntConsumer) state -> statesInBottomComponents.put(state, component));
    components.add(component);
    exitProbability.keySet().removeAll(states);
    componentReachability.keySet().removeAll(states);
    return component;
  }

  private void sample(int initialState) {
    assert !statesInBottomComponents.containsKey(initialState);

    Distribution.WeightFunction samplingWeight;
    Set<BottomComponent> updateComponentReachability = new HashSet<>(2);

    if (mode == Mode.SAMPLE_TARGET) {
      List<SamplingTarget> samplingTargets = new ArrayList<>(2 * components.size());
      double exitProbability = getExitProbabilityUpperBound(initialState);
      for (BottomComponent component : components) {
        Bound reachabilityBounds = getComponentReachabilityBounds(component, initialState);
        if (component.error() > 0.0) {
          samplingTargets.add(new SamplingTarget(reachabilityBounds.upper() * component.error(), component,
              (s, p) -> p * getComponentReachabilityUpperBound(component, s)));
        }
        double reachabilityScore = reachabilityBounds.difference();
        if (reachabilityScore > 0.0) {
          samplingTargets.add(new SamplingTarget(reachabilityScore, component,
              (s, p) -> p * getComponentReachabilityBounds(component, s).difference()));
        }
      }
      samplingTargets.add(new SamplingTarget(exitProbability, null, (s, p) -> p * getExitProbabilityUpperBound(s)));
      SamplingTarget target = Sample.sampleWeighted(samplingTargets, SamplingTarget::weight).orElseThrow();
      if (target.component() != null) {
        updateComponentReachability.add(target.component());
      }
      samplingWeight = target.weights();
    } else if (mode == Mode.SAMPLE_NAIVE) {
      samplingWeight = (s, p) -> p;
    } else {
      throw new AssertionError();
    }

    int exploreDuringSample = 0;
    int currentState = initialState;

    boolean checkForComponents = false;

    IntSet visitedStateSet = new IntOpenHashSet();
    IntList visitedStates = new IntArrayList();
    IntStack visitStack = (IntStack) visitedStates;

    while (true) {
      assert explorer.isExploredState(currentState);

      BottomComponent component = statesInBottomComponents.get(currentState);
      if (component != null) {
        updateComponentReachability.add(component);
        computedComponentUpdates += 1L;
        component.update(currentState);
        component.countVisit();
        break;
      }

      assert !visitedStateSet.contains(currentState);
      visitedStateSet.add(currentState);
      visitStack.push(currentState);

      // Sample the successor
      Distribution transitions = model.getTransitions(currentState);
      if (transitions == null || transitions.isEmpty()) { // Deadlock state
        updateComponentReachability.add(createComponent(NatBitSets.singleton(currentState)));
        visitedStateSet.remove(currentState);
        newStatesSinceComponentSearch.remove(currentState);
        break;
      }

      int nextState = transitions.sampleWeightedFiltered(samplingWeight, s -> !visitedStateSet.contains(s));
      if (nextState == -1) {
        checkForComponents = true;
        abortedSamples += 1;
        break;
      }
      assert !visitedStateSet.contains(nextState);

      sampledStatesCount += 1L;
      if (!explorer.isExploredState(nextState)) {
        if (exploreDuringSample == MAX_EXPLORES_PER_SAMPLE) {
          break;
        }

        exploreDuringSample += 1;
        explorer.exploreState(nextState);
        newStatesSinceComponentSearch.add(nextState);
      }
      currentState = nextState;
    }

    if (checkForComponents) {
      loopCount += 1;
      if (loopCount > loopStopsUntilCollapse) {
        loopCount = 0;

        if (!newStatesSinceComponentSearch.isEmpty()) {
          assert newStatesSinceComponentSearch.intStream().noneMatch(this.statesInBottomComponents::containsKey);
          List<NatBitSet> components = SccDecomposition.computeSccs(this.model::getSuccessors, newStatesSinceComponentSearch,
                  s -> this.explorer.isExploredState(s) && !statesInBottomComponents.containsKey(s), false)
              .stream().filter(component -> SccDecomposition.isBscc(this.model::getSuccessors, component)).toList();
          newStatesSinceComponentSearch.clear();
          if (logger.isLoggable(Level.FINE)) {
            logger.log(Level.FINE, String.format("Found %d BSCCs with %d states",
                components.size(), components.stream().mapToInt(NatBitSet::size).sum()));
          }
          for (NatBitSet component : components) {
            createComponent(component);
            visitedStateSet.removeAll(component);
          }
        }
        loopStopsUntilCollapse = explorer.exploredStates().size();
      }
    }

    // Propagate values backwards along the path
    boolean updateExitProbability = !preExplore;
    while (!visitStack.isEmpty() && (updateExitProbability || !updateComponentReachability.isEmpty())) {
      int state = visitStack.popInt();

      if (!visitedStateSet.contains(state)) {
        continue;
      }
      assert !statesInBottomComponents.containsKey(state);
      Distribution transitions = model.getTransitions(state);

      if (updateExitProbability) {
        double exitProbability = transitions.sumWeightedExceptJacobi(this::getExitProbabilityUpperBound, state);
        if (Double.isNaN(exitProbability) || Util.isOne(exitProbability)) {
          updateExitProbability = false;
        } else {
          assert Util.lessOrEqual(0.0, exitProbability);
          double previous = this.exitProbability.put(state, exitProbability);
          assert Double.isNaN(previous) || Util.lessOrEqual(exitProbability, previous) :
              "Updating exit bound of %d from %.5g to %.5g".formatted(state, previous, exitProbability);
        }
      }

      if (!updateComponentReachability.isEmpty()) {
        Iterator<BottomComponent> iterator = updateComponentReachability.iterator();
        while (iterator.hasNext()) {
          computedTransientUpdates += 1L;

          BottomComponent component = iterator.next();
          double lowerBound = transitions.sumWeightedExceptJacobi(s -> getComponentReachabilityLowerBound(component, s), state);
          double upperBound = transitions.sumWeightedExceptJacobi(s -> getComponentReachabilityUpperBound(component, s), state);
          assert Double.isNaN(lowerBound) == Double.isNaN(upperBound);
          if (Double.isNaN(lowerBound) || (Util.isZero(lowerBound) && Util.isOne(upperBound))) {
            iterator.remove();
          } else {
            updateComponentReachabilityBound(component, state, Bound.of(lowerBound, upperBound));
          }
        }
      }
    }
  }

  public enum Mode {
    SAMPLE_NAIVE, SAMPLE_TARGET
  }
}
