package stationary;

import de.tum.in.naturals.map.Nat2ObjectDenseArrayMap;
import de.tum.in.probmodels.explorer.DefaultExplorer;
import de.tum.in.probmodels.explorer.SelfLoopHandling;
import de.tum.in.probmodels.generator.Generator;
import de.tum.in.probmodels.graph.Component;
import de.tum.in.probmodels.graph.SccDecomposition;
import de.tum.in.probmodels.model.MutableDenseSystem;
import de.tum.in.probmodels.model.distribution.Distribution;
import de.tum.in.probmodels.model.impl.DynamicQuotient;
import de.tum.in.probmodels.util.Sample;
import de.tum.in.probmodels.util.Util;
import de.tum.in.probmodels.util.Util.KahanSum;
import de.tum.in.probmodels.values.Bounds;
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
import java.util.Collection;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.function.BiConsumer;
import java.util.function.IntConsumer;
import java.util.function.IntFunction;
import java.util.function.IntToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import javax.annotation.Nullable;
import parser.State;
import stationary.component.AbsorbingComponent;
import stationary.component.BottomComponent;
import stationary.component.NontrivialApproximatingComponent;
import stationary.component.NontrivialApproximationSolvingComponent;
import stationary.component.NontrivialSolvingComponent;
import stationary.util.FrequencyRecord;

public final class StationaryDistributionEstimator {
  private enum ComponentMode {
    APPROXIMATION, SOLVING, FULL_APPROXIMATION_REACH, FULL_APPROXIMATION
  }

  private record SamplingTarget(double weight, @Nullable BottomComponent component, Distribution.WeightFunction weights) {}

  private static final Logger logger = Logger.getLogger(StationaryDistributionEstimator.class.getName());
  private static final int MAX_EXPLORES_PER_SAMPLE = 10;

  private final double precision;
  private final Mode mode;

  private final DefaultExplorer<State, MutableDenseSystem> explorer;
  private final DynamicQuotient<MutableDenseSystem> quotient;
  private final int initialState;
  private final List<BottomComponent> components = new ArrayList<>();
  private final IntSet newStatesSinceComponentSearch = new IntOpenHashSet();
  private final Int2ObjectMap<BottomComponent> statesInBottomComponents = new Nat2ObjectDenseArrayMap<>(1024);
  // private final Int2ObjectMap<Int2ObjectMap<Bound>> componentReachability = new Int2ObjectOpenHashMap<>();
  private final List<Int2ObjectMap<Bounds>> componentStateReachability = new ArrayList<>();
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
  private ComponentMode componentMode;

  public StationaryDistributionEstimator(Generator<State> generator, double precision, Mode mode,
      boolean preExplore, boolean solveComponents) {
    this.precision = precision;
    this.mode = mode;
    this.preExplore = preExplore;
    this.componentMode = solveComponents ? ComponentMode.SOLVING : ComponentMode.APPROXIMATION;

    this.explorer = DefaultExplorer.of(generator, SelfLoopHandling.KEEP);
    initialState = explorer.onlyInitialStateId();
    loopStopsUntilCollapse = 10;

    // Fail-fast if these are accessed for a non-transient state
    exitProbability.defaultReturnValue(Double.NaN);
    quotient = new DynamicQuotient<>(explorer.partialSystem(), SelfLoopHandling.INLINE);
  }

  public State getState(int s) {
    return explorer.getState(s);
  }

  private void preExplore() {
    explorer.exploreReachable(explorer.initialStateIds());
    var components = quotient.updateComponents(explorer.initialStateIds(), s -> true);

    if (logger.isLoggable(Level.FINE)) {
      logger.log(Level.FINE, String.format("Found %d BSCCs with %d states",
          components.size(), components.values().stream().mapToInt(Component::size).sum()));
    }
    if (components.size() == 1) {
      if (components.get(0).size() == explorer.exploredStateCount()) {
        logger.log(Level.INFO, "Single SCC model, solving fully");
        componentMode = ComponentMode.FULL_APPROXIMATION;
      } else {
        logger.log(Level.INFO, "Model has a single SCC");
        componentMode = ComponentMode.FULL_APPROXIMATION_REACH;
      }
    }

    components.values().forEach(this::createComponent);
  }

  public Int2ObjectMap<FrequencyRecord> solve() {
    lastProgressUpdate = System.currentTimeMillis();
    newStatesSinceComponentSearch.addAll(explorer.initialStateIds());

    if (preExplore) {
      preExplore();
    }
    while (computeTotalError(initialState) > precision) {
      mainLoopCount += 1L;
      BottomComponent initialComponent = statesInBottomComponents.get(initialState);
      if (initialComponent == null) {
        sample(quotient.onlyInitialState());
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
      Bounds reachability = getComponentReachabilityBounds(component).apply(initialState);
      double error = reachability.upperBound() * component.error();
      component.states().intStream().forEach((int s) ->
          frequency.put(s, new FrequencyRecord(component.frequency(s).lowerBound() * reachability.lowerBound(), error)));
    }

    assert Util.lessOrEqual((1.0 - getAbsorptionLowerBound(initialState)) + components.stream().mapToDouble(component ->
        getComponentReachabilityUpperBound(component).applyAsDouble(initialState) * component.error()).max().orElseThrow(), precision);

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

      int initialState = this.initialState;

      String boundsString;
      if (components.isEmpty()) {
        boundsString = "None";
      } else {
        boundsString = components.stream()
            .sorted(Comparator.comparingDouble(component ->
                -getComponentReachabilityLowerBound(component).applyAsDouble(initialState) * component.error()))
            .limit(10L)
            .map(component -> "%s@[%.5f](%d)".formatted(getComponentReachabilityBounds(component).apply(initialState), component.error(),
                component.getVisitCount()))
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
    double bound = 0.0;
    for (Int2ObjectMap<Bounds> map : componentStateReachability) {
      bound += map.getOrDefault(state, Bounds.reachUnknown()).lowerBound();
    }
    return bound;
  }

  private double computeTotalError(int state) {
    double error;
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      KahanSum exploredWeightedError = new KahanSum();
      KahanSum absorptionProbability = new KahanSum();
      var iterator = componentStateReachability.listIterator();
      while (iterator.hasNext()) {
        var component = components.get(iterator.nextIndex());
        var stateReachability = iterator.next();
        double componentError = component.error();
        Bounds bound = stateReachability.get(state);
        if (bound != null) {
          if (componentError != 0.0) {
            exploredWeightedError.add(componentError * bound.lowerBound());
          }
          absorptionProbability.add(bound.lowerBound());
        }
      }
      assert Util.lessOrEqual(exploredWeightedError.get(), 1.0) : "Got unexpected error: %s".formatted(exploredWeightedError);
      assert Util.lessOrEqual(absorptionProbability.get(), 1.0) : "Got unexpected absorption: %s".formatted(absorptionProbability);
      error = exploredWeightedError.get() + (1.0 - absorptionProbability.get());
    } else {
      error = stateComponent.error();
    }
    assert Util.lessOrEqual(0.0, error) && Util.lessOrEqual(error, 1.0) : "Got unexpected frequency error estimate %.5f".formatted(error);
    return error;
  }

  private IntFunction<Bounds> getComponentReachabilityBounds(BottomComponent component) {
    var reachability = componentStateReachability.get(component.index);
    return state -> {
      BottomComponent stateComponent = statesInBottomComponents.get(state);
      if (stateComponent == null) {
        return reachability.getOrDefault(state, Bounds.reachUnknown());
      }
      return stateComponent.equals(component) ? Bounds.reachOne() : Bounds.reachZero();
    };
  }

  private void updateFromLowerBounds() {
    Int2ObjectMap<KahanSum> stateAbsorptionProbabilities = new Nat2ObjectDenseArrayMap<>(explorer.exploredStateCount());
    for (Int2ObjectMap<Bounds> bounds : componentStateReachability) {
      for (Int2ObjectMap.Entry<Bounds> entry : bounds.int2ObjectEntrySet()) {
        stateAbsorptionProbabilities.computeIfAbsent(entry.getIntKey(), k -> new KahanSum()).add(entry.getValue().lowerBound());
      }
    }

    for (Int2ObjectMap<Bounds> bounds : componentStateReachability) {
      for (Int2ObjectMap.Entry<KahanSum> entry : stateAbsorptionProbabilities.int2ObjectEntrySet()) {
        int state = entry.getIntKey();
        Bounds currentBounds = bounds.getOrDefault(state, Bounds.reachUnknown());
        KahanSum stateAbsorption = entry.getValue();
        double globalUpper = KahanSum.of(1.0, -stateAbsorption.get(), currentBounds.lowerBound());
        if (currentBounds.upperBound() > globalUpper) {
          bounds.put(state, currentBounds.withUpper(globalUpper));
        }
      }
    }
  }

  private IntToDoubleFunction getComponentReachabilityLowerBound(BottomComponent component) {
    var bounds = getComponentReachabilityBounds(component);
    return state -> bounds.apply(state).lowerBound();
  }

  private IntToDoubleFunction getComponentReachabilityUpperBound(BottomComponent component) {
    var bounds = getComponentReachabilityBounds(component);
    return state -> bounds.apply(state).upperBound();
  }

  private BiConsumer<Integer, Bounds> updateComponentReachabilityBound(BottomComponent component) {
    var reachability = componentStateReachability.get(component.index);
    return (state, bound) -> {
      assert !statesInBottomComponents.containsKey(state.intValue());
      reachability.merge(state.intValue(), bound, Bounds::shrink);
      // reachability.put(state.intValue(), bound);
      // assert (oldBounds == null ? Bounds.reachUnknown() : oldBounds).contains(bound, Util.WEAK_EPS)
      //    : "Updating reachability of %d from %s to %s".formatted(state, oldBounds, bound);
    };
  }

  private BottomComponent createComponent(Component component) {
    assert component.stateStream().noneMatch(statesInBottomComponents::containsKey) : "States %s already in component".formatted(component);
    assert SccDecomposition.isBscc(explorer.partialSystem()::successorsIterator, component.states());
    int index = components.size();
    BottomComponent wrapper;
    if (component.size() == 1) {
      wrapper = new AbsorbingComponent(index, component.states().iterator().nextInt());
    } else {
      wrapper = switch (componentMode) {
        case APPROXIMATION -> new NontrivialApproximatingComponent(index, component);
        case SOLVING -> new NontrivialSolvingComponent(index, component);
        case FULL_APPROXIMATION_REACH -> new NontrivialApproximationSolvingComponent(index, component, precision / 2);
        case FULL_APPROXIMATION -> new NontrivialApproximationSolvingComponent(index, component, precision);
      };
    }

    component.states().forEach((IntConsumer) state -> statesInBottomComponents.put(state, wrapper));
    components.add(wrapper);
    exitProbability.keySet().removeAll(component.states());
    for (Int2ObjectMap<Bounds> bounds : componentStateReachability) {
      bounds.keySet().removeAll(component.states());
    }
    componentStateReachability.add(new Nat2ObjectDenseArrayMap<>(1024));
    return wrapper;
  }

  private void sample(int initialState) {
    assert !statesInBottomComponents.containsKey(initialState);

    Distribution.WeightFunction samplingWeight;
    Set<BottomComponent> updateComponentReachability = new HashSet<>();

    if (mode == Mode.SAMPLE_TARGET) {
      List<SamplingTarget> samplingTargets = new ArrayList<>(1 + components.size());

      double exitProbability = getExitProbabilityUpperBound(initialState);
      samplingTargets.add(new SamplingTarget(exitProbability, null, (s, p) -> p * getExitProbabilityUpperBound(s)));

      for (BottomComponent component : components) {
        var componentBounds = getComponentReachabilityBounds(component);
        Bounds reachabilityBounds = componentBounds.apply(initialState);
        double score = reachabilityBounds.upperBound() * component.error() + reachabilityBounds.difference();
        if (!Util.isZero(score)) {
          samplingTargets.add(new SamplingTarget(score, component, (s, p) -> p * componentBounds.apply(s).upperBound()));
        }

        // Heuristic
        if (Sample.random.nextDouble() * reachabilityBounds.difference() > precision) {
          updateComponentReachability.add(component);
        }
      }
      SamplingTarget target = Sample.sampleWeighted(samplingTargets, SamplingTarget::weight).orElseThrow();

      assert samplingTargets.stream().mapToDouble(SamplingTarget::weight).max().orElseThrow() > precision :
          samplingTargets.stream().mapToDouble(SamplingTarget::weight).toString();
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
        if (!Util.isZero(component.error())) {
          component.update(currentState);
        }
        component.countVisit();
        break;
      }

      assert !visitedStateSet.contains(currentState);
      visitedStateSet.add(currentState);
      visitStack.push(currentState);

      // Sample the successor
      Distribution transitions = quotient.onlyDistribution(currentState);

      int nextState = transitions.sampleWeightedExcept(samplingWeight, visitedStateSet::contains);
      if (nextState == -1) {
        checkForComponents = true;
        abortedSamples += 1;
        break;
      }
      assert !visitedStateSet.contains(nextState) : transitions;

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

    if (!preExplore && checkForComponents) {
      loopCount += 1;
      if (loopCount > loopStopsUntilCollapse) {
        loopCount = 0;

        if (!newStatesSinceComponentSearch.isEmpty()) {
          assert newStatesSinceComponentSearch.intStream().noneMatch(this.statesInBottomComponents::containsKey);
          var components = quotient.updateComponents(quotient.initialStates(), explorer::isExploredState);
          newStatesSinceComponentSearch.clear();

          if (logger.isLoggable(Level.FINE)) {
            logger.log(Level.FINE, String.format("Found %d BSCCs with %d states",
                components.size(), components.values().stream().mapToInt(Component::size).sum()));
          }

          Int2ObjectMap<KahanSum> stateAbsorptionProbabilities = new Int2ObjectOpenHashMap<>(explorer.exploredStateCount());
          for (Int2ObjectMap<Bounds> bounds : componentStateReachability) {
            for (Int2ObjectMap.Entry<Bounds> entry : bounds.int2ObjectEntrySet()) {
              double lowerBound = entry.getValue().lowerBound();
              if (lowerBound > 0.0) {
                stateAbsorptionProbabilities.computeIfAbsent(entry.getIntKey(), k -> new KahanSum()).add(lowerBound);
              }
            }
          }

          for (Component component : components.values()) {
            BottomComponent bottomComponent = createComponent(component);
            visitedStateSet.removeAll(component.states());

            var componentBounds = componentStateReachability.get(bottomComponent.index);
            for (Int2ObjectMap.Entry<KahanSum> entry : stateAbsorptionProbabilities.int2ObjectEntrySet()) {
              componentBounds.put(entry.getIntKey(), Bounds.reach(0.0, KahanSum.of(1.0, -entry.getValue().get())));
            }
          }
        }
        loopStopsUntilCollapse = explorer.exploredStates().size();
      }
    }

    // Propagate values backwards along the path
    boolean updateExitProbability = !preExplore;
    Collection<ComponentUpdate> componentReachabilityUpdates = new ArrayList<>(updateComponentReachability.size());
    for (BottomComponent component : updateComponentReachability) {
      componentReachabilityUpdates.add(new ComponentUpdate(component,
          getComponentReachabilityBounds(component),
          updateComponentReachabilityBound(component)));
    }

    while (!visitStack.isEmpty() && (updateExitProbability || !componentReachabilityUpdates.isEmpty())) {
      int state = visitStack.popInt();

      if (!visitedStateSet.contains(state)) {
        continue;
      }
      assert !statesInBottomComponents.containsKey(state);
      Distribution transitions = quotient.onlyDistribution(state);

      if (updateExitProbability) {
        double exitProbability = transitions.sumWeighted(this::getExitProbabilityUpperBound);
        if (Util.isOne(exitProbability)) {
          updateExitProbability = false;
        } else {
          assert Util.lessOrEqual(0.0, exitProbability) : exitProbability;
          double previous = this.exitProbability.put(state, exitProbability);
          assert Double.isNaN(previous) || Util.lessOrEqual(exitProbability, previous) :
              "Updating exit bound of %d from %.5g to %.5g".formatted(state, previous, exitProbability);
        }
      }

      Iterator<ComponentUpdate> iterator = componentReachabilityUpdates.iterator();
      while (iterator.hasNext()) {
        computedTransientUpdates += 1L;

        ComponentUpdate component = iterator.next();
        Bounds bounds = transitions.sumWeightedBounds(component.reachabilityBounds);
        assert !bounds.isNaN();
        if (bounds.equalsUpTo(Bounds.reachUnknown())) {
          iterator.remove();
        } else {
          component.updateFunction.accept(state, bounds);
        }
      }
    }
  }

  public enum Mode {
    SAMPLE_NAIVE, SAMPLE_TARGET
  }

  private record ComponentUpdate(BottomComponent component, IntFunction<Bounds> reachabilityBounds,
                                 BiConsumer<Integer, Bounds> updateFunction) {}
}
