package stationary;

import de.tum.in.naturals.Indices;
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
import it.unimi.dsi.fastutil.HashCommon;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMaps;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayFIFOQueue;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntHeapPriorityQueue;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntPriorityQueue;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntSets;
import it.unimi.dsi.fastutil.ints.IntStack;
import java.time.Duration;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;
import java.util.Locale;
import java.util.Set;
import java.util.concurrent.TimeUnit;
import java.util.function.IntConsumer;
import java.util.function.IntFunction;
import java.util.function.IntUnaryOperator;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import jeigen.DenseMatrix;
import parser.State;

public final class StationaryDistributionEstimator {
  private static final Logger logger = Logger.getLogger(StationaryDistributionEstimator.class.getName());
  private static final int MAX_EXPLORES_PER_SAMPLE = 10;

  private static final Bound UNKNOWN = new Bound(0.0, 1.0);
  private static final Bound ONE = new Bound(1.0, 1.0);
  private static final Bound ZERO = new Bound(0.0, 0.0);
  private static final double FREQUENCY_PRECISION_CHECK = 1.0e-10;

  private final double precision;
  private final MarkovChain model;
  private final DefaultExplorer<State, MarkovChain> explorer;
  private final List<BottomComponent> components = new ArrayList<>();
  private final IntSet newStatesSinceComponentSearch = new IntOpenHashSet();
  private final Int2ObjectMap<BottomComponent> statesInBottomComponents = new Int2ObjectOpenHashMap<>();
  private final Int2ObjectMap<Int2ObjectMap<Bound>> componentReachability = new Int2ObjectOpenHashMap<>();
  private final Int2DoubleMap propagatedErrorBound = new Int2DoubleOpenHashMap();
  private final Int2DoubleMap exitProbability = new Int2DoubleOpenHashMap();
  private int loopCount = 0;
  private int loopStopsUntilCollapse;
  private long sampleRunCount = 0L;
  private long sampledStatesCount = 0L;
  private long computedComponentUpdates = 0L;
  private long computedTransientUpdates = 0L;
  private long lastProgressUpdate = 0L;
  private long lastProgressUpdateSamples = 0L;
  private SamplingMode samplingMode = SamplingMode.TARGET;

  public StationaryDistributionEstimator(DtmcGenerator generator, double precision) {
    this.precision = precision;
    this.model = new MarkovChain();
    this.explorer = DefaultExplorer.of(model, generator, SelfLoopHandling.KEEP);
    loopStopsUntilCollapse = 10;

    // Fail-fast if these are accessed for a non-transient state
    exitProbability.defaultReturnValue(Double.NaN);
    propagatedErrorBound.defaultReturnValue(Double.NaN);
  }

  private static boolean checkFrequency(Int2DoubleMap frequency, IntFunction<Distribution> successors, double precision) {
    Int2DoubleMap incoming = new Int2DoubleOpenHashMap();
    frequency.int2DoubleEntrySet().forEach(entry -> {
      Distribution distribution = successors.apply(entry.getIntKey());
      double freq = entry.getDoubleValue();
      distribution.forEach((t, p) -> incoming.mergeDouble(t, p * freq, Double::sum));
    });
    for (Int2DoubleMap.Entry entry : incoming.int2DoubleEntrySet()) {
      int state = entry.getIntKey();
      double stateFrequency = frequency.get(state);
      double incomingFrequency = entry.getDoubleValue();
      assert Util.lessOrEqual(Math.abs(incomingFrequency - stateFrequency), precision) :
          "Frequency mismatch in state %d: Got %.5g, expected %.5g".formatted(state, incomingFrequency, stateFrequency);
    }
    return true;
  }

  private static Bound bound(double lower, double upper) {
    assert Util.lessOrEqual(lower, upper);
    return new Bound(lower, upper);
  }

  private void preExplore() {
    int initialState = model.getFirstInitialState();
    IntPriorityQueue queue = new IntArrayFIFOQueue();
    queue.enqueue(initialState);
    while (!queue.isEmpty()) {
      int state = queue.dequeueInt();
      for (Distribution choice : explorer.getChoices(state)) {
        for (Integer s : choice.support()) {
          if (!explorer.isExploredState(s)) {
            explorer.exploreState(s);
            queue.enqueue(s.intValue());
          }
        }
      }
    }
    List<NatBitSet> components = SccDecomposition.computeSccs(this.model::getSuccessors, IntSet.of(initialState),
        s -> this.explorer.isExploredState(s) && !statesInBottomComponents.containsKey(s), false);
    if (logger.isLoggable(Level.FINE)) {
      logger.log(Level.FINE, String.format("Found %d components with %d states",
          components.size(), components.stream().mapToInt(NatBitSet::size).sum()));
    }
    components.stream().filter(component -> SccDecomposition.isBscc(this.model::getSuccessors, component))
        .forEach(this::createComponent);
  }

  public void solve() {
    int initialState = model.getFirstInitialState();
    lastProgressUpdate = System.currentTimeMillis();
    newStatesSinceComponentSearch.add(initialState);

    while (computeTotalError(initialState) > precision) {
      sampleRunCount += 1L;
      sample(initialState);
      logProgress(false);
    }
    logProgress(true);
    /*
    logger.log(Level.INFO, "Switching to reachability iteration");
    solveReachabilityWithIteration();
    logProgress(true);
     */

    Int2ObjectMap<FrequencyRecord> frequency = new Int2ObjectOpenHashMap<>(statesInBottomComponents.size());
    for (BottomComponent component : components) {
      Int2DoubleMap componentFrequency = component.frequency();
      double reachability = getComponentReachabilityLowerBound(component, initialState);
      componentFrequency.int2DoubleEntrySet().forEach(entry ->
          frequency.put(entry.getIntKey(), new FrequencyRecord(entry.getDoubleValue() * reachability, component.error())));
    }

    assert Util.lessOrEqual((1.0 - getAbsorptionLowerBound(initialState)) + components.stream().mapToDouble(component ->
        getComponentReachabilityUpperBound(component, initialState) * component.error()).max().orElseThrow(), precision);

    frequency.int2ObjectEntrySet().stream().sorted(Comparator
            .comparingDouble((Int2ObjectMap.Entry<FrequencyRecord> e) -> e.getValue().frequency()).reversed()
            .thenComparing(Int2ObjectMap.Entry::getIntKey))
        .forEach(entry -> System.out.printf("%d: %.6g (+ %.5g)%n", entry.getIntKey(),
            entry.getValue().frequency(), entry.getValue().error()));
    System.out.printf("%nRemaining: %.5g%n", 1.0 - frequency.values().stream().mapToDouble(FrequencyRecord::frequency).sum());
  }

  private void solveReachabilityWithIteration() {
    int initialState = model.getFirstInitialState();

    List<BottomComponent> orderedComponents = new ArrayList<>(this.components);
    /* orderedComponents.sort(Comparator.comparingDouble((BottomComponent c) ->
        getComponentReachabilityLowerBound(c, initialState))); */
    orderedComponents.sort(Comparator.comparingInt(BottomComponent::getVisitCount).reversed()
        .thenComparingDouble((BottomComponent c) -> getComponentReachabilityLowerBound(c, initialState)));

    int states = model.getNumStates();

    IntSet[] predecessors = new IntSet[states];
    Arrays.setAll(predecessors, k -> new IntOpenHashSet());
    explorer.exploredStates().forEach((int s) -> {
      if (!statesInBottomComponents.containsKey(s)) {
        model.getTransitions(s).support().forEach((int t) -> predecessors[t].add(s));
      }
    });
    for (int i = 0; i < predecessors.length; i++) {
      if (predecessors[i].isEmpty()) {
        predecessors[i] = IntSet.of();
      }
    }

    IntStack maybeOutsideDfs = new IntArrayList();
    IntStream.range(0, states).filter(s -> !explorer.isExploredState(s)).forEach(maybeOutsideDfs::push);
    IntSet mayReachOutsideStates = new IntOpenHashSet();
    while (!maybeOutsideDfs.isEmpty()) {
      predecessors[maybeOutsideDfs.popInt()].forEach((int s) -> {
        if (mayReachOutsideStates.add(s)) {
          maybeOutsideDfs.push(s);
        }
      });
    }
    assert mayReachOutsideStates.intStream().mapToObj(s -> predecessors[s]).allMatch(mayReachOutsideStates::containsAll);

    double[] lower = new double[states];
    double[] newLower = new double[states];
    double[] upper = new double[states];
    double[] newUpper = new double[states];
    int iteratedComponents = 0;

    for (BottomComponent component : orderedComponents) {
      iteratedComponents += 1;

      IntSet maybeStates = new IntOpenHashSet(mayReachOutsideStates);
      IntStack maybeComponentDfs = new IntArrayList(component.getStates());
      while (!maybeComponentDfs.isEmpty()) {
        predecessors[maybeComponentDfs.popInt()].forEach((int t) -> {
          if (maybeStates.add(t)) {
            maybeComponentDfs.push(t);
          }
        });
      }

      IntSet yesStates = new IntOpenHashSet(component.getStates());
      IntSet yesStatesBfs = new IntOpenHashSet(yesStates);
      while (!yesStatesBfs.isEmpty()) {
        IntSet nextLayer = new IntOpenHashSet(yesStatesBfs.size());
        yesStatesBfs.forEach((int s) -> {
          assert explorer.isExploredState(s);
          assert component.contains(s) || maybeStates.contains(s);
          if (!yesStates.contains(s) && yesStates.containsAll(model.getTransitions(s).support())) {
            if (yesStates.add(s)) {
              nextLayer.addAll(predecessors[s]);
            }
          }
        });
        yesStatesBfs.clear();
        yesStatesBfs.addAll(nextLayer);
      }
      maybeStates.removeAll(yesStates);

      assert yesStates.intStream().allMatch(s -> component.contains(s) || maybeStates.contains(s));
      assert maybeStates.contains(initialState);

      if (logger.isLoggable(Level.FINE)) {
        logger.log(Level.FINE, "Iterating component %s (%d of %d): %d yes, %d maybe (of %d), remaining error %.5g, lower bound %.4g"
            .formatted(component, iteratedComponents, orderedComponents.size(), yesStates.size(), maybeStates.size(),
                explorer.exploredStateCount(), computeTotalError(initialState),
                getComponentReachabilityLowerBound(component, initialState)));
      }

      for (int s = 0; s < states; s++) {
        Bound bound;
        if (yesStates.contains(s)) {
          bound = ONE;
          if (!statesInBottomComponents.containsKey(s)) {
            updateComponentReachabilityBound(component, s, ONE);
          }
        } else if (maybeStates.contains(s)) {
          bound = getComponentReachabilityBounds(component, s);
        } else if (explorer.isExploredState(s)) {
          bound = ZERO;
          assert getComponentReachabilityLowerBound(component, s) == 0.0;
        } else {
          bound = UNKNOWN;
        }
        lower[s] = bound.lower();
        upper[s] = bound.upper();
        newLower[s] = bound.lower();
        newUpper[s] = bound.upper();
      }

      // TODO Adaptive precision bounds
      double componentPrecision = (precision / 2.0 + (1.0 - component.error()) * precision / 2.0) / components.size();
      while (upper[initialState] - lower[initialState] > componentPrecision) {
        double[] currentLower = lower;
        double[] currentUpper = upper;
        double[] nextLower = newLower;
        double[] nextUpper = newUpper;
        explorer.exploredStates().forEach((int s) -> {
          if (statesInBottomComponents.containsKey(s)) {
            return;
          }
          double[] result = model.getTransitions(s).sumWeightedExceptJacobi(new double[][] {currentLower, currentUpper}, s);
          if (result == null) {
            nextLower[s] = currentLower[s];
            nextUpper[s] = currentUpper[s];
          } else {
            assert !Double.isNaN(result[0]) && !Double.isNaN(result[1])
                : "Got NaN result %s %s in state %d".formatted(result[0], result[1], s);
            assert Util.lessOrEqual(currentLower[s], result[0]) : "Lower bound from %.5g to %.5g".formatted(currentLower[s], result[0]);
            assert Util.lessOrEqual(result[1], currentUpper[s]) : "Upper bound from %.5g to %.5g".formatted(currentUpper[s], result[1]);
            nextLower[s] = result[0];
            nextUpper[s] = result[1];
          }
        });
        double[] lowerSwap = lower;
        lower = newLower;
        newLower = lowerSwap;

        double[] upperSwap = upper;
        upper = newUpper;
        newUpper = upperSwap;
      }
      assert lower[initialState] < upper[initialState] && upper[initialState] - lower[initialState] < precision;

      double[] computedLowerBound = lower;
      double[] computedUpperBound = upper;
      maybeStates.forEach((int s) -> updateComponentReachabilityBound(component, s, bound(computedLowerBound[s], computedUpperBound[s])));

      if (logger.isLoggable(Level.FINE)) {
        logger.log(Level.FINE, "Finished component, remaining error %.5g, bound %s"
            .formatted(computeTotalError(initialState), getComponentReachabilityBounds(component, initialState)));
      }
      if (computeTotalError(initialState) < precision) {
        break;
      }
    }
    assert Util.lessOrEqual(getPropagatedErrorBound(initialState), getExitProbabilityUpperBound(initialState));
    assert Util.lessOrEqual(computeTotalError(initialState), precision);
    logger.log(Level.FINE, String.format("Iterated %d of %d components", iteratedComponents, orderedComponents.size()));
  }

  private void logProgress(boolean force) {
    if (logger.isLoggable(Level.INFO)) {
      long now = System.currentTimeMillis();
      long timeDelta = now - this.lastProgressUpdate;
      if (!force && timeDelta < TimeUnit.SECONDS.toMillis(5L)) {
        return;
      }
      long samples = sampledStatesCount - this.lastProgressUpdateSamples;
      this.lastProgressUpdate = now;
      this.lastProgressUpdateSamples = sampledStatesCount;

      int initialState = model.getFirstInitialState();

      logger.log(Level.INFO, String.format("Progress Report:%n"
              + "    %d sample runs, %d samples (%d/second)%n"
              + "    loop condition %.5g, exit probability %.5g in initial state%n"
              + "    %d explored states, %d states in BSSCs (%d absorbing)%n"
              + "    %d component updates, %d transient state updates%n"
              + "    Bounds: %s",
          sampleRunCount, sampledStatesCount, samples * TimeUnit.SECONDS.toMillis(1L) / timeDelta,
          computeTotalError(initialState), getExitProbabilityUpperBound(initialState),
          explorer.exploredStateCount(), statesInBottomComponents.size(),
          statesInBottomComponents.values().stream().map(BottomComponent::getStates).mapToInt(IntSet::size).filter(i -> i == 1).count(),
          computedComponentUpdates, computedTransientUpdates,
          components.stream()
              .sorted(Comparator.comparingDouble(component -> getComponentReachabilityLowerBound(component, initialState)))
              .limit(10L)
              .map(component -> getComponentReachabilityBounds(component, initialState).toString() + "@" + component.error())
              .collect(Collectors.joining(", ")) + (components.size() > 10 ? " ... (" + components.size() + ")" : "")
      ));
    }
  }

  private double getExitProbabilityUpperBound(int state) {
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

  private double getPropagatedErrorBound(int state) {
    BottomComponent component = statesInBottomComponents.get(state);
    assert !propagatedErrorBound.containsKey(state) || component == null;
    double bound = component == null
        ? propagatedErrorBound.getOrDefault(state, 1.0)
        : component.error();
    assert Util.lessOrEqual(0.0, bound) && Util.lessOrEqual(bound, 1.0);
    //noinspection FloatingPointEquality
    assert explorer.isExploredState(state) || bound == 1.0;
    return bound;
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
    return stateComponent.equals(component) ? ONE : ZERO;
  }

  private Bound computeComponentReachabilityBoundsTransient(BottomComponent component, int state) {
    assert !statesInBottomComponents.containsKey(state);
    Int2ObjectMap<Bound> reachability = componentReachability.get(state);
    return reachability == null ? UNKNOWN : reachability.getOrDefault(component.index, UNKNOWN);
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
    assert reachability.getOrDefault(component.index, UNKNOWN).contains(bound)
        : "Updating reachability of %d from %s to %s".formatted(state, reachability.get(component.index), bound);
    reachability.put(component.index, bound);
  }

  private double getComponentWeightedError(BottomComponent component, int state) {
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      return getComponentWeightedErrorTransient(component, state);
    }
    return stateComponent.equals(component) ? component.error() : 0.0;
  }

  private double getComponentWeightedErrorTransient(BottomComponent component, int state) {
    assert !statesInBottomComponents.containsKey(state);
    Bound bound = computeComponentReachabilityBoundsTransient(component, state);
    return bound.upper() * component.error() + bound.difference();
  }

  private void createComponent(NatBitSet states) {
    assert states.intStream().noneMatch(statesInBottomComponents::containsKey) : "States %s already in component".formatted(states);
    assert SccDecomposition.isBscc(this.model::getSuccessors, states);
    int index = components.size();
    BottomComponent component = states.size() == 1
        ? new AbsorbingComponent(index, states.firstInt())
        : new NontrivialApproximatingComponent(index, states, model::getTransitions);
    states.forEach((IntConsumer) state -> statesInBottomComponents.put(state, component));
    components.add(component);
    exitProbability.keySet().removeAll(states);
    propagatedErrorBound.keySet().removeAll(states);
    componentReachability.keySet().removeAll(states);
  }

  private void sample(int initialState) {
    if (statesInBottomComponents.containsKey(initialState)) {
      BottomComponent component = statesInBottomComponents.get(initialState);
      component.update(initialState);
      return;
    }

    Distribution.WeightFunction samplingWeight;
    Set<BottomComponent> updateComponentReachability = new HashSet<>(2);

    if (samplingMode == SamplingMode.TARGET) {
      if (components.isEmpty()) {
        samplingWeight = (s, p) -> p * getExitProbabilityUpperBound(s);
      } else {
        double[] selectionScores = new double[2 * components.size()];
        int index = 0;
        for (BottomComponent component : components) {
          Bound componentBounds = getComponentReachabilityBounds(component, initialState);
          selectionScores[index] = componentBounds.upper() * component.error();
          selectionScores[index + 1] = componentBounds.difference();
          index += 2;
        }
        int selection = Sample.sample(selectionScores);
        BottomComponent target = components.get(selection / 2);
        updateComponentReachability.add(target);
        if (selection % 2 == 0) {
          samplingWeight = (s, p) -> p * getComponentReachabilityUpperBound(target, s);
        } else {
          samplingWeight = (s, p) -> p * getComponentReachabilityBounds(target, s).difference();
        }
      }
    } else if (samplingMode == SamplingMode.NAIVE) {
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

      if (visitedStateSet.add(currentState)) {
        visitStack.push(currentState);
      } else {
        break;
      }

      // Sample the successor
      Distribution transitions = model.getTransitions(currentState);
      if (transitions == null) { // Deadlock state
        createComponent(NatBitSets.singleton(currentState));
        visitedStateSet.remove(currentState);
        newStatesSinceComponentSearch.remove(currentState);
        break;
      }

      int nextState = transitions.sampleWeightedFiltered(samplingWeight, s -> !visitedStates.contains(s));
      if (nextState == -1) {
        checkForComponents = !transitions.isEmpty();
        break;
      }

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
        loopStopsUntilCollapse += explorer.exploredStates().size();
      }
    }

    // Propagate values backwards along the path
    boolean updateExitProbability = true;
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
          if ((Double.isNaN(lowerBound) && Double.isNaN(upperBound)) || (Util.isZero(lowerBound) && Util.isOne(upperBound))) {
            iterator.remove();
          } else {
            updateComponentReachabilityBound(component, state, bound(lowerBound, upperBound));
          }
        }
      }
    }
  }

  public enum SamplingMode {
    NAIVE, TARGET
  }

  private abstract static class BottomComponent {
    public final int index;
    private int visits = 0;

    public BottomComponent(int index) {
      this.index = index;
    }

    public abstract IntSet getStates();

    public abstract void update(int state);

    public abstract double error();

    public abstract Int2DoubleMap frequency();

    public boolean contains(int state) {
      return getStates().contains(state);
    }

    public void countVisit() {
      visits += 1;
    }

    public int getVisitCount() {
      return visits;
    }

    public int size() {
      return getStates().size();
    }

    @Override
    public boolean equals(Object o) {
      if (this == o) {
        return true;
      }
      if (!(o instanceof BottomComponent)) {
        return false;
      }
      return index == ((BottomComponent) o).index;
    }

    @Override
    public int hashCode() {
      return HashCommon.murmurHash3(index);
    }
  }

  private static final class AbsorbingComponent extends BottomComponent {
    private final int state;

    public AbsorbingComponent(int index, int state) {
      super(index);
      this.state = state;
    }

    @Override
    public IntSet getStates() {
      return IntSet.of(state);
    }

    @Override
    public boolean contains(int state) {
      return state == this.state;
    }

    @Override
    public void update(int state) {
      assert state == this.state;
      // empty
    }

    @Override
    public double error() {
      return 0.0;
    }

    @Override
    public Int2DoubleMap frequency() {
      return Int2DoubleMaps.singleton(state, 1.0);
    }

    @Override
    public int size() {
      return 1;
    }

    @Override
    public String toString() {
      return "{%d}[%d]".formatted(state, index);
    }
  }

  private abstract static class NontrivialComponent extends BottomComponent {
    final NatBitSet states;
    final IntFunction<Distribution> successors;

    private long updateTime = 0L;

    public NontrivialComponent(int index, NatBitSet states, IntFunction<Distribution> successors) {
      super(index);
      this.states = states;
      this.successors = successors;
    }

    @Override
    public IntSet getStates() {
      return IntSets.unmodifiable(states);
    }

    @Override
    public boolean contains(int state) {
      return states.contains(state);
    }

    @Override
    public final void update(int state) {
      long start = System.currentTimeMillis();
      doUpdate(state);
      updateTime += System.currentTimeMillis() - start;
    }

    protected abstract void doUpdate(int state);

    @Override
    public int size() {
      return states.size();
    }

    public Duration updateTime() {
      return Duration.ofMillis(updateTime);
    }

    @Override
    public String toString() {
      return "%s[%d]".formatted(states.size() < 10 ? states.toString() : "<%d>".formatted(states.size()), index);
    }
  }

  private static final class NontrivialSolvingComponent extends NontrivialComponent {
    private final Int2DoubleMap frequency;

    public NontrivialSolvingComponent(int index, NatBitSet component, IntFunction<Distribution> successors) {
      super(index, component, successors);

      frequency = new Int2DoubleOpenHashMap(component.size());
    }

    @Override
    protected void doUpdate(int state) {
      if (frequency.isEmpty()) {
        int size = states.size();
        DenseMatrix matrix = new DenseMatrix(size, size);
        {
          IntUnaryOperator stateToIndexMap = Indices.elementToIndexMap(states);
          IntIterator iterator = states.iterator();
          int index = 0;
          while (iterator.hasNext()) {
            int s = iterator.nextInt();
            Distribution distribution = successors.apply(s);
            int sIndex = index;
            distribution.forEach((t, p) -> matrix.set(stateToIndexMap.applyAsInt(t), sIndex, p));
            index += 1;
          }
        }

        DenseMatrix.EigenResult eig = matrix.eig();
        DenseMatrix real = eig.values.real();
        int maximumIndex = 0;
        double maximum = 0.0;
        for (int i = 0; i < size; i++) {
          double v = real.get(i, 0);
          if (maximum < v) {
            maximum = v;
            maximumIndex = i;
          }
        }
        assert Util.isOne(maximum) : "Expected maximal vector 1, got %5g".formatted(maximum);

        double[] steadyState = eig.vectors.abs().col(maximumIndex).getValues();
        double sum = Util.kahanSum(steadyState);

        IntIterator iterator = states.iterator();
        for (int index = 0; index < size; index++) {
          frequency.put(iterator.nextInt(), steadyState[index] / sum);
        }

        assert checkFrequency(frequency, successors, FREQUENCY_PRECISION_CHECK);
      }
    }

    @Override
    public double error() {
      return frequency.isEmpty() ? 1.0 : 0.0;
    }

    @Override
    public Int2DoubleMap frequency() {
      return Int2DoubleMaps.unmodifiable(frequency);
    }
  }

  private static final class NontrivialApproximatingComponent extends NontrivialComponent {
    private boolean runSampling;

    private final Int2IntMap stateSamplingCounts;
    private int maximalSampleLength;
    private int samplesPerUpdateRun;

    private int iterationsPerUpdate;

    private final Int2DoubleMap stateError;
    private final IntPriorityQueue queue;
    private final Int2ObjectMap<Int2DoubleMap> iteration;
    private final Int2DoubleMap frequencyLowerBounds;

    private double precisionGuide = 1.0;

    private long totalSampledStates = 0L;
    private long totalIterations = 0L;

    public NontrivialApproximatingComponent(int index, NatBitSet component, IntFunction<Distribution> successors) {
      super(index, component, successors);

      int size = component.size();
      this.stateSamplingCounts = new Int2IntOpenHashMap(size);
      component.forEach((int s) -> stateSamplingCounts.put(s, 0));

      iteration = new Int2ObjectOpenHashMap<>(size);
      component.forEach((int s) -> iteration.put(s, new Int2DoubleOpenHashMap()));

      stateError = new Int2DoubleOpenHashMap(size);
      component.forEach((int s) -> stateError.put(s, 1.0));

      frequencyLowerBounds = new Int2DoubleOpenHashMap(size);
      component.forEach((int s) -> frequencyLowerBounds.put(s, 0.0));

      queue = new IntHeapPriorityQueue((int a, int b) -> Double.compare(stateError.get(b), stateError.get(a)));
      component.forEach(queue::enqueue);

      maximalSampleLength = 2 * size;
      samplesPerUpdateRun = 10;
      iterationsPerUpdate = 10;
      runSampling = false;
    }

    @Override
    protected void doUpdate(int initialState) {
      int size = states.size();

      if (runSampling) {
        for (int i = 0; i < samplesPerUpdateRun; i++) {
          int s = Sample.sampleWeighted(states, (int k) -> (double) stateSamplingCounts.get(k) + 1.0);

          stateSamplingCounts.mergeInt(s, 1, Integer::sum);
          totalSampledStates += 1L;
          int sampleLength = Sample.random.nextInt(maximalSampleLength);
          for (int j = 0; j < sampleLength; j++) {
            s = successors.apply(s).sample();
            stateSamplingCounts.mergeInt(s, 1, Integer::sum);
            totalSampledStates += 1L;
          }
        }
        if (totalSampledStates > 10L * (long) size) {
          runSampling = false;
        }
        assert stateSamplingCounts.values().intStream().asLongStream().sum() == totalSampledStates;
        return;
      }

      // Perform iteration

      for (int i = 0; i < 10; i++) {
        if (queue.isEmpty()) {
          break;
        }
        int s = queue.dequeueInt();

        Int2DoubleMap map = iteration.get(s);
        map.defaultReturnValue(0.0);
        if (map.isEmpty() && !stateSamplingCounts.isEmpty()) {
          states.forEach((int t) -> map.put(t, stateSamplingCounts.get(t)));
        }
        Int2DoubleMap current = map;
        Int2DoubleMap next = new Int2DoubleOpenHashMap(size);
        next.defaultReturnValue(0.0);

        double error = Double.NaN;
        double difference = Double.NaN;
        for (int step = 0; step < iterationsPerUpdate; step++) {
          double minimalDifference = Double.MAX_VALUE;
          double maximalDifference = Double.MIN_VALUE;

          IntIterator inner = states.iterator();
          while (inner.hasNext()) {
            int t = inner.nextInt();

            double value = successors.apply(t).sumWeighted(current);
            if (s == t) {
              value += 1.0;
            }
            assert value >= current.get(t);

            double stepDifference;
            if (value > 0.0) {
              next.put(t, value);
              stepDifference = value - current.get(t);
            } else {
              assert current.get(t) == 0.0;
              stepDifference = 0.0;
            }
            if (stepDifference > maximalDifference) {
              maximalDifference = stepDifference;
            }
            if (stepDifference < minimalDifference) {
              minimalDifference = stepDifference;
            }
          }

          difference = minimalDifference;
          error = maximalDifference - minimalDifference;
          totalIterations += 1L;

          Int2DoubleMap swap = next;
          next = current;
          current = swap;
          //noinspection ObjectEquality
          if (error < precisionGuide && current == map) {
            break;
          }
        }

        frequencyLowerBounds.put(s, difference);
        stateError.put(s, error);
        if (!Util.isZero(error)) {
          queue.enqueue(s);
        }
      }
      if (!queue.isEmpty() && stateError.get(queue.firstInt()) < precisionGuide) {
        increasePrecision();
      }
    }

    public void increasePrecision() {
      this.precisionGuide /= 2.0;
    }

    public void setMinimalPrecision(double precision) {
      precisionGuide = Math.min(precisionGuide, precision);
    }

    @Override
    public double error() {
      double error = this.stateError.values().doubleStream().max().orElseThrow();
      assert Util.lessOrEqual(0.0, error) && Util.lessOrEqual(error, 1.0);
      return error;
    }

    @Override
    public Int2DoubleMap frequency() {
      assert checkFrequency(frequencyLowerBounds, successors, error());
      return Int2DoubleMaps.unmodifiable(frequencyLowerBounds);
    }

    public String toInfoString() {
      return ("Component %s: %d samples, %d iterations, %s, %d non-converged, %d below precision %.3g, total update time %s").formatted(
          this, totalSampledStates, totalIterations, runSampling ? "sampling" : "iteration",
          queue.size(), stateError.values().doubleStream().filter(d -> d < precisionGuide).count(), precisionGuide, updateTime()
      );
    }
  }

  private record Bound(double lower, double upper) {
    public double difference() {
      return upper - lower;
    }

    public double average() {
      return (upper + lower) / 2.0;
    }

    @Override
    public String toString() {
      return String.format(Locale.ENGLISH, "[%.5f+%.5f]", lower, difference());
    }

    public boolean contains(Bound bound) {
      return lower() <= bound.lower() && bound.upper() <= upper();
    }
  }

  private record FrequencyRecord(double frequency, double error) {}
}
