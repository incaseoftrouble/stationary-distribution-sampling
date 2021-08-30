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
import java.util.List;
import java.util.Locale;
import java.util.function.IntConsumer;
import java.util.function.IntFunction;
import java.util.function.IntUnaryOperator;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import java.util.stream.IntStream;
import javax.annotation.Nullable;
import jeigen.DenseMatrix;
import parser.State;

public final class StationaryDistributionEstimator {
  private static final int MAX_EXPLORES_PER_SAMPLE = 10;
  private static final Logger logger = Logger.getLogger(StationaryDistributionEstimator.class.getName());

  private static final Bound UNKNOWN = new Bound(0.0, 1.0);
  private static final Bound ONE = new Bound(1.0, 1.0);
  private static final Bound ZERO = new Bound(0.0, 0.0);

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

    public abstract int size();

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

  private static boolean checkFrequency(Int2DoubleMap frequency, IntFunction<Distribution> successors, double precision) {
    Int2DoubleMap incoming = new Int2DoubleOpenHashMap();
    frequency.int2DoubleEntrySet().forEach(entry -> {
      Distribution distribution = successors.apply(entry.getIntKey());
      double f = entry.getDoubleValue();
      distribution.forEach((t, p) -> incoming.mergeDouble(t, p * f, Double::sum));
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

  private abstract static class NontrivialComponent extends BottomComponent {
    final NatBitSet component;
    final IntFunction<Distribution> successors;

    private long updateTime = 0;

    public NontrivialComponent(int index, NatBitSet component, IntFunction<Distribution> successors) {
      super(index);
      this.component = component;
      this.successors = successors;
    }

    @Override
    public IntSet getStates() {
      return IntSets.unmodifiable(component);
    }

    @Override
    public boolean contains(int state) {
      return component.contains(state);
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
      return component.size();
    }

    public Duration updateTime() {
      return Duration.ofMillis(updateTime);
    }

    @Override
    public String toString() {
      return "%s[%d]".formatted(component.size() < 10 ? component.toString() : "<%d>".formatted(component.size()), index);
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
        int size = component.size();
        DenseMatrix matrix = new DenseMatrix(size, size);
        {
          IntUnaryOperator stateToIndexMap = Indices.elementToIndexMap(component);
          IntIterator iterator = component.iterator();
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

        IntIterator iterator = component.iterator();
        for (int index = 0; index < size; index++) {
          frequency.put(iterator.nextInt(), steadyState[index] / sum);
        }

        assert checkFrequency(frequency, successors, 1e-10);
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

    private long totalSampledStates = 0;
    private long totalIterations = 0;

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
      int size = component.size();

      if (runSampling) {
        for (int i = 0; i < samplesPerUpdateRun; i++) {
          int s = Sample.sampleWeighted(component, k -> stateSamplingCounts.get(k) + 1.0);

          stateSamplingCounts.mergeInt(s, 1, Integer::sum);
          totalSampledStates += 1;
          int sampleLength = Sample.random.nextInt(maximalSampleLength);
          for (int j = 0; j < sampleLength; j++) {
            s = successors.apply(s).sample();
            stateSamplingCounts.mergeInt(s, 1, Integer::sum);
            totalSampledStates += 1;
          }
        }
        if (totalSampledStates > 10L * size) {
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
          component.forEach((int t) -> map.put(t, stateSamplingCounts.get(t)));
        }
        Int2DoubleMap current = map;
        Int2DoubleMap next = new Int2DoubleOpenHashMap(size);
        next.defaultReturnValue(0.0);

        double error = Double.NaN;
        double difference = Double.NaN;
        for (int step = 0; step < iterationsPerUpdate; step++) {
          double minimalDifference = Double.MAX_VALUE;
          double maximalDifference = Double.MIN_VALUE;

          IntIterator inner = component.iterator();
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
          totalIterations += 1;

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
      this.precisionGuide /= 2;
    }

    public void setMinimalPrecision(double precision) {
      precisionGuide = Math.min(precisionGuide, precision);
    }

    @Override
    public double error() {
      return this.stateError.values().doubleStream().max().orElseThrow();
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

  private static Bound bound(double lower, double upper) {
    assert Util.lessOrEqual(lower, upper);
    return new Bound(lower, upper);
  }

  private record Bound(double lower, double upper) {
    public double difference() {
      return upper - lower;
    }

    public double average() {
      return (upper + lower) / 2;
    }

    @Override
    public String toString() {
      return String.format(Locale.ENGLISH, "[%.5f+%.3g]", lower, difference());
    }
  }

  private final double precision;

  private final MarkovChain model;
  private final DefaultExplorer<State, MarkovChain> explorer;

  private final IntSet newStatesSinceComponentSearch = new IntOpenHashSet();
  private final Int2ObjectMap<BottomComponent> statesInBottomComponents = new Int2ObjectOpenHashMap<>();
  private final Int2ObjectMap<Int2DoubleMap> transientStateComponentReachability = new Int2ObjectOpenHashMap<>();

  private final Int2DoubleMap transientWeightedErrorBound = new Int2DoubleOpenHashMap();
  private final List<BottomComponent> components = new ArrayList<>();

  private int loopCount = 0;
  private int loopStopsUntilCollapse;

  private long sampleRunCount = 0;
  private long sampledStatesCount = 0;
  private long computedComponentUpdates = 0;
  private long computedTransientUpdates = 0;

  private long lastProgressUpdate = 0;

  public StationaryDistributionEstimator(DtmcGenerator generator, double precision) {
    this.precision = precision;
    this.model = new MarkovChain();
    this.explorer = DefaultExplorer.of(model, generator, SelfLoopHandling.KEEP);
    loopStopsUntilCollapse = 10;

    // Fail-fast if these are accessed for a non-transient state
    transientWeightedErrorBound.defaultReturnValue(Double.NaN);
  }

  private void preExplore() {
    int initialState = model.getFirstInitialState();
    IntArrayFIFOQueue queue = new IntArrayFIFOQueue();
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

  private static record FrequencyRecord(double frequency, double error) {

  }

  public void solve() {
    int initialState = model.getFirstInitialState();
    lastProgressUpdate = System.currentTimeMillis();
    newStatesSinceComponentSearch.add(initialState);

    while (getWeightedErrorBound(initialState) > precision / 8) {
      sample(initialState);
      logProgress(false);
    }
    logProgress(true);
    logger.log(Level.INFO, "Switching to reachability iteration");
    solveReachabilityWithIteration();
    logProgress(true);

    Int2ObjectMap<FrequencyRecord> frequency = new Int2ObjectOpenHashMap<>(statesInBottomComponents.size());

    for (BottomComponent component : components) {
      Int2DoubleMap componentFrequency = component.frequency();
      double reachability = getComponentReachabilityLowerBound(component, initialState);
      componentFrequency.int2DoubleEntrySet().forEach(entry ->
          frequency.put(entry.getIntKey(), new FrequencyRecord(entry.getDoubleValue() * reachability, component.error())));
    }

    frequency.int2ObjectEntrySet().stream().sorted(Comparator
            .comparingDouble((Int2ObjectMap.Entry<FrequencyRecord> e) -> e.getValue().frequency()).reversed()
            .thenComparing(Int2ObjectMap.Entry::getIntKey))
        .forEach(entry -> System.out.printf("%d: %.6g (+ %.3g)%n", entry.getIntKey(),
            entry.getValue().frequency(), entry.getValue().error()));
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
                explorer.exploredStateCount(), getFrequencyErrorBound(initialState),
                getComponentReachabilityLowerBound(component, initialState)));
      }

      for (int s = 0; s < states; s++) {
        Bound bound;
        if (yesStates.contains(s)) {
          bound = ONE;
          if (!statesInBottomComponents.containsKey(s)) {
            updateTransientComponentReachabilityLowerBound(component, s, 1.0);
          }
        } else if (maybeStates.contains(s)) {
          bound = computeComponentReachabilityBounds(component, s);
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
      double componentPrecision = (precision / 2 + (1 - component.error()) * precision / 2) / components.size();
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
      maybeStates.forEach((int s) -> updateTransientComponentReachabilityLowerBound(component, s, computedLowerBound[s]));

      if (logger.isLoggable(Level.FINE)) {
        logger.log(Level.FINE, "Finished component, remaining error %.5g, bound %s"
            .formatted(getFrequencyErrorBound(initialState), computeComponentReachabilityBounds(component, initialState)));
      }
      if (getFrequencyErrorBound(initialState) < precision) {
        break;
      }
    }
    assert Util.lessOrEqual(getWeightedErrorBound(initialState), 1.0 - getAbsorptionLowerBound(initialState));
    assert Util.lessOrEqual(getFrequencyErrorBound(initialState), precision);
    logger.log(Level.FINE, String.format("Iterated %d of %d components", iteratedComponents, orderedComponents.size()));
  }

  private void logProgress(boolean force) {
    if (logger.isLoggable(Level.INFO)) {
      long now = System.currentTimeMillis();
      if (!force && now - this.lastProgressUpdate < 5000) {
        return;
      }
      this.lastProgressUpdate = now;

      int initialState = model.getFirstInitialState();
      logger.log(Level.INFO, String.format("Progress Report:%n"
              + "    %d sample runs, %d samples%n"
              + "    weighted error heuristic %.5g, frequency error %.5g, absorption %.5f in initial state%n"
              + "    %d explored states, %d states in BSSCs (%d absorbing)%n"
              + "    %d component updates, %d transient state updates%n"
              + "    Reachability bounds: %s, Error bounds: %s",
          sampleRunCount, sampledStatesCount,
          getWeightedErrorBound(initialState), getFrequencyErrorBound(initialState), getAbsorptionLowerBound(initialState),
          explorer.exploredStateCount(), statesInBottomComponents.size(),
          statesInBottomComponents.values().stream().map(BottomComponent::getStates).mapToInt(IntSet::size).filter(i -> i == 1).count(),
          computedComponentUpdates, computedTransientUpdates,
          components.size() < 10
              ? components.stream().map(component -> computeComponentReachabilityBounds(component, initialState))
              .sorted(Comparator.comparingDouble(Bound::lower).reversed())
              .map(Bound::toString)
              .collect(Collectors.joining(";", "[", "]"))
              : "[...]",
          components.size() < 10
              ? components.stream().mapToDouble(BottomComponent::error).mapToObj("%.5g"::formatted)
              .collect(Collectors.joining(";", "[", "]"))
              : "[...]"
      ));
    }
  }

  private double getAbsorptionLowerBound(int state) {
    return statesInBottomComponents.containsKey(state) ? 1.0 : getTransientAbsorptionLowerBound(state);
  }

  private double getTransientAbsorptionLowerBound(int state) {
    assert !statesInBottomComponents.containsKey(state);
    Int2DoubleMap reachability = transientStateComponentReachability.get(state);
    return reachability == null ? 0.0 : reachability.values().doubleStream().sum();
  }

  private double getWeightedErrorBound(int state) {
    BottomComponent component = statesInBottomComponents.get(state);
    assert !transientWeightedErrorBound.containsKey(state) || component == null;
    double bound = component == null
        ? transientWeightedErrorBound.getOrDefault(state, 1.0)
        : component.error();
    assert Util.lessOrEqual(0.0, bound) && Util.lessOrEqual(bound, 1.0);
    //noinspection FloatingPointEquality
    assert explorer.isExploredState(state) || bound == 1.0;
    return bound;
  }

  private double getFrequencyErrorBound(int state) {
    double error;
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      Int2DoubleMap reachabilityMap = transientStateComponentReachability.get(state);
      if (reachabilityMap == null) {
        error = 1.0;
      } else {
        KahanSum exploredWeightedError = new KahanSum();
        KahanSum absorptionProbability = new KahanSum();
        for (BottomComponent component : components) {
          double componentError = component.error();
          double reachabilityLower = reachabilityMap.getOrDefault(component.index, 0.0);
          exploredWeightedError.add(componentError * reachabilityLower);
          absorptionProbability.add(reachabilityLower);
        }
        assert Util.lessOrEqual(exploredWeightedError.get(), 1.0) && Util.lessOrEqual(absorptionProbability.get(), 1.0)
            : "Got unexpected error / absorption: %s, %s".formatted(exploredWeightedError, absorptionProbability);
        error = exploredWeightedError.get() + (1.0 - absorptionProbability.get());
      }
    } else {
      error = stateComponent.error();
    }
    assert Util.lessOrEqual(0.0, error) && Util.lessOrEqual(error, 1.0)
        : "Got unexpected frequency error estimate %.5f".formatted(error);
    return error;
  }

  private Bound computeComponentReachabilityBounds(BottomComponent component, int state) {
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      return computeTransientComponentReachabilityBounds(component, state);
    }
    return stateComponent.equals(component) ? ONE : ZERO;
  }

  private Bound computeTransientComponentReachabilityBounds(BottomComponent component, int state) {
    assert !statesInBottomComponents.containsKey(state);
    Int2DoubleMap reachability = transientStateComponentReachability.get(state);
    if (reachability == null) {
      return UNKNOWN;
    }
    double lowerBound = reachability.get(component.index);
    double absorption = reachability.values().doubleStream().sum();
    assert Util.lessOrEqual(absorption, 1.0);
    double upperBound = lowerBound + 1.0 - absorption;
    assert Util.lessOrEqual(lowerBound, upperBound) && Util.lessOrEqual(upperBound, 1.0);
    return bound(lowerBound, upperBound);
  }

  private double getComponentReachabilityLowerBound(BottomComponent component, int state) {
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      return getTransientComponentReachabilityLowerBound(component, state);
    }
    return stateComponent.equals(component) ? 1.0 : 0.0;
  }

  private double getTransientComponentReachabilityLowerBound(BottomComponent component, int state) {
    assert !statesInBottomComponents.containsKey(state);
    Int2DoubleMap reachability = transientStateComponentReachability.get(state);
    return reachability == null ? 0.0 : reachability.getOrDefault(component.index, 0.0);
  }

  private void updateTransientComponentReachabilityLowerBound(BottomComponent component, int state, double lowerBound) {
    assert !statesInBottomComponents.containsKey(state);
    if (Util.isOne(lowerBound)) {
      return;
    }
    Int2DoubleMap reachability = transientStateComponentReachability.computeIfAbsent(state, k -> new Int2DoubleOpenHashMap());
    assert !reachability.containsKey(component.index) || Util.lessOrEqual(reachability.get(component.index), lowerBound)
        : "Updating reachability of %d from %.5g to %.5g".formatted(state, reachability.get(component.index), lowerBound);
    reachability.put(component.index, lowerBound);
  }

  private void createComponent(NatBitSet states) {
    assert states.intStream().noneMatch(statesInBottomComponents::containsKey)
        : "States %s already in component".formatted(states);
    assert SccDecomposition.isBscc(this.model::getSuccessors, states);
    int index = components.size();
    BottomComponent component = states.size() == 1
        ? new AbsorbingComponent(index, states.firstInt())
        : new NontrivialApproximatingComponent(index, states, model::getTransitions);
    states.forEach((IntConsumer) state -> statesInBottomComponents.put(state, component));
    components.add(component);
    transientWeightedErrorBound.keySet().removeAll(states);
    transientStateComponentReachability.keySet().removeAll(states);
  }

  private void sample(int initialState) {
    sampleRunCount += 1;

    @Nullable
    Distribution.WeightFunction samplingWeight;

    IntList visitedStates = new IntArrayList();
    IntStack visitStack = (IntStack) visitedStates;
    IntSet visitedStateSet = new IntOpenHashSet();

    if (components.isEmpty() || Sample.random.nextDouble() >= getWeightedErrorBound(initialState)) {
      if (Sample.random.nextBoolean()) {
        samplingWeight = (s, p) -> visitedStateSet.contains(s) ? 0.0 : getWeightedErrorBound(s);
      } else {
        samplingWeight = (s, p) -> visitedStateSet.contains(s) ? 0.0 : p * getWeightedErrorBound(s);
      }
    } else {
      samplingWeight = (s, p) -> visitedStateSet.contains(s) ? 0.0 : p;
    }

    int exploreDuringSample = 0;
    int currentState = initialState;

    @Nullable
    BottomComponent visitedComponent = null;
    boolean checkForComponents = false;
    while (true) {
      assert explorer.isExploredState(currentState);

      BottomComponent component = statesInBottomComponents.get(currentState);
      if (component != null) {
        visitedComponent = component;
        computedComponentUpdates += 1;
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

      int nextState = samplingWeight == null
          ? transitions.sampleExcept(currentState)
          : transitions.sampleWeighted(samplingWeight);

      if (nextState == -1) {
        checkForComponents = !transitions.isEmpty();
        break;
      }

      sampledStatesCount += 1;
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
    boolean updateErrorBound = true;
    boolean updateComponentReachability = visitedComponent != null;
    while (!visitStack.isEmpty() && (updateErrorBound || updateComponentReachability)) {
      int state = visitStack.popInt();

      if (!visitedStateSet.contains(state)) {
        continue;
      }
      assert !statesInBottomComponents.containsKey(state);
      Distribution transitions = model.getTransitions(state);

      if (updateErrorBound) {
        double errorBound = transitions.sumWeightedExceptJacobi(this::getWeightedErrorBound, state);
        if (Double.isNaN(errorBound) || Util.isOne(errorBound)) {
          updateErrorBound = false;
        } else {
          assert Util.lessOrEqual(0.0, errorBound);
          computedTransientUpdates += 1;
          double previous = transientWeightedErrorBound.put(state, errorBound);
          assert Double.isNaN(previous) || Util.lessOrEqual(errorBound, previous) :
              "Updating weighted bound of %d from %.5g to %.5g".formatted(state, previous, errorBound);
        }
      }

      if (updateComponentReachability) {
        BottomComponent finalVisitedComponent = visitedComponent;
        double lowerBound = transitions.sumWeightedExceptJacobi(s -> getComponentReachabilityLowerBound(finalVisitedComponent, s), state);
        if (Double.isNaN(lowerBound) || Util.isZero(lowerBound)) {
          updateComponentReachability = false;
        } else {
          updateTransientComponentReachabilityLowerBound(visitedComponent, state, lowerBound);
        }
      }
    }
  }
}
