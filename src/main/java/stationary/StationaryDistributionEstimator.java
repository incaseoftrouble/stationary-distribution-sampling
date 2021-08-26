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
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntStack;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Locale;
import java.util.function.IntConsumer;
import java.util.function.IntFunction;
import java.util.function.IntToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.stream.Collectors;
import javax.annotation.Nullable;
import parser.State;

public final class StationaryDistributionEstimator {
  private static final int MAX_EXPLORES_PER_SAMPLE = 10;
  private static final Logger logger = Logger.getLogger(StationaryDistributionEstimator.class.getName());

  private static final Bound UNKNOWN = new Bound(0.0, 1.0);
  private static final Bound ONE = new Bound(1.0, 1.0);
  private static final Bound ZERO = new Bound(0.0, 0.0);

  private abstract static class BottomComponent {
    public final int index;
    private final Int2DoubleMap transientStatesReachabilityLowerBound = new Int2DoubleOpenHashMap();
    private int visitCount = 0;

    public BottomComponent(int index) {
      this.index = index;
      transientStatesReachabilityLowerBound.defaultReturnValue(Double.NaN);
    }

    public void visit() {
      visitCount += 1;
    }

    public int visitCount() {
      return visitCount;
    }

    public abstract IntSet getStates();

    public abstract void update(int state);

    public abstract double error();

    public abstract Int2DoubleMap frequencyLowerBounds();

    public boolean contains(int state) {
      return getStates().contains(state);
    }

    public double getReachabilityLower(int state) {
      return transientStatesReachabilityLowerBound.getOrDefault(state, 0.0);
    }

    public void updateReachabilityLower(int state, double value) {
      double previous = transientStatesReachabilityLowerBound.put(state, value);
      assert Double.isNaN(previous) || Util.lessOrEqual(previous, value) : "%f -> %f".formatted(previous, value);
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
    public void update(int state) {
      assert state == this.state;
      // empty
    }

    @Override
    public double error() {
      return 0.0;
    }

    @Override
    public Int2DoubleMap frequencyLowerBounds() {
      return Int2DoubleMaps.singleton(state, 1.0);
    }

    @Override
    public int size() {
      return 1;
    }

    @Override
    public String toString() {
      return "Component {%d}[%d]".formatted(state, index);
    }
  }

  private static final class NontrivialComponent extends BottomComponent {
    private final NatBitSet component;
    private final IntFunction<Distribution> successors;
    private final Int2IntMap stateSamplingCounts;

    private int maximalSampleLength;
    private int samplesPerUpdateRun;
    private long totalSampledStates = 0;
    private long totalIterationSteps = 0;
    private boolean runSampling = true;
    private boolean useErrorEstimate = true;
    private double validatedError = 1.0;
    private int iterationsPerUpdate = 10;

    private Int2DoubleMap iteration;
    private Int2DoubleMap nextIteration;

    public NontrivialComponent(int index, NatBitSet component, IntFunction<Distribution> successors) {
      super(index);
      this.component = component;
      this.successors = successors;

      int size = component.size();
      this.stateSamplingCounts = new Int2IntOpenHashMap(size);
      component.forEach((int s) -> stateSamplingCounts.put(s, 0));

      iteration = new Int2DoubleOpenHashMap(size);
      nextIteration = new Int2DoubleOpenHashMap(size);

      this.maximalSampleLength = 2 * size;
      this.samplesPerUpdateRun = 10;
    }

    @Override
    public IntSet getStates() {
      return component;
    }

    @Override
    public void update(int initialState) {
      int size = component.size();
      if (runSampling) {
        for (int i = 0; i < samplesPerUpdateRun; i++) {
          int s = Sample.sampleWeighted(component, k -> iteration.get(k) + 1.0);

          stateSamplingCounts.mergeInt(s, 1, Integer::sum);
          totalSampledStates += 1;
          int sampleLength = Sample.random.nextInt(maximalSampleLength);
          for (int j = 0; j < sampleLength; j++) {
            s = successors.apply(s).sample();
            stateSamplingCounts.mergeInt(s, 1, Integer::sum);
            totalSampledStates += 1;
          }
        }
        assert stateSamplingCounts.values().intStream().asLongStream().sum() == totalSampledStates;
      }

      if (iteration.isEmpty()) {
        // Initialize iteration - make sure that every state has non-zero frequency
        IntToDoubleFunction initialValue = stateSamplingCounts.isEmpty()
            ? s -> 1.0 / component.size()
            : s -> (1.0 + stateSamplingCounts.getOrDefault(s, 0)) / (totalSampledStates + component.size());
        component.forEach((int s) -> iteration.put(s, initialValue.applyAsDouble(s)));
      }

      // Perform iteration

      for (int step = 0; step < iterationsPerUpdate; step++) {
        totalIterationSteps += 1;

        double errorEstimate = 0.0;
        IntIterator iterator = component.iterator();
        while (iterator.hasNext()) {
          int s = iterator.nextInt();

          double nextValue = successors.apply(s).sumWeighted(iteration);
          nextIteration.put(s, nextValue);
          if (useErrorEstimate) {
            double currentValue = iteration.get(s);
            errorEstimate = Math.max(errorEstimate, Math.abs(currentValue - nextValue));
          }
        }
        if (useErrorEstimate && errorEstimate < 1e-8) {
          logger.log(Level.INFO, "Switching to correct error estimation");
          useErrorEstimate = false;
          iterationsPerUpdate *= 2;
          break;
        }
        Int2DoubleMap swap = nextIteration;
        nextIteration = iteration;
        iteration = swap;
      }

      if (useErrorEstimate) {
        iterationsPerUpdate += 10;
      } else {
        double minimalFrequency = iteration.values().doubleStream().min().orElseThrow();
        assert minimalFrequency > 0.0;

        IntIterator iterator = component.iterator();
        double maximalError = Double.MIN_VALUE;
        while (iterator.hasNext()) {
          int s = iterator.nextInt();

          IntToDoubleFunction estimate = t -> stateSamplingCounts.get(t);

          double minimalDifference = Double.MAX_VALUE;
          double maximalDifference = Double.MIN_VALUE;

          IntIterator inner = component.iterator();
          while (inner.hasNext()) {
            int t = inner.nextInt();

            double stepDifference = successors.apply(t).sumWeighted(estimate) - estimate.applyAsDouble(t);
            if (s == t) {
              stepDifference += 1.0;
            }
            if (stepDifference > maximalDifference) {
              maximalDifference = stepDifference;
            }
            if (stepDifference < minimalDifference) {
              minimalDifference = stepDifference;
            }
          }
          double errorBound = maximalDifference - minimalDifference;
          if (errorBound > maximalError) {
            maximalError = errorBound;
          }
          if (Util.isOne(maximalError)) {
            break;
          }
        }
        validatedError = maximalError;
        System.out.println(validatedError);
      }
    }

    @Override
    public double error() {
      return validatedError;
    }

    @Override
    public Int2DoubleMap frequencyLowerBounds() {
      Int2DoubleMap frequencies = new Int2DoubleOpenHashMap(component.size());
      // component.forEach(state -> frequencies.put(state, this.stateFrequencyLowerBound(state)));
      return frequencies;
    }

    @Override
    public int size() {
      return component.size();
    }

    @Override
    public String toString() {
      return "Component %s[%d]: %d samples, %d iterations, %s update".formatted(
          component.size() < 10 ? component.toString() : component.size(), index,
          totalSampledStates, totalIterationSteps, runSampling ? "sampling" : "iteration"
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

  private final Int2DoubleMap absorptionLowerBound = new Int2DoubleOpenHashMap();
  private final Int2DoubleMap transientErrorBounds = new Int2DoubleOpenHashMap();
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
    absorptionLowerBound.defaultReturnValue(Double.NaN);
    transientErrorBounds.defaultReturnValue(Double.NaN);
  }

  public void solve() {
    int initialState = model.getFirstInitialState();
    lastProgressUpdate = System.currentTimeMillis();

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

    while (getWeightedErrorBound(initialState) > precision / 4) {
      sample(initialState);
      logProgress(false);
    }
    logProgress(true);
    logger.log(Level.INFO, "Switching to reachability iteration");
    solveWithIteration();
    logProgress(true);
  }

  private void solveWithIteration() {
    int initialState = model.getFirstInitialState();

    List<BottomComponent> components = new ArrayList<>(this.components);
    components.sort(Comparator.comparing(BottomComponent::visitCount).reversed());

    int states = model.getNumStates();
    IntSet transientStates = new IntOpenHashSet(explorer.exploredStateCount() - statesInBottomComponents.size());
    explorer.exploredStates().forEach((int s) -> {
      if (!statesInBottomComponents.containsKey(s)) {
        transientStates.add(s);
      }
    });

    double[] lower = new double[states];
    double[] newLower = new double[states];
    double[] upper = new double[states];
    double[] newUpper = new double[states];
    int iteratedComponents = 0;
    for (BottomComponent component : components) {
      iteratedComponents += 1;

      for (int s = 0; s < states; s++) {
        if (explorer.isExploredState(s)) {
          Bound bounds = getReachabilityBounds(component, s);
          lower[s] = bounds.lower();
          upper[s] = bounds.upper();
          newLower[s] = bounds.lower();
          newUpper[s] = bounds.upper();
        } else {
          lower[s] = 0.0;
          upper[s] = 1.0;
          newLower[s] = 0.0;
          newUpper[s] = 1.0;
        }
      }

      while (upper[initialState] - lower[initialState] > precision / 2) {
        double[] currentLower = lower;
        double[] currentUpper = upper;
        double[] nextLower = newLower;
        double[] nextUpper = newUpper;
        transientStates.forEach((int s) -> {
          double[] result = model.getTransitions(s).sumWeightedExceptJacobi(new double[][] {currentLower, currentUpper}, s);
          if (result == null) {
            nextLower[s] = currentLower[s];
            nextUpper[s] = currentUpper[s];
          } else {
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

      double[] computedLowerBound = lower;
      transientStates.forEach((int s) -> {
        component.updateReachabilityLower(s, computedLowerBound[s]);
        double absorption = 0.0;
        for (BottomComponent bottomComponent : components) {
          absorption += bottomComponent.getReachabilityLower(s);
        }
        absorptionLowerBound.put(s, absorption);
      });
      if (logger.isLoggable(Level.FINER)) {
        logger.log(Level.FINER, "Finished reachability of %s (%d of %d), remaining error %f".formatted(
            component, iteratedComponents, components.size(), getFrequencyError(initialState)));
      }
      if (getFrequencyError(initialState) < precision) {
        break;
      }
    }
    assert Util.lessOrEqual(getFrequencyError(initialState), precision);
    logger.log(Level.FINE, String.format("Iterated %d of %d components", iteratedComponents, components.size()));
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
              + "    weighted error bound %.5f, frequency error %.5f, absorption %.5f in initial state%n"
              + "    %d explored states, %d states in BSSCs (%d absorbing)%n"
              + "    %d component updates, %d transient state updates%n"
              + "    Reachability bounds: %s",
          sampleRunCount, sampledStatesCount,
          getWeightedErrorBound(initialState), getFrequencyError(initialState), getAbsorptionLowerBound(initialState),
          explorer.exploredStateCount(), statesInBottomComponents.size(),
          statesInBottomComponents.values().stream().map(BottomComponent::getStates).mapToInt(IntSet::size).filter(i -> i == 1).count(),
          computedComponentUpdates, computedTransientUpdates,
          components.stream().mapToInt(BottomComponent::size).sum() < 20
              ? components.stream().map(component -> getReachabilityBounds(component, initialState))
              .sorted(Comparator.comparingDouble(Bound::lower).reversed())
              .map(Bound::toString).collect(Collectors.joining(";", "[", "]"))
              : "[...]"
      ));
    }
  }

  private double getAbsorptionLowerBound(int state) {
    return statesInBottomComponents.containsKey(state)
        ? 1.0
        : absorptionLowerBound.getOrDefault(state, 0.0);
  }

  private double getWeightedErrorBound(int state) {
    BottomComponent component = statesInBottomComponents.get(state);
    double bound = component == null
        ? transientErrorBounds.getOrDefault(state, 1.0)
        : component.error();
    assert bound >= 0.0;
    return bound;
  }

  private double getFrequencyError(int state) {
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      double exploredMinimalError = 0.0;
      double maximalError = 0.0;
      for (BottomComponent component : components) {
        double error = component.error();
        exploredMinimalError += error * component.getReachabilityLower(state);
        if (maximalError < error) {
          maximalError = error;
        }
      }
      double absorptionProbability = absorptionLowerBound.getOrDefault(state, 0.0);
      return exploredMinimalError + (1.0 - absorptionProbability) * (1.0 + maximalError);
    }
    return stateComponent.error();
  }

  private Bound getReachabilityBounds(BottomComponent component, int state) {
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      double lowerBound = component.getReachabilityLower(state);
      double upperBound = lowerBound + 1.0 - absorptionLowerBound.getOrDefault(state, 0.0);
      assert Util.isEqual(absorptionLowerBound.getOrDefault(state, 0.0),
          components.stream().mapToDouble(c -> c.getReachabilityLower(state)).sum());
      return bound(lowerBound, upperBound);
    }
    return stateComponent.equals(component) ? ONE : ZERO;
  }

  private double getReachabilityLowerBounds(BottomComponent component, int state) {
    BottomComponent stateComponent = statesInBottomComponents.get(state);
    if (stateComponent == null) {
      return component.getReachabilityLower(state);
    }
    return stateComponent.equals(component) ? 1.0 : 0.0;
  }

  private void createComponent(NatBitSet states) {
    assert states.intStream().noneMatch(statesInBottomComponents::containsKey)
        : "States %s already in component".formatted(states);
    assert SccDecomposition.isBscc(this.model::getSuccessors, states);
    int index = components.size();
    BottomComponent component = states.size() == 1
        ? new AbsorbingComponent(index, states.firstInt())
        : new NontrivialComponent(index, states, model::getTransitions);
    states.forEach((IntConsumer) state -> statesInBottomComponents.put(state, component));
    components.add(component);
    transientErrorBounds.keySet().removeAll(states);
    absorptionLowerBound.keySet().removeAll(states);
  }

  private void sample(int initialState) {
    sampleRunCount += 1;

    @Nullable
    Distribution.WeightFunction samplingWeight;
    boolean countVisits;
    if (components.isEmpty() || Sample.random.nextDouble() >= getWeightedErrorBound(initialState)) {
      samplingWeight = (s, p) -> getWeightedErrorBound(s);
      countVisits = false;
    } else {
      samplingWeight = null;
      countVisits = true;
    }

    IntList visitedStates = new IntArrayList();
    IntStack visitStack = (IntStack) visitedStates;
    IntSet visitedStateSet = new IntOpenHashSet();

    int exploreDuringSample = 0;
    int currentState = initialState;
    int stateRevisit = 0;

    @Nullable
    BottomComponent visitedComponent = null;
    while (stateRevisit < 10) {
      assert explorer.isExploredState(currentState);

      BottomComponent component = statesInBottomComponents.get(currentState);
      if (component != null) {
        visitedComponent = component;
        if (countVisits) {
          component.visit();
        }
        computedComponentUpdates += 1;
        component.update(currentState);
        break;
      }

      visitStack.push(currentState);
      if (!visitedStateSet.add(currentState)) {
        stateRevisit += 1;
      }

      // Sample the successor
      Distribution transitions = model.getTransitions(currentState);
      if (transitions == null) { // Deadlock state
        createComponent(NatBitSets.singleton(currentState));
        visitedStateSet.remove(currentState);
        break;
      }

      int nextState = samplingWeight == null
          ? transitions.sample()
          : transitions.sampleWeighted(samplingWeight);

      if (nextState == -1) {
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

    if (stateRevisit > 5) {
      loopCount += 1;
      if (loopCount > loopStopsUntilCollapse) {
        loopCount = 0;

        if (newStatesSinceComponentSearch.isEmpty()) {
          logger.fine("Skipping SCC computation");
        } else {
          assert newStatesSinceComponentSearch.intStream().noneMatch(this.statesInBottomComponents::containsKey);
          List<NatBitSet> components = SccDecomposition.computeSccs(this.model::getSuccessors, newStatesSinceComponentSearch,
              s -> this.explorer.isExploredState(s) && !statesInBottomComponents.containsKey(s), false);
          newStatesSinceComponentSearch.clear();
          if (logger.isLoggable(Level.FINE)) {
            logger.log(Level.FINE, String.format("Found %d components with %d states",
                components.size(), components.stream().mapToInt(NatBitSet::size).sum()));
          }

          components.stream().filter(component -> SccDecomposition.isBscc(this.model::getSuccessors, component))
              .forEach(component -> {
                createComponent(component);
                visitedStateSet.removeAll(component);
              });
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
        double lowerBound = transitions.sumWeightedExceptJacobi(this::getWeightedErrorBound, state);
        computedTransientUpdates += 1;
        if (Double.isNaN(lowerBound) || Util.isOne(lowerBound)) {
          updateErrorBound = false;
        } else {
          assert lowerBound >= 0.0;
          transientErrorBounds.put(state, lowerBound);
        }
      }

      if (updateComponentReachability) {
        BottomComponent finalVisitedComponent = visitedComponent;
        double lowerBound = transitions.sumWeightedExceptJacobi(s -> getReachabilityLowerBounds(finalVisitedComponent, s), state);
        if (Double.isNaN(lowerBound) || Util.isZero(lowerBound)) {
          updateComponentReachability = false;
        } else {
          visitedComponent.updateReachabilityLower(state, lowerBound);

          double absorption = 0.0;
          for (BottomComponent component : components) {
            absorption += component.getReachabilityLower(state);
          }
          absorptionLowerBound.put(state, absorption);
        }
      }
    }
  }
}
