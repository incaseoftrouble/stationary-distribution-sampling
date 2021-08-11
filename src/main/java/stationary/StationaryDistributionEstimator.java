package stationary;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.probmodels.SelfLoopHandling;
import de.tum.in.probmodels.explorer.DefaultExplorer;
import de.tum.in.probmodels.generator.DtmcGenerator;
import de.tum.in.probmodels.graph.SccDecomposition;
import de.tum.in.probmodels.model.Distribution;
import de.tum.in.probmodels.model.MarkovChain;
import de.tum.in.probmodels.util.Util;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleMaps;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2IntMap;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntList;
import it.unimi.dsi.fastutil.ints.IntOpenHashSet;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntStack;
import java.util.List;
import java.util.Random;
import java.util.function.IntConsumer;
import java.util.function.IntFunction;
import java.util.function.IntToDoubleFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.annotation.Nullable;
import parser.State;

public final class StationaryDistributionEstimator {
  private static final Random random = new Random();

  private abstract static class BottomComponent {
    public final int index;

    public BottomComponent(int index) {
      this.index = index;
    }

    public abstract IntSet getStates();

    public abstract void update();

    public abstract double precisionLowerBound();

    public abstract Int2DoubleMap computeFrequencies();
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
    public void update() {
      // empty
    }

    @Override
    public double precisionLowerBound() {
      return 1.0;
    }

    @Override
    public Int2DoubleMap computeFrequencies() {
      return Int2DoubleMaps.singleton(state, 1.0);
    }
  }

  private static final class NontrivialComponent extends BottomComponent {
    private final NatBitSet component;
    private final IntFunction<Distribution> successors;
    private double computedErrorBound = 0.0;
    private final Int2IntMap stateSamplingCounts;

    private int maximalSampleLength;
    private int samplesPerUpdateRun;
    private long totalSampledStates;

    public NontrivialComponent(int index, NatBitSet component, IntFunction<Distribution> successors) {
      super(index);
      this.component = component;
      this.successors = successors;

      int size = component.size();
      this.stateSamplingCounts = new Int2IntOpenHashMap(size);
      this.stateSamplingCounts.defaultReturnValue(0);

      this.maximalSampleLength = 2 * size;
      this.samplesPerUpdateRun = size;
    }

    @Override
    public IntSet getStates() {
      return component;
    }

    @Override
    public void update() {
      for (int i = 0; i < samplesPerUpdateRun; i++) {
        int initialIndex = random.nextInt(this.component.size());
        IntIterator iterator = this.component.iterator();
        iterator.skip(initialIndex);
        int state = iterator.nextInt();

        stateSamplingCounts.mergeInt(state, 1, Integer::sum);
        totalSampledStates += 1;
        for (int j = 0; j < this.maximalSampleLength; j++) {
          state = successors.apply(state).sample();
          stateSamplingCounts.mergeInt(state, 1, Integer::sum);
          totalSampledStates += 1;
        }
      }
      samplesPerUpdateRun *= 2;

      assert stateSamplingCounts.values().intStream().asLongStream().sum() == totalSampledStates;

      IntIterator iterator = component.iterator();
      double lowerBound = 0.0;

      while (iterator.hasNext()) {
        int state = iterator.nextInt();
        // apply one step of value iteration for mean payoff

        double minimalDifference = Double.POSITIVE_INFINITY;
        IntIterator computationIterator = component.iterator();
        while (computationIterator.hasNext()) {
          int s = computationIterator.nextInt();

          double value = successors.apply(s).sumWeighted(stateSamplingCounts::get) + (s == state ? 1.0 : 0.0);
          double difference = value - stateSamplingCounts.get(s);
          if (difference < minimalDifference) {
            minimalDifference = difference;
          }
        }
        if (minimalDifference > 0.0) {
          lowerBound += minimalDifference;
        }
      }
      assert lowerBound <= 1.0;
      double newBound = lowerBound;
      if (Util.lessOrEqual(newBound, computedErrorBound)) {
        maximalSampleLength *= 2;
      }
      computedErrorBound = newBound;
    }

    @Override
    public double precisionLowerBound() {
      return computedErrorBound;
    }

    @Override
    public Int2DoubleMap computeFrequencies() {
      Int2DoubleMap frequencies = new Int2DoubleOpenHashMap(component.size());
      IntIterator iterator = component.iterator();

      while (iterator.hasNext()) {
        int state = iterator.nextInt();
        // apply one step of value iteration for mean payoff

        double minimalDifference = Double.POSITIVE_INFINITY;
        double maximalDifference = Double.NEGATIVE_INFINITY;
        IntIterator computationIterator = component.iterator();
        while (computationIterator.hasNext()) {
          int s = computationIterator.nextInt();

          double value = successors.apply(s).sumWeighted(stateSamplingCounts::get) + (s == state ? 1.0 : 0.0);
          double difference = value - stateSamplingCounts.get(s);
          if (difference < minimalDifference) {
            minimalDifference = difference;
          }
          if (difference > maximalDifference) {
            maximalDifference = difference;
          }
        }
        double value = (maximalDifference + minimalDifference) / 2;
        frequencies.put(state, value);
      }
      return frequencies;
    }
  }

  private static final int MAX_EXPLORES_PER_SAMPLE = 10;
  private final double precision;

  private final Logger logger = Logger.getLogger(StationaryDistributionEstimator.class.getName());

  private final MarkovChain model;
  private final DefaultExplorer<State, MarkovChain> explorer;

  private final IntSet newStatesSinceComponentSearch = new IntOpenHashSet();
  private int newComponentIndex = 0;
  private final Int2ObjectMap<BottomComponent> statesInBottomComponents = new Int2ObjectOpenHashMap<>();
  private final Int2DoubleMap transientStateLowerBounds = new Int2DoubleOpenHashMap();

  private final Int2ObjectMap<Int2DoubleMap> componentReachabilityLowerBound = new Int2ObjectOpenHashMap<>();

  private int loopCount = 0;
  private int loopStopsUntilCollapse;

  private double skipUpdateErrorBound;
  private long skipUpdateCount = 0;
  private long skipUpdatePrecisionLimit;

  private long sampleRunCount = 0;
  private long sampledStatesCount = 0;
  private long computedComponentUpdates = 0;
  private long computedTransientUpdates = 0;

  public StationaryDistributionEstimator(DtmcGenerator generator, double precision) {
    this.precision = precision;
    skipUpdateErrorBound = precision / 2;

    this.model = new MarkovChain();
    this.explorer = DefaultExplorer.of(model, generator, SelfLoopHandling.KEEP);

    loopStopsUntilCollapse = 10;
    skipUpdatePrecisionLimit = 1024;
  }

  public void solve() {
    int progress = 0;
    int initialState = model.getFirstInitialState();
    while (getLowerBound(initialState) < 1.0 - precision) {
      sample(initialState);
      progress += 1;
      if (progress % 10000 == 0) {
        logProgress();
      }
    }
    logProgress();
  }

  private void logProgress() {
    if (logger.isLoggable(Level.INFO)) {
      logger.log(Level.INFO, String.format("Progress Report:%n"
              + "    %d sample runs, %d samples, precision %.5f in initial state%n"
              + "    %d explored states, %d states in BSSCs (%d absorbing)%n"
              + "    %d component updates, %d transient state updates",
          sampleRunCount, sampledStatesCount, getLowerBound(model.getFirstInitialState()),
          explorer.exploredStateCount(), statesInBottomComponents.size(),
          statesInBottomComponents.values().stream().map(BottomComponent::getStates).mapToInt(IntSet::size).filter(i -> i == 1).count(),
          computedComponentUpdates, computedTransientUpdates));
    }
  }

  private double getLowerBound(int state) {
    BottomComponent component = statesInBottomComponents.get(state);
    return component == null
        ? transientStateLowerBounds.getOrDefault(state, 0.0)
        : component.precisionLowerBound();
  }

  private void sample(int initialState) {
    sampleRunCount += 1;

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
        if (component.precisionLowerBound() > skipUpdateErrorBound) {
          computedComponentUpdates += 1;
          component.update();
        } else if (component instanceof NontrivialComponent) {
          skipUpdateCount += 1;
        }
        break;
      }

      visitStack.push(currentState);
      if (!visitedStateSet.add(currentState)) {
        stateRevisit += 1;
      }

      // Sample the successor
      Distribution transitions = model.getTransitions(currentState);
      if (transitions == null) { // Deadlock state
        statesInBottomComponents.put(currentState, new AbsorbingComponent(newComponentIndex, currentState));
        visitedStateSet.remove(currentState);
        break;
      }

      // Bias towards uncertain areas!
      int nextState = transitions.sampleWeighted((s, p) -> (1 - getLowerBound(s)) * p);
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
      } else if (getLowerBound(nextState) > 1 - skipUpdateErrorBound) {
        skipUpdateCount += 1;
        break;
      }
      currentState = nextState;
    }

    if (skipUpdateCount > skipUpdatePrecisionLimit) {
      skipUpdatePrecisionLimit *= 2;
      skipUpdateErrorBound /= 2;
      skipUpdateCount = 0;
    }

    if (stateRevisit > 5) {
      loopCount += 1;
      if (loopCount > loopStopsUntilCollapse) {
        loopCount = 0;

        if (newStatesSinceComponentSearch.isEmpty()) {
          logger.fine("Skipping SCC computation");
        } else {
          List<NatBitSet> components = SccDecomposition.computeSccs(this.model::getSuccessors,
              newStatesSinceComponentSearch, this.explorer::isExploredState, false);
          newStatesSinceComponentSearch.clear();
          if (logger.isLoggable(Level.FINE)) {
            logger.log(Level.FINE, String.format("Found %d components with %d states",
                components.size(), components.stream().mapToInt(NatBitSet::size).sum()));
          }
          for (NatBitSet bottomComponent : components) {
            assert bottomComponent.intStream().noneMatch(statesInBottomComponents::containsKey);
            if (!SccDecomposition.isBscc(this.model::getSuccessors, bottomComponent)) {
              continue;
            }
            BottomComponent component = bottomComponent.size() == 1
                ? new AbsorbingComponent(newComponentIndex, bottomComponent.firstInt())
                : new NontrivialComponent(newComponentIndex, bottomComponent, model::getTransitions);
            bottomComponent.forEach((IntConsumer) state -> statesInBottomComponents.put(state, component));
            newComponentIndex += 1;

            visitedStateSet.removeAll(bottomComponent);
            transientStateLowerBounds.keySet().removeAll(bottomComponent);
          }
        }
        loopStopsUntilCollapse += explorer.exploredStates().size();
      }
    }

    @Nullable
    IntToDoubleFunction visitedComponentBounds;
    @Nullable
    Int2DoubleMap visitedComponentMap;
    if (visitedComponent == null) {
      visitedComponentBounds = null;
      visitedComponentMap = null;
    } else {
      int index = visitedComponent.index;
      visitedComponentMap = componentReachabilityLowerBound.computeIfAbsent(index, k -> new Int2DoubleOpenHashMap());
      visitedComponentBounds = s -> {
        BottomComponent component = statesInBottomComponents.get(s);
        if (component == null) {
          return visitedComponentMap.getOrDefault(s, 0.0);
        }
        return component.index == index ? component.precisionLowerBound() : 0.0;
      };
    }

    // Propagate values backwards along the path
    boolean updateReachability = true;
    boolean updateComponentReachability = visitedComponentBounds != null;
    while (!visitStack.isEmpty() && (updateReachability || updateComponentReachability)) {
      int state = visitStack.popInt();
      if (visitedStateSet.contains(state)) {
        Distribution transitions = model.getTransitions(state);
        if (updateReachability) {
          double lowerBound = transitions.sumWeightedExceptJacobi(this::getLowerBound, state);
          if (Util.isZero(lowerBound)) {
            updateReachability = false;
          } else {
            computedTransientUpdates += 1;
            transientStateLowerBounds.put(state, lowerBound);
          }
        }

        if (updateComponentReachability) {
          double successorsReachability = transitions.sumWeightedExceptJacobi(visitedComponentBounds, state);
          if (successorsReachability > 0.0) {
            visitedComponentMap.put(state, successorsReachability);
          } else {
            updateComponentReachability = false;
          }
        }
      }
    }
  }
}
