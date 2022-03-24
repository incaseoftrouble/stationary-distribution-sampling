package stationary.component;

import de.tum.in.probmodels.graph.Component;
import de.tum.in.probmodels.util.Util;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2LongLinkedOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2LongMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntHeapPriorityQueue;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntPriorityQueue;
import stationary.util.Bound;
import stationary.util.Check;

public final class NontrivialApproximatingComponent extends NontrivialComponent {
  private final IntPriorityQueue queue;
  private final Int2ObjectMap<Int2DoubleMap> iteration;
  private final Int2ObjectMap<Bound> frequencyBounds;
  private final Int2LongMap samplingCounts;

  private double lowerBoundCache = 0.0;
  private double precisionGuide = 1.0;

  public NontrivialApproximatingComponent(int index, Component component) {
    super(index, component);

    int size = component.size();

    iteration = new Int2ObjectOpenHashMap<>(size);
    component.states().forEach((int s) -> iteration.put(s, new Int2DoubleOpenHashMap()));

    frequencyBounds = new Int2ObjectOpenHashMap<>(size);
    component.states().forEach((int s) -> frequencyBounds.put(s, Bound.UNKNOWN));

    samplingCounts = new Int2LongLinkedOpenHashMap(size);
    component.states().forEach((int s) -> {
      int current = s;
      for (int i = 0; i < size; i++) {
        current = component.onlyChoice(current).distribution().sample();
      }
      for (int i = 0; i < size; i++) {
        current = component.onlyChoice(current).distribution().sample();
        samplingCounts.mergeLong(current, 1L, Long::sum);
      }
    });

    queue = new IntHeapPriorityQueue((int a, int b) -> {
      double aError = stateErrorBound(a);
      double bError = stateErrorBound(b);
      if (Util.isEqual(aError, bError)) {
        return Long.compare(samplingCounts.get(b), samplingCounts.get(a));
      }
      return aError < bError ? 1 : -1;
    });
    component.states().forEach(queue::enqueue);
  }

  private double stateErrorBound(int s) {
    Bound bound = frequencyBounds.get(s);
    double upperBound = 1.0 - lowerBoundCache + bound.lower();
    return Math.min(bound.difference(), upperBound);
  }

  @Override
  protected void doUpdate(int initialState) {
    int size = component.size();

    for (int i = 0; i < 10; i++) {
      if (queue.isEmpty()) {
        break;
      }
      int s = queue.dequeueInt();

      Int2DoubleMap map = iteration.get(s);
      map.defaultReturnValue(0.0);
      Int2DoubleMap current = map;
      Int2DoubleMap next = new Int2DoubleOpenHashMap(size);
      next.defaultReturnValue(0.0);

      Bound bound = null;
      for (int step = 0; step < size; step++) {
        double minimalDifference = Double.MAX_VALUE;
        double maximalDifference = Double.MIN_VALUE;

        IntIterator iterator = component.states().iterator();
        while (iterator.hasNext()) {
          int t = iterator.nextInt();

          double value = component.onlyChoice(t).distribution().sumWeighted(current);
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

        bound = Bound.of(minimalDifference, maximalDifference);
        Int2DoubleMap swap = next;
        next = current;
        current = swap;
        if (bound.difference() < precisionGuide) {
          break;
        }
      }
      //noinspection ObjectEquality
      if (current != map) {
        iteration.put(s, current);
      }

      assert bound != null;
      frequencyBounds.put(s, bound);
      if (!Util.isZero(bound.difference())) {
        queue.enqueue(s);
      }
    }
    lowerBoundCache = frequencyBounds.values().stream().mapToDouble(Bound::lower).sum();
    if (!queue.isEmpty() && frequencyBounds.get(queue.firstInt()).difference() < precisionGuide) {
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
    // double error = this.stateError.values().doubleStream().max().orElseThrow();
    // assert Util.lessOrEqual(0.0, error) && Util.lessOrEqual(error, 1.0);
    return 1.0 - lowerBoundCache;
  }

  @Override
  public Bound frequency(int state) {
    assert Check.checkFrequency(component.states(),
        s -> frequencyBounds.get(s).lower(),
        s -> component.onlyChoice(s).distribution(), error());
    return frequencyBounds.get(state);
  }
}
