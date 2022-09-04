package stationary.component;

import de.tum.in.probmodels.graph.Component;
import de.tum.in.probmodels.model.distribution.Distribution;
import de.tum.in.probmodels.util.Util;
import de.tum.in.probmodels.values.Bounds;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntHeapPriorityQueue;
import it.unimi.dsi.fastutil.ints.IntPriorityQueue;
import java.util.Arrays;
import stationary.util.Check;

public final class NontrivialApproximatingComponent extends NontrivialComponent {

  private double lowerBoundCache = 0.0;
  private double maximalErrorCache = 1.0;
  private double precisionGuide = 0.1;

  private boolean enableRegularization = false;

  private final Distribution[] distributions;
  private final IntPriorityQueue queue;
  private final double[][] iteration;
  private double[] next;
  private final Bounds[] frequencyBounds;
  private final Int2IntOpenHashMap renumbering;

  public NontrivialApproximatingComponent(int index, Component component) {
    super(index, component);

    int size = component.size();

    renumbering = new Int2IntOpenHashMap(size);

    int n = 0;
    var renumberIterator = component.states().iterator();
    while (renumberIterator.hasNext()) {
      int s = renumberIterator.nextInt();
      renumbering.put(s, n);
      n+= 1;
    }

    n = 0;
    distributions = new Distribution[size];
    var remapIterator = component.states().iterator();
    while (remapIterator.hasNext()) {
      distributions[n] = component.onlyChoice(remapIterator.nextInt()).distribution().map(renumbering).build();
      n+= 1;
    }

    iteration = new double[size][size];
    next = new double[size];

    frequencyBounds = new Bounds[size];
    Arrays.fill(frequencyBounds, Bounds.reachUnknown());

    long[] samplingCounts = new long[size];
    for (int s = 0; s < size; s++) {
      int current = s;
      for (int i = 0; i < size; i++) {
        current = distributions[current].sample();
      }
      for (int i = 0; i < size; i++) {
        current = distributions[current].sample();
        samplingCounts[current] += 1;
      }
    }

    // queue = new IntArrayFIFOQueue();
    queue = new IntHeapPriorityQueue((int a, int b) -> {
      Bounds aBounds = frequencyBounds[a];
      double aError = aBounds.difference() * aBounds.upperBound();
      Bounds bBounds = frequencyBounds[b];
      double bError = bBounds.difference() * bBounds.upperBound();
      if (Util.isEqual(aError, bError)) {
        return Long.compare(samplingCounts[b], samplingCounts[a]);
      }
      return aError < bError ? 1 : -1;
    });
    for (int i = 0; i < size; i++) {
      queue.enqueue(i);
    }
  }

  @Override
  protected void doUpdate(int initialState) {
    int size = component.size();

    if (queue.isEmpty()) {
      increasePrecision();
    }
    if (lowerBoundCache == 0.0 && getVisitCount() > size) {
      enableRegularization = true;
    }

    var next = this.next;
    for (int i = 0; i < Math.max(size / 10, 10); i++) {
      if (queue.isEmpty()) {
        break;
      }
      int s = queue.dequeueInt();

      double[] current = iteration[s];

      Bounds bound = null;
      for (int step = 0; step < size; step++) {
        double minimalDifference = Double.MAX_VALUE;
        double maximalDifference = Double.MIN_VALUE;

        for (int t = 0; t < size; t++) {
          double currentValue = current[t];
          double value = distributions[t].sumWeighted(current);
          if (enableRegularization) {
            value = 0.01 * currentValue + 0.99 * value;
          }
          if (s == t) {
            value += 1.0;
          }
          assert value >= currentValue;

          double stepDifference;
          next[t] = value;
          if (value > 0.0) {
            stepDifference = value - currentValue;
          } else {
            assert currentValue == 0.0;
            stepDifference = 0.0;
          }
          if (stepDifference > maximalDifference) {
            maximalDifference = stepDifference;
          }
          if (stepDifference < minimalDifference) {
            minimalDifference = stepDifference;
          }
        }

        bound = Bounds.reach(minimalDifference, maximalDifference);
        double[] swap = next;
        next = current;
        current = swap;
        if (bound.difference() < precisionGuide) {
          break;
        }
      }
      //noinspection ObjectEquality
      if (current != next) {
        iteration[s] = next;
        next = current;
      }

      assert bound != null;
      frequencyBounds[s] =  bound;
      if (bound.difference() > precisionGuide) {
        queue.enqueue(s);
      }
    }
    this.next = next;
    lowerBoundCache = Arrays.stream(frequencyBounds).mapToDouble(Bounds::lowerBound).sum();
    maximalErrorCache = Arrays.stream(frequencyBounds).mapToDouble(Bounds::difference).max().orElseThrow();
  }

  public void increasePrecision() {
    this.precisionGuide /= 2.0;
    for (int s = 0; s < frequencyBounds.length; s++) {
      if (frequencyBounds[s].difference() > precisionGuide) {
        queue.enqueue(s);
      }
    }
  }

  public void setMinimalPrecision(double precision) {
    precisionGuide = Math.min(precisionGuide, precision);
  }

  @Override
  public double error() {
    // double error = this.stateError.values().doubleStream().max().orElseThrow();
    // assert Util.lessOrEqual(0.0, error) && Util.lessOrEqual(error, 1.0);
    return Math.max(Math.min(1.0 - lowerBoundCache, maximalErrorCache), 0.0);
  }

  @Override
  public Bounds frequency(int state) {
    assert Check.checkFrequency(component.states(),
        s -> frequencyBounds[renumbering.get(s)].lowerBound(),
        s -> component.onlyChoice(s).distribution(), error());
    Bounds bound = frequencyBounds[renumbering.get(state)];
    double upperBound = 1.0 - lowerBoundCache + bound.lowerBound();
    return bound.upperBound() <= upperBound ? bound : bound.withUpper(upperBound);
  }
}
