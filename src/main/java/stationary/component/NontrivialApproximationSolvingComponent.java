package stationary.component;

import de.tum.in.probmodels.graph.Component;
import de.tum.in.probmodels.model.distribution.Distribution;
import de.tum.in.probmodels.util.Util;
import de.tum.in.probmodels.values.Bounds;
import it.unimi.dsi.fastutil.ints.Int2IntOpenHashMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import java.util.logging.Level;
import java.util.logging.Logger;
import stationary.util.Check;

public final class NontrivialApproximationSolvingComponent extends NontrivialComponent {
  private static final Logger logger = Logger.getLogger(NontrivialApproximationSolvingComponent.class.getName());

  private final Int2ObjectMap<Bounds> frequencyBounds;
  private double error = 1.0;
  private final double precision;

  public NontrivialApproximationSolvingComponent(int index, Component component, double precision) {
    super(index, component);
    this.precision = precision;

    int size = component.size();
    frequencyBounds = new Int2ObjectOpenHashMap<>(size);
  }

  @Override
  protected void doUpdate(int initialState) {
    if (!frequencyBounds.isEmpty()) {
      return;
    }

    int size = component.size();
    boolean enableRegularization = false;
    double lastReport = System.currentTimeMillis();

    int[] reverse = new int[size];
    Int2IntOpenHashMap renumbering = new Int2IntOpenHashMap(size);

    int index = 0;
    var renumberIterator = component.states().iterator();
    while (renumberIterator.hasNext()) {
      int s = renumberIterator.nextInt();
      renumbering.put(s, index);
      reverse[index] = s;
      index+= 1;
    }

    index = 0;
    Distribution[] distributions = new Distribution[size];
    var remapIterator = component.states().iterator();
    while (remapIterator.hasNext()) {
      distributions[index] = component.onlyChoice(remapIterator.nextInt()).distribution().map(renumbering).build();
      index+= 1;
    }

    double[] next = new double[size];

    for (int s = 0; s < size; s++) {
      double[] current = new double[size];

      int steps = 0;
      while (true) {
        double minimalDifference = Double.MAX_VALUE;
        double maximalDifference = Double.MIN_VALUE;

        for (int t = 0; t < size; t++) {
          double currentValue = current[t];
          double value = distributions[t].sumWeighted(current);
          if (enableRegularization) {
            value = 0.05 * currentValue + 0.95 * value;
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

        Bounds bound = Bounds.reach(minimalDifference, maximalDifference);
        double[] swap = next;
        next = current;
        current = swap;
        if (bound.difference() < precision) {
          frequencyBounds.put(reverse[s], bound);
          break;
        }
        steps += 1;
        if (steps == size && Util.isOne(bound.difference())) {
          enableRegularization = true;
        }
      }

      if (System.currentTimeMillis() - lastReport > 5000) {
        lastReport = System.currentTimeMillis();
        logger.log(Level.INFO, "Solved %d / %d (%.1f %%) states".formatted(s, component.size(), 100.0 * s / component.size()));
      }
    }
    error = precision;
  }

  @Override
  public double error() {
    return error;
  }

  @Override
  public Bounds frequency(int state) {
    assert Check.checkFrequency(component.states(),
        s -> frequencyBounds.get(s).lowerBound(),
        s -> component.onlyChoice(s).distribution(), error());
    return frequencyBounds.getOrDefault(state, Bounds.unknownReach());
  }
}
