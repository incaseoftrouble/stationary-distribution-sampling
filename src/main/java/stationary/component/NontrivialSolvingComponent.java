package stationary.component;

import de.tum.in.naturals.Indices;
import de.tum.in.probmodels.graph.Component;
import de.tum.in.probmodels.model.distribution.Distribution;
import de.tum.in.probmodels.util.Util;
import de.tum.in.probmodels.values.Bounds;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntIterator;
import java.util.function.IntUnaryOperator;
import jeigen.DenseMatrix;
import stationary.util.Check;

public final class NontrivialSolvingComponent extends NontrivialComponent {
  private static final double FREQUENCY_PRECISION_CHECK = 1.0e-10;

  private final Int2DoubleMap frequency;

  public NontrivialSolvingComponent(int index, Component component) {
    super(index, component);

    frequency = new Int2DoubleOpenHashMap(component.size());
  }

  @Override
  protected void doUpdate(int state) {
    if (frequency.isEmpty()) {
      int size = component.size();
      DenseMatrix matrix = new DenseMatrix(size, size);
      {
        IntUnaryOperator stateToIndexMap = Indices.elementToIndexMap(component.states());
        IntIterator iterator = component.states().iterator();
        int index = 0;
        while (iterator.hasNext()) {
          int s = iterator.nextInt();
          Distribution distribution = component.onlyChoice(s).distribution();
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
      assert Util.isOne(maximum) : "Expected maximal vector 1, got %.5f".formatted(maximum);

      double[] steadyState = eig.vectors.abs().col(maximumIndex).getValues();
      double sum = Util.kahanSum(steadyState);

      IntIterator iterator = component.states().iterator();
      for (int index = 0; index < size; index++) {
        frequency.put(iterator.nextInt(), steadyState[index] / sum);
      }

      assert Check.checkFrequency(component.states(), frequency::get,
          s -> component.onlyChoice(s).distribution(), FREQUENCY_PRECISION_CHECK);
    }
  }

  @Override
  public double error() {
    return frequency.isEmpty() ? 1.0 : 0.0;
  }

  @Override
  public Bounds frequency(int state) {
    return Bounds.reach(frequency.get(state));
  }
}
