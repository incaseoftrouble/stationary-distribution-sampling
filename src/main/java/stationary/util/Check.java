package stationary.util;

import de.tum.in.probmodels.model.Distribution;
import de.tum.in.probmodels.util.Util;
import it.unimi.dsi.fastutil.ints.Int2DoubleMap;
import it.unimi.dsi.fastutil.ints.Int2DoubleOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntSet;
import java.util.function.IntFunction;
import java.util.function.IntToDoubleFunction;

public final class Check {
  private Check() {}

  public static boolean checkFrequency(IntSet states, IntToDoubleFunction frequency, IntFunction<Distribution> successors, double precision) {
    Int2DoubleMap incoming = new Int2DoubleOpenHashMap();
    states.forEach((int s) -> {
      Distribution distribution = successors.apply(s);
      double freq = frequency.applyAsDouble(s);
      distribution.forEach((t, p) -> incoming.mergeDouble(t, p * freq, Double::sum));
    });
    for (Int2DoubleMap.Entry entry : incoming.int2DoubleEntrySet()) {
      int state = entry.getIntKey();
      double stateFrequency = frequency.applyAsDouble(state);
      double incomingFrequency = entry.getDoubleValue();
      assert Util.lessOrEqual(Math.abs(incomingFrequency - stateFrequency), precision) :
          "Frequency mismatch in state %d: Got %.5g, expected %.5g".formatted(state, incomingFrequency, stateFrequency);
    }
    return true;
  }
}
