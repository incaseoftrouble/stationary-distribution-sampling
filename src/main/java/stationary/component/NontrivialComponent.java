package stationary.component;

import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.probmodels.model.Distribution;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntSets;
import java.time.Duration;
import java.util.function.IntFunction;

public abstract class NontrivialComponent extends BottomComponent {
  final NatBitSet states;
  final IntFunction<Distribution> successors;

  private long updateTime = 0L;

  public NontrivialComponent(int index, NatBitSet states, IntFunction<Distribution> successors) {
    super(index);
    this.states = states;
    this.successors = successors;
  }

  @Override
  public IntSet states() {
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
