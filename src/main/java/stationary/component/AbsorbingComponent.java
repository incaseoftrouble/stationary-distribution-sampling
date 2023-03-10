package stationary.component;

import de.tum.in.probmodels.values.Bounds;
import it.unimi.dsi.fastutil.ints.IntSet;

public final class AbsorbingComponent extends BottomComponent {
  private final int state;

  public AbsorbingComponent(int index, int state) {
    super(index);
    this.state = state;
  }

  @Override
  public IntSet states() {
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
  public Bounds frequency(int s) {
    return Bounds.one();
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
