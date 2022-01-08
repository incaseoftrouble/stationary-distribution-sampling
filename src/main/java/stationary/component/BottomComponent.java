package stationary.component;

import it.unimi.dsi.fastutil.HashCommon;
import it.unimi.dsi.fastutil.ints.IntSet;
import stationary.util.Bound;

public abstract class BottomComponent {
  public final int index;
  private int visits = 0;

  public BottomComponent(int index) {
    this.index = index;
  }

  public abstract IntSet states();

  public abstract void update(int state);

  public abstract double error();

  public abstract Bound frequency(int state);

  public boolean contains(int state) {
    return states().contains(state);
  }

  public void countVisit() {
    visits += 1;
  }

  public int getVisitCount() {
    return visits;
  }

  public int size() {
    return states().size();
  }

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
