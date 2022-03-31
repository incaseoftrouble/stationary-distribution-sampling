package stationary.component;

import de.tum.in.probmodels.graph.Component;
import it.unimi.dsi.fastutil.ints.IntSet;
import it.unimi.dsi.fastutil.ints.IntSets;
import java.time.Duration;

public abstract class NontrivialComponent extends BottomComponent {
  private long updateTime = 0L;
  final Component component;

  public NontrivialComponent(int index, Component component) {
    super(index);
    this.component = component;
  }

  @Override
  public IntSet states() {
    return IntSets.unmodifiable(component.states());
  }

  @Override
  public boolean contains(int state) {
    return component.contains(state);
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
    return component.size();
  }

  public Duration updateTime() {
    return Duration.ofMillis(updateTime);
  }

  @Override
  public String toString() {
    return "%s[%d]".formatted(component.size() < 10 ? component.states().toString() : "<%d>".formatted(component.size()), index);
  }
}
