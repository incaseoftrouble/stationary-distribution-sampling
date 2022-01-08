package stationary.util;

import de.tum.in.probmodels.explorer.Explorer;
import de.tum.in.probmodels.model.Distribution;
import de.tum.in.probmodels.model.Model;
import it.unimi.dsi.fastutil.ints.IntArrayList;
import it.unimi.dsi.fastutil.ints.IntIterator;
import it.unimi.dsi.fastutil.ints.IntStack;
import prism.PrismException;

public final class Explore {
  private Explore() {}

  public static <S, M extends Model> void explore(Explorer<S, M> explorer, int initialState) throws PrismException {
    IntStack search = new IntArrayList();
    search.push(initialState);
    while (!search.isEmpty()) {
      int state = search.popInt();
      for (Distribution choice : explorer.getChoices(state)) {
        IntIterator support = choice.support().intIterator();
        while (support.hasNext()) {
          int s = support.nextInt();
          if (!explorer.isExploredState(s)) {
            explorer.exploreState(s);
            search.push(s);
          }
        }
      }
    }
  }
}
