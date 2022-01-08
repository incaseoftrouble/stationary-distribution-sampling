package stationary;

import de.tum.in.naturals.bitset.BitSets;
import de.tum.in.naturals.set.NatBitSet;
import de.tum.in.probmodels.SelfLoopHandling;
import de.tum.in.probmodels.explorer.DefaultExplorer;
import de.tum.in.probmodels.generator.DtmcGenerator;
import de.tum.in.probmodels.graph.SccDecomposition;
import de.tum.in.probmodels.model.MarkovChain;
import explicit.DTMCModelChecker;
import explicit.ModelCheckerResult;
import explicit.ProbModelChecker;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
import it.unimi.dsi.fastutil.ints.IntSet;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import parser.State;
import prism.PrismDevNullLog;
import prism.PrismException;
import prism.PrismSettings;
import stationary.component.AbsorbingComponent;
import stationary.component.BottomComponent;
import stationary.component.NontrivialSolvingComponent;
import stationary.util.Explore;
import stationary.util.FrequencyRecord;

public final class StationaryDistributionSolver {
  private static final Logger logger = Logger.getLogger(StationaryDistributionSolver.class.getName());

  private final MarkovChain model;
  private final DefaultExplorer<State, MarkovChain> explorer;
  private final double precision;

  public StationaryDistributionSolver(DtmcGenerator generator, double precision) {
    this.model = new MarkovChain();
    this.explorer = DefaultExplorer.of(model, generator, SelfLoopHandling.KEEP);
    this.precision = precision;
  }

  public State getState(int s) {
    return explorer.getState(s);
  }

  public Int2ObjectMap<FrequencyRecord> solve() throws PrismException {
    int initialState = model.getFirstInitialState();
    Explore.explore(explorer, initialState);
    model.findDeadlocks(true);

    List<NatBitSet> bsccs = SccDecomposition.computeSccs(this.model::getSuccessors, IntSet.of(initialState), s -> true, false)
        .stream().filter(component -> SccDecomposition.isBscc(this.model::getSuccessors, component))
        .toList();
    if (logger.isLoggable(Level.FINE)) {
      logger.log(Level.FINE, String.format("Found %d BSCCs with %d states",
          bsccs.size(), bsccs.stream().mapToInt(NatBitSet::size).sum()));
    }
    List<BottomComponent> components = new ArrayList<>(bsccs.size());
    for (NatBitSet states : bsccs) {
      int index = components.size();
      BottomComponent component = states.size() == 1
          ? new AbsorbingComponent(index, states.firstInt())
          : new NontrivialSolvingComponent(index, states, model::getTransitions);
      components.add(component);
    }

    DTMCModelChecker mc = new DTMCModelChecker(null);
    mc.setSettings(new PrismSettings());
    mc.setLog(new PrismDevNullLog());
    mc.setTermCrit(ProbModelChecker.TermCrit.ABSOLUTE);
    mc.setTermCritParam(precision);
    mc.setDoIntervalIteration(true);

    double[] reachability = new double[components.size()];
    for (BottomComponent component : components) {
      component.update(component.states().intIterator().nextInt());
      ModelCheckerResult result = mc.computeReachProbs(model, BitSets.of(component.states()));
      reachability[component.index] = result.soln[initialState];
    }

    Int2ObjectMap<FrequencyRecord> frequency = new Int2ObjectOpenHashMap<>(components.stream().mapToInt(BottomComponent::size).sum());
    for (BottomComponent component : components) {
      double componentReachability = reachability[component.index];
      component.states().intStream().forEach((int s) ->
          frequency.put(s, new FrequencyRecord(component.frequency(s).lower() * componentReachability, 1.0e-6)));
    }

    return frequency;
  }
}
