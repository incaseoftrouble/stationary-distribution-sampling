package stationary;

import de.tum.in.naturals.bitset.BitSets;
import de.tum.in.probmodels.SelfLoopHandling;
import de.tum.in.probmodels.explorer.DefaultExplorer;
import de.tum.in.probmodels.generator.Generator;
import de.tum.in.probmodels.graph.BsccComponentAnalyser;
import de.tum.in.probmodels.graph.Component;
import de.tum.in.probmodels.model.impl.DenseDeterministicStochasticSystem;
import de.tum.in.probmodels.prism.model.DenseMarkovChainView;
import explicit.DTMCModelChecker;
import explicit.ModelCheckerResult;
import explicit.ProbModelChecker;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import it.unimi.dsi.fastutil.ints.Int2ObjectOpenHashMap;
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
import stationary.util.FrequencyRecord;

public final class StationaryDistributionSolver {
  private static final Logger logger = Logger.getLogger(StationaryDistributionSolver.class.getName());

  private final DefaultExplorer<State, DenseDeterministicStochasticSystem> explorer;
  private final double precision;
  private final DenseDeterministicStochasticSystem model;
  private final BsccComponentAnalyser analyser = new BsccComponentAnalyser();

  public StationaryDistributionSolver(Generator<State> generator, double precision) {
    this.model = new DenseDeterministicStochasticSystem();
    this.explorer = DefaultExplorer.of(this.model, generator, SelfLoopHandling.KEEP);
    this.precision = precision;
  }

  public State getState(int s) {
    return explorer.getState(s);
  }

  public Int2ObjectMap<FrequencyRecord> solve() throws PrismException {
    int initialState = model.onlyInitialState();
    explorer.exploreReachable();

    var bottomComponents = analyser.findComponents(model);
    if (logger.isLoggable(Level.FINE)) {
      logger.log(Level.FINE, String.format("Found %d BSCCs with %d states",
          bottomComponents.size(), bottomComponents.stream().mapToInt(Component::size).sum()));
    }
    List<BottomComponent> components = new ArrayList<>(bottomComponents.size());
    for (Component component : bottomComponents) {
      int index = components.size();
      BottomComponent wrapper = component.size() == 1
          ? new AbsorbingComponent(index, component.states().iterator().nextInt())
          : new NontrivialSolvingComponent(index, component);
      components.add(wrapper);
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
      ModelCheckerResult result = mc.computeReachProbs(new DenseMarkovChainView(model), BitSets.of(component.states()));
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
