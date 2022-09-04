package stationary;

import static picocli.CommandLine.Command;
import static picocli.CommandLine.Option;
import static picocli.CommandLine.Parameters;
import static picocli.CommandLine.ParentCommand;

import de.tum.in.naturals.set.NatBitSets;
import de.tum.in.naturals.set.RoaringNatBitSetFactory;
import de.tum.in.probmodels.explorer.DefaultExplorer;
import de.tum.in.probmodels.explorer.SelfLoopHandling;
import de.tum.in.probmodels.generator.Action;
import de.tum.in.probmodels.generator.Generator;
import de.tum.in.probmodels.graph.BsccComponentAnalyser;
import de.tum.in.probmodels.impl.prism.PrismProblemInstance;
import de.tum.in.probmodels.impl.prism.generator.NondeterministicStrategyGenerator;
import de.tum.in.probmodels.impl.prism.model.DenseMarkovChainView;
import de.tum.in.probmodels.model.distribution.ArrayDistributionBuilder;
import de.tum.in.probmodels.model.distribution.Distribution;
import de.tum.in.probmodels.model.distribution.DistributionBuilder;
import explicit.CTMC;
import explicit.DTMCModelChecker;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.function.BiFunction;
import java.util.function.IntFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import jdd.JDDNode;
import parser.State;
import parser.ast.ModulesFile;
import picocli.CommandLine;
import prism.ECComputer;
import prism.Model;
import prism.NondetModel;
import prism.Prism;
import prism.PrismDevNullLog;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismPrintStreamLog;
import prism.PrismSettings;
import prism.SCCComputer;
import prism.StateListMTBDD;
import stationary.util.FrequencyRecord;

@Command(name = "dist-est", mixinStandardHelpOptions = true, version = "0.1",
         description = "Computes/approximates the stationary distribution of a given Markov chain in various ways",
         subcommands = {Main.Approximate.class, Main.Solve.class, Main.PrismCall.class, Main.UniformizationConstant.class,
             Main.ModelSize.class, Main.ModelComponents.class})
public class Main {
  private static final Logger logger = Logger.getLogger(Main.class.getName());

  private static final double gamma = 0.9;
  private static final BiFunction<State, List<Action<State>>, Distribution> strategy = (s, a) -> {
    DistributionBuilder builder = new ArrayDistributionBuilder();
    double p = 1.0;
    for (int i = 0; i < a.size(); i++) {
      p *= gamma;
      builder.add(i, 1.0 / a.size());
    }
    return builder.scaled();
  };

  @Parameters(index = "0", description = "The PRISM MC file")
  private Path path;

  @Option(names = {"-c", "--const"}, description = "Constants", split = ",", arity = "0..*")
  private List<String> constants = List.of();

  @Option(names = {"--unif", "--uniformization"}, description = "Uniformization constant for CTMC")
  private Double uniformizationConstant;

  @Command(name = "approximate")
  public static class Approximate implements Callable<Void> {
    @SuppressWarnings("InstanceVariableMayNotBeInitialized")
    @ParentCommand
    private Main main;

    @Option(names = "--sampling", description = "Sampling mode to choose")
    private StationaryDistributionEstimator.Mode samplingMode = StationaryDistributionEstimator.Mode.SAMPLE_TARGET;

    @Option(names = "--precision", description = "Desired precision")
    private double precision = 1.0e-6;

    @Option(names = "--explore", description = "Eagerly build the complete model")
    private boolean preExplore = false;

    @Option(names = "--solve-bsccs", description = "Use equation solving for BSCCs")
    private boolean solveComponents = false;

    @Override
    public Void call() throws IOException {
      var estimator = new StationaryDistributionEstimator(main.obtainGenerator(), precision, samplingMode, preExplore, solveComponents);
      Int2ObjectMap<FrequencyRecord> frequency = estimator.solve();
      Main.printResult(frequency, estimator::getState);
      return null;
    }
  }

  @Command(name = "solve")
  public static class Solve implements Callable<Void> {
    @SuppressWarnings("InstanceVariableMayNotBeInitialized")
    @ParentCommand
    private Main main;

    @Option(names = "--precision", description = "Desired precision")
    private double precision = 1.0e-6;

    @Override
    public Void call() throws PrismException, IOException {
      StationaryDistributionSolver solver = new StationaryDistributionSolver(main.obtainGenerator(), precision);
      Int2ObjectMap<FrequencyRecord> frequency = solver.solve();
      Main.printResult(frequency, solver::getState);
      return null;
    }
  }

  @Command(name = "prism")
  public static class PrismCall implements Callable<Void> {
    @ParentCommand
    private Main main;

    @Option(names = "--precision", description = "Desired precision")
    private double precision = 1.0e-6;

    @Option(names = "--explicit", description = "Use explicit engine")
    private boolean explicit = false;

    @Override
    public Void call() throws PrismException, IOException {
      PrismLog log = new PrismPrintStreamLog(System.out);
      Prism prism = new Prism(log);
      PrismSettings settings = new PrismSettings();
      prism.setSettings(settings);
      settings.set(PrismSettings.PRISM_TERM_CRIT, "Absolute");
      settings.set(PrismSettings.PRISM_TERM_CRIT_PARAM, precision);
      if (explicit) {
        settings.set(PrismSettings.PRISM_ENGINE, "Explicit");
      }
      settings.set(PrismSettings.PRISM_MAX_ITERS, 10 * 1000 * 1000);
      prism.initialise();

      PrismProblemInstance parse = PrismProblemInstance.of(prism, main.path, null, main.constants, main.uniformizationConstant);
      if (explicit) {
        switch (parse.generator().getModelType()) {
          case DTMC, CTMC -> {
            prism.loadPRISMModel(parse.modulesFile());
            prism.doSteadyState();
          }
          case MDP -> {
            var explorer = DefaultExplorer.of(new NondeterministicStrategyGenerator(parse.generator(),
                strategy, true), SelfLoopHandling.KEEP);
            explorer.exploreReachable(explorer.initialStateIds());
            new DTMCModelChecker(prism).doSteadyState(new DenseMarkovChainView(explorer.partialSystem()));
          }
          default -> throw new AssertionError(parse.generator().getModelType());
        }
      } else {
        // var model = new Modules2MTBDD(prism, parse.modulesFile()).translate();
        switch (parse.generator().getModelType()) {
          case DTMC, CTMC -> {
            prism.loadPRISMModel(parse.modulesFile());
            prism.doSteadyState();
          }
          case MDP -> throw new IllegalArgumentException("MDP only works with explicit engine");
          default -> throw new AssertionError(parse.generator().getModelType());
        }
      }
      return null;
    }
  }

  @Command(name = "uniformization")
  public static class UniformizationConstant implements Callable<Void> {
    @ParentCommand
    private Main main;

    @Override
    public Void call() throws Exception {
      Prism prism = new Prism(new PrismDevNullLog());
      PrismSettings settings = new PrismSettings();
      settings.set(PrismSettings.PRISM_ENGINE, "Explicit");
      prism.setSettings(settings);
      prism.initialise();

      PrismProblemInstance instance = PrismProblemInstance.of(prism, main.path, null, main.constants, null);
      prism.loadPRISMModel(instance.modulesFile());
      prism.buildModel();
      System.out.println(((CTMC) prism.getBuiltModelExplicit()).getDefaultUniformisationRate());
      return null;
    }
  }

  @Command(name = "stats")
  public static class ModelSize implements Callable<Void> {
    @ParentCommand
    private Main main;

    @Option(names = "--components", description = "Compute components")
    private boolean components;

    @Override
    public Void call() throws Exception {
      Prism prism = new Prism(new PrismDevNullLog());
      prism.setSettings(new PrismSettings());
      prism.initialise();

      PrismProblemInstance instance = PrismProblemInstance.of(prism, main.path, null, main.constants, main.uniformizationConstant);
      ModulesFile modulesFile = instance.modulesFile();
      prism.loadPRISMModel(modulesFile);
      prism.buildModel();
      Model model = prism.getBuiltModel();
      System.out.println(model.getNumStates());
      if (components) {
        List<JDDNode> components;
        if (model instanceof NondetModel nondetModel) {
          ECComputer ecComputer = ECComputer.createECComputer(prism, nondetModel);
          ecComputer.computeMECStates();
          components = ecComputer.getMECStates();
        } else {
          SCCComputer sccComputer = SCCComputer.createSCCComputer(prism, model);
          sccComputer.computeBSCCs();
          components = sccComputer.getBSCCs();
        }
        System.out.println(components.size());
        System.out.println(components.stream()
            .map(c -> new StateListMTBDD(c, model).sizeString())
            .mapToLong(c -> c.isEmpty() ? Long.MAX_VALUE : Long.parseLong(c))
            .max().orElseThrow());
      }
      return null;
    }
  }

  @Command(name = "stats-explicit")
  public static class ModelComponents implements Callable<Void> {
    @ParentCommand
    private Main main;

    @Override
    public Void call() throws Exception {
      var explorer = DefaultExplorer.of(main.obtainGenerator(), SelfLoopHandling.KEEP);
      explorer.exploreReachable(explorer.initialStateIds());
      System.out.println(explorer.partialSystem().stateCount());
      System.out.println(new BsccComponentAnalyser().findComponents(explorer.partialSystem()).size());
      return null;
    }
  }

  public static void main(String[] args) {
    NatBitSets.setFactory(new RoaringNatBitSetFactory());
    logger.log(Level.INFO, "Command line: {0}", String.join(" ", args));
    System.exit(new CommandLine(new Main()).execute(args));
  }


  protected Generator<State> obtainGenerator() throws IOException {
    Prism prism = new Prism(new PrismDevNullLog());
    PrismSettings settings = new PrismSettings();
    prism.setSettings(settings);
    PrismProblemInstance instance = PrismProblemInstance.of(prism, path, null, constants, uniformizationConstant);
    Generator<State> model = instance.model();
    if (model.isDeterministic()) {
      return model;
    }
    return new NondeterministicStrategyGenerator(instance.generator(), strategy, true);
  }

  protected static <V> void printResult(Int2ObjectMap<FrequencyRecord> frequency, IntFunction<V> stateFunction) {
    frequency.int2ObjectEntrySet().stream().sorted(Comparator
            .comparingDouble((Int2ObjectMap.Entry<FrequencyRecord> e) -> e.getValue().frequency()).reversed()
            .thenComparing(Int2ObjectMap.Entry::getIntKey))
        .forEach(entry -> System.out.printf("%s: %.6g (+ %.5g)%n", stateFunction.apply(entry.getIntKey()),
            entry.getValue().frequency(), entry.getValue().error()));
  }
}
