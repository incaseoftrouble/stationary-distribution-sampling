package stationary;

import static java.util.Objects.requireNonNull;
import static picocli.CommandLine.Command;
import static picocli.CommandLine.Option;
import static picocli.CommandLine.Parameters;
import static picocli.CommandLine.ParentCommand;

import de.tum.in.naturals.set.NatBitSets;
import de.tum.in.naturals.set.RoaringNatBitSetFactory;
import de.tum.in.probmodels.SelfLoopHandling;
import de.tum.in.probmodels.explorer.DefaultExplorer;
import de.tum.in.probmodels.generator.Generator;
import de.tum.in.probmodels.model.impl.DenseDeterministicStochasticSystem;
import de.tum.in.probmodels.prism.generator.CtmcUniformizingGenerator;
import de.tum.in.probmodels.prism.generator.NondeterministicUniformGenerator;
import de.tum.in.probmodels.prism.generator.StochasticGenerator;
import de.tum.in.probmodels.prism.model.DenseMarkovChainView;
import de.tum.in.probmodels.util.PrismHelper;
import explicit.CTMC;
import explicit.DTMCModelChecker;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.function.IntFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import parser.State;
import parser.ast.ModulesFile;
import picocli.CommandLine;
import prism.ModelType;
import prism.Modules2MTBDD;
import prism.Prism;
import prism.PrismDevNullLog;
import prism.PrismException;
import prism.PrismLog;
import prism.PrismPrintStreamLog;
import prism.PrismSettings;
import prism.ProbModelChecker;
import prism.StochModelChecker;
import simulator.ModulesFileModelGenerator;
import stationary.util.FrequencyRecord;

@Command(name = "dist-est", mixinStandardHelpOptions = true, version = "0.1",
         description = "Computes/approximates the stationary distribution of a given Markov chain in various ways",
         subcommands = {Main.Approximate.class, Main.Solve.class, Main.PrismCall.class, Main.UniformizationConstant.class,
             Main.ModelStatistics.class})
public class Main {
  private static final Logger logger = Logger.getLogger(Main.class.getName());

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
    public Void call() throws PrismException, IOException {
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

      PrismHelper.PrismParseResult parse = PrismHelper.parse(main.path, null, String.join(",", main.constants));
      ModulesFile modulesFile = parse.modulesFile();
      if (explicit) {
        ModelType modelType = modulesFile.getModelType();
        switch (modelType) {
          case DTMC, CTMC -> {
            prism.loadPRISMModel(modulesFile);
            prism.doSteadyState();
          }
          case MDP -> {
            var generator = new NondeterministicUniformGenerator(new ModulesFileModelGenerator(modulesFile, prism), true);
            DenseDeterministicStochasticSystem system = new DenseDeterministicStochasticSystem();
            DefaultExplorer.of(system, generator, SelfLoopHandling.KEEP).exploreReachable();
            new DTMCModelChecker(prism).doSteadyState(new DenseMarkovChainView(system));
          }
          default -> throw new AssertionError(modelType);
        }
      } else {
        var model = new Modules2MTBDD(prism, modulesFile).translate();
        switch (model.getModelType()) {
          case DTMC -> new ProbModelChecker(prism, model, null).doSteadyState();
          case MDP -> throw new IllegalArgumentException("MDP only works with explicit engine");
          case CTMC -> new StochModelChecker(prism, model, null).doSteadyState();
          default -> throw new IllegalArgumentException("Unsupported model type " + model.getModelType());
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
      PrismHelper.PrismParseResult parse = PrismHelper.parse(main.path, null, String.join(",", main.constants));
      ModulesFile modulesFile = parse.modulesFile();
      Prism prism = new Prism(new PrismDevNullLog());
      PrismSettings settings = new PrismSettings();
      settings.set(PrismSettings.PRISM_ENGINE, "Explicit");
      prism.setSettings(settings);
      prism.initialise();
      prism.loadPRISMModel(modulesFile);
      prism.buildModel();
      System.out.println(((CTMC) prism.getBuiltModelExplicit()).getDefaultUniformisationRate());
      return null;
    }
  }

  @Command(name = "stats")
  public static class ModelStatistics implements Callable<Void> {
    @ParentCommand
    private Main main;

    @Override
    public Void call() throws Exception {
      PrismHelper.PrismParseResult parse = PrismHelper.parse(main.path, null, String.join(",", main.constants));
      ModulesFile modulesFile = parse.modulesFile();
      Prism prism = new Prism(new PrismDevNullLog());
      prism.setSettings(new PrismSettings());
      prism.initialise();
      prism.loadPRISMModel(modulesFile);
      prism.buildModel();
      System.out.println(prism.getBuiltModel().getNumStates());
      return null;
    }
  }

  public static void main(String[] args) {
    NatBitSets.setFactory(new RoaringNatBitSetFactory());
    logger.log(Level.INFO, "Command line: {0}", String.join(" ", args));
    System.exit(new CommandLine(new Main()).execute(args));
  }


  protected Generator<State> obtainGenerator() throws PrismException, IOException {
    PrismHelper.PrismParseResult parse = PrismHelper.parse(path, null, String.join(",", constants));
    ModulesFile modulesFile = parse.modulesFile();
    Prism prism = new Prism(new PrismDevNullLog());
    PrismSettings settings = new PrismSettings();
    prism.setSettings(settings);

    var generator = new ModulesFileModelGenerator(modulesFile, prism);
    return switch (generator.getModelType()) {
      case CTMC -> new CtmcUniformizingGenerator(generator, requireNonNull(uniformizationConstant));
      case DTMC -> new StochasticGenerator(generator, true);
      case MDP -> new NondeterministicUniformGenerator(generator, true);
      default -> throw new IllegalArgumentException("Illegal model type " + generator.getModelType());
    };
  }

  protected static <V> void printResult(Int2ObjectMap<FrequencyRecord> frequency, IntFunction<V> stateFunction) {
    frequency.int2ObjectEntrySet().stream().sorted(Comparator
            .comparingDouble((Int2ObjectMap.Entry<FrequencyRecord> e) -> e.getValue().frequency()).reversed()
            .thenComparing(Int2ObjectMap.Entry::getIntKey))
        .forEach(entry -> System.out.printf("%s: %.6g (+ %.5g)%n", stateFunction.apply(entry.getIntKey()),
            entry.getValue().frequency(), entry.getValue().error()));
    System.out.printf("%nRemaining: %.5g%n", 1.0 - frequency.values().stream().mapToDouble(FrequencyRecord::frequency).sum());
  }
}
