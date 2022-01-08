package stationary;

import de.tum.in.naturals.set.NatBitSets;
import de.tum.in.naturals.set.RoaringNatBitSetFactory;
import de.tum.in.probmodels.generator.DtmcGenerator;
import de.tum.in.probmodels.util.PrismHelper;
import it.unimi.dsi.fastutil.ints.Int2ObjectMap;
import java.io.IOException;
import java.nio.file.Path;
import java.util.Comparator;
import java.util.List;
import java.util.concurrent.Callable;
import java.util.function.IntFunction;
import java.util.logging.Level;
import java.util.logging.Logger;
import parser.ast.ModulesFile;
import picocli.CommandLine;
import prism.ModelGenerator;
import prism.Prism;
import prism.PrismDevNullLog;
import prism.PrismException;
import simulator.ModulesFileModelGenerator;
import stationary.util.FrequencyRecord;

@SuppressWarnings({"FieldMayBeFinal", "ProhibitedExceptionDeclared"})
@CommandLine.Command(name = "dist-est", mixinStandardHelpOptions = true, version = "0.1",
                     description = "Computes/approximates the stationary distribution of a given Markov chain in various ways",
                     subcommands = {Main.Approximate.class, Main.Solve.class})
public class Main {
  private static final Logger logger = Logger.getLogger(Main.class.getName());

  @CommandLine.Parameters(index = "0", description = "The PRISM MC file")
  private Path path = null;

  @CommandLine.Option(names = {"-c", "--const"}, description = "Constants", split = ",", arity = "0..*")
  private List<String> constants = List.of();

  @CommandLine.Command(name = "approximate")
  public static class Approximate implements Callable<Void> {
    @SuppressWarnings("InstanceVariableMayNotBeInitialized")
    @CommandLine.ParentCommand
    private Main main;

    @CommandLine.Option(names = "--sampling", description = "Sampling mode to choose")
    private StationaryDistributionEstimator.Mode samplingMode = StationaryDistributionEstimator.Mode.SAMPLE_TARGET;

    @CommandLine.Option(names = "--precision", description = "Desired precision")
    private double precision = 1.0e-6;

    @CommandLine.Option(names = "--explore", description = "Eagerly build the complete model")
    private boolean preExplore = false;

    @Override
    public Void call() throws Exception {
      StationaryDistributionEstimator estimator = new StationaryDistributionEstimator(
          new DtmcGenerator(main.obtainGenerator()), precision, samplingMode, preExplore);
      Int2ObjectMap<FrequencyRecord> frequency = estimator.solve();
      Main.printResult(frequency, estimator::getState);
      return null;
    }
  }

  @CommandLine.Command(name = "solve")
  public static class Solve implements Callable<Void> {
    @SuppressWarnings("InstanceVariableMayNotBeInitialized")
    @CommandLine.ParentCommand
    private Main main;

    @CommandLine.Option(names = "--precision", description = "Desired precision")
    private double precision = 1.0e-6;

    @Override
    public Void call() throws Exception {
      StationaryDistributionSolver solver = new StationaryDistributionSolver(new DtmcGenerator(main.obtainGenerator()), precision);
      Int2ObjectMap<FrequencyRecord> frequency = solver.solve();
      Main.printResult(frequency, solver::getState);
      return null;
    }
  }

  public static void main(String[] args) {
    logger.log(Level.INFO, "Command line: {0}", String.join(" ", args));
    System.exit(new CommandLine(new Main()).execute(args));
  }

  protected ModelGenerator obtainGenerator() throws PrismException, IOException {
    NatBitSets.setFactory(new RoaringNatBitSetFactory());
    PrismHelper.PrismParseResult parse = PrismHelper.parse(path, null, String.join(",", constants));
    ModulesFile modulesFile = parse.modulesFile();
    Prism prism = new Prism(new PrismDevNullLog());
    return new ModulesFileModelGenerator(modulesFile, prism);
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
