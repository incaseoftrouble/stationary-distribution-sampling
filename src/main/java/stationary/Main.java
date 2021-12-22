package stationary;

import de.tum.in.naturals.set.NatBitSets;
import de.tum.in.naturals.set.RoaringNatBitSetFactory;
import de.tum.in.probmodels.generator.DtmcGenerator;
import de.tum.in.probmodels.util.PrismHelper;
import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import parser.ast.ModulesFile;
import prism.ModelGenerator;
import prism.Prism;
import prism.PrismDevNullLog;
import prism.PrismException;
import simulator.ModulesFileModelGenerator;

public class Main {
  private static final Logger logger = Logger.getLogger(Main.class.getName());

  public static void main(String[] args) throws PrismException, IOException {
    logger.log(Level.INFO, "Command line: {0}", String.join(" ", args));
    String model = args[0];
    String constants = args[1];
    double precision = Double.parseDouble(args[2]);
    NatBitSets.setFactory(new RoaringNatBitSetFactory());

    PrismHelper.PrismParseResult parse = PrismHelper.parse(model, null, constants);
    ModulesFile modulesFile = parse.modulesFile();
    Prism prism = new Prism(new PrismDevNullLog());
    ModelGenerator generator = new ModulesFileModelGenerator(modulesFile, prism);

    new StationaryDistributionEstimator(new DtmcGenerator(generator), precision).solve();
  }
}
