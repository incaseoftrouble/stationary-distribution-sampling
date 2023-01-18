# Stationary Distribution Sampling

## Building

* Extract sources if necessary
* Run `docker build -t sds .`

## Running

* To run models: Execute `docker run --rm --cpus=1 sds timeout 10m time -v ...` with
    * PRISM: `prism -steadystate <model path> -const <model constants>`
    * Classic: `sds <model path> --const <model constants> solve`
    * Naive: `sds <model path> --const <model constants> approximate --precision 1e-4 --sampling SAMPLE_NAIVE`
    * Sample: `sds <model path> --const <model constants> approximate --precision 1e-4 --sampling SAMPLE_TARGET`

* Models from the paper:
    * `brp`: `models/dtmc/brp/brp.pm`, `N=64,MAX=5`
    * `nand`: `models/dtmc/nand/nand.pm`, `N=15,K=2`
    * `zeroconf_dl`: `models/mdp/zeroconf_dl/zeroconf_dl.nm`, `reset=false,deadline=40,N=1000,K=1`
    * `phil4`: `models/mdp/phil/phil4.nm` (no constants) (add `--explore` to the SDS call)
    * `rabin3`: `models/mdp/rabin/rabin3.nm` (no constants) (add `--explore` to the SDS call)
    * `branch`: `models/dtmc/crafted_branch/branch.pm`, `loops=1000,steps=100,p=0.008,levels=1000`

* Some further interesting models:
    * `brp`: `models/dtmc/brp/brp.pm`, `N=100,MAX=100`
    * `crowds`: `models/dtmc/crowds/crowds.pm`, `TotalRuns=3,CrowdSize=10`
    * `nand`: `models/dtmc/nand/nand.pm`, `N=10,K=5`
    * `loop`: `models/dtmc/crafted_loop/loop.pm`, `size=1000` (add `-explicit` to PRISM call)
    * `error`: `models/dtmc/crafted_error/error.pm`, `p=0.5,e=0.0000001` (add `-maxiters 1000000000` to PRISM the call to increase the maximal iterations and `-power` or `-jor` to see the wrong results under the respective iteration schemes)
      Note: PRISM's hybrid engine does eventually converge to the correct result, see below.
    * `error_hybrid`: `models/dtmc/crafted_error/error_hybrid.pm`, `e=0.0000001` (a variant model where the default hybrid scheme immediately terminates and yields a wrong solution)

* To execute the full benchmark suite run:
  `docker run --cpus=1 -t -v $(pwd)/results:/opt/results:rw sds python3 run.py --results /opt/results/data.json --timelimit 60 <further args>`
  (or `--timelimit 300` to fully replicate the experiments in the paper)
* This can also be configured to exclude some approaches or models, `--help` explains further arguments
* To view results, run `docker run -v $(pwd)/results:/opt/results:ro sds python3 run.py /opt/results/data.json`).
  Append `--only single_scc` or `--only one_state_scc` to apply the respective filters.
  Add `--scatter approx,solve` to obtain the data for scatter plots.
