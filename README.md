# Stationary Distribution Sampling

## Building

* Run `docker build -t sds .`

## Running

* To run models: Execute `docker run --rm --cpus=1 sds timeout 10m time -v ...` with
  * PRISM: `prism -steadystate <model path> -const <model constants>`
  * Classic: `sds <model path> --const <model constants> solve`
  * Naive: `sds <model path> --const <model constants> approximate --precision 1e-6 --sampling SAMPLE_NAIVE`
  * Sample: `sds <model path> --const <model constants> approximate --precision 1e-6 --sampling SAMPLE_TARGET`

* Some used models:
  * `brp`: `models/brp/brp.pm`, `N=100,MAX=100`
  * `crowds`: `models/crowds/crowds.pm`, `TotalRuns=3,CrowdSize=10`
  * `nand`: `models/nand/nand.pm`, `N=10,K=5`
  * `loop`: `models/crafted/loop.pm`, `size=1000` (add `-explicit` to PRISM call)
  * `branch`: `models/crafted/branch.pm`, `loops=500,steps=100,p=0.008,levels=1000`
  * `error`: `models/crafted/error.pm`, `p=0.5,e=0.0000001` (add `-explicit` to PRISM call to see the wrong results)

* To execute the full suite (can also be configured to exclude some approaches or models):
  `docker run --cpus=1 -v $(pwd)/results:/opt/results:rw sds python3 run.py --results /opt/results/data.json --timelimit 60 <further args>`

Note: PRISM's hybrid engine does eventually converge to the correct result on `error`, run `models/crafted/error_hybrid.pm` with `p=0.1,e=0.0000001` for a model where the hybrid engine does not converge.