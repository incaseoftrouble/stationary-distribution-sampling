# Stationary Distribution Sampling

## Building

* Run `docker build -t sds .`

## Running

* To run models: Execute `docker run --rm --cpus=1 sds timeout 10m time -v ...` with
  * PRISM: `prism -steadystate <model path> -const <model constants>`
  * Classic: `sds <model path> --const <model constants> solve`
  * Naive: `sds <model path> --const <model constants> approximate --precision 1e-6 --sampling SAMPLE_NAIVE`
  * Sample: `sds <model path> --const <model constants> approximate --precision 1e-6 --sampling SAMPLE_TARGET`

* Used models:
  * `brp`: `models/brp/brp.pm`, `N=100,MAX=100`
  * `crowds`: `models/crowds/crowds.pm`, `TotalRuns=3,CrowdSize=10`
  * `nand`: `models/nand/nand.pm`, `N=10,K=5`
  * `loop`: 
  * `branch`: `models/crafted/branch.pm`, `loops=500,steps=100,p=0.001,levels=1000`