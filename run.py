import abc
import argparse
import csv
import dataclasses
import json
import pathlib
import re
import subprocess
import time
from typing import List, Optional


@dataclasses.dataclass
class ModelInstance(object):
    size: Optional[int]
    path: pathlib.Path
    constants: str

    def __str__(self):
        return f"{self.path.name}/{self.constants}" if self.constants else self.path.name

    def __lt__(self, other):
        size = self.size if self.size else 0
        other_size = other.size if other.size else 0
        return (size, self.path.name, self.constants) < (other_size, other.path.name, other.constants)


@dataclasses.dataclass
class Model(object):
    model_type: str
    model_name: str
    instances: List[ModelInstance]


@dataclasses.dataclass
class Result(abc.ABC):
    output: str
    error: str
    time: float

    @abc.abstractmethod
    def to_json(self):
        pass


@dataclasses.dataclass
class Timeout(Result):
    def to_json(self):
        return {
            "type": "timeout"
        }

    def __str__(self):
        return f"Timeout"


@dataclasses.dataclass
class Error(Result):
    code: int
    error_type: str

    def to_json(self):
        return {
            "type": "error",
            "code": self.code,
            "time": self.time,
            "error": self.error_type
        }

    def __str__(self):
        return f"Error({self.error_type},{self.code})"


@dataclasses.dataclass
class Success(Result):
    def to_json(self):
        return {
            "type": "success",
            "time": self.time
        }


class Approach(abc.ABC):
    @abc.abstractmethod
    def supports(self, model, instance):
        pass

    @abc.abstractmethod
    def create_args(self, instance, uniformization_constant=None):
        pass


class SdsApproach(Approach):
    def __init__(self, name, executable, arguments, supported_models=None):
        self.name = name
        self.executable = executable
        self.arguments = arguments
        self.supported_models = supported_models

    def supports(self, model, instance):
        return self.supported_models is None or model.model_type in self.supported_models

    def create_args(self, instance, uniformization_constant=None):
        return [self.executable, instance.path, "--const", instance.constants,
                "--uniformization", "1.0" if uniformization_constant is None else str(uniformization_constant), *self.arguments]

    def __hash__(self):
        return hash(self.name)

    def __eq__(self, other):
        return self.name == other.name

    def __str__(self):
        return self.name


def run(arguments, timelimit):
    start = time.time()
    process = subprocess.Popen(arguments, stderr=subprocess.PIPE, stdout=subprocess.PIPE, encoding="utf-8", text=True,
                               env={"JAVA_OPTS": "-Xmx10G"})
    try:
        stdout, stderr = process.communicate(timeout=timelimit)
        duration = time.time() - start
        if process.returncode:
            if "java.lang.StackOverflowError" in stderr:
                error_type = "stackoverflow"
            elif "java.lang.OutOfMemoryError" in stderr:
                error_type = "memory"
            elif "Can not create DoubleVector" in stderr:
                error_type = "bounds"
            elif "Iterative method did not converge" in stderr:
                error_type = "convergence"
            else:
                print("Unrecognized error:")
                print(stderr)
                error_type = "generic"
            return Error(stdout, stderr, duration, process.returncode, error_type)
        return Success(stdout, stderr, duration)
    except subprocess.TimeoutExpired:
        process.kill()
        stdout = "\n".join(process.stdout.readlines())
        stderr = "\n".join(process.stderr.readlines())
        return Timeout(stdout, stderr, timelimit)


def execute(args):
    models_path = pathlib.Path(__file__).parent / "models"
    models = list()
    for model_type_folder in models_path.iterdir():
        if not model_type_folder.is_dir():
            continue

        model_type = model_type_folder.name
        for family_folder in model_type_folder.iterdir():
            if not family_folder.is_dir():
                continue
            models_file = family_folder / "models.csv"
            if not models_file.is_file():
                continue
            family_name = family_folder.name
            family_models = []

            updated_csv = []
            changed = False
            with models_file.open(mode="rt") as f:
                reader = csv.reader(f)
                updated_csv.append(next(reader))
                for row in reader:
                    name, constants, _, states = row[:4]
                    model_path = family_folder / name
                    if not states:
                        print(f"Computing states of {name}/{constants}")
                        result = run([args.executable, model_path, "--const", constants, "stats"], 3600)
                        if isinstance(result, Success):
                            states = int(result.output)
                            row[3] = str(states)
                            row[4] = str(result.time)
                            changed = True
                    updated_csv.append(row)
                    family_models.append(ModelInstance(int(states) if states else None, model_path, constants))
            if changed:
                with (models_file).open(mode="wt") as f:
                    csv.writer(f).writerows(updated_csv)
            models.append(Model(model_type, family_name, sorted(family_models)))

    if args.results is not None and args.results.exists():
        with args.results.open(mode="rt") as f:
            results = json.load(f)
    else:
        results = {}

    if args.approach:
        approach_names = set(args.approach)
    else:
        approach_names = {"solve", "approximate", "approximate-naive", "approximate-solve", "prism", "prism-explicit"}
    if args.exclude_approach:
        approach_names -= set(args.exclude_approach)
    approaches = set()
    if "solve" in approach_names:
        approaches.add(SdsApproach("solve", args.executable, ["solve"]))
    if "approximate" in approach_names:
        approaches.add(SdsApproach("approx", args.executable, ["approximate", "--precision", "1e-6", "--sampling", "SAMPLE_TARGET"]))
    if "approximate-naive" in approach_names:
        approaches.add(SdsApproach("approx-naive", args.executable, ["approximate", "--precision", "1e-6", "--sampling", "SAMPLE_NAIVE"]))
    if "approximate-solve" in approach_names:
        approaches.add(SdsApproach("approx-solve", args.executable, ["approximate", "--precision", "1e-6", "--sampling", "SAMPLE_TARGET", "--solve-bsccs"]))
    if "prism" in approach_names:
        approaches.add(SdsApproach("prism", args.executable, ["prism", "--precision", "1e-6"], {"ctmc", "dtmc"}))
    if "prism-explicit" in approach_names:
        approaches.add(SdsApproach("prism-explicit", args.executable, ["prism", "--precision", "1e-6", "--explicit"]))
    timeout = args.timelimit if args.timelimit else 60

    for model in models:
        if args.type and model.model_type not in args.type:
            continue
        if args.name and not any(pattern.matches(model.model_name) for pattern in args.name):
            continue

        if model.model_name not in results:
            results[model.model_name] = {}
        model_data = results[model.model_name]
        model_data["type"] = model.model_type
        if "results" not in model_data:
            model_data["results"] = {}
        model_results = model_data["results"]

        print(model.model_name)
        instances = model.instances
        skipped_approaches = set()
        if args.familylimit and len(instances) > args.familylimit:
            instances = instances[:args.familylimit]
        for instance in instances:
            instance_key = str(instance)
            print("  " + instance_key)
            if instance_key not in model_results:
                model_results[instance_key] = {}
            instance_data = model_results[instance_key]
            if instance.size:
                instance_data["size"] = instance.size

            if model.model_type == "ctmc":
                if "uniformization" in instance_data:
                    if instance_data["uniformization"] == "error":
                        break
                    uniformization_constant = instance_data["uniformization"]
                else:
                    result = run([args.executable, instance.path, "--const", instance.constants, "uniformization"], timeout)
                    if isinstance(result, Success):
                        uniformization_constant = float(result.output)
                        instance_data["uniformization"] = uniformization_constant
                    else:
                        print("    Finding uniformization constant failed")
                        instance_data["uniformization"] = "error"
                        break
            else:
                uniformization_constant = None

            if "results" not in instance_data:
                instance_data["results"] = {}
            instance_results = instance_data["results"]
            for approach in approaches:
                if approach.name in instance_results and instance_results[approach.name]["type"] == "skipped":
                    skipped_approaches.add(approach)
                    continue

                if args.rerun_error:
                    if approach.name in instance_results and (instance_results[approach.name]["type"] != "error"
                                                              or instance_results[approach.name].get("error", "") != "generic"):
                        continue
                else:
                    if approach.name in instance_results:
                        continue

                if approach in skipped_approaches or not approach.supports(model, instance):
                    instance_results[approach.name] = {"type": "skipped"}
                else:
                    print(f"    {approach}")
                    execution = approach.create_args(instance, uniformization_constant)
                    result = run(execution, timeout)
                    instance_results[approach.name] = result.to_json()
                    if not isinstance(result, Success):
                        print(f"      {approach} failed ({result}), skipping")
                        skipped_approaches.add(approach)
            if args.results is not None:
                with args.results.open(mode="wt") as f:
                    json.dump(results, f)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Run a (sub-)set of experiments')

    parser.add_argument('--type', help="Filter by model type (DTMC, CTMC, MDP)", nargs="*")
    parser.add_argument('--name', help="Filter by name", nargs="*", type=re.Pattern)
    parser.add_argument('--approach', help="Filter by approach (solve, approximate, prism, prism-explicit)", nargs="*")
    parser.add_argument('--exclude-approach', help="Filter by approach (solve, approximate, prism, prism-explicit)", nargs="*")
    parser.add_argument('--timelimit', help="Execution time limit", type=float)
    parser.add_argument('--familylimit', help="Limit number of executions per model family (same model with different constants)",
                        type=int)
    parser.add_argument('--results', help="Results file to load partial results", type=pathlib.Path)
    parser.add_argument('--executable', help="Path to SDS executable", type=str, default="sds")
    parser.add_argument('--rerun-error', help="Rerun failed experiments", default=False, action='store_true')
    execute(parser.parse_args())
