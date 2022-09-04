import abc
import argparse
import csv
import dataclasses
import json
import pathlib
import re
import subprocess
import time
from typing import List, Optional, Dict

ALL_APPROACHES = {"solve", "approximate", "approximate-naive", "approximate-solve", "prism", "prism-explicit"}


@dataclasses.dataclass
class ModelInstance(object):
    data: Dict
    path: pathlib.Path
    constants: str

    def __str__(self):
        return f"{self.path.name}/{self.constants}" if self.constants else self.path.name

    def __lt__(self, other):
        size = self.data["states"] if self.data.get("states", None) else 0
        other_size = other.data["states"] if other.data.get("states", None) else 0
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
            elif "ArrayIndexOutOfBounds" in stderr:
                error_type = "memory"
            elif "NegativeArraySizeException" in stderr:
                error_type = "memory"
            elif "IllegalArgumentException" in stderr:
                error_type = "internal"
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

        model_type = model_type_folder.name.lower()
        for family_folder in model_type_folder.iterdir():
            family_folder: pathlib.Path
            if not family_folder.is_dir():
                continue
            models_file = family_folder / "models.csv"
            if not models_file.is_file():
                continue
            models_data_file = family_folder / "data.json"
            if models_data_file.exists():
                try:
                    with models_data_file.open(mode="rt") as f:
                        models_data = json.load(f)
                except json.DecodeError:
                    models_data_file.unlink()
            else:
                models_data = {}

            family_name = family_folder.name
            family_models = []

            updated_csv = []
            changed = False
            with models_file.open(mode="rt") as f:
                reader = csv.reader(f)
                updated_csv.append(next(reader))

                for row in reader:
                    name, constants, _ = row[:3]
                    model_path = family_folder / name

                    identifier = f"{name}/{constants}"

                    if identifier not in models_data:
                        models_data[identifier] = {}
                    data = models_data[identifier]

                    if not row[3]:
                        print(f"Computing states of {identifier} ... ", end="")
                        result = run([args.executable, model_path, "--const", constants, "stats"], 300)
                        if isinstance(result, Success):
                            print(result.output.strip())
                            row[3] = str(int(result.output))
                            row[4] = f"{result.time:.2f}"
                            changed = True
                        else:
                            print("failed")
                    updated_csv.append(row)

                    data["states"] = int(row[3]) if row[3] else "error"
                    if "components" not in data:
                        print(f"Computing components of {name}/{constants} ... ", end="")
                        result = run([args.executable, model_path, "--const", constants, "stats", "--components"], 300)
                        if isinstance(result, Success):
                            print(",".join(result.output.split()))
                            _, components, largest_component = map(int, result.output.split())
                            data["components"] = components
                            data["largest_component"] = largest_component
                        else:
                            print("failed")
                            data["components"] = "error"
                            data["largest_component"] = "error"

                    if model_type.lower() == "ctmc" and "uniformization" not in data:
                        print(f"Computing uniformization constant of {name}/{constants} ... ", end="")
                        result = run([args.executable, model_path, "--const", constants, "uniformization"], 300)
                        if isinstance(result, Success):
                            data["uniformization"] = float(result.output)
                            print("done")
                        else:
                            data["uniformization"] = "error"
                            print("failed")

                    family_models.append(ModelInstance(data, model_path, constants))
            with models_data_file.open(mode="wt") as f:
                json.dump(models_data, f)
            if changed:
                with models_file.open(mode="wt") as f:
                    csv.writer(f).writerows(updated_csv)
            models.append(Model(model_type, family_name, sorted(family_models)))

    if args.results is not None and args.results.exists():
        with args.results.open(mode="rt") as f:
            results = json.load(f)
    else:
        results = {}

    approach_names = set(args.approach)
    if args.exclude_approach:
        approach_names -= set(args.exclude_approach)
    approaches = set()
    if "solve" in approach_names:
        approaches.add(SdsApproach("solve", args.executable, ["solve"]))
    if "approximate" in approach_names:
        approaches.add(SdsApproach("approx", args.executable,
                                   ["approximate", "--precision", "1e-4", "--sampling", "SAMPLE_TARGET"]))
    if "approximate-naive" in approach_names:
        approaches.add(SdsApproach("approx-naive", args.executable,
                                   ["approximate", "--precision", "1e-4", "--sampling", "SAMPLE_NAIVE"]))
    if "approximate-solve" in approach_names:
        approaches.add(SdsApproach("approx-solve", args.executable,
                                   ["approximate", "--precision", "1e-6", "--sampling", "SAMPLE_TARGET", "--solve-bsccs"]))
    if "prism" in approach_names:
        approaches.add(SdsApproach("prism", args.executable, ["prism", "--precision", "1e-4"], {"ctmc", "dtmc"}))
    if "prism-explicit" in approach_names:
        approaches.add(SdsApproach("prism-explicit", args.executable, ["prism", "--precision", "1e-4", "--explicit"]))
    timeout = args.timelimit

    for model in models:
        if args.type and model.model_type not in args.type:
            continue
        if args.name and not any(re.search(pattern, model.model_name) for pattern in args.name):
            continue

        if model.model_name not in results:
            results[model.model_name] = {}
        model_data = results[model.model_name]
        model_data["type"] = model.model_type
        if "results" not in model_data:
            model_data["results"] = {}
        model_results = model_data["results"]

        print(model.model_name)

        instances: List[ModelInstance] = model.instances
        skipped_approaches = set()
        if args.familylimit and len(instances) > args.familylimit:
            instances = instances[:args.familylimit]
        for instance in instances:
            instance: ModelInstance
            instance_key = str(instance)

            single_scc = instance.data["states"] == instance.data["largest_component"]
            if single_scc:
                print(f"  {instance_key} (single SCC)")
            else:
                print(f"  {instance_key}")
            if instance_key not in model_results:
                model_results[instance_key] = {}
            instance_data = model_results[instance_key]
            instance_data["data"] = instance.data

            if "results" not in instance_data:
                instance_data["results"] = {}
            instance_results = instance_data["results"]

            if "uniformization" in instance.data:
                if instance.data["uniformization"] == "error":
                    print("    Missing uniformization constant")
                    continue
                uniformization_constant = float(instance.data["uniformization"])
            else:
                uniformization_constant = None

            instance_approaches = set(approaches)
            if single_scc:
                non_approx_approaches = {approach for approach in instance_approaches
                                         if not isinstance(approach, SdsApproach) or approach.arguments[0] != "approximate"}
                if non_approx_approaches != instance_approaches:
                    non_approx_approaches.add(SdsApproach("approx", args.executable, ["approximate", "--precision", "1e-4",
                                                                                      "--sampling", "SAMPLE_TARGET", "--explore"]))
                    instance_approaches = non_approx_approaches

            for approach in instance_approaches:
                if args.rerun and approach.name in instance_results:
                    del instance_results[approach.name]

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

                if approach in skipped_approaches:
                    instance_results[approach.name] = {"type": "skipped"}
                    continue

                if not approach.supports(model, instance):
                    continue
                execution = approach.create_args(instance, uniformization_constant)
                print(f"    {approach}")
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
    parser.add_argument('--name', help="Filter by name", nargs="*")
    parser.add_argument('--approach', help=f"Filter by approach ({', '.join(ALL_APPROACHES)})", nargs="*",
                        default=list(ALL_APPROACHES))
    parser.add_argument('--exclude-approach', help="Exclude approaches", nargs="*")
    parser.add_argument('--timelimit', help="Execution time limit in seconds", type=float, default=60)
    parser.add_argument('--familylimit', help="Limit number of executions per model family (same model with different constants)",
                        type=int)
    parser.add_argument('--results', help="Results file to load partial results and write results", type=pathlib.Path)
    parser.add_argument('--executable', help="Path to SDS executable", type=str, default="sds")
    parser.add_argument('--rerun-error', help="Rerun failed experiments", default=False, action='store_true')
    parser.add_argument('--rerun', help="Rerun all experiments", default=False, action='store_true')
    execute(parser.parse_args())
