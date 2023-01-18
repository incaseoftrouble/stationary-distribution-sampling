import json
import pathlib
import sys
import csv
from collections import defaultdict
import argparse

import tabulate

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Evaluate results")

    parser.add_argument("data", help="Results file", type=pathlib.Path)
    parser.add_argument("--only", help="Filter by type (single_scc, one_state_scc)")
    parser.add_argument("--latex", help="LaTeX output", action="store_true")
    parser.add_argument(
        "--scatter", help="Instead, write the scatter plot data", type=str
    )
    parser.add_argument(
        "--scatter-to",
        help="Value to replace timeout with in scatter plot",
        type=int,
        default=600,
    )
    args = parser.parse_args()

    with args.data.open(mode="rt") as f:
        results = json.load(f)

    approaches = set()
    data = defaultdict(dict)
    for family, family_results in results.items():
        for instance, instance_data in family_results.get("results", {}).items():
            d = instance_data["data"]
            size = d.get("states", 0)
            if size == "error":
                size = -1
            components = d.get("components", 0)
            if components == "error":
                components = -1
            largest_component = d.get("largest_component", 0)
            if largest_component == "error":
                largest_component = -1
            if args.only == "one_state_scc":
                if largest_component != 1:
                    continue
            elif args.only == "single_scc":
                if largest_component != size:
                    continue

            for approach, approach_result in instance_data.get("results", {}).items():
                approaches.add(approach)
                data[
                    (
                        family_results["type"].lower(),
                        family,
                        size,
                        components,
                        largest_component,
                        instance,
                    )
                ][approach] = approach_result
    approaches = list(sorted(approaches))

    if args.scatter:
        if len(args.scatter.split(",")) != 2:
            sys.exit("Give scatter as <left>,<right>")
        left_name, right_name = args.scatter.split(",")
        if left_name not in approaches:
            sys.exit(f"{left_name} not known")
        if right_name not in approaches:
            sys.exit(f"{right_name} not known")

        writer = csv.writer(sys.stdout)
        writer.writerow([left_name, right_name])
        for key in data.keys():
            model_type, family, size, components, largest_component, instance = key

            left, right = data[key].get(left_name, None), data[key].get(
                right_name, None
            )
            if (
                left is None
                or left.get("type", None) not in {"success", "timeout", "skipped"}
                or right is None
                or right.get("type", None) not in {"success", "timeout", "skipped"}
            ):
                continue
            if left["type"] != "success" and right["type"] != "success":
                continue

            writer.writerow(
                [
                    round(left["time"], 3)
                    if left["type"] == "success"
                    else args.scatter_to,
                    round(right["time"], 3)
                    if right["type"] == "success"
                    else args.scatter_to,
                ]
            )
        sys.exit()

    header = [
        "type",
        "instance",
        "size",
        "components",
        "largest_component",
    ] + approaches
    rows = []
    for key in sorted(data.keys()):
        model_type, family, size, components, largest_component, instance = key
        row = [model_type, f"{family}/{instance}", size, components, largest_component]
        for approach in approaches:
            result = data[key].get(approach, None)
            if result is None or "type" not in result:
                row.append("?")
            else:
                result_type = result["type"]
                if result_type == "error":
                    error_type = result.get("error", "generic")
                    if error_type == "stackoverflow":
                        row.append("Stack Overflow")
                    elif error_type == "memory" or result["code"] == -9:
                        row.append("M/O")
                    elif error_type == "bounds":
                        row.append("Internal error")
                    elif error_type == "convergence":
                        row.append("Convergence error")
                    else:
                        row.append(f"E({result['code']})")
                elif result_type == "timeout":
                    row.append("T/O")
                elif result_type == "skipped":
                    row.append("-")
                elif result_type == "success":
                    row.append(f"{result['time']:.2f}")
                else:
                    row.append(f"?{result_type}?")
        rows.append(row)
    print(tabulate.tabulate(rows, header, tablefmt="latex" if args.latex else "simple"))
