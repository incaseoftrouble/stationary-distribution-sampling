import json
import sys
from collections import defaultdict

import tabulate

if __name__ == "__main__":
    with open(sys.argv[1]) as f:
        results = json.load(f)

    approaches = set()
    data = defaultdict(dict)
    for family, family_results in results.items():
        for instance, instance_data in family_results.get("results", {}).items():
            size = instance_data.get("size", 0)
            for approach, approach_result in instance_data.get("results", {}).items():
                approaches.add(approach)
                data[(family_results.get("type", "?"), family, size, instance)][approach] = approach_result
    approaches = list(sorted(approaches))
    header = ["model", "size"] + approaches
    rows = []
    for key in sorted(data.keys()):
        model_type, family, size, instance = key
        row = [f"{model_type}/{family}/{instance}", size]
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
                    elif error_type == "memory" or result['code'] == -9:
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
    print(tabulate.tabulate(rows, header))
