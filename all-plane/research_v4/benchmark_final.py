#!/usr/bin/env python3
"""Same-binary new-filter comparison plus the frozen external v3 executable."""
import hashlib
import json
import pathlib
import statistics
import subprocess
from datetime import datetime, timezone

ROOT = pathlib.Path(__file__).resolve().parent.parent
# Preserve the completed release measurements in final_benchmark.json.
OUT = ROOT / "research_v4" / "final_benchmark_rerun.json"
names = ["all_plane_gap_v4.cpp", "all_plane_gap_v3.cpp", "all_plane_sieved.cpp",
         "research_v4_local/lifted_projection.hpp", "research_v4_math/critical_extension.hpp"]
binary_names = ["all_plane_gap_v3", "all_plane_gap_v4"]
def hashes(paths):
    return {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest() for name in paths}

record = {
    "description": "Alternating single-process timing; exact integers and no CSV output in timed runs.",
    "compiler": subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
    "started_utc": datetime.now(timezone.utc).isoformat(),
    "source_sha256": hashes(names),
    "binary_sha256": hashes(binary_names),
    "mark_tuning": [], "runs": [],
}

def run(binary, bound, flags):
    command = [str(ROOT / binary), str(bound), *flags]
    completed = subprocess.run(command, text=True, capture_output=True, check=True)
    return command, json.loads(completed.stdout.strip().splitlines()[-1])

tuning = {"fast_mark": [], "direct_mark": ["--slow-five-adic-mark"],
          "late_fifth": ["--five-adic-late"]}
for round_number, order in enumerate([
    ["fast_mark", "direct_mark", "late_fifth"],
    ["late_fifth", "direct_mark", "fast_mark"],
    ["direct_mark", "fast_mark", "late_fifth"],
], 1):
    for mode in order:
        command, result = run("all_plane_gap_v4", 100000, tuning[mode])
        record["mark_tuning"].append({"round": round_number, "mode": mode,
                                      "command": command, "result": result})
        print("mark", round_number, mode, result["seconds_search"], flush=True)
        OUT.write_text(json.dumps(record, indent=2) + "\n")
assert len({r["result"]["exact_parameter_pairs"] for r in record["mark_tuning"]}) == 1
record["mark_tuning_medians"] = {
    mode: statistics.median(r["result"]["seconds_search"] for r in record["mark_tuning"] if r["mode"] == mode)
    for mode in tuning
}

modes = {
    "external_v3": ("all_plane_gap_v3", []),
    "v4_no_new_filters": ("all_plane_gap_v4", ["--no-new-filters"]),
    "v4_default": ("all_plane_gap_v4", []),
}
for bound in [100000, 200000]:
    for round_number, order in enumerate([
        ["external_v3", "v4_no_new_filters", "v4_default"],
        ["v4_default", "v4_no_new_filters", "external_v3"],
        ["v4_no_new_filters", "external_v3", "v4_default"],
    ], 1):
        for mode in order:
            binary, flags = modes[mode]
            command, result = run(binary, bound, flags)
            record["runs"].append({"bound": bound, "round": round_number, "mode": mode,
                                   "command": command, "result": result})
            OUT.write_text(json.dumps(record, indent=2) + "\n")
            print(bound, round_number, mode, result["seconds_search"],
                  result["exact_parameter_pairs"], flush=True)
    rows = [r for r in record["runs"] if r["bound"] == bound]
    for mode in modes:
        for field in ["exact_parameter_pairs", "primitive_three_candidates", "solutions"]:
            assert len({r["result"][field] for r in rows if r["mode"] == mode}) == 1
    assert {r["result"]["exact_parameter_pairs"] for r in rows if r["mode"] != "v4_default"}.__len__() == 1

record["summary"] = {}
for bound in [100000, 200000]:
    rows = [r for r in record["runs"] if r["bound"] == bound]
    summary = {}
    for mode in modes:
        times = [next(r["result"]["seconds_search"] for r in rows if r["mode"] == mode and r["round"] == k)
                 for k in range(1, 4)]
        summary[mode] = {"seconds": times, "median_seconds": statistics.median(times)}
    for baseline in ["external_v3", "v4_no_new_filters"]:
        ratios = [a / b for a, b in zip(summary[baseline]["seconds"], summary["v4_default"]["seconds"])]
        summary["v4_default"]["paired_speedup_vs_" + baseline] = ratios
        summary["v4_default"]["median_paired_speedup_vs_" + baseline] = statistics.median(ratios)
    record["summary"][str(bound)] = summary
record["hash_verification"] = {
    "verified_utc": datetime.now(timezone.utc).isoformat(),
    "verification_timing": "at_end_of_benchmark",
    "source_sha256": hashes(names),
    "binary_sha256": hashes(binary_names),
}
record["hash_verification"]["sources_match_recorded_start"] = (
    record["hash_verification"]["source_sha256"] == record["source_sha256"])
record["hash_verification"]["binaries_match_recorded_start"] = (
    record["hash_verification"]["binary_sha256"] == record["binary_sha256"])
OUT.write_text(json.dumps(record, indent=2) + "\n")
assert record["hash_verification"]["sources_match_recorded_start"], "source changed during benchmark"
assert record["hash_verification"]["binaries_match_recorded_start"], "binary changed during benchmark"
print(json.dumps({"mark_tuning_medians": record["mark_tuning_medians"], "summary": record["summary"]}, indent=2), flush=True)
