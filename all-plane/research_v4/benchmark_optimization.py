#!/usr/bin/env python3
"""Replay the archived prototype2 optimization experiment; no new sieves.

The current v4 has additional mathematical filters. Recompile the archived
source so rerunning this script neither changes that meaning nor overwrites the
historical optimization_benchmark.json.
"""
import hashlib
import json
import pathlib
import statistics
import subprocess

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "research_v4" / "optimization_replay.json"
ARCHIVE = ROOT / "research_v4" / "prototype2_source.txt"
BINARY = ROOT / "research_v4" / "prototype2_replay"
archive_record = json.loads((ROOT / "research_v4" / "prototype2_record.json").read_text())
assert hashlib.sha256(ARCHIVE.read_bytes()).hexdigest() == archive_record["source_sha256"]
assert hashlib.sha256((ROOT / "all_plane_gap_v3.cpp").read_bytes()).hexdigest() == archive_record["frozen_v3_sha256"]
subprocess.run(["g++", "-O3", "-std=c++17", "-x", "c++", "-I", str(ROOT),
                str(ARCHIVE), "-o", str(BINARY)], check=True)
MODES = {
    "v3_reference": ["--reference-v3"],
    "v4_late": ["--gamma-stage", "late"],
    "v4_generation": ["--gamma-stage", "generation"],
}
record = {
    "description": "Same binary, same mathematical predicates; three alternating rounds per bound, one process at a time.",
    "compiler": subprocess.check_output(["g++", "--version"], text=True).splitlines()[0],
    "source_sha256": {name: hashlib.sha256((ROOT / name).read_bytes()).hexdigest()
                      for name in ["research_v4/prototype2_source.txt", "all_plane_gap_v3.cpp", "all_plane_sieved.cpp"]},
    "runs": [],
}
for bound in [100000, 200000]:
    for round_number, modes in enumerate([
        ["v3_reference", "v4_late", "v4_generation"],
        ["v4_generation", "v4_late", "v3_reference"],
        ["v4_late", "v3_reference", "v4_generation"],
    ], 1):
        for mode in modes:
            command = [str(BINARY), str(bound), *MODES[mode]]
            completed = subprocess.run(command, text=True, capture_output=True, check=True)
            result = json.loads(completed.stdout.strip().splitlines()[-1])
            record["runs"].append({"bound": bound, "round": round_number, "mode": mode,
                                   "command": command, "result": result})
            OUT.write_text(json.dumps(record, indent=2) + "\n")
            print(bound, round_number, mode, result["seconds_search"],
                  result["exact_parameter_pairs"], flush=True)
    rows = [r for r in record["runs"] if r["bound"] == bound]
    for field in ["exact_parameter_pairs", "primitive_three_candidates", "solutions"]:
        assert len({r["result"][field] for r in rows}) == 1, (bound, field)

record["summary"] = {}
for bound in [100000, 200000]:
    rows = [r for r in record["runs"] if r["bound"] == bound]
    record["summary"][str(bound)] = {}
    for mode in MODES:
        times = [r["result"]["seconds_search"] for r in rows if r["mode"] == mode]
        reference = [next(r["result"]["seconds_search"] for r in rows
                          if r["mode"] == "v3_reference" and r["round"] == k)
                     for k in range(1, 4)]
        record["summary"][str(bound)][mode] = {
            "seconds": times,
            "median_seconds": statistics.median(times),
            "paired_speedup_ratios": [a / b for a, b in zip(reference, times)],
            "median_paired_speedup": statistics.median(a / b for a, b in zip(reference, times)),
        }
OUT.write_text(json.dumps(record, indent=2) + "\n")
print(json.dumps(record["summary"], indent=2), flush=True)
