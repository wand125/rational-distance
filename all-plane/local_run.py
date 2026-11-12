#!/usr/bin/env python3
"""Small, bounded local launcher; the frozen search sources are never edited."""
import argparse
import csv
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import resource
import shlex
import shutil
import signal
import subprocess
import sys
import tempfile
import time

ROOT = Path(__file__).resolve().parent
SOURCES = ("all_plane_gap_v2.cpp", "all_plane_gap_v3.cpp", "all_plane_gap_v4.cpp",
           "all_plane_sieved.cpp", "research_v4_local/lifted_projection.hpp",
           "research_v4_math/critical_extension.hpp", "run_sharded_v2.py",
           "run_sharded_v3.py", "run_sharded_v4.py", "search.cpp")
TARGETS = (("all_plane_gap_v2_local", "all_plane_gap_v2.cpp", "c++17"),
           ("all_plane_gap_v3_local", "all_plane_gap_v3.cpp", "c++17"),
           ("all_plane_gap_v4_local", "all_plane_gap_v4.cpp", "c++17"),
           ("square_search_local", "search.cpp", "c++20"))
STATE = ROOT / "local_build.json"
SMOKE_BOUND = 100000
SMOKE_REFERENCE_JSON = "research_v4/final_100000.json"
SMOKE_REFERENCE_CSV = "research_v4/final_three_100000.csv"


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def source_hashes():
    return {name: digest(ROOT / name) for name in SOURCES}


def write_json(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")


def output_dir(label):
    base = ROOT / "local_runs"
    base.mkdir(exist_ok=True)
    stamp = datetime.datetime.now().astimezone().strftime("%Y%m%dT%H%M%S%z")
    return Path(tempfile.mkdtemp(prefix=label + "_" + stamp + "_", dir=str(base)))


def check_manifest():
    count = 0
    for line in (ROOT / "SHA256SUMS").read_text(encoding="utf-8").splitlines():
        expected, name = line.split("  ", 1)
        path = (ROOT / name).resolve()
        path.relative_to(ROOT)  # Reject a manifest path outside this project.
        if not path.is_file() or digest(path) != expected:
            raise RuntimeError("SHA256 mismatch or missing file: " + name)
        count += 1
    if not count:
        raise RuntimeError("Empty SHA256SUMS")
    print("SHA256 OK: {} packaged files".format(count), flush=True)
    return count


def child_limits():
    # The driver and its <=2 workers inherit these conservative limits.
    # RLIMIT_AS is intentionally Linux-only: macOS virtual address reservations
    # are not a useful upper bound on physical RAM.
    limits = [(resource.RLIMIT_CPU, 30), (resource.RLIMIT_FSIZE, 16 * 1024 * 1024)]
    if sys.platform.startswith("linux"):
        limits.append((resource.RLIMIT_AS, 512 * 1024 * 1024))
    for kind, cap in limits:
        _, hard = resource.getrlimit(kind)
        value = cap if hard == resource.RLIM_INFINITY else min(cap, hard)
        resource.setrlimit(kind, (value, value))


def execute(command, folder, name, deadline):
    remaining = deadline - time.monotonic()
    if remaining <= 0:
        raise RuntimeError("50-second local batch deadline reached; run is incomplete")
    print("Running: " + shlex.join([str(x) for x in command]), flush=True)
    with (folder / (name + ".stdout.txt")).open("x") as stdout, \
            (folder / (name + ".stderr.txt")).open("x") as stderr:
        process = subprocess.Popen(command, cwd=str(ROOT), stdout=stdout, stderr=stderr,
                                   start_new_session=True, preexec_fn=child_limits)
        try:
            code = process.wait(timeout=remaining)
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            process.wait()
            raise
    if code:
        raise RuntimeError("Command failed ({}); see {}".format(code, folder / (name + ".stderr.txt")))


def compiler_command():
    if os.environ.get("CXX"):
        command = shlex.split(os.environ["CXX"])
        if command and shutil.which(command[0]):
            return command
        raise RuntimeError("CXX does not name an available compiler")
    for name in ("clang++", "g++", "c++"):
        if shutil.which(name):
            return [shutil.which(name)]
    raise RuntimeError("No C++ compiler found; install a Clang/GCC toolchain first")


def build(deadline):
    folder, compiler, before = output_dir("build"), compiler_command(), source_hashes()
    record = {"system": platform.system(), "machine": platform.machine(),
              "python": platform.python_version(), "source_sha256": before, "targets": {}}
    for name, source, standard in TARGETS:
        fd, temporary = tempfile.mkstemp(prefix=name + "_", dir=str(ROOT))
        os.close(fd)
        try:
            command = compiler + ["-O3", "-std=" + standard, str(ROOT / source), "-o", temporary]
            execute(command, folder, name, deadline)
            path = ROOT / name
            os.chmod(temporary, 0o755)
            os.replace(temporary, path)
            record["targets"][name] = {"sha256": digest(path), "command": command}
        finally:
            Path(temporary).unlink(missing_ok=True)
    if source_hashes() != before:
        raise RuntimeError("Source changed during build")
    write_json(folder / "build.json", record)
    write_json(STATE, record)
    print("Build complete: " + str(folder), flush=True)


def require_build():
    if not STATE.exists():
        raise RuntimeError("Run 'python3 local_run.py build' first")
    record = json.loads(STATE.read_text(encoding="utf-8"))
    if record["source_sha256"] != source_hashes() or \
            (record["system"], record["machine"]) != (platform.system(), platform.machine()):
        raise RuntimeError("Source or platform changed; rebuild first")
    for name, _, _ in TARGETS:
        path = ROOT / name
        if not path.exists() or digest(path) != record["targets"][name]["sha256"]:
            raise RuntimeError("Local executable changed or missing; rebuild first: " + name)


def sharded(bound, workers, folder, deadline, engine="v4"):
    require_build()
    result = folder / "shards"
    command = [sys.executable, "-I", str(ROOT / ("run_sharded_" + engine + ".py")),
               "--bound", str(bound), "--workers", str(workers),
               "--binary", str(ROOT / ("all_plane_gap_" + engine + "_local")), "--output", str(result)]
    execute(command, folder, "driver", deadline)
    summary = json.loads((result / "summary.json").read_text())
    if not summary.get("all_intervals_completed") or not summary.get("source_and_binary_unchanged_during_run"):
        raise RuntimeError("Driver did not confirm a completed, stable run")
    if engine in ("v3", "v4") and (summary.get("generator_version") != engine or
                                       not summary.get("generator_configuration_checked")):
        raise RuntimeError("Driver did not verify the requested generator version and configuration")
    return result, summary


def points(path):
    with Path(path).open(newline="") as source:
        rows = [tuple(int(row[k]) for k in ("N", "p", "q", "X", "Y", "A"))
                for row in csv.DictReader(source)]
    if len(rows) != len(set(rows)):
        raise RuntimeError("Duplicate three-distance candidate")
    return set(rows)


def smoke(deadline):
    require_build()
    reference_paths = (SMOKE_REFERENCE_JSON, SMOKE_REFERENCE_CSV)
    reference_hashes = {name: digest(ROOT / name) for name in reference_paths}
    expected_stats = json.loads((ROOT / SMOKE_REFERENCE_JSON).read_text(encoding="utf-8"))
    expected = points(ROOT / SMOKE_REFERENCE_CSV)
    if (expected_stats.get("generator") != "gap_CRT_v4" or
            expected_stats.get("scope") != "whole_plane" or
            expected_stats.get("min_denominator") != 1 or
            expected_stats.get("max_minimal_coordinate_denominator") != SMOKE_BOUND or
            expected_stats.get("primitive_three_candidates") != len(expected) or not expected):
        raise RuntimeError("The canonical v4 smoke reference is incomplete or inconsistent")
    folder = output_dir("smoke")
    execute([str(ROOT / "all_plane_gap_v4_local"), str(SMOKE_BOUND), "--validate"],
            folder, "generator_validation", deadline)
    result, summary = sharded(SMOKE_BOUND, 2, folder, deadline)
    actual = points(result / "three_candidates.csv")
    if actual != expected:
        raise RuntimeError("Nonempty three-distance reference set did not match")
    for key in ("exact_parameter_pairs", "primitive_three_candidates", "solutions"):
        if summary["counts"][key] != expected_stats[key]:
            raise RuntimeError("v4 canonical reference count did not match: " + key)
    for n, p, q, x, y, a in actual:
        if not (math.gcd(math.gcd(x, y), n) == 1 and x*x+y*y == a*a and
                (n-x)**2+y*y == (a+p)**2 and x*x+(n-y)**2 == (a+q)**2):
            raise RuntimeError("Exact integer witness check failed")
    execute([str(ROOT / "square_search_local"), "--limit", str(SMOKE_BOUND)],
            folder, "interior", deadline)
    interior = json.loads((folder / "interior.stdout.txt").read_text())
    reference = json.loads((ROOT / "result_100000.json").read_text())
    for key, value in reference.items():
        if key.startswith("seconds_") or key == "peak_rss_kib":
            continue
        if interior.get(key) != value:
            raise RuntimeError("Interior reference mismatch: " + key)
    require_build()
    if {name: digest(ROOT / name) for name in reference_paths} != reference_hashes:
        raise RuntimeError("Canonical v4 reference changed during smoke verification")
    report = {"status": "passed", "system": platform.system(), "machine": platform.machine(),
              "engine": "v4", "bound": SMOKE_BOUND, "workers": 2,
              "exact_parameter_pairs": summary["counts"]["exact_parameter_pairs"],
              "matched_three_distance_candidates": len(actual),
              "four_distance_solutions": summary["counts"]["solutions"],
              "nonempty_reference_set_matches": True, "interior_reference_matches": True,
              "canonical_reference_sha256": reference_hashes,
              "source_sha256": source_hashes(), "user_mac_verified": platform.system() == "Darwin"}
    write_json(folder / "smoke.json", report)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    print("Smoke check complete: " + str(folder), flush=True)


def main():
    if platform.system() not in ("Darwin", "Linux"):
        raise RuntimeError("Use macOS or Linux; native Windows is not supported by this POSIX launcher")
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    for name in ("check", "build", "smoke"):
        commands.add_parser(name)
    run = commands.add_parser("run")
    run.add_argument("--bound", type=int, default=100000)
    run.add_argument("--workers", type=int, choices=(1, 2), default=2)
    run.add_argument("--engine", choices=("v2", "v3", "v4"), default="v4")
    args = parser.parse_args()
    deadline = time.monotonic() + 50
    if args.command == "check":
        check_manifest()
    elif args.command == "build":
        build(deadline)
    elif args.command == "smoke":
        smoke(deadline)
    else:
        if not 24 <= args.bound <= 2000000:
            parser.error("bound must be 24..2000000")
        folder = output_dir("run")
        result, summary = sharded(args.bound, args.workers, folder, deadline, args.engine)
        print("Completed: " + str(result / "summary.json"), flush=True)
        print("Four-distance solutions: " + str(summary["counts"]["solutions"]), flush=True)


if __name__ == "__main__":
    try:
        main()
    except (RuntimeError, ValueError, OSError, subprocess.TimeoutExpired) as error:
        print("ERROR: " + str(error), file=sys.stderr)
        sys.exit(1)
