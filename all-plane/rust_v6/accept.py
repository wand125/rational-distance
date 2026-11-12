"""Sequential reference comparisons, with a <60 second whole-command deadline."""
import json
import pathlib
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
ROOT = HERE.parent
ACCEPT = HERE / "accept"
ACCEPT.mkdir(exist_ok=True)
DEADLINE = time.monotonic() + 54
CASES = {
    "100000": (["100000"], [1]),
    "200000_min150001": (["200000", "--min-denominator", "150001"], [1]),
    "100000_no_new_filters": (["100000", "--no-new-filters"], [1]),
    "200000_threads": (["200000"], [1, 4]),
}


def run(command, **kwargs):
    remaining = DEADLINE - time.monotonic()
    if remaining <= 0:
        raise RuntimeError("54-second harness deadline; rerun individual cases")
    return subprocess.run(command, cwd=ROOT, timeout=remaining, check=True, **kwargs)


def search(prefix, program, arguments):
    paths = {kind: ACCEPT / f"{prefix}.{kind}" for kind in ("three.csv", "parameters.csv", "stdout", "stderr")}
    # Every tested search bound is at most 200000.
    assert 1 <= int(arguments[0]) <= 200000
    with paths["stdout"].open("wb") as out, paths["stderr"].open("wb") as err:
        run([str(program), *arguments, "--audit-three", str(paths["three.csv"]),
             "--audit-parameters", str(paths["parameters.csv"])], stdout=out, stderr=err)
    rows = [json.loads(line) for line in paths["stdout"].read_text().splitlines()]
    return paths, rows


def compare(reference, rust, threads):
    ref_paths, ref_rows = reference
    rust_paths, rust_rows = rust
    for kind in ("three.csv", "parameters.csv"):
        run(["cmp", str(ref_paths[kind]), str(rust_paths[kind])], stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    assert ref_rows[:-1] == rust_rows[:-1], "solution JSON/order mismatch"
    ref, got = ref_rows[-1], rust_rows[-1]
    assert got["threads"] == threads and type(got["threads"]) is int
    assert ref["generator"] == "gap_CRT_v5" and got["generator"] == "gap_CRT_v6_rust"
    assert ref["reconstruction"] == "exact_u256_v5" and got["reconstruction"] == "exact_u256_v6_rust"
    assert [k for k in got if k != "threads"] == list(ref), "JSON keys/order mismatch"
    for key, expected in ref.items():
        if key in ("generator", "reconstruction") or key.startswith("seconds_"):
            continue
        assert type(got[key]) is type(expected) and got[key] == expected, (key, expected, got[key])


def main():
    requested = sys.argv[1:] or list(CASES)
    assert all(name in CASES for name in requested), f"cases: {', '.join(CASES)}"
    # Each invocation is a command lasting <60 seconds. A subset can be run
    # independently if a slower machine cannot finish all cases in that time.
    with (ACCEPT / "ACCEPT_LOG.txt").open("a", buffering=1) as log:
        def report(message):
            print(message, flush=True)
            print(message, file=log, flush=True)
        for name in requested:
            arguments, threads = CASES[name]
            try:
                reference = search(f"{name}.cpp", ROOT / "all_plane_gap_v5_local", arguments)
                rust_runs = []
                for count in threads:
                    rust = search(f"{name}.rust_t{count}", HERE / "target/release/rsq6",
                                  [*arguments, "--threads", str(count)])
                    compare(reference, rust, count)
                    rust_runs.append(rust)
                    report(f"PASS {name} threads={count}: both CSVs byte-identical; solutions and all non-timing JSON values/keys/order match")
                    report(f"TIMING {name} threads={count}: C++ seconds_search={reference[1][-1]['seconds_search']}; rsq6 seconds_search={rust[1][-1]['seconds_search']}")
                if len(rust_runs) == 2:
                    for kind in ("three.csv", "parameters.csv"):
                        run(["cmp", str(rust_runs[0][0][kind]), str(rust_runs[1][0][kind])])
                    report("PASS 200000 threads=1 vs threads=4: both CSVs byte-identical (both match C++)")
            except Exception as error:
                report(f"FAIL {name}: {error}")
                return 1
        report("PASS all requested acceptance cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
