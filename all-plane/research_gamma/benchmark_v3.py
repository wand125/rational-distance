#!/usr/bin/env python3
"""Bounded sequential paired benchmark of the frozen final v3 binary."""
import hashlib
import json
from pathlib import Path
import statistics
import subprocess
import sys

ROOT = Path(__file__).resolve().parent.parent
OUT = Path(__file__).resolve().parent
COMMANDS = {
    'baseline': ['--no-gamma', '--no-critical-content', '--no-critical-projections',
                 '--no-field-wheel19'],
    'refined': [],
}


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    binary = ROOT / 'all_plane_gap_v3'
    paths = [ROOT / 'all_plane_gap_v3.cpp', ROOT / 'all_plane_sieved.cpp', binary]
    before = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    timings = {name: [] for name in COMMANDS}
    counts = {}
    for repetition in range(3):
        order = list(COMMANDS) if repetition % 2 == 0 else list(reversed(COMMANDS))
        for name in order:
            command = [str(binary), '100000'] + COMMANDS[name]
            result = subprocess.run(command, capture_output=True, text=True, timeout=30, check=True)
            record = json.loads(result.stdout.splitlines()[-1])
            (OUT / ('v3_paired_{}_{}.json'.format(repetition, name))).write_text(
                json.dumps(record, indent=2) + '\n')
            value = {k: v for k, v in record.items() if isinstance(v, int) and not isinstance(v, bool)}
            if name in counts and counts[name] != value:
                raise RuntimeError('Nondeterministic integer counts')
            counts[name] = value
            timings[name].append(record['seconds_total_excluding_validation'])
    after = {str(p.relative_to(ROOT)): sha(p) for p in paths}
    if before != after:
        raise RuntimeError('Source or binary changed during paired benchmark')
    medians = {name: statistics.median(values) for name, values in timings.items()}
    report = {
        'bound': 100000, 'workers': 1, 'paired_repetitions': 3,
        'scope': 'whole_plane_minimal_coordinate_denominator',
        'commands': COMMANDS, 'hashes': before, 'sources_and_binary_stable': True,
        'seconds_total_excluding_validation': timings, 'median_seconds': medians,
        'median_speed_ratio_baseline_over_refined': medians['baseline'] / medians['refined'],
        'counts': counts,
        'generated_pair_reduction_fraction': 1 - counts['refined']['generated_pairs'] / counts['baseline']['generated_pairs'],
        'exact_parameter_reduction_fraction': 1 - counts['refined']['exact_parameter_pairs'] / counts['baseline']['exact_parameter_pairs'],
        'caveat': 'Shared Linux environment; these are measured times, not a speed guarantee.',
    }
    (OUT / 'v3_comparison.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
