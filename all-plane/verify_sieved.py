#!/usr/bin/env python3
"""Check disjoint interval coverage against a complete nonempty run.

The C++ --validate mode independently checks CRT candidate sets, finite-field
tables, and primitive three-distance points against the original evaluator.
This script additionally checks the min-denominator sharding interface.
"""
import csv
import json
from pathlib import Path
import subprocess
import tempfile


def main():
    binary = str(Path(__file__).with_name('all_plane_sieved').resolve())
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        sets = []
        reports = []
        for index, (low, high) in enumerate([(1, 1000), (1, 501), (502, 1000)]):
            audit = tmp / f'points_{index}.csv'
            args = [binary, str(high), '--min-denominator', str(low),
                    '--no-local', '--no-difference', '--audit-three', str(audit)]
            report = json.loads(subprocess.check_output(args, text=True))
            reports.append(report)
            with audit.open() as f:
                rows = [tuple(map(int, row)) for row in list(csv.reader(f))[1:]]
            assert len(rows) == len(set(rows)), 'duplicate three-distance witnesses'
            assert all(low <= row[0] <= high for row in rows), 'witness outside interval'
            sets.append(set(rows))
        assert sets[0], 'positive comparison set is empty'
        assert not sets[1].intersection(sets[2]), 'overlapping shards'
        assert sets[0] == sets[1].union(sets[2]), 'shards do not cover complete run'
        for key in ['dense_even_pairs', 'u_rows', 'generated_pairs',
                    'exact_evaluations', 'primitive_three_candidates', 'solutions']:
            assert reports[0][key] == reports[1][key] + reports[2][key], key
        result = {'passed': True, 'bound': 1000, 'split': [[1, 501], [502, 1000]],
                  'full_primitive_three_points': len(sets[0]),
                  'shard_primitive_three_points': [len(sets[1]), len(sets[2])],
                  'candidate_sets_equal': True, 'selected_counters_additive': True}
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
