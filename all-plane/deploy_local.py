#!/usr/bin/env python3
"""Install this verified source bundle and run a small local acceptance check.

No downloads, package installation, Git operations, or long searches are run.
The destination must be new or empty. Python 3.9+; macOS/Linux.
"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path, PurePosixPath
import platform
import shutil
import subprocess
import sys


def manifest_files(source):
    manifest = source / 'SHA256SUMS'
    entries = []
    seen = set()
    for line in manifest.read_text(encoding='utf-8').splitlines():
        digest, relative = line.split('  ', 1)
        part = PurePosixPath(relative)
        if (len(digest) != 64 or any(c not in '0123456789abcdef' for c in digest)
                or part.is_absolute() or '..' in part.parts or relative in seen):
            raise RuntimeError('Invalid manifest entry: ' + relative)
        path = source.joinpath(*part.parts)
        if not path.is_file() or path.is_symlink():
            raise RuntimeError('Missing or nonregular bundle file: ' + relative)
        if hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise RuntimeError('SHA-256 mismatch: ' + relative)
        entries.append((relative, path))
        seen.add(relative)
    required = {'local_run.py', 'all_plane_gap_v2.cpp', 'all_plane_gap_v3.cpp',
                'all_plane_gap_v4.cpp', 'all_plane_sieved.cpp', 'search.cpp',
                'run_sharded_v2.py', 'run_sharded_v3.py', 'run_sharded_v4.py',
                'research_v4_local/lifted_projection.hpp',
                'research_v4_math/critical_extension.hpp',
                'research_v4/final_100000.json', 'research_v4/final_three_100000.csv',
                'result_100000.json', 'LOCAL_RUN.md', 'LOCAL_HANDOFF.md'}
    if not required.issubset(seen):
        raise RuntimeError('Incomplete bundle: ' + ', '.join(sorted(required - seen)))
    return entries, hashlib.sha256(manifest.read_bytes()).hexdigest()


def main():
    default = (Path.home() / 'Files/03_ThinkAndPlay/WandsBox/Mathematics'
               / 'RationalDistance/rational_square_search')
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--destination', type=Path, default=default,
                        help='New or empty installation directory; default: %(default)s')
    parser.add_argument('--copy-only', action='store_true',
                        help='Copy verified sources without building or running')
    args = parser.parse_args()
    if platform.system() not in ('Darwin', 'Linux'):
        raise RuntimeError('Use macOS, Linux, or a Linux environment in WSL.')
    source = Path(__file__).resolve().parent
    destination = args.destination.expanduser().absolute()
    if destination.is_symlink():
        raise RuntimeError('Destination must not be a symbolic link.')
    destination = destination.resolve()
    if destination == source or source in destination.parents:
        raise RuntimeError('Choose a destination outside the extracted source folder.')
    if destination.exists() and (not destination.is_dir() or any(destination.iterdir())):
        raise RuntimeError('Destination is not empty; existing files were not changed: '
                           + str(destination))
    entries, manifest_hash = manifest_files(source)
    destination.mkdir(parents=True, exist_ok=True)
    for relative, src in entries:
        dst = destination / relative
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
    shutil.copy2(source / 'SHA256SUMS', destination / 'SHA256SUMS')
    receipt = {
        'created_at_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
        'destination': str(destination),
        'system': platform.system(), 'machine': platform.machine(),
        'python_version': platform.python_version(),
        'manifest_sha256': manifest_hash, 'source_file_count': len(entries),
        'copied_and_source_hashes_verified': True,
        'commands': [], 'status': 'copied_only',
        'full_search_started': False,
    }
    receipt_path = destination / 'deployment_receipt.json'
    log_path = destination / 'deployment.log'
    try:
        if not args.copy_only:
            with log_path.open('w', encoding='utf-8') as log:
                for action in ('check', 'build', 'smoke'):
                    command = [sys.executable, str(destination / 'local_run.py'), action]
                    print('Running:', action, flush=True)
                    log.write('\nCOMMAND: ' + json.dumps(command, ensure_ascii=False) + '\n')
                    log.flush()
                    result = subprocess.run(command, cwd=destination, stdout=log, stderr=log)
                    receipt['commands'].append({'action': action, 'returncode': result.returncode})
                    if result.returncode:
                        raise RuntimeError(action + ' failed; see ' + str(log_path))
            receipt['status'] = 'installed_and_smoke_verified'
    except (OSError, RuntimeError, KeyboardInterrupt) as exc:
        receipt['status'] = 'verification_incomplete'
        receipt['error'] = str(exc)
        raise
    finally:
        receipt_path.write_text(json.dumps(receipt, ensure_ascii=False, indent=2) + '\n',
                                encoding='utf-8')
    print(json.dumps(receipt, ensure_ascii=False, indent=2))


if __name__ == '__main__':
    try:
        main()
    except (OSError, ValueError, RuntimeError, KeyboardInterrupt) as error:
        print('Deployment did not complete: ' + str(error), file=sys.stderr)
        sys.exit(1)
