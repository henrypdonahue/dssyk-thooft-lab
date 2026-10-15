#!/usr/bin/env python3
"""Run the Makefile's validation commands in a preserved, logged snapshot.

Only the existing `cd MODULE && $(PY) ARGS` recipe form is accepted.
No shell interpolation is used. Continue after failures to expose the full
result; return nonzero if any command fails. Run with the pinned environment.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import platform
import shlex
import shutil
import subprocess
import sys
import time
from datetime import datetime, timezone

TARGETS = ('test', 'test-all', 'reproduce', 'reproduce-slow')


def recipes(root, target):
    active = False
    out = []
    for line in (root / 'Makefile').read_text().splitlines():
        if line and not line[0].isspace() and not line.startswith('#'):
            active = line.startswith(target + ':')
        if active and line.startswith('\t'):
            tokens = shlex.split(line.strip())
            if len(tokens) < 5 or tokens[0] != 'cd' or tokens[2:4] != ['&&', '$(PY)']:
                raise ValueError('Unsupported recipe; review runner before executing: ' + line)
            module = Path(tokens[1])
            if module.is_absolute() or '..' in module.parts:
                raise ValueError('Recipe directory must stay within repository')
            out.append((str(module), tokens[4:]))
    if not out:
        raise ValueError('No recipes found for ' + target)
    return out


def git(root, *args):
    p = subprocess.run(['git', '-C', str(root), *args], capture_output=True, text=True)
    return p.stdout.strip() if p.returncode == 0 else None


# Fields that record wall-clock cost rather than results.
RUNTIME_KEYS = ('runtime_s', 'runtime', 'elapsed_s', 'seconds')

# A value below this fraction of the largest magnitude in its own record
# (the nearest enclosing JSON object) is a roundoff zero and is not
# compared: e.g. the dynamical histogram of an exactly conserved channel,
# 1e-22 against a total weight of 1.
ROUNDOFF_FLOOR = 1e-14


def _is_number(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _magnitudes(value):
    if _is_number(value):
        return [abs(value)] if value == value else []   # skip NaN
    if isinstance(value, dict):
        value = list(value.values())
    if isinstance(value, list):
        return [m for v in value for m in _magnitudes(v)]
    return []


def _record(out, rel, path):
    if rel > out['max_rel']:
        out['max_rel'], out['max_rel_at'] = rel, path


def _walk(old, new, path, out, floor):
    """Collect the numeric relative differences and every non-numeric
    change between two parsed JSON values."""
    if isinstance(old, dict) and isinstance(new, dict):
        floor = ROUNDOFF_FLOOR * max(_magnitudes(old) + _magnitudes(new) + [0.0])
        for key in sorted(set(old) | set(new)):
            if key in RUNTIME_KEYS:
                continue
            if key not in old or key not in new:
                out['structure'].append(f'{path}/{key}')
            else:
                _walk(old[key], new[key], f'{path}/{key}', out, floor)
    elif isinstance(old, list) and isinstance(new, list):
        if len(old) != len(new):
            out['structure'].append(f'{path} (length {len(old)} -> {len(new)})')
            return
        if old and all(_is_number(v) for v in old + new):
            # a numeric array (a spectrum, a histogram, a series): measure each
            # difference against the array's largest entry
            scale = max(abs(v) for v in old + new)
            if scale > floor:
                for i, (a, b) in enumerate(zip(old, new)):
                    _record(out, abs(a - b) / scale, f'{path}[{i}]')
            return
        for i, (a, b) in enumerate(zip(old, new)):
            _walk(a, b, f'{path}[{i}]', out, floor)
    elif _is_number(old) and _is_number(new):
        scale = max(abs(old), abs(new))
        if old != new and scale > floor:
            _record(out, abs(old - new) / scale, path)
    elif old != new:
        out['text'].append(path)


def _parse_csv(text):
    rows = []
    for line in text.splitlines():
        cells = []
        for cell in line.split(','):
            try:
                cells.append(float(cell))
            except ValueError:
                cells.append(cell)
        rows.append(cells)
    return rows


def compare_artifact(old_path, new_path):
    """Compare one committed artifact with its regenerated copy.  JSON and
    CSV are compared value by value: numbers by relative difference
    (arrays against their largest entry, roundoff zeros skipped), runtime
    fields ignored.  Anything else is reported as changed bytes."""
    suffix = old_path.suffix.lower()
    if suffix not in ('.json', '.csv'):
        return dict(kind='binary', changed=True)
    try:
        if suffix == '.json':
            old, new = json.loads(old_path.read_text()), json.loads(new_path.read_text())
        else:
            old, new = _parse_csv(old_path.read_text()), _parse_csv(new_path.read_text())
    except (ValueError, UnicodeDecodeError) as exc:
        return dict(kind=suffix[1:], error=f'unparseable: {exc}')
    out = dict(max_rel=0.0, max_rel_at=None, structure=[], text=[])
    top_floor = ROUNDOFF_FLOOR * max(_magnitudes(old) + _magnitudes(new) + [0.0])
    _walk(old, new, '', out, top_floor)
    return dict(kind=suffix[1:], max_rel=out['max_rel'], max_rel_at=out['max_rel_at'],
                structure_changes=out['structure'][:20], text_changes=out['text'][:20],
                n_structure_changes=len(out['structure']), n_text_changes=len(out['text']))


def write_json(path, value):
    temp = path.with_suffix('.tmp')
    temp.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')
    temp.replace(path)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('target', choices=TARGETS)
    ap.add_argument('--list', action='store_true', help='show commands without running')
    ap.add_argument('--output', type=Path, help='new output directory; default validation-runs/UTC timestamp')
    args = ap.parse_args()
    root = Path(__file__).resolve().parent.parent
    commands = recipes(root, args.target)
    if args.list:
        for mod, argv in commands:
            print(mod + ': ' + shlex.join([sys.executable, *argv]))
        return 0
    stamp = datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')
    out = (args.output or root / 'validation-runs' / (stamp + '-' + args.target)).resolve()
    out.mkdir(parents=True, exist_ok=False)
    snapshot = out / 'snapshot'
    # Preserve precisely the tracked + nonignored files, including dirty work.
    names = subprocess.check_output(['git', '-C', str(root), 'ls-files', '-z', '--cached', '--others', '--exclude-standard']).decode().split('\0')
    hashes = {}
    for name in sorted(set(names) - {''}):
        src = root / name
        if not src.is_file() or src.is_symlink() or out == src or out in src.parents:
            continue
        dst = snapshot / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        hashes[name] = hashlib.sha256(src.read_bytes()).hexdigest()
    versions = {}
    for package in ('numpy', 'scipy', 'mpmath', 'sympy', 'pytest', 'matplotlib', 'setuptools'):
        try:
            versions[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            versions[package] = None
    meta = dict(started_utc=stamp, target=args.target, status='running',
                source_commit=git(root, 'rev-parse', 'HEAD'),
                source_status=git(root, 'status', '--short'), source_sha256=hashes,
                python=sys.version, executable=sys.executable, platform=platform.platform(),
                packages=versions, results=[])
    env = dict(os.environ, OMP_NUM_THREADS='2', OPENBLAS_NUM_THREADS='2', MKL_NUM_THREADS='2',
               MPLCONFIGDIR=str(out / 'mpl-cache'), XDG_CACHE_HOME=str(out / 'cache'),
               PYTHONUNBUFFERED='1', PYTHONHASHSEED='0')
    meta['execution_environment'] = {k: env[k] for k in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'MKL_NUM_THREADS', 'PYTHONHASHSEED')}
    write_json(out / 'manifest.json', meta)
    try:
        for index, (module, argv) in enumerate(commands, 1):
            command = [sys.executable, *argv]
            log = out / f'{index:02d}.log'
            print(f'[{index}/{len(commands)}] {module}: {shlex.join(command)}', flush=True)
            start = time.monotonic()
            with log.open('w') as fh:
                p = subprocess.run(command, cwd=snapshot / module, env=env, stdout=fh, stderr=subprocess.STDOUT)
            meta['results'].append(dict(module=module, command=command, exit_code=p.returncode,
                                        seconds=time.monotonic()-start, log=log.name))
            write_json(out / 'manifest.json', meta)
            print(f'  exit {p.returncode}', flush=True)
    except KeyboardInterrupt:
        meta['status'] = 'interrupted'
        write_json(out / 'manifest.json', meta)
        return 130
    meta['status'] = 'failed' if any(r['exit_code'] for r in meta['results']) else 'passed'
    generated = {}
    for name, original in hashes.items():
        path = snapshot / name
        if path.is_file():
            digest = hashlib.sha256(path.read_bytes()).hexdigest()
            if digest != original:
                generated[name] = digest
    for path in snapshot.rglob('*'):
        name = str(path.relative_to(snapshot))
        if path.is_file() and name not in hashes and '__pycache__' not in path.parts and '.pytest_cache' not in path.parts:
            generated[name] = hashlib.sha256(path.read_bytes()).hexdigest()
    meta['changed_or_new_snapshot_files_sha256'] = generated
    # Value-level comparison of every committed artifact the run rewrote.
    # The baseline is the working-tree file, valid only if it is still the
    # content the snapshot started from.
    comparison = {}
    for name in sorted(generated):
        if name not in hashes:
            continue
        src = root / name
        if not src.is_file() or hashlib.sha256(src.read_bytes()).hexdigest() != hashes[name]:
            comparison[name] = dict(kind='skipped', error='source changed during the run')
        else:
            comparison[name] = compare_artifact(src, snapshot / name)
    meta['artifact_comparison'] = comparison
    write_json(out / 'manifest.json', meta)
    if comparison:
        print('regenerated committed artifacts (runtime fields ignored):')
        for name, c in comparison.items():
            if c['kind'] == 'binary':
                note = 'bytes differ (not compared)'
            elif 'error' in c:
                note = c['error']
            else:
                note = f"max rel diff {c['max_rel']:.1e}"
                if c['n_structure_changes'] or c['n_text_changes']:
                    note += (f", {c['n_structure_changes']} structural and "
                             f"{c['n_text_changes']} text changes")
            print(f'  {name}: {note}')
    print(f'{meta["status"]}: {out / "manifest.json"}')
    return int(meta['status'] != 'passed')


if __name__ == '__main__':
    sys.exit(main())
