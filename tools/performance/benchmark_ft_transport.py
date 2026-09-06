#!/usr/bin/env python3
"""Matched, sequential FT timings with existing exact caches and fresh numerics."""
import argparse, decimal, hashlib, json, os, pathlib, signal, statistics, subprocess, time


def run(command, output, timeout):
    start = time.monotonic()
    with output.with_suffix('.json').open('w') as out, output.with_suffix('.log').open('w') as log:
        process = subprocess.Popen(list(map(str, command)), stdout=out, stderr=log, start_new_session=True)
        try:
            code = process.wait(timeout=timeout)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid, signal.SIGTERM)
            try:
                process.wait(timeout=3)
            except subprocess.TimeoutExpired:
                os.killpg(process.pid, signal.SIGKILL)
                process.wait()
            code = 124
    return {'exit_code': code, 'wall_seconds': time.monotonic() - start}


def compare_coefficients(left, right):
    """Conservative L1 difference of printed balls, normalized by boundary scale."""
    D = decimal.Decimal
    def ball(text):
        text = text.strip('[] ')
        if '+/-' in text:
            midpoint, radius = text.split('+/-')
            return D(midpoint.strip() or '0'), D(radius.strip())
        return D(text), D(0)
    def values(result):
        return {(row, order + result['epsilon_low']): tuple(ball(value[key]) for key in ('real', 'imaginary'))
                for row, coefficients in enumerate(result['coefficients']) for order, value in enumerate(coefficients)}
    a, b = values(left), values(right)
    if a.keys() != b.keys():
        raise ValueError('Matched coefficient windows differ')
    with decimal.localcontext() as context:
        context.prec = 50 + max(len(str(value)) for result in (a, b) for pair in result.values() for component in pair for value in component)
        worst = D(0)
        for key in a:
            error = sum(abs(x[0] - y[0]) + x[1] + y[1] for x, y in zip(a[key], b[key]))
            scale = 1 + max(D(0), *(abs(x[0]) - x[1] for x in a[key]))
            worst = max(worst, error / scale)
        return {'coefficients': len(a), 'normalized_difference_bound': str(worst), 'pass': worst < D('1e-25')}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--executable', type=pathlib.Path, required=True)
    parser.add_argument('--fire', type=pathlib.Path, required=True)
    parser.add_argument('--checker', type=pathlib.Path, required=True)
    parser.add_argument('--case', action='append', nargs=4, metavar=('LABEL', 'REFERENCE', 'CONFIG', 'CACHE'), required=True)
    parser.add_argument('--output', type=pathlib.Path, required=True)
    parser.add_argument('--repeats', type=int, default=3)
    parser.add_argument('--timeout', type=int, default=600)
    args = parser.parse_args()
    if not 1 <= args.timeout <= 600 or not 1 <= args.repeats <= 10:
        parser.error('timeout must be 1..600 seconds; repeats must be 1..10')
    if args.output.exists() and any(args.output.iterdir()):
        parser.error('Use an empty output directory')
    if len({case[0] for case in args.case}) != len(args.case):
        parser.error('Case labels must be unique')
    for label, _, config, cache in args.case:
        if not label or pathlib.Path(label).name != label or label in ('.', '..'):
            parser.error('Case labels must be plain directory names')
        if not pathlib.Path(config).is_file() or not pathlib.Path(cache).is_dir():
            parser.error('Each case requires a configuration and an existing exact cache')
    args.output.mkdir(parents=True, exist_ok=True)
    summary = {'schema': 'DiffExp.FTTransportBenchmark/v1', 'repeats': args.repeats,
               'cache_policy': 'existing exact reductions; no numerical checkpoint reads or writes', 'cases': []}
    passed = True
    for label, reference, config, cache in args.case:
        folder = args.output / label
        folder.mkdir()
        source = pathlib.Path(config).read_bytes()
        copied = folder / 'family.json'
        copied.write_bytes(source)
        record = {'label': label, 'reference': reference, 'configuration_sha256': hashlib.sha256(source).hexdigest(), 'runs': []}
        for repeat in range(args.repeats):
            # Alternate order to avoid always favoring the second run.
            for method in (['taylor', 'auto'] if repeat % 2 == 0 else ['auto', 'taylor']):
                output = folder / f'{method}-{repeat}'
                command = [args.executable, 'ft', copied, '--cache', cache, '--fire', args.fire,
                           '--fire-threads', '1', '--fire-simplifier-threads', '1',
                           '--no-numerical-cache', '--ft-transport', method, '--json']
                item = dict(method=method, repeat=repeat, **run(command, output, args.timeout))
                if item['exit_code'] == 0:
                    result = json.loads(output.with_suffix('.json').read_text())
                    item['result'] = result
                    if result['systems_built'] != 0:
                        raise RuntimeError('Exact preparation cache miss: this is not a matched warm-preparation run')
                    check_path = folder / f'{method}-{repeat}-reference'
                    item['reference_run'] = run([args.checker, reference, copied, output.with_suffix('.json')], check_path, 120)
                    if item['reference_run']['exit_code'] == 0:
                        item['validation'] = json.loads(check_path.with_suffix('.json').read_text())
                passed &= item.get('validation', {}).get('status') == 'pass'
                record['runs'].append(item)
                print(label, method, repeat, item['wall_seconds'], item.get('validation', {}).get('status', 'FAILED'), flush=True)
        record['medians'] = {}
        for method in ['taylor', 'auto']:
            rows = [item for item in record['runs'] if item['method'] == method and item.get('validation', {}).get('status') == 'pass']
            if len(rows) == args.repeats:
                record['medians'][method] = {'wall_seconds': statistics.median(item['wall_seconds'] for item in rows),
                    **{key: statistics.median(item['result']['timings'][key] for item in rows)
                       for key in ['total_seconds', 'preparation_seconds', 'numerical_seconds', 'ordinary_seconds']}}
        record['matched_values'] = []
        for repeat in range(args.repeats):
            pair = {item['method']: item for item in record['runs'] if item['repeat'] == repeat}
            if all('result' in pair[method] for method in ('taylor', 'auto')):
                comparison = compare_coefficients(pair['taylor']['result'], pair['auto']['result'])
                passed &= comparison['pass']
                record['matched_values'].append(comparison)
        summary['cases'].append(record)
        summary['pass'] = passed
        (args.output / 'summary.json').write_text(json.dumps(summary, indent=2) + '\n')
    if not passed:
        raise SystemExit(1)


if __name__ == '__main__':
    main()
