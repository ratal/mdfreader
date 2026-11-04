"""
benchmark.py — mdfreader performance benchmark
Reproduces the methodology from asammdf's benchmark page and produces an
HTML report with Chart.js bar charts.

Usage:
    python benchmark.py [--runs N] [--output FILE]

Options:
    --runs N        Number of timed runs per measurement (default: 3)
    --output FILE   HTML output file (default: benchmark_results.html)
"""

import argparse
import gc
import os
import platform
import sys
import tempfile
import threading
import time
import json
from datetime import datetime

# ── File paths ─────────────────────────────────────────────────────────────
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)
MDF3_FILE = os.path.join(_REPO_ROOT, 'mdfreader', 'tests', 'mdf3', 'test.mdf')
MDF4_FILE = os.path.join(_REPO_ROOT, 'mdfreader', 'tests', 'MDF4', 'MDF4.3', 'Base_Standard', 'Examples', 'Simple', 'test.mf4')

# ── Reference numbers: mdfreader 4.2, Windows 11, AMD64 (from asammdf docs) ─
REFERENCE = {
    'open_default':     {'mdf3': (1785, 523),  'mdf4': (3513, 544)},
    'open_nodata':      {'mdf3': (644,  252),  'mdf4': (2348, 288)},
    'open_compress':    {'mdf3': (1369, 389),  'mdf4': (3326, 413)},
    'save':             {'mdf3': (6241, 548),  'mdf4': (10259, 562)},
    'get_all_channels': {'mdf3': (35,   523),  'mdf4': (45,   543)},
    'merge_3':          {'mdf3': (20639, 1140), 'mdf4': (27437, 1170)},
}

OPERATION_LABELS = {
    'open_default':     'Open (default)',
    'open_nodata':      'Open (no_data_loading=True)',
    'open_compress':    'Open (compression=blosc)',
    'save':             'Save',
    'get_all_channels': 'Get all channels',
    'merge_3':          'Merge 3 files',
}

# ── Memory polling via /proc/self/status ────────────────────────────────────

class _PeakMemoryPoller:
    """Background thread that samples VmRSS every 50ms and tracks peak."""

    def __init__(self):
        self._peak = 0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _sample(self):
        try:
            with open('/proc/self/status') as f:
                for line in f:
                    if line.startswith('VmRSS:'):
                        return int(line.split()[1])  # kB
        except Exception:
            pass
        return 0

    def _run(self):
        while not self._stop.is_set():
            val = self._sample()
            if val > self._peak:
                self._peak = val
            self._stop.wait(0.05)

    def start(self):
        self._peak = self._sample()
        self._thread.start()

    def stop_and_get(self):
        self._stop.set()
        self._thread.join(timeout=1)
        return self._peak / 1024  # MB


def _peak_ram_fallback_mb():
    """Fallback when /proc is unavailable (non-Linux): returns 0."""
    return 0


# ── Measurement harness ─────────────────────────────────────────────────────

def measure(fn, n_runs=3):
    """Run fn() n_runs times. Return (median_ms, peak_ram_mb)."""
    times = []
    peak_mb = 0

    for i in range(n_runs):
        gc.collect()
        poller = _PeakMemoryPoller()
        try:
            poller.start()
            t0 = time.perf_counter()
            fn()
            elapsed_ms = (time.perf_counter() - t0) * 1000
            peak_mb = poller.stop_and_get()
        except Exception:
            poller.stop_and_get()
            raise
        times.append(elapsed_ms)

    times.sort()
    return times[n_runs // 2], peak_mb


# ── Benchmark functions ─────────────────────────────────────────────────────

def _bench_open_default(path):
    import mdfreader
    mdfreader.Mdf(path)


def _bench_open_nodata(path):
    import mdfreader
    mdfreader.Mdf(path, no_data_loading=True)


def _bench_open_compress(path):
    import mdfreader
    mdfreader.Mdf(path, compression='blosc')


def _bench_save(path):
    import mdfreader
    mdf = mdfreader.Mdf(path)
    with tempfile.NamedTemporaryFile(suffix='.mdf', delete=False) as f:
        tmp = f.name
    try:
        mdf.write(tmp)
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass


def _bench_get_all_channels(path):
    import mdfreader
    mdf = mdfreader.Mdf(path)
    for ch in mdf.keys():
        mdf.get_channel_data(ch)


def _bench_merge_3(path):
    import mdfreader
    mdf = mdfreader.Mdf(path)
    mdf.merge_mdf(mdfreader.Mdf(path))
    mdf.merge_mdf(mdfreader.Mdf(path))


# ── Optional: asammdf comparison ───────────────────────────────────────────

def _asammdf_bench_open_default(path):
    import asammdf
    asammdf.MDF(path)


def _asammdf_bench_open_nodata(path):
    import asammdf
    asammdf.MDF(path, use_display_names=False)


def _asammdf_bench_save(path):
    import asammdf
    mdf = asammdf.MDF(path)
    with tempfile.NamedTemporaryFile(suffix='.mf4', delete=False) as f:
        tmp = f.name
    try:
        mdf.save(tmp, overwrite=True)
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass


def _asammdf_bench_get_all_channels(path):
    import asammdf
    mdf = asammdf.MDF(path)
    for ch in mdf.channels_db:
        try:
            mdf.get(ch)
        except Exception:
            pass


def _asammdf_bench_merge_3(path):
    import asammdf
    asammdf.MDF.concatenate([path, path, path])


# ── Run all benchmarks ──────────────────────────────────────────────────────

def run_benchmarks(n_runs=3, with_compress=True, with_asammdf=True):
    import mdfreader

    has_blosc = False
    if with_compress:
        try:
            import blosc  # noqa: F401
            has_blosc = True
            print('  blosc available — compression mode enabled')
        except ImportError:
            print('  blosc not installed — skipping compression mode')

    has_asammdf = False
    if with_asammdf:
        try:
            import asammdf  # noqa: F401
            has_asammdf = True
            print(f'  asammdf {asammdf.__version__} available — side-by-side comparison enabled')
        except ImportError:
            print('  asammdf not installed — skipping side-by-side')

    files = {'mdf3': MDF3_FILE, 'mdf4': MDF4_FILE}
    results = {}  # results[op_key][lib_key][fmt] = (ms, mb) or None

    ops = [
        ('open_default',     _bench_open_default,     None),
        ('open_nodata',      _bench_open_nodata,      None),
        ('open_compress',    _bench_open_compress,    'blosc'),
        ('save',             _bench_save,             None),
        ('get_all_channels', _bench_get_all_channels, None),
        ('merge_3',          _bench_merge_3,          None),
    ]

    asammdf_ops = {
        'open_default':     _asammdf_bench_open_default,
        'open_nodata':      _asammdf_bench_open_nodata,
        'save':             _asammdf_bench_save,
        'get_all_channels': _asammdf_bench_get_all_channels,
        'merge_3':          _asammdf_bench_merge_3,
    }

    total_ops = sum(
        len(files) for (op_key, fn, req) in ops
        if req is None or (req == 'blosc' and has_blosc)
    )
    done = 0

    for op_key, fn, req in ops:
        results[op_key] = {}

        if req == 'blosc' and not has_blosc:
            results[op_key]['mdfreader'] = {'mdf3': None, 'mdf4': None}
            if has_asammdf and op_key in asammdf_ops:
                results[op_key]['asammdf'] = {'mdf3': None, 'mdf4': None}
            continue

        results[op_key]['mdfreader'] = {}
        for fmt, path in files.items():
            done += 1
            label = OPERATION_LABELS[op_key]
            print(f'  [{done}/{total_ops}] mdfreader {label} ({fmt.upper()}) ...', end=' ', flush=True)
            try:
                ms, mb = measure(lambda p=path: fn(p), n_runs)
                results[op_key]['mdfreader'][fmt] = (round(ms), round(mb))
                print(f'{ms:.0f} ms  {mb:.0f} MB')
            except Exception as e:
                print(f'ERROR: {e}')
                results[op_key]['mdfreader'][fmt] = None

        if has_asammdf and op_key in asammdf_ops:
            results[op_key]['asammdf'] = {}
            afn = asammdf_ops[op_key]
            for fmt, path in files.items():
                print(f'       asammdf  {OPERATION_LABELS[op_key]} ({fmt.upper()}) ...', end=' ', flush=True)
                try:
                    ms, mb = measure(lambda p=path: afn(p), n_runs)
                    results[op_key]['asammdf'][fmt] = (round(ms), round(mb))
                    print(f'{ms:.0f} ms  {mb:.0f} MB')
                except Exception as e:
                    print(f'ERROR: {e}')
                    results[op_key]['asammdf'][fmt] = None
        elif op_key in asammdf_ops:
            results[op_key]['asammdf'] = {'mdf3': None, 'mdf4': None}

    return results, has_asammdf


# ── HTML generation ─────────────────────────────────────────────────────────

_HTML_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>mdfreader benchmark results</title>
<script src="https://cdn.jsdelivr.net/npm/chart.js@4/dist/chart.umd.min.js"></script>
<style>
  body {{ font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, sans-serif;
          max-width: 1100px; margin: 40px auto; padding: 0 20px; color: #333; }}
  h1 {{ border-bottom: 2px solid #4a90e2; padding-bottom: 8px; }}
  h2 {{ color: #4a90e2; margin-top: 40px; }}
  .meta {{ background: #f5f7fa; border-left: 4px solid #4a90e2; padding: 10px 16px;
           border-radius: 4px; margin: 20px 0; font-size: 0.9em; line-height: 1.7; }}
  table {{ border-collapse: collapse; width: 100%; margin: 16px 0 32px; font-size: 0.88em; }}
  th {{ background: #4a90e2; color: #fff; padding: 8px 12px; text-align: left; }}
  td {{ padding: 7px 12px; border-bottom: 1px solid #ddd; }}
  tr:nth-child(even) {{ background: #f9f9f9; }}
  tr.ref {{ color: #888; font-style: italic; }}
  tr.ref td {{ background: #fafafa; }}
  .better {{ color: #27ae60; font-weight: bold; }}
  .worse  {{ color: #e74c3c; }}
  .na     {{ color: #aaa; }}
  .chart-row {{ display: flex; gap: 24px; margin-bottom: 40px; flex-wrap: wrap; }}
  .chart-box {{ flex: 1 1 420px; background: #fff; border: 1px solid #ddd;
                border-radius: 6px; padding: 16px; }}
  canvas {{ max-width: 100%; }}
  .note {{ font-size: 0.8em; color: #777; margin-top: 4px; }}
</style>
</head>
<body>
<h1>mdfreader benchmark results</h1>

<div class="meta">
  <strong>Date:</strong> {date}<br>
  <strong>Platform:</strong> {platform}<br>
  <strong>Python:</strong> {python}<br>
  <strong>mdfreader:</strong> {mdfreader_version}<br>
  {asammdf_line}
  <strong>Runs per measurement:</strong> {n_runs} (median reported)<br>
  <strong>Memory:</strong> Peak VmRSS sampled every 50 ms (Linux /proc/self/status)<br>
  <strong>Test files:</strong> test.mdf (MDF3, 176 MB, 183 DGs, 36 424 channels) &amp;
    test.mf4 (MDF4, 192 MB, 183 DGs, 36 424 channels)
</div>

<p class="note">
  Reference row (<em>4.2 ref</em>) shows mdfreader 4.2 numbers from the
  <a href="https://asammdf.readthedocs.io/en/latest/benchmarks.html" target="_blank">asammdf benchmark page</a>
  (Windows 11, AMD64). Speed-up column compares current mdfreader vs. the 4.2 reference.
</p>

{sections}

<script>
{charts_js}
</script>
</body>
</html>
"""

_SECTION_TEMPLATE = """\
<h2>{title}</h2>
<table>
  <thead>
    <tr>
      <th>Library / mode</th>
      <th>MDF3 time [ms]</th><th>MDF3 RAM [MB]</th>
      <th>MDF4 time [ms]</th><th>MDF4 RAM [MB]</th>
      <th>Speed-up vs 4.2 (MDF4)</th>
    </tr>
  </thead>
  <tbody>
{rows}
  </tbody>
</table>
<div class="chart-row">
  <div class="chart-box"><canvas id="{chart_id}_time"></canvas></div>
  <div class="chart-box"><canvas id="{chart_id}_ram"></canvas></div>
</div>
"""


def _fmt(val, ref=None, is_time=True):
    if val is None:
        return '<span class="na">N/A</span>'
    return str(val)


def _speedup(cur_ms, ref_ms):
    if cur_ms is None or ref_ms is None:
        return '<span class="na">—</span>'
    ratio = ref_ms / cur_ms
    cls = 'better' if ratio >= 1.0 else 'worse'
    return f'<span class="{cls}">{ratio:.2f}×</span>'


def _row(label, mdf3, mdf4, ref_op, css_class=''):
    ref = REFERENCE.get(ref_op, {})
    ref_mdf4 = ref.get('mdf4')
    m3t = _fmt(mdf3[0] if mdf3 else None)
    m3r = _fmt(mdf3[1] if mdf3 else None)
    m4t = _fmt(mdf4[0] if mdf4 else None)
    m4r = _fmt(mdf4[1] if mdf4 else None)
    su = _speedup(mdf4[0] if mdf4 else None, ref_mdf4[0] if ref_mdf4 else None)
    tr_class = f' class="{css_class}"' if css_class else ''
    return (f'    <tr{tr_class}>'
            f'<td>{label}</td>'
            f'<td>{m3t}</td><td>{m3r}</td>'
            f'<td>{m4t}</td><td>{m4r}</td>'
            f'<td>{su}</td>'
            f'</tr>')


def build_html(results, has_asammdf, n_runs):
    import mdfreader

    try:
        import asammdf as _asammdf
        asammdf_ver = _asammdf.__version__
    except ImportError:
        asammdf_ver = None

    sections_html = []
    charts_js_parts = []

    # Group operations the same way asammdf docs do
    groups = [
        ('Open', ['open_default', 'open_nodata', 'open_compress']),
        ('Save', ['save']),
        ('Get all channels', ['get_all_channels']),
        ('Merge 3 files', ['merge_3']),
    ]

    for group_title, op_keys in groups:
        rows = []
        chart_data_time = {'labels': [], 'mdf3': [], 'mdf4': [], 'ref_mdf3': [], 'ref_mdf4': []}
        chart_data_ram  = {'labels': [], 'mdf3': [], 'mdf4': [], 'ref_mdf3': [], 'ref_mdf4': []}
        chart_id = group_title.lower().replace(' ', '_').replace('/', '_')

        for op_key in op_keys:
            label = OPERATION_LABELS[op_key]
            res = results.get(op_key, {})
            mr = res.get('mdfreader', {})
            mdf3_cur = mr.get('mdf3')
            mdf4_cur = mr.get('mdf4')

            # mdfreader current
            rows.append(_row(f'mdfreader 4.3 — {label}', mdf3_cur, mdf4_cur, op_key))

            # asammdf if available
            if has_asammdf and 'asammdf' in res:
                ar = res.get('asammdf', {})
                adf3 = ar.get('mdf3')
                adf4 = ar.get('mdf4')
                rows.append(_row(f'asammdf — {label}', adf3, adf4, None))

            # reference row
            ref = REFERENCE.get(op_key, {})
            ref3 = ref.get('mdf3')
            ref4 = ref.get('mdf4')
            ref3_tuple = (ref3[0], ref3[1]) if ref3 else None
            ref4_tuple = (ref4[0], ref4[1]) if ref4 else None
            rows.append(_row(f'mdfreader 4.2 ref — {label}', ref3_tuple, ref4_tuple, None, css_class='ref'))

            # chart data
            chart_data_time['labels'].append(label)
            chart_data_time['mdf3'].append(mdf3_cur[0] if mdf3_cur else None)
            chart_data_time['mdf4'].append(mdf4_cur[0] if mdf4_cur else None)
            chart_data_time['ref_mdf3'].append(ref3[0] if ref3 else None)
            chart_data_time['ref_mdf4'].append(ref4[0] if ref4 else None)

            chart_data_ram['labels'].append(label)
            chart_data_ram['mdf3'].append(mdf3_cur[1] if mdf3_cur else None)
            chart_data_ram['mdf4'].append(mdf4_cur[1] if mdf4_cur else None)
            chart_data_ram['ref_mdf3'].append(ref3[1] if ref3 else None)
            chart_data_ram['ref_mdf4'].append(ref4[1] if ref4 else None)

        rows_html = '\n'.join(rows)
        sections_html.append(_SECTION_TEMPLATE.format(
            title=group_title,
            rows=rows_html,
            chart_id=chart_id,
        ))

        # Build Chart.js JS for time and RAM
        for metric, cdata, unit, suffix in [
            ('Time', chart_data_time, 'ms', 'time'),
            ('RAM', chart_data_ram, 'MB', 'ram'),
        ]:
            canvas_id = f'{chart_id}_{suffix}'
            datasets = [
                {'label': 'mdfreader 4.3 MDF3', 'data': cdata['mdf3'],
                 'backgroundColor': 'rgba(74,144,226,0.75)'},
                {'label': 'mdfreader 4.3 MDF4', 'data': cdata['mdf4'],
                 'backgroundColor': 'rgba(39,174,96,0.75)'},
                {'label': 'mdfreader 4.2 ref MDF3', 'data': cdata['ref_mdf3'],
                 'backgroundColor': 'rgba(74,144,226,0.25)',
                 'borderColor': 'rgba(74,144,226,0.9)', 'borderWidth': 1},
                {'label': 'mdfreader 4.2 ref MDF4', 'data': cdata['ref_mdf4'],
                 'backgroundColor': 'rgba(39,174,96,0.25)',
                 'borderColor': 'rgba(39,174,96,0.9)', 'borderWidth': 1},
            ]
            charts_js_parts.append(f"""\
new Chart(document.getElementById('{canvas_id}'), {{
  type: 'bar',
  data: {{
    labels: {json.dumps(cdata['labels'])},
    datasets: {json.dumps(datasets)}
  }},
  options: {{
    plugins: {{ title: {{ display: true, text: '{group_title} — {metric} [{unit}]' }} }},
    scales: {{ y: {{ beginAtZero: true, title: {{ display: true, text: '{unit}' }} }} }},
    responsive: true
  }}
}});""")

    asammdf_line = (f'<strong>asammdf:</strong> {asammdf_ver}<br>' if asammdf_ver else
                    '<strong>asammdf:</strong> not installed<br>')

    html = _HTML_TEMPLATE.format(
        date=datetime.now().strftime('%Y-%m-%d %H:%M'),
        platform=platform.platform(),
        python=sys.version.split()[0],
        mdfreader_version=mdfreader.__version__,
        asammdf_line=asammdf_line,
        n_runs=n_runs,
        sections='\n'.join(sections_html),
        charts_js='\n'.join(charts_js_parts),
    )
    return html


# ── Entry point ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument('--runs', type=int, default=3,
                        help='Number of timed runs per measurement (default: 3)')
    parser.add_argument('--output', default='benchmark_results.html',
                        help='Output HTML file (default: benchmark_results.html)')
    parser.add_argument('--no-compress', action='store_true',
                        help='Skip blosc compression mode')
    parser.add_argument('--no-asammdf', action='store_true',
                        help='Skip asammdf side-by-side comparison')
    args = parser.parse_args()

    # Validate test files exist
    for label, path in [('MDF3', MDF3_FILE), ('MDF4', MDF4_FILE)]:
        if not os.path.exists(path):
            sys.exit(f'ERROR: {label} test file not found: {path}\n'
                     f'Run from the repo root or adjust the path constants in benchmark.py.')

    import mdfreader
    print(f'mdfreader {mdfreader.__version__} benchmark')
    print(f'  MDF3: {MDF3_FILE} ({os.path.getsize(MDF3_FILE) // 1024 // 1024} MB)')
    print(f'  MDF4: {MDF4_FILE} ({os.path.getsize(MDF4_FILE) // 1024 // 1024} MB)')
    print(f'  Runs: {args.runs} (median)\n')

    t_start = time.perf_counter()
    results, has_asammdf = run_benchmarks(
        n_runs=args.runs,
        with_compress=not args.no_compress,
        with_asammdf=not args.no_asammdf,
    )
    elapsed = time.perf_counter() - t_start
    print(f'\nTotal benchmark time: {elapsed:.1f}s')

    html = build_html(results, has_asammdf, n_runs=args.runs)
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f'Report written to: {args.output}')


if __name__ == '__main__':
    main()
