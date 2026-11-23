"""
benchmark.py — mdfreader performance benchmark
Reproduces the methodology from asammdf's benchmark page, compares
mdfreader side-by-side against mdfr (and optionally asammdf) on core
operations as well as sorted and unsorted MDF files, and produces an
HTML report with Chart.js bar charts.

Usage:
    python benchmark.py [--runs N] [--output FILE] [--no-compress]
                        [--no-mdfr] [--no-asammdf] [--include-large]

Options:
    --runs N          Number of timed runs per measurement (default: 3)
    --output FILE     HTML output file (default: benchmark_results.html)
    --no-compress     Skip blosc compression mode
    --no-mdfr         Skip mdfr side-by-side comparison
    --no-asammdf      Skip asammdf side-by-side comparison
    --include-large   Include >1 GB sorted files in the layout benchmark
"""

import argparse
from datetime import datetime
import gc
import json
import os
import platform
import subprocess
import sys
import tempfile
import threading
import time

# ── File paths ─────────────────────────────────────────────────────────────
_REPO_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if _REPO_ROOT not in sys.path:
    sys.path.insert(0, _REPO_ROOT)

MDF3_FILE = os.path.join(_REPO_ROOT, 'mdfreader', 'tests', 'mdf3', 'test.mdf')
MDF4_FILE = os.path.join(
    _REPO_ROOT, 'mdfreader', 'tests', 'MDF4', 'MDF4.3',
    'Base_Standard', 'Examples', 'Simple', 'test.mf4',
)

# (category, key, short_name, layout_desc, rel_path, is_large)
LAYOUT_FILES = [
    # Sorted files
    ('sorted', 'test_mdf3', 'test.mdf',
     'MDF3 sorted row-wise (183 DGs, 36 424 ch)',
     'mdfreader/tests/mdf3/test.mdf', False),
    ('sorted', 'test_mf4', 'test.mf4',
     'MDF4 sorted row-wise (183 DGs, 36 424 ch)',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/Simple/test.mf4', False),
    ('sorted', 'fichier_mf4', 'FichierMF4.mf4',
     'MDF4 sorted row-wise (multi-DG)',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/Simple/FichierMF4.mf4', False),
    ('sorted', 'measure_mf4', 'Measure.mf4',
     'MDF4 sorted row-wise',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/Simple/Measure.mf4', False),
    ('sorted', 'error_col', 'error_col.mf4',
     'MDF4 sorted column-oriented ##LD (36 424 ch)',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/Simple/error_col.mf4', False),
    ('sorted', 'zeitschrieb', 'Zeitschrieb_9022.mf4',
     'MDF4 sorted compressed ##HL/##DZ',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/Simple/Zeitschrieb_9022.mf4', False),
    ('sorted', 'nedc_row', 'T3_121121_000_6NEDC.mf4',
     'MDF4 large sorted row-wise',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/Simple/T3_121121_000_6NEDC.mf4', True),
    ('sorted', 'nedc_col', 'T3_121121_000_6NEDC_col.mf4',
     'MDF4 large sorted column-oriented',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/Simple/T3_121121_000_6NEDC_col.mf4', True),
    # Unsorted files
    ('unsorted', 'kra_4010', 'KRA_4010_..._log_23.mf4',
     'MDF4 unsorted multi-CG',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/Simple/KRA_4010_2022_08_25_065833_log_23.mf4', False),
    ('unsorted', 'vec_vlsd', 'Vector_Unsorted_VLSD.MF4',
     'MDF4 unsorted VLSD',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/UnsortedData/VLSD/Vector_Unsorted_VLSD.MF4', False),
    ('unsorted', 'vec_dz', 'Vector_SingleDZ_Unsorted.MF4',
     'MDF4 unsorted compressed ##DZ',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/CompressedData/Unsorted/Vector_SingleDZ_Unsorted.MF4', False),
    ('unsorted', 'rac_vlsd', 'RAC_MDF430_Unsorted_VLSD_Compact.mf4',
     'MDF4 unsorted compact VLSD',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/UnsortedData/VLSD/RAC_MDF430_Unsorted_VLSD_Compact_Structure.mf4', False),
    ('unsorted', 'rac_vlsc', 'RAC_MDF430_Unsorted_VLSC_Compact.mf4',
     'MDF4 unsorted compact VLSC',
     'mdfreader/tests/MDF4/MDF4.3/Base_Standard/Examples/UnsortedData/VLSC/RAC_MDF430_Unsorted_VLSC_Compact_Structure.mf4', False),
]

# ── Reference numbers: mdfreader 4.2, Windows 11, AMD64 (from asammdf docs) ─
REFERENCE = {
    'open_default':     {'mdf3': (1785, 523),   'mdf4': (3513, 544)},
    'open_nodata':      {'mdf3': (644,  252),   'mdf4': (2348, 288)},
    'open_compress':    {'mdf3': (1369, 389),   'mdf4': (3326, 413)},
    'save':             {'mdf3': (6241, 548),   'mdf4': (10259, 562)},
    'get_all_channels': {'mdf3': (35,   523),   'mdf4': (45,   543)},
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


# ── Memory polling / /proc/self/status helpers ─────────────────────────────

def _read_proc_status_kb(field):
    """Return a memory counter (in kB) from /proc/self/status, or 0."""
    prefix = f'{field}:'
    try:
        with open('/proc/self/status', encoding='utf-8') as f:
            for line in f:
                if line.startswith(prefix):
                    return int(line.split()[1])
    except Exception:
        pass
    return 0


class _PeakMemoryPoller:
    """Background thread that samples VmRSS every 50ms and tracks peak."""

    def __init__(self):
        self._peak = 0
        self._stop = threading.Event()
        self._thread = threading.Thread(target=self._run, daemon=True)

    def _sample(self):
        return _read_proc_status_kb('VmRSS')

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


# ── Measurement harness ─────────────────────────────────────────────────────

def measure(fn, n_runs=3):
    """Run fn() n_runs times in-process. Return (median_ms, peak_ram_mb)."""
    times = []
    peak_mb = 0.0

    for _ in range(n_runs):
        gc.collect()
        poller = _PeakMemoryPoller()
        try:
            poller.start()
            t0 = time.perf_counter()
            fn()
            elapsed_ms = (time.perf_counter() - t0) * 1000
            peak_mb = max(peak_mb, poller.stop_and_get())
        except Exception:
            poller.stop_and_get()
            raise
        times.append(elapsed_ms)

    times.sort()
    return times[n_runs // 2], peak_mb


# ── Benchmark operations per library ────────────────────────────────────────

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
    suffix = os.path.splitext(path)[1] or '.mdf'
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as f:
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


def _mdfr_bench_open_default(path):
    import mdfr
    r = mdfr.Mdfr(path)
    r.load_all_channels_data_in_memory()


def _mdfr_bench_open_nodata(path):
    import mdfr
    mdfr.Mdfr(path)


def _mdfr_bench_save(path):
    import mdfr
    r = mdfr.Mdfr(path)
    r.load_all_channels_data_in_memory()
    suffix = os.path.splitext(path)[1] or '.mf4'
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as f:
        tmp = f.name
    try:
        r.write(tmp)
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass


def _mdfr_bench_get_all_channels(path):
    import mdfr
    r = mdfr.Mdfr(path)
    r.load_all_channels_data_in_memory()
    for ch in r.get_channel_names_set():
        r.get_channel_data(ch)


def _mdfr_bench_merge_3(path):
    import mdfr
    r1 = mdfr.Mdfr(path)
    r1.load_all_channels_data_in_memory()
    r2 = mdfr.Mdfr(path)
    r2.load_all_channels_data_in_memory()
    r3 = mdfr.Mdfr(path)
    r3.load_all_channels_data_in_memory()
    r1.concat_mdf(r2)
    r1.concat_mdf(r3)


def _asammdf_bench_open_default(path):
    import asammdf
    asammdf.MDF(path)


def _asammdf_bench_open_nodata(path):
    import asammdf
    asammdf.MDF(path, use_display_names=False)


def _asammdf_bench_save(path):
    import asammdf
    mdf = asammdf.MDF(path)
    suffix = os.path.splitext(path)[1] or '.mf4'
    with tempfile.NamedTemporaryFile(suffix=suffix, delete=False) as f:
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


_LIB_OPS = {
    'mdfreader': {
        'open_default':     _bench_open_default,
        'open_nodata':      _bench_open_nodata,
        'open_compress':    _bench_open_compress,
        'save':             _bench_save,
        'get_all_channels': _bench_get_all_channels,
        'merge_3':          _bench_merge_3,
    },
    'mdfr': {
        'open_default':     _mdfr_bench_open_default,
        'open_nodata':      _mdfr_bench_open_nodata,
        'save':             _mdfr_bench_save,
        'get_all_channels': _mdfr_bench_get_all_channels,
        'merge_3':          _mdfr_bench_merge_3,
    },
    'asammdf': {
        'open_default':     _asammdf_bench_open_default,
        'open_nodata':      _asammdf_bench_open_nodata,
        'save':             _asammdf_bench_save,
        'get_all_channels': _asammdf_bench_get_all_channels,
        'merge_3':          _asammdf_bench_merge_3,
    },
}


# ── Isolated subprocess measurement ─────────────────────────────────────────

def _run_worker(lib_name, op_key, path, n_runs):
    """Import target library first, run n_runs iterations, print JSON stats."""
    __import__(lib_name)
    fn = _LIB_OPS[lib_name][op_key]
    poller = _PeakMemoryPoller() if _read_proc_status_kb('VmHWM') == 0 else None
    if poller is not None:
        poller.start()

    times = []
    try:
        for _ in range(n_runs):
            gc.collect()
            t0 = time.perf_counter()
            fn(path)
            times.append((time.perf_counter() - t0) * 1000)
    finally:
        poller_mb = poller.stop_and_get() if poller is not None else 0.0

    hwm_kb = _read_proc_status_kb('VmHWM')
    peak_mb = hwm_kb / 1024 if hwm_kb > 0 else poller_mb
    times.sort()
    payload = {'ms': times[n_runs // 2], 'mb': peak_mb}
    print(json.dumps(payload))


def measure_isolated(lib_name, op_key, path, n_runs=3):
    """Measure (median_ms, peak_ram_mb) in a clean subprocess.

    Importing the target library happens inside the subprocess before timing
    starts, while Linux /proc/self/status VmHWM captures exact process peak RSS
    without cross-library allocator bleed.
    """
    cmd = [
        sys.executable, os.path.abspath(__file__),
        '--_worker', lib_name, op_key, path, '--runs', str(n_runs),
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        err = (proc.stderr or proc.stdout).strip().splitlines()
        raise RuntimeError(err[-1] if err else f'worker exited with {proc.returncode}')
    lines = [line.strip() for line in proc.stdout.splitlines() if line.strip()]
    data = json.loads(lines[-1])
    return data['ms'], data['mb']


def _round_ms(ms):
    return round(ms, 1) if ms < 20.0 else round(ms)


def _format_ms(val):
    if val is None:
        return 'N/A'
    if isinstance(val, float) and val < 20.0:
        return f'{val:.1f}'
    return str(round(val))


# ── Library discovery ───────────────────────────────────────────────────────

def _detect_libraries(with_compress=True, with_mdfr=True, with_asammdf=True):
    """Probe optional libraries and print availability."""
    has_blosc = False
    if with_compress:
        try:
            import blosc  # noqa: F401
            has_blosc = True
            print('  blosc available — compression mode enabled')
        except ImportError:
            print('  blosc not installed — skipping compression mode')

    mdfr_ver = None
    if with_mdfr:
        try:
            import mdfr  # noqa: F401
            from importlib import metadata
            try:
                mdfr_ver = metadata.version('mdfr')
            except Exception:
                mdfr_ver = getattr(mdfr, '__version__', 'installed')
            print(f'  mdfr {mdfr_ver} available — side-by-side comparison enabled')
        except ImportError:
            print('  mdfr not installed — skipping mdfr comparison')

    asammdf_ver = None
    if with_asammdf:
        try:
            import asammdf
            asammdf_ver = asammdf.__version__
            print(f'  asammdf {asammdf_ver} available — side-by-side comparison enabled')
        except ImportError:
            print('  asammdf not installed — skipping asammdf comparison')

    return has_blosc, mdfr_ver, asammdf_ver


# ── Run core benchmarks ─────────────────────────────────────────────────────

def run_benchmarks(n_runs=3, has_blosc=False, has_mdfr=False, has_asammdf=False):
    """Run core operations suite (Open, Save, Get all channels, Merge 3)."""
    files = {'mdf3': MDF3_FILE, 'mdf4': MDF4_FILE}
    results = {}  # results[op_key][lib_key][fmt] = (ms, mb) or None

    ops = [
        ('open_default',     None),
        ('open_nodata',      None),
        ('open_compress',    'blosc'),
        ('save',             None),
        ('get_all_channels', None),
        ('merge_3',          None),
    ]

    libs = ['mdfreader']
    if has_mdfr:
        libs.append('mdfr')
    if has_asammdf:
        libs.append('asammdf')

    total_steps = sum(
        len(files) for (op_key, req) in ops
        if req is None or (req == 'blosc' and has_blosc)
    )
    done = 0

    for op_key, req in ops:
        results[op_key] = {}
        label = OPERATION_LABELS[op_key]

        if req == 'blosc' and not has_blosc:
            for lib in libs:
                if op_key in _LIB_OPS[lib]:
                    results[op_key][lib] = {'mdf3': None, 'mdf4': None}
            continue

        for lib in libs:
            if op_key not in _LIB_OPS[lib]:
                continue
            results[op_key][lib] = {}
            for fmt, path in files.items():
                if lib == 'mdfreader':
                    done += 1
                    prefix = f'  [{done}/{total_steps}]'
                else:
                    prefix = '        '
                print(f'{prefix} {lib:<9s} {label} ({fmt.upper()}) ...', end=' ', flush=True)
                try:
                    ms, mb = measure_isolated(lib, op_key, path, n_runs)
                    results[op_key][lib][fmt] = (_round_ms(ms), round(mb))
                    print(f'{_format_ms(ms)} ms  {mb:.0f} MB')
                except Exception as e:
                    print(f'ERROR: {e}')
                    results[op_key][lib][fmt] = None

    return results


# ── Run sorted vs unsorted layout benchmarks ────────────────────────────────

def run_layout_benchmarks(n_runs=3, has_mdfr=False, has_asammdf=False, include_large=False):
    """Compare libraries on curated sorted and unsorted MDF files (open_default)."""
    libs = ['mdfreader']
    if has_mdfr:
        libs.append('mdfr')
    if has_asammdf:
        libs.append('asammdf')

    layout_results = []
    candidates = [
        entry for entry in LAYOUT_FILES
        if (include_large or not entry[5])
        and os.path.exists(os.path.join(_REPO_ROOT, entry[4]))
    ]

    total = len(candidates)
    print(f'\nSorted & Unsorted open_default matrix ({total} files):')
    for idx, (category, key, short_name, desc, rel_path, _is_large) in enumerate(candidates, 1):
        full_path = os.path.join(_REPO_ROOT, rel_path)
        size_mb = os.path.getsize(full_path) / (1024 * 1024)
        entry_res = {
            'category': category,
            'key': key,
            'name': short_name,
            'desc': desc,
            'rel_path': rel_path,
            'size_mb': size_mb,
            'libs': {},
        }
        print(f'  [{idx}/{total}] [{category}] {short_name} ({size_mb:.2f} MB)')
        for lib in libs:
            print(f'      {lib:<9s} ...', end=' ', flush=True)
            try:
                ms, mb = measure_isolated(lib, 'open_default', full_path, n_runs)
                entry_res['libs'][lib] = (_round_ms(ms), round(mb))
                print(f'{_format_ms(ms)} ms  {mb:.0f} MB')
            except Exception as e:
                print(f'ERROR: {e}')
                entry_res['libs'][lib] = None
        layout_results.append(entry_res)

    return layout_results


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
          max-width: 1180px; margin: 40px auto; padding: 0 20px; color: #333; }}
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
  {mdfr_line}
  {asammdf_line}
  <strong>Runs per measurement:</strong> {n_runs} (median reported, library imported prior to timer)<br>
  <strong>Memory:</strong> Peak process RSS (VmHWM) in an isolated subprocess per (library, test)<br>
  <strong>Core test files:</strong> test.mdf (MDF3, 168 MB, 183 DGs, 36 424 channels) &amp;
    test.mf4 (MDF4, 184 MB, 183 DGs, 36 424 channels)
</div>

<p class="note">
  Reference row (<em>4.2 ref</em>) shows mdfreader 4.2 numbers from the
  <a href="https://asammdf.readthedocs.io/en/latest/benchmarks.html" target="_blank">asammdf benchmark page</a>
  (Windows 11, AMD64). Speed-up columns show <code>ref_ms / cur_ms</code> vs 4.2 (for core ops)
  and <code>mdfr_ms / mdfreader_ms</code> (for sorted/unsorted files: &gt;1&times; means mdfreader is faster;
  &lt;1&times; means mdfr is faster).
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


def _fmt(val):
    if val is None:
        return '<span class="na">N/A</span>'
    return _format_ms(val)


def _speedup(cur_ms, ref_ms):
    if cur_ms is None or ref_ms is None or cur_ms <= 0:
        return '<span class="na">—</span>'
    ratio = ref_ms / cur_ms
    cls = 'better' if ratio >= 1.0 else 'worse'
    return f'<span class="{cls}">{ratio:.2f}&times;</span>'


def _row(label, mdf3, mdf4, ref_op, css_class=''):
    ref = REFERENCE.get(ref_op, {}) if ref_op else {}
    ref_mdf4 = ref.get('mdf4')
    m3t = _fmt(mdf3[0] if mdf3 else None)
    m3r = _fmt(mdf3[1] if mdf3 else None)
    m4t = _fmt(mdf4[0] if mdf4 else None)
    m4r = _fmt(mdf4[1] if mdf4 else None)
    su = _speedup(mdf4[0] if mdf4 else None, ref_mdf4[0] if ref_mdf4 else None)
    tr_class = f' class="{css_class}"' if css_class else ''
    return (
        f'    <tr{tr_class}>'
        f'<td>{label}</td>'
        f'<td>{m3t}</td><td>{m3r}</td>'
        f'<td>{m4t}</td><td>{m4r}</td>'
        f'<td>{su}</td>'
        f'</tr>'
    )


def _build_chart_js(canvas_id, title, unit, labels, datasets):
    return f"""\
new Chart(document.getElementById('{canvas_id}'), {{
  type: 'bar',
  data: {{
    labels: {json.dumps(labels)},
    datasets: {json.dumps(datasets)}
  }},
  options: {{
    plugins: {{ title: {{ display: true, text: {json.dumps(f'{title} — {unit}')} }} }},
    scales: {{ y: {{ beginAtZero: true, title: {{ display: true, text: {json.dumps(unit)} }} }} }},
    responsive: true
  }}
}});"""


def _build_layout_section(category, title, entries, mdfreader_ver, mdfr_ver, asammdf_ver):
    if not entries:
        return '', []

    rows = []
    labels = []
    time_series = {'mdfreader': [], 'mdfr': [], 'asammdf': []}
    ram_series = {'mdfreader': [], 'mdfr': [], 'asammdf': []}

    for item in entries:
        libs = item['libs']
        mr = libs.get('mdfreader')
        mf = libs.get('mdfr')
        am = libs.get('asammdf')
        size_str = f"{item['size_mb']:.2f}" if item['size_mb'] < 10 else f"{item['size_mb']:.1f}"

        cells = [
            f"<td><strong>{item['name']}</strong></td>",
            f"<td>{size_str}</td>",
            f"<td>{item['desc']}</td>",
            f"<td>{_fmt(mr[0] if mr else None)}</td>",
            f"<td>{_fmt(mr[1] if mr else None)}</td>",
        ]
        if mdfr_ver is not None:
            cells.extend([
                f"<td>{_fmt(mf[0] if mf else None)}</td>",
                f"<td>{_fmt(mf[1] if mf else None)}</td>",
            ])
        if asammdf_ver is not None:
            cells.extend([
                f"<td>{_fmt(am[0] if am else None)}</td>",
                f"<td>{_fmt(am[1] if am else None)}</td>",
            ])
        if mdfr_ver is not None:
            ratio_cell = _speedup(mr[0] if mr else None, mf[0] if mf else None)
            cells.append(f"<td>{ratio_cell}</td>")

        rows.append(f"    <tr>{''.join(cells)}</tr>")

        labels.append(item['name'])
        for lib_key in ('mdfreader', 'mdfr', 'asammdf'):
            val = libs.get(lib_key)
            time_series[lib_key].append(val[0] if val else None)
            ram_series[lib_key].append(val[1] if val else None)

    headers = [
        '<th>File</th>',
        '<th>Size [MB]</th>',
        '<th>Layout / Type</th>',
        f'<th>mdfreader {mdfreader_ver} time [ms]</th>',
        f'<th>mdfreader {mdfreader_ver} RAM [MB]</th>',
    ]
    if mdfr_ver is not None:
        headers.extend([
            f'<th>mdfr {mdfr_ver} time [ms]</th>',
            f'<th>mdfr {mdfr_ver} RAM [MB]</th>',
        ])
    if asammdf_ver is not None:
        headers.extend([
            f'<th>asammdf {asammdf_ver} time [ms]</th>',
            f'<th>asammdf {asammdf_ver} RAM [MB]</th>',
        ])
    if mdfr_ver is not None:
        headers.append('<th>Speed-up vs mdfr</th>')

    chart_id = f'layout_{category}'
    section_html = f"""\
<h2>{title}</h2>
<table>
  <thead>
    <tr>{''.join(headers)}</tr>
  </thead>
  <tbody>
{chr(10).join(rows)}
  </tbody>
</table>
<div class="chart-row">
  <div class="chart-box"><canvas id="{chart_id}_time"></canvas></div>
  <div class="chart-box"><canvas id="{chart_id}_ram"></canvas></div>
</div>"""

    charts = []
    for suffix, unit, series_map in [
        ('time', 'Time [ms]', time_series),
        ('ram', 'Peak RAM [MB]', ram_series),
    ]:
        datasets = [
            {'label': f'mdfreader {mdfreader_ver}', 'data': series_map['mdfreader'],
             'backgroundColor': 'rgba(74,144,226,0.8)'},
        ]
        if mdfr_ver is not None:
            datasets.append({
                'label': f'mdfr {mdfr_ver}', 'data': series_map['mdfr'],
                'backgroundColor': 'rgba(230,126,34,0.8)',
            })
        if asammdf_ver is not None:
            datasets.append({
                'label': f'asammdf {asammdf_ver}', 'data': series_map['asammdf'],
                'backgroundColor': 'rgba(155,89,182,0.8)',
            })
        charts.append(_build_chart_js(f'{chart_id}_{suffix}', title, unit, labels, datasets))

    return section_html, charts


def build_html(results, layout_results, mdfr_ver, asammdf_ver, n_runs):
    import mdfreader
    mdfreader_ver = mdfreader.__version__

    sections_html = []
    charts_js_parts = []

    # 1. Sorted & Unsorted file open_default matrices
    for category, title in [
        ('sorted', 'Sorted files — Open (default)'),
        ('unsorted', 'Unsorted files — Open (default)'),
    ]:
        cat_entries = [e for e in layout_results if e['category'] == category]
        sec_html, sec_charts = _build_layout_section(
            category, title, cat_entries, mdfreader_ver, mdfr_ver, asammdf_ver
        )
        if sec_html:
            sections_html.append(sec_html)
            charts_js_parts.extend(sec_charts)

    # 2. Core operations grouped the same way as asammdf docs
    groups = [
        ('Core: Open', ['open_default', 'open_nodata', 'open_compress']),
        ('Core: Save', ['save']),
        ('Core: Get all channels', ['get_all_channels']),
        ('Core: Merge 3 files', ['merge_3']),
    ]

    for group_title, op_keys in groups:
        rows = []
        chart_data_time = {
            'labels': [], 'mdf3': [], 'mdf4': [],
            'mdfr_mdf3': [], 'mdfr_mdf4': [],
            'ref_mdf3': [], 'ref_mdf4': [],
        }
        chart_data_ram = {
            'labels': [], 'mdf3': [], 'mdf4': [],
            'mdfr_mdf3': [], 'mdfr_mdf4': [],
            'ref_mdf3': [], 'ref_mdf4': [],
        }
        chart_id = group_title.lower().replace(' ', '_').replace(':', '').replace('/', '_')

        for op_key in op_keys:
            label = OPERATION_LABELS[op_key]
            res = results.get(op_key, {})
            mr = res.get('mdfreader', {})
            mdf3_cur = mr.get('mdf3')
            mdf4_cur = mr.get('mdf4')

            rows.append(_row(f'mdfreader {mdfreader_ver} — {label}', mdf3_cur, mdf4_cur, op_key))

            mf = res.get('mdfr', {})
            if mdfr_ver is not None and 'mdfr' in res:
                rows.append(_row(f'mdfr {mdfr_ver} — {label}', mf.get('mdf3'), mf.get('mdf4'), op_key))

            if asammdf_ver is not None and 'asammdf' in res:
                ar = res.get('asammdf', {})
                rows.append(_row(f'asammdf {asammdf_ver} — {label}', ar.get('mdf3'), ar.get('mdf4'), op_key))

            ref = REFERENCE.get(op_key, {})
            ref3 = ref.get('mdf3')
            ref4 = ref.get('mdf4')
            ref3_tuple = (ref3[0], ref3[1]) if ref3 else None
            ref4_tuple = (ref4[0], ref4[1]) if ref4 else None
            rows.append(_row(f'mdfreader 4.2 ref — {label}', ref3_tuple, ref4_tuple, None, css_class='ref'))

            chart_data_time['labels'].append(label)
            chart_data_time['mdf3'].append(mdf3_cur[0] if mdf3_cur else None)
            chart_data_time['mdf4'].append(mdf4_cur[0] if mdf4_cur else None)
            chart_data_time['mdfr_mdf3'].append(mf['mdf3'][0] if mf.get('mdf3') else None)
            chart_data_time['mdfr_mdf4'].append(mf['mdf4'][0] if mf.get('mdf4') else None)
            chart_data_time['ref_mdf3'].append(ref3[0] if ref3 else None)
            chart_data_time['ref_mdf4'].append(ref4[0] if ref4 else None)

            chart_data_ram['labels'].append(label)
            chart_data_ram['mdf3'].append(mdf3_cur[1] if mdf3_cur else None)
            chart_data_ram['mdf4'].append(mdf4_cur[1] if mdf4_cur else None)
            chart_data_ram['mdfr_mdf3'].append(mf['mdf3'][1] if mf.get('mdf3') else None)
            chart_data_ram['mdfr_mdf4'].append(mf['mdf4'][1] if mf.get('mdf4') else None)
            chart_data_ram['ref_mdf3'].append(ref3[1] if ref3 else None)
            chart_data_ram['ref_mdf4'].append(ref4[1] if ref4 else None)

        sections_html.append(_SECTION_TEMPLATE.format(
            title=group_title,
            rows='\n'.join(rows),
            chart_id=chart_id,
        ))

        for metric, cdata, unit, suffix in [
            ('Time', chart_data_time, 'ms', 'time'),
            ('RAM', chart_data_ram, 'MB', 'ram'),
        ]:
            datasets = [
                {'label': f'mdfreader {mdfreader_ver} MDF3', 'data': cdata['mdf3'],
                 'backgroundColor': 'rgba(74,144,226,0.8)'},
                {'label': f'mdfreader {mdfreader_ver} MDF4', 'data': cdata['mdf4'],
                 'backgroundColor': 'rgba(39,174,96,0.8)'},
            ]
            if mdfr_ver is not None:
                datasets.extend([
                    {'label': f'mdfr {mdfr_ver} MDF3', 'data': cdata['mdfr_mdf3'],
                     'backgroundColor': 'rgba(230,126,34,0.8)'},
                    {'label': f'mdfr {mdfr_ver} MDF4', 'data': cdata['mdfr_mdf4'],
                     'backgroundColor': 'rgba(211,84,0,0.8)'},
                ])
            datasets.extend([
                {'label': 'mdfreader 4.2 ref MDF3', 'data': cdata['ref_mdf3'],
                 'backgroundColor': 'rgba(74,144,226,0.2)',
                 'borderColor': 'rgba(74,144,226,0.9)', 'borderWidth': 1},
                {'label': 'mdfreader 4.2 ref MDF4', 'data': cdata['ref_mdf4'],
                 'backgroundColor': 'rgba(39,174,96,0.2)',
                 'borderColor': 'rgba(39,174,96,0.9)', 'borderWidth': 1},
            ])
            charts_js_parts.append(
                _build_chart_js(f'{chart_id}_{suffix}', f'{group_title} — {metric}', unit, cdata['labels'], datasets)
            )

    mdfr_line = (f'<strong>mdfr:</strong> {mdfr_ver}<br>' if mdfr_ver else
                 '<strong>mdfr:</strong> not installed / skipped<br>')
    asammdf_line = (f'<strong>asammdf:</strong> {asammdf_ver}<br>' if asammdf_ver else
                    '<strong>asammdf:</strong> not installed / skipped<br>')

    return _HTML_TEMPLATE.format(
        date=datetime.now().strftime('%Y-%m-%d %H:%M'),
        platform=platform.platform(),
        python=sys.version.split()[0],
        mdfreader_version=mdfreader_ver,
        mdfr_line=mdfr_line,
        asammdf_line=asammdf_line,
        n_runs=n_runs,
        sections='\n'.join(sections_html),
        charts_js='\n'.join(charts_js_parts),
    )


# ── Entry point ─────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument('--runs', type=int, default=3,
                        help='Number of timed runs per measurement (default: 3)')
    parser.add_argument('--output', default='benchmark_results.html',
                        help='Output HTML file (default: benchmark_results.html)')
    parser.add_argument('--no-compress', action='store_true',
                        help='Skip blosc compression mode')
    parser.add_argument('--no-mdfr', action='store_true',
                        help='Skip mdfr side-by-side comparison')
    parser.add_argument('--no-asammdf', action='store_true',
                        help='Skip asammdf side-by-side comparison')
    parser.add_argument('--include-large', action='store_true',
                        help='Include >1 GB sorted files in the layout benchmark')
    parser.add_argument('--_worker', nargs=3, metavar=('LIB', 'OP', 'PATH'),
                        help=argparse.SUPPRESS)
    args = parser.parse_args()

    if args._worker:
        lib_name, op_key, path = args._worker
        _run_worker(lib_name, op_key, path, args.runs)
        return

    for label, path in [('MDF3', MDF3_FILE), ('MDF4', MDF4_FILE)]:
        if not os.path.exists(path):
            sys.exit(
                f'ERROR: {label} test file not found: {path}\n'
                f'Run from the repo root or adjust the path constants in benchmark.py.'
            )

    import mdfreader
    print(f'mdfreader {mdfreader.__version__} benchmark')
    print(f'  MDF3: {MDF3_FILE} ({os.path.getsize(MDF3_FILE) // 1024 // 1024} MB)')
    print(f'  MDF4: {MDF4_FILE} ({os.path.getsize(MDF4_FILE) // 1024 // 1024} MB)')
    print(f'  Runs: {args.runs} (median)\n')

    has_blosc, mdfr_ver, asammdf_ver = _detect_libraries(
        with_compress=not args.no_compress,
        with_mdfr=not args.no_mdfr,
        with_asammdf=not args.no_asammdf,
    )

    t_start = time.perf_counter()
    results = run_benchmarks(
        n_runs=args.runs,
        has_blosc=has_blosc,
        has_mdfr=mdfr_ver is not None,
        has_asammdf=asammdf_ver is not None,
    )
    layout_results = run_layout_benchmarks(
        n_runs=args.runs,
        has_mdfr=mdfr_ver is not None,
        has_asammdf=asammdf_ver is not None,
        include_large=args.include_large,
    )
    elapsed = time.perf_counter() - t_start
    print(f'\nTotal benchmark time: {elapsed:.1f}s')

    html = build_html(
        results, layout_results,
        mdfr_ver=mdfr_ver,
        asammdf_ver=asammdf_ver,
        n_runs=args.runs,
    )
    with open(args.output, 'w', encoding='utf-8') as f:
        f.write(html)
    print(f'Report written to: {args.output}')


if __name__ == '__main__':
    main()
