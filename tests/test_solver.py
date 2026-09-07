import json
import math
import os
import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope='session')
def executable():
  path = Path(os.environ.get('PURINECFD_BIN', ROOT / 'build/purinecfd')).resolve()
  assert path.is_file(), f'Build the solver first or set PURINECFD_BIN: {path}'
  return path


@pytest.fixture(scope='session')
def mesh(tmp_path_factory):
  directory = tmp_path_factory.mktemp('cylinder-mesh')
  subprocess.run(
    [sys.executable, str(ROOT / 'cases/cylinder/generate_case.py'), str(directory)],
    check=True,
    capture_output=True,
  )
  return directory / 'mesh.txt'


def settings(mesh):
  config = json.loads((ROOT / 'config.json').read_text())
  config['io'] = {'mesh': str(mesh), 'log': 'run.log', 'field': 'field'}
  config['solver']['max_steps'] = 3
  return config


def run_case(executable, directory, config, threads=1, affinity=None):
  directory.mkdir(parents=True, exist_ok=True)
  path = directory / 'config.json'
  path.write_text(config if isinstance(config, str) else json.dumps(config))
  environment = dict(os.environ)
  if threads is None:
    environment.pop('OMP_NUM_THREADS', None)
  else:
    environment['OMP_NUM_THREADS'] = str(threads)
  return subprocess.run(
    [str(executable), str(path)],
    cwd=directory.parent,
    env=environment,
    preexec_fn=(lambda: os.sched_setaffinity(0, affinity)) if affinity else None,
    capture_output=True,
    text=True,
    timeout=30,
  )


def read_field(path):
  return [[float(value) for value in row.split()] for row in path.read_text().splitlines()[2:]]


@pytest.mark.parametrize('model', ['sa', 'laminar'])
def test_model_regression_and_thread_consistency(executable, mesh, tmp_path, model):
  config = settings(mesh)
  config['solver']['max_steps'] = 200
  config['solver']['model'] = model
  for threads in [1, 4]:
    directory = tmp_path / str(threads)
    result = run_case(executable, directory, config, threads)
    assert result.returncode == 0, result.stderr
  first = tmp_path / '1/field/step_000200.dat'
  assert first.read_bytes() == (tmp_path / '4/field/step_000200.dat').read_bytes()
  stats = []
  for threads in [1, 4]:
    lines = (tmp_path / str(threads) / 'run.log').read_text().splitlines()
    stats.append([line for line in lines if line.lstrip()[:1].isdigit()])
  assert len(stats[0]) == 199 and stats[0] == stats[1]
  data = read_field(first)
  reference = 'reference_200.json' if model == 'sa' else 'reference_laminar_200.json'
  expected = json.loads((ROOT / 'tests' / reference).read_text())
  if model == 'laminar':
    assert all(row[8] == 0 for row in data)
  else:
    assert any(row[8] > 0 for row in data)
  assert len(data) == expected['cells']
  for row in data:
    assert len(row) == 9 and all(math.isfinite(value) for value in row)
    assert all(row[i] > 0 for i in [2, 5, 6]) and row[8] >= 0
    assert row[7] == pytest.approx(
      math.hypot(row[3], row[4]) / math.sqrt(1.4 * 287.05 * row[5]), rel=1e-7
    )
  for index, row in expected['probes'].items():
    assert data[int(index)] == pytest.approx(row, rel=2e-7, abs=1e-10)
  for label, function in [('minimum', min), ('maximum', max)]:
    actual = [function(row[i] for row in data) for i in range(9)]
    assert actual == pytest.approx(expected[label], rel=2e-7, abs=1e-10)
  log = (tmp_path / '1/run.log').read_text()
  label = 'Steady SA-RANS' if model == 'sa' else 'Steady laminar Navier-Stokes'
  assert label in log and 'convergence criterion not satisfied' in log
  assert f'model={model}' in first.read_text().splitlines()[0]
  assert not (tmp_path / '1/field/checkpoint.dat').exists()


@pytest.mark.parametrize(
  'section,key,value',
  [
    ('solver', 'model', 'sst'),
    ('solver', 'model', True),
    ('solver', 'global_dt', False),
    ('solver', 'rk_coeff', [1.0]),
    ('solver', 'cfl', 0),
    ('solver', 'cfl', True),
    ('solver', 'cfl', '1.0'),
    ('solver', 'max_steps', 1.5),
    ('solver', 'max_steps', 2**40),
    ('solver', 'dump_interval', -1),
    ('solver', 'convergence_interval', 1),
    ('farfield', 'Ma', 1.2),
    ('farfield', 'T', 0),
    ('farfield', 'p', -1),
    ('farfield', 'AOA_deg', 0),
    ('io', 'mesh', ''),
    ('io', 'mesh', 'bad\x00path'),
  ],
)
def test_invalid_settings_fail(executable, mesh, tmp_path, section, key, value):
  config = settings(mesh)
  config[section][key] = value
  result = run_case(executable, tmp_path / 'case', config)
  assert result.returncode != 0 and 'Error:' in result.stderr


@pytest.mark.parametrize(
  'text',
  [
    '{',
    '{"solver": [}',
    '{"io":{},"io":{}}',
    '{"io":{},}',
    '{} trailing',
    '{"solver":{"cfl":1e999}}',
  ],
)
def test_malformed_json_fails_without_hanging(executable, tmp_path, text):
  result = run_case(executable, tmp_path / 'case', text)
  assert result.returncode != 0 and 'Error:' in result.stderr


def test_missing_and_legacy_keys_rejected(executable, mesh, tmp_path):
  for name, config in [('missing', settings(mesh)), ('legacy', settings(mesh))]:
    if name == 'missing':
      del config['farfield']['T']
    else:
      config['viscous'] = {'ifviscous': True, 'model': 'SA'}
    result = run_case(executable, tmp_path / name, config)
    assert result.returncode != 0


@pytest.mark.parametrize(
  'damage', ['truncated', 'node_index', 'boundary_type', 'face_index', 'cell_shape', 'zero_length']
)
def test_invalid_mesh_fails(executable, mesh, tmp_path, damage):
  lines = mesh.read_text().splitlines()
  node_start = lines.index('(node)') + 1
  edge_start = lines.index('(edge)') + 2
  cell_start = lines.index('(cell)') + 1
  if damage == 'truncated':
    lines = lines[: node_start + 2]
  elif damage == 'node_index':
    lines[node_start] = '0 0.5 0'
  elif damage == 'boundary_type':
    lines[2] = 'cylinder = VIL'
  elif damage == 'face_index':
    parts = lines[edge_start].split()
    parts[1] = '99999999'
    lines[edge_start] = ' '.join(parts)
  elif damage == 'cell_shape':
    lines[cell_start] += ' 7'
  elif damage == 'zero_length':
    a = lines[node_start].split()
    b = lines[node_start + 1].split()
    lines[node_start + 1] = ' '.join([b[0], *a[1:]])
  path = tmp_path / 'broken.txt'
  path.write_text('\n'.join(lines) + '\n')
  result = run_case(executable, tmp_path / 'case', settings(path))
  assert result.returncode != 0 and 'Error:' in result.stderr


def test_missing_input_and_output_failure(executable, mesh, tmp_path):
  result = run_case(executable, tmp_path / 'missing', settings(tmp_path / 'absent.mesh'))
  assert result.returncode != 0 and 'Cannot open mesh' in result.stderr
  bad_output = tmp_path / 'not-a-directory'
  bad_output.write_text('keep')
  config = settings(mesh)
  config['io']['field'] = str(bad_output)
  result = run_case(executable, tmp_path / 'bad-output', config)
  assert result.returncode != 0 and bad_output.read_text() == 'keep'


def test_log_cannot_overwrite_mesh(executable, mesh, tmp_path):
  config = settings(mesh)
  config['io']['log'] = str(mesh)
  before = mesh.read_bytes()
  result = run_case(executable, tmp_path / 'case', config)
  assert result.returncode != 0 and 'must not overwrite' in result.stderr
  assert mesh.read_bytes() == before


def test_auto_threads_respect_affinity(executable, mesh, tmp_path):
  cpu = min(os.sched_getaffinity(0))
  result = run_case(executable, tmp_path / 'one-core', settings(mesh), threads=None, affinity={cpu})
  assert result.returncode == 0, result.stderr
  log = (tmp_path / 'one-core/run.log').read_text()
  assert 'OpenMP threads=1 | thread_policy=available-physical-cores' in log


def test_auto_threads_do_not_count_smt_twice(executable, mesh, tmp_path):
  allowed = os.sched_getaffinity(0)
  pair = None
  for cpu in sorted(allowed):
    text = Path(f'/sys/devices/system/cpu/cpu{cpu}/topology/thread_siblings_list').read_text()
    siblings = set()
    for part in text.strip().split(','):
      lo, _, hi = part.partition('-')
      siblings.update(range(int(lo), int(hi or lo) + 1))
    candidates = sorted(siblings & allowed)
    if len(candidates) >= 2:
      pair = set(candidates[:2])
      break
  if pair is None:
    pytest.skip('No accessible SMT sibling pair on this machine')
  result = run_case(executable, tmp_path / 'smt-pair', settings(mesh), threads=None, affinity=pair)
  assert result.returncode == 0, result.stderr
  assert 'OpenMP threads=1 ' in (tmp_path / 'smt-pair/run.log').read_text()


def test_explicit_thread_count_overrides_auto(executable, mesh, tmp_path):
  result = run_case(executable, tmp_path / 'explicit', settings(mesh), threads=2)
  assert result.returncode == 0, result.stderr
  assert (
    'OpenMP threads=2 | thread_policy=OMP_NUM_THREADS'
    in (tmp_path / 'explicit/run.log').read_text()
  )


def test_invalid_flow_reports_same_first_cell_across_threads(executable, mesh, tmp_path):
  config = settings(mesh)
  config['solver']['cfl'] = 1e6
  errors = []
  for threads in [1, 4]:
    result = run_case(executable, tmp_path / str(threads), config, threads)
    assert result.returncode != 0
    assert 'Invalid flow state at step' in result.stderr
    errors.append(result.stderr)
  assert errors[0] == errors[1]
