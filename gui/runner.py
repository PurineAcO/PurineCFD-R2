"""运行编排: run 目录编号、临时配置生成、命令行组装、日志增量读取."""

import datetime
import glob
import json
import os
import re
import shutil

from gui import farfield

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = ('cylinder', 'naca0012', 'oat15a')
PANELS = ('u', 'v', 'rho', 'miubl', 'T', 'Ma', 'Cp', 'schlieren', 'vorticity')
DEFAULT_BIN = os.path.join(ROOT, 'build', 'purinecfd')
POST_PREFIX = ['uv', 'run', '--with', 'numpy', '--with', 'matplotlib', '--with', 'scipy', 'python']

# 内置兜底模板: 算例目录下没有可用 json 时使用, 键名取自 src/io.hpp 白名单
TEMPLATE_BUILTIN = '（内置默认）'
DEFAULT_L = {'cylinder': 1.0, 'naca0012': 1.0, 'oat15a': 0.23}
DEFAULT_MA, DEFAULT_T, DEFAULT_RE = 0.2, 300.0, 1000.0
DEFAULT_SOLVER = {
  'end_step': 800,
  'cfl': 100.0,
  'dump_interval': 25,
  'convergence_interval': 10,
  'threads': 32,
  'urans': {
    'dt': 1e-3,
    'inner': 40,
    'inner_tol': 1e-3,
    'sweeps': 4,
    'steady_iters': 1500,
    'steady_cfl': 50.0,
    'wall_interval': 5,
    'seed': 0.02,
  },
}

BENCH_ITERS = 25
BENCH_DIR = '_bench'


def case_dir(case):
  return os.path.join(ROOT, 'cases', case)


def list_templates(case):
  d = case_dir(case)
  if not os.path.isdir(d):
    return []
  out = []
  for name in sorted(os.listdir(d)):
    if not name.endswith('.json') or name.startswith('_'):
      continue
    try:
      with open(os.path.join(d, name), encoding='utf-8') as f:
        obj = json.load(f)
    except (OSError, json.JSONDecodeError):
      continue
    if isinstance(obj, dict) and {'io', 'solver', 'farfield'} <= set(obj):
      out.append(name)
  return out


def builtin_template(case):
  """没有模板文件时的默认配置; p 由 Ma/Re/T/L 反推, 保证界面里 Re 能原样显示."""
  d = farfield.from_re(DEFAULT_MA, DEFAULT_RE, DEFAULT_T, DEFAULT_L.get(case, 1.0))
  return {
    'io': {'mesh': '', 'log': '', 'field': ''},
    'solver': json.loads(json.dumps(DEFAULT_SOLVER)),
    'farfield': {'Ma': DEFAULT_MA, 'T': DEFAULT_T, 'p': d['p'], 'alpha': 0.0},
  }


def load_template(case, name):
  if not name or name == TEMPLATE_BUILTIN:
    return builtin_template(case)
  with open(os.path.join(case_dir(case), name), encoding='utf-8') as f:
    return json.load(f)


def run_root():
  return os.path.join(ROOT, 'run')


def next_run_number():
  n = 0
  root = run_root()
  if os.path.isdir(root):
    for name in os.listdir(root):
      m = re.match(r'(\d+)-', name)
      if m:
        n = max(n, int(m.group(1)))
  return n + 1


def tag_of(case, alpha):
  if case == 'cylinder':
    return 'cyc'
  stem = 'naca' if case == 'naca0012' else 'oat'
  return f'{stem}-a{alpha:g}'.replace('.', '')


def make_run_dir(case, alpha):
  number = next_run_number()
  stamp = datetime.datetime.now().strftime('%Y%m%d-%H%M')
  path = os.path.join(run_root(), f'{number}-{tag_of(case, alpha)}-{stamp}')
  os.makedirs(path, exist_ok=True)
  return number, path


def write_config(case, template, number, run_dir, mesh, structured, solver, far, restart=''):
  cfg = json.loads(json.dumps(template))
  base = case_dir(case)
  cfg['io'] = dict(cfg.get('io', {}))
  cfg['io']['mesh'] = os.path.relpath(mesh, base)
  cfg['io']['field'] = os.path.relpath(run_dir, base)
  cfg['io']['log'] = os.path.relpath(os.path.join(run_dir, 'run.log'), base)
  if structured:
    cfg['io']['structured'] = os.path.relpath(structured, base)
  else:
    cfg['io'].pop('structured', None)
  if restart:
    cfg['io']['restart'] = os.path.relpath(restart, base)
  else:
    cfg['io'].pop('restart', None)
  cfg['solver'] = solver
  cfg['farfield'] = dict(Ma=far['Ma'], T=far['T'], p=far['p'], alpha=far['alpha'])
  path = os.path.join(base, f'_run_{number:02d}.json')
  with open(path, 'w', encoding='utf-8') as f:
    json.dump(cfg, f, indent=1, ensure_ascii=False)
    f.write('\n')
  return path


def build_post_cmd(cfg_path, case, run_dir, opt):
  cmd = list(POST_PREFIX) + [
    os.path.join(ROOT, 'tests', 'post.py'),
    cfg_path,
    '--case',
    case,
    '--out',
    run_dir,
    '--probe',
    f'{opt["probe_pos"]:g}',
  ]
  cmd += ['--panels', ','.join(opt['panels']) if opt['panels'] else 'none']
  if opt['readme']:
    cmd.append('--readme')
  if opt['wall'] is not None:
    cmd += ['--wall', f'{opt["wall"]:.1f}']
  if opt['u'] == 'pos':
    cmd.append('--u-positive')
  elif opt['u'] == 'sym':
    cmd.append('--u-symmetric')
  if opt['phase'] != 'auto':
    cmd += ['--phase', opt['phase']]
  for key in ('history', 'residual', 'cp_dist', 'cp_rms', 'probe', 'sequence'):
    if not opt[key]:
      cmd.append('--no-' + key.replace('_', '-'))
  return cmd


class LogTail:
  """按增量读取日志文件, 处理文件尚未创建 / 被截断的情况."""

  def __init__(self, path):
    self.path = path
    self.pos = 0

  def reset(self):
    self.pos = 0

  def poll(self):
    try:
      size = os.path.getsize(self.path)
    except OSError:
      return ''
    if size < self.pos:
      self.pos = 0
    if size == self.pos:
      return ''
    with open(self.path, 'rb') as f:
      f.seek(self.pos)
      data = f.read()
      self.pos = f.tell()
    return data.decode('utf-8', 'replace')


# ---------------------------------------------------------------- 线程测速
def bench_dir():
  return os.path.join(run_root(), BENCH_DIR)


def thread_candidates():
  """候选线程数, 不超过本机核心数."""
  cpu = os.cpu_count() or 32
  out = sorted(n for n in (4, 8, 16, 32, 64, 96, 128) if n <= cpu)
  return out or [cpu]


def restart_step(path):
  """从 step_*.dat 的 TITLE 行读步号, 读不出来返回 None."""
  try:
    with open(path, encoding='utf-8', errors='replace') as f:
      head = f.readline()
  except OSError:
    return None
  hit = re.search(r'step\s+(\d+)', head)
  return int(hit.group(1)) if hit else None


def latest_step(run_dir):
  """run 目录里步号最大的场文件, 没有就返回 None."""
  best, best_n = None, -1
  for path in glob.glob(os.path.join(run_dir, 'step_*.dat')):
    try:
      n = int(os.path.basename(path)[5:11])
    except ValueError:
      continue
    if n >= best_n:
      best, best_n = path, n
  return best


def write_bench_config(case, template, mesh, structured, solver, far, threads, out_dir, restart=''):
  """测速配置: 只跑少量定常预热步 + 1 个物理步, 输出全部落在 out_dir 内."""
  cfg = json.loads(json.dumps(template))
  base = case_dir(case)
  cfg['io'] = dict(cfg.get('io', {}))
  cfg['io']['mesh'] = os.path.relpath(mesh, base)
  cfg['io']['field'] = os.path.relpath(os.path.join(out_dir, f'field_{threads}'), base)
  cfg['io']['log'] = os.path.relpath(os.path.join(out_dir, f'{threads}.log'), base)
  if structured:
    cfg['io']['structured'] = os.path.relpath(structured, base)
  else:
    cfg['io'].pop('structured', None)
  if restart:
    cfg['io']['restart'] = os.path.relpath(restart, base)
  else:
    cfg['io'].pop('restart', None)
  s = json.loads(json.dumps(solver))
  s['max_steps'] = 1
  s.pop('end_step', None)  # max_steps 与 end_step 互斥
  s['dump_interval'] = 1000
  s['threads'] = threads
  u = s.setdefault('urans', json.loads(json.dumps(DEFAULT_SOLVER['urans'])))
  for key, value in DEFAULT_SOLVER['urans'].items():
    u.setdefault(key, value)
  u['steady_iters'] = BENCH_ITERS
  u['steady_cfl'] = 50.0
  u.pop('wall_interval', None)  # 缺省表示不输出壁面文件, 写成 0 会被 io.hpp 判为非法
  cfg['solver'] = s
  cfg['farfield'] = dict(Ma=far['Ma'], T=far['T'], p=far['p'], alpha=far['alpha'])
  path = os.path.join(base, f'_run_bench_{threads}.json')
  with open(path, 'w', encoding='utf-8') as f:
    json.dump(cfg, f, indent=1, ensure_ascii=False)
    f.write('\n')
  return path


def bench_log(threads):
  return os.path.join(bench_dir(), f'{threads}.log')


def last_steady_res(log_path):
  """日志里最后一条定常迭代的残差, 用来比较不同线程数下的收敛能力."""
  res = None
  try:
    with open(log_path, encoding='utf-8', errors='replace') as f:
      for line in f:
        parts = line.split()
        if len(parts) < 3 or parts[0] != 'steady' or 'res' not in parts:
          continue
        try:
          res = float(parts[parts.index('res') + 1])
        except (IndexError, ValueError):
          pass
  except OSError:
    pass
  return res


def clean_bench(case):
  """清掉测速产生的临时日志、临时场输出与临时配置."""
  shutil.rmtree(bench_dir(), ignore_errors=True)
  for path in glob.glob(os.path.join(case_dir(case), '_run_bench_*.json')):
    try:
      os.remove(path)
    except OSError:
      pass
