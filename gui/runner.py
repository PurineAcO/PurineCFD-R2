"""运行编排: run 目录编号、临时配置生成、命令行组装、日志增量读取."""

import datetime
import json
import os
import re

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CASES = ('cylinder', 'naca0012', 'oat15a')
PANELS = ('u', 'v', 'rho', 'miubl', 'T', 'Ma', 'Cp', 'schlieren', 'vorticity')
DEFAULT_BIN = os.path.join(ROOT, 'build', 'purinecfd')
POST_PREFIX = ['uv', 'run', '--with', 'numpy', '--with', 'matplotlib',
               '--with', 'scipy', 'python']


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


def load_template(case, name):
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


def write_config(case, template, number, run_dir, mesh, structured, solver, far):
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
    cfg['solver'] = solver
    cfg['farfield'] = dict(Ma=far['Ma'], T=far['T'], p=far['p'], alpha=far['alpha'])
    path = os.path.join(base, f'_run_{number:02d}.json')
    with open(path, 'w', encoding='utf-8') as f:
        json.dump(cfg, f, indent=1, ensure_ascii=False)
        f.write('\n')
    return path


def build_post_cmd(cfg_path, case, run_dir, opt):
    cmd = list(POST_PREFIX) + [os.path.join(ROOT, 'tests', 'post.py'), cfg_path,
                               '--case', case, '--out', run_dir,
                               '--probe', f"{opt['probe_pos']:g}"]
    cmd += ['--panels', ','.join(opt['panels']) if opt['panels'] else 'none']
    if opt['readme']:
        cmd.append('--readme')
    if opt['wall'] is not None:
        cmd += ['--wall', f"{opt['wall']:.1f}"]
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
