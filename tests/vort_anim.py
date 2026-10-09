#!/usr/bin/env python3
"""把 run 目录里的 step_*.dat 序列做成动画 (默认涡量, GIF).

用法:
  uv run --with matplotlib --with numpy --with pillow python tests/vort_anim.py <config> --out <run 目录>
可选项: --field vorticity|schlieren|Ma|Cp   --fps 8   --level 4   --stride 1   --size 8
结构化网格(圆柱/NACA)走 reshape; 非结构(OAT)走三角网。
"""
import argparse
import glob
import json
import math
import os
from collections import defaultdict

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.animation import PillowWriter

GAMMA, R = 1.4, 287.05

# 与 post.py 的 CASES 保持一致, 保证云图与动画用同一视图
VIEWS = {
    'cylinder': (-5.0, 15.0, -10.0, 10.0),
    'naca0012': (-0.7, 0.8, -0.35, 0.55),
    'oat15a': (-0.5, 1.5, -0.6, 0.6),
}

ap = argparse.ArgumentParser()
ap.add_argument('config')
ap.add_argument('--out', default=None, help='run 目录(默认取 config 的 io.field)')
ap.add_argument('--field', default='vorticity',
                choices=['vorticity', 'schlieren', 'Ma', 'Cp', 'u', 'v'])
ap.add_argument('--fps', type=float, default=8.0)
ap.add_argument('--level', type=float, default=0.0,
                help='涡量色标 ±LEVEL (L/U 归一); <=0 表示按首帧 99%% 分位自动定')
ap.add_argument('--stride', type=int, default=1, help='每隔几帧取一帧')
ap.add_argument('--size', type=float, default=8.0, help='图宽英寸')
ap.add_argument('--case', default=None, choices=list(VIEWS),
                help='用该算例的默认视图(与 post.py 的云图一致)')
ap.add_argument('--view', default=None,
                help='x0,x1,y0,y1 覆盖视图; 逗号列表以负号开头时必须写成 --view=a,b,c,d')
a = ap.parse_args()

base = os.path.dirname(os.path.abspath(a.config))
cfg = json.load(open(a.config))
out = a.out or os.path.normpath(os.path.join(base, cfg['io']['field']))
far = cfg['farfield']
U = far['Ma'] * math.sqrt(GAMMA * R * far['T'])
q_inf = 0.5 * GAMMA * far['p'] * far['Ma'] ** 2
LREF = 1.0  # 圆柱 D / 翼型弦长归一

# ---- 网格 ----
if cfg['io'].get('structured'):
    S = int(open(os.path.join(base, cfg['io']['structured'])).readline().split()[0])
else:
    S = None
    _L = open(os.path.join(base, cfg['io']['mesh'])).read().split('\n')
    N, NF, CEL = map(int, _L[0].split()[:3])
    i0, i1, ic = _L.index('(node)'), _L.index('(edge)'), _L.index('(cell)')
    pts = np.array([[float(v) for v in _L[i0 + 1 + k].split()[1:3]] for k in range(N)])
    faces, cells, cur = {}, {}, None
    for k in range(i1 + 1, ic):
        t = _L[k]
        if t in ('(end)', ''):
            continue
        vv = t.split()
        if len(vv) == 1:
            cur = t
        elif len(vv) >= 4:
            faces[int(vv[0])] = (int(vv[1]), int(vv[2]))
    for k in range(ic + 1, ic + 1 + CEL):
        vv = _L[k].split()
        cells[int(vv[0])] = [int(x) for x in vv[1:5]]

    def cell_ring(fs):
        n0, n1 = faces[fs[0]]
        o, used = [n0, n1], {fs[0]}
        while len(o) < 4:
            hit = None
            for fc in fs:
                if fc in used:
                    continue
                p, q = faces[fc]
                if p == o[-1]:
                    hit = q
                elif q == o[-1]:
                    hit = p
                if hit is not None:
                    o.append(hit)
                    used.add(fc)
                    break
            if hit is None:
                return None
        return o

    quads = [q for q in (cell_ring(cells[c]) for c in range(1, CEL + 1)) if q]
    tri = np.array([[q[0] - 1, q[1] - 1, q[2] - 1] for q in quads]
                   + [[q[0] - 1, q[2] - 1, q[3] - 1] for q in quads])

def col_names(fn):
    """从 VARIABLES 行读列名, 兼容不同版本的 dump_field(9 列 / 12 列)"""
    import re
    with open(fn) as f:
        f.readline()
        return re.findall(r'"([^"]+)"', f.readline())


KEYS = []  # 由首帧的实际列名决定


def load(fn):
    df = np.loadtxt(fn, skiprows=2)
    if S:
        raw = {k: df[:, i].reshape(-1, S) for i, k in enumerate(KEYS)}

        def ds(v):
            return 0.5 * (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1))

        def dn(v):
            return np.gradient(v, axis=0)

        jac = ds(raw['x']) * dn(raw['y']) - dn(raw['x']) * ds(raw['y'])
        vx = (ds(raw['v']) * dn(raw['y']) - dn(raw['v']) * ds(raw['y'])) / jac
        uy = (dn(raw['u']) * ds(raw['x']) - ds(raw['u']) * dn(raw['x'])) / jac
        raw['vorticity'] = (vx - uy) * LREF / U
        # 收尾闭合, 免得 contourf 在接缝处缺一条
        for k in list(raw):
            raw[k] = np.hstack([raw[k], raw[k][:, :1]])
        return raw
    raw = {k: np.zeros(len(pts)) for k in KEYS}
    for i, k in enumerate(KEYS):
        raw[k] = df[:, i]
    raw['Cp'] = (raw['p'] - far['p']) / q_inf

    def ng(f, tag):
        t0, t1, t2 = tri[:, 0], tri[:, 1], tri[:, 2]
        x1, y1 = pts[t0, 0], pts[t0, 1]
        x2, y2 = pts[t1, 0], pts[t1, 1]
        x3, y3 = pts[t2, 0], pts[t2, 1]
        det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
        f1, f2, f3 = f[t0], f[t1], f[t2]
        gx = ((f2 - f1) * (y3 - y1) - (f3 - f1) * (y2 - y1)) / det
        gy = ((x2 - x1) * (f3 - f1) - (x3 - x1) * (f2 - f1)) / det
        sx, sy, w = np.zeros(len(pts)), np.zeros(len(pts)), np.zeros(len(pts))
        for t in (t0, t1, t2):
            np.add.at(sx, t, gx)
            np.add.at(sy, t, gy)
            np.add.at(w, t, 1.0)
        return sx / w, sy / w

    ux, uy = ng(raw['u'], 'u')
    vx, vy = ng(raw['v'], 'v')
    raw['vorticity'] = (vx - uy) * LREF / U
    gx, gy = ng(raw['rho'], 'rho')
    grad = np.hypot(gx, gy)
    raw['schlieren'] = np.exp(-15 * grad / np.percentile(grad, 99.5))
    return raw


files = sorted(glob.glob(os.path.join(out, 'step_*.dat')))[::a.stride]
if not files:
    raise SystemExit(f'{out} 里没有 step_*.dat')
KEYS = col_names(files[0])
print(f'共 {len(files)} 帧, 字段 {a.field}, 列 {KEYS}')

# ---- 视图与色标 ----
d0 = load(files[0])
if a.view:
    x0, x1, y0, y1 = (float(v) for v in a.view.split(','))
elif a.case:
    x0, x1, y0, y1 = VIEWS[a.case]
else:
    x0, x1 = float(d0['x'].min()), float(d0['x'].max())
    y0, y1 = float(d0['y'].min()), float(d0['y'].max())
if a.field == 'vorticity':
    # 色标自动定: 排除壁面/远场的极端值, 取 |ω| 的 99 分位; 全套帧共用同一色标
    if a.level > 0:
        L = a.level
    else:
        L = float(np.nanpercentile(np.abs(d0['vorticity']), 95.0))
        L = round(L, 1) or 1.0
        print(f'涡量色标自动: ±{L:.1f}')
    levels = np.linspace(-L, L, 41)
    cmap, extend = 'RdBu_r', 'both'
    bar = r'$\omega D/U$'
else:
    allv = np.concatenate([load(f)[a.field].ravel() for f in files[::max(1, len(files)//8)]])
    lo, hi = float(np.nanpercentile(allv, 0.5)), float(np.nanpercentile(allv, 99.5))
    if a.field == 'schlieren':
        lo, hi = 0.0, 1.0
    levels = np.linspace(lo, hi, 41)
    cmap, extend = ('gray', 'neither') if a.field == 'schlieren' else ('jet', 'both')
    bar = a.field

fig, ax = plt.subplots(figsize=(a.size, a.size * (y1 - y0) / (x1 - x0)))
path = os.path.join(out, f'anim_{a.field}.gif')
w = PillowWriter(fps=a.fps)
w.setup(fig, path, dpi=110)

times = {}
hp = os.path.join(out, 'history.csv')
if os.path.exists(hp):
    h = np.genfromtxt(hp, delimiter=',', names=True, dtype=None, encoding=None)
    if 'time' in (h.dtype.names or ()):
        times = dict(zip(h['step'], h['time']))

for i, fn in enumerate(files):
    d = load(fn)
    st = int(os.path.basename(fn)[5:11])
    ax.clear()
    if S:
        c = ax.contourf(d['x'], d['y'], d[a.field], levels=levels, cmap=cmap, extend=extend,
                        antialiased=False)
    else:
        c = ax.tricontourf(pts[:, 0], pts[:, 1], tri, d[a.field], levels=levels, cmap=cmap,
                           extend=extend, antialiased=False)
    ax.fill(np.cos(np.linspace(0, 2 * np.pi, 200)) * 0.5,
            np.sin(np.linspace(0, 2 * np.pi, 200)) * 0.5, color='w', edgecolor='k', lw=0.8,
            zorder=5)
    ax.set_xlim(x0, x1)
    ax.set_ylim(y0, y1)
    ax.set_aspect('equal')
    tt = times.get(st)
    ax.set_title(f'step {st}' + (f'   t = {tt:.3f} s' if tt is not None else ''), fontsize=11)
    if i == 0:
        fig.colorbar(c, ax=ax, fraction=0.03, pad=0.02).set_label(bar, fontsize=10)
    w.grab_frame()
    if (i + 1) % 8 == 0:
        print(f'  {i + 1}/{len(files)}')

w.finish()
plt.close(fig)
print(f'写出 {path}  ({os.path.getsize(path) / 1e6:.1f} MB)')
