#!/usr/bin/env python3
"""三款算例统一后处理: 圆柱 / NACA0012 / OAT15A.

用法:
  python tests/post.py cases/cylinder/urans_d1.json          --case cylinder
  python tests/post.py cases/naca0012/urans_m07a5.json       --case naca0012
  python tests/post.py cases/oat15a/rans_m073_re33e6.json    --case oat15a
可加 --readme --wall 170 写 README.md(含壁钟耗时); --u-positive/--u-symmetric 覆盖默认色标.

三款算例的差异集中在 CASES 表: 网格是否结构化、有无翼型、定常或 URANS、单位与视图、
u 色标默认值; 其余(力系数历程/主频、壁面 Cp 上下表面剥离、激波位置、云图面板)全部共用.
非结构网格的壁面剥离不能按 y 正负粗分: OAT 网格的壁面环里除翼型表面还缠着钝尾缘小面
和一条 x/c 0.75~1.0 的尾缘楔形内层曲线, 故统一走"网格真实连接 -> 从 LE 取 x 单调段".
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
from matplotlib.ticker import MultipleLocator

GAMMA, R = 1.4, 287.05

CASES = {
  'cylinder': dict(
    structured=True, airfoil=False, phase='U', kind='angle',
    xlabel='$x/D$', ylabel='$y/D$', view=(-5.0, 15.0, -10.0, 10.0), u_positive=False,
    note='圆柱绕流',
  ),
  'naca0012': dict(
    structured=True, airfoil=True, phase='U', kind='surface',
    xlabel='$x/c$', ylabel='$y/c$', view=(-0.7, 0.8, -0.35, 0.55), u_positive=True,
    note='NACA0012 URANS',
  ),
  'oat15a': dict(
    structured=False, airfoil=True, phase='S', kind='surface',
    xlabel='$x/c$', ylabel='$y/c$', view=(-0.5, 1.5, -0.6, 0.6), u_positive=True,
    note='OAT15A 定常 RANS',
  ),
}

ap = argparse.ArgumentParser()
ap.add_argument('config')
ap.add_argument('--case', required=True, choices=list(CASES))
ap.add_argument('--u-positive', dest='u_pos', action='store_true', default=None,
                help='u 色标只覆盖正值(负值走 under)')
ap.add_argument('--u-symmetric', dest='u_pos', action='store_false', default=None,
                help='u 色标对称')
ap.add_argument('--panels', default='auto', help="云图面板, 'auto' 或逗号分隔的名字, 'none' 跳过")
ap.add_argument('--readme', action='store_true', help='生成 README.md')
ap.add_argument('--phase', default='auto', choices=['auto', 'U', 'S'],
                help="用哪一段历史: auto 按 U 行数自动判别(定常算例只有 1 行 U)")
ap.add_argument('--out', default=None, help='覆盖输出目录(默认取 config 的 io.field)')
ap.add_argument('--wall', type=float, default=float('nan'), help='壁钟耗时 [s], 写入 README')
ap.add_argument('--probe', type=float, default=0.45,
                help='壁面压力探针位置 x/c (输出 Cp(t) 与 PSD, 默认 0.45)')
a = ap.parse_args()
P = CASES[a.case]
if a.u_pos is None:
  a.u_pos = P['u_positive']

base = os.path.dirname(os.path.abspath(a.config))
cfg = json.load(open(a.config))
out = a.out or os.path.normpath(os.path.join(base, cfg['io']['field']))
far = cfg['farfield']
U = far['Ma'] * math.sqrt(GAMMA * R * far['T'])
q_inf = 0.5 * GAMMA * far['p'] * far['Ma'] ** 2
print(f'算例 {a.case}: {out}')

# ============================ 网格 ============================
L = open(os.path.join(base, cfg['io']['mesh'])).read().split('\n')
N, NF, CEL = map(int, L[0].split()[:3])
i0, i1, ic = L.index('(node)'), L.index('(edge)'), L.index('(cell)')
# 头部名字->类型映射: 圆柱 'cylinder = WALL' / 翼型 'airfoil = WALL' / OAT 'wall = WALL',
# 组名不统一, 故按类型挑出壁面组
gmap = {}
for ln in L[1:i0]:
  if '=' in ln:
    nm, tp = (s.strip() for s in ln.split('=', 1))
    gmap[nm] = tp
WALLG = {nm for nm, tp in gmap.items() if tp.upper().startswith('WALL')}
pts = np.array([[float(v) for v in L[i0 + 1 + k].split()[1:3]] for k in range(N)])
faces, cells, cur, wallf = {}, {}, None, []
for k in range(i1 + 1, ic):
  t = L[k]
  if t in ('(end)', ''):
    continue
  v = t.split()
  if len(v) == 1:
    cur = t
  elif len(v) >= 4:
    faces[int(v[0])] = (int(v[1]), int(v[2]))
    if cur in WALLG:
      wallf.append((int(v[1]), int(v[2])))
for k in range(ic + 1, ic + 1 + CEL):
  v = L[k].split()
  cells[int(v[0])] = [int(x) for x in v[1:5]]
xw = pts[[n - 1 for pq in wallf for n in pq], 0]
X0, CHORD = xw.min(), xw.max() - xw.min()
LREF = 1.0 if a.case != 'oat15a' else CHORD  # 求解器力系数参考长度硬编码为 1

# 壁面节点环(所有壁面节点度数为 2): 从 LE 起顺序走一圈
adj = defaultdict(list)
for p, q in wallf:
  adj[p].append(q)
  adj[q].append(p)
start = min(adj, key=lambda n: pts[n - 1, 0])
order, prev, cur = [start], None, start
while True:
  nxt = [n for n in adj[cur] if n != prev]
  if not nxt or nxt[0] == start:
    break
  order.append(nxt[0])
  prev, cur = cur, nxt[0]


def monotone(seq):
  """从 LE 沿环走 x 单调不减的一段, 到 TE(x 首次减小)为止."""
  o = [seq[0]]
  for n in seq[1:]:
    if pts[n - 1, 0] < pts[o[-1] - 1, 0] - 1e-12:
      break
    o.append(n)
  return o


if P['kind'] == 'surface':
  s1, s2 = monotone(order), monotone([order[0]] + order[:0:-1])
  if np.mean([pts[n - 1, 1] for n in s1]) < np.mean([pts[n - 1, 1] for n in s2]):
    s1, s2 = s2, s1
  sides = [('upper', np.array(s1)), ('lower', np.array(s2))]
else:
  ang = np.degrees(np.arctan2([pts[n - 1, 1] for n in order], [pts[n - 1, 0] for n in order]))
  sides = [('surface', np.array(order))]

# ============================ 力系数历程 ============================
h = np.genfromtxt(os.path.join(out, 'history.csv'), delimiter=',', names=True, dtype=None,
                  encoding=None)
if a.phase == 'auto':
  # 定常算例的 URANS 段只有 max_steps 次(通常 1 行), URANS 算例则有几百行
  P['phase'] = 'U' if (h['phase'] == 'U').sum() > 20 else 'S'
  print(f"自动判别 phase = {P['phase']}")
else:
  P['phase'] = a.phase
u = h[h['phase'] == P['phase']]
ok = np.isfinite(u['res']) & np.isfinite(u['Cl']) & np.isfinite(u['Cd'])
div_step = None
if not ok.all():                      # 中途发散: 截断到最后一个有效步, 后续统计照常出
  div_step = int(u['step'][ok][-1])
  print(f'警告: {int((~ok).sum())} 个历史行非有限(求解器发散), 已截断到 step {div_step}')
  u = u[ok]
lines = []


def peak_freq(y, dt, fmin=0.0):
  """Hann 窗 + 补零定位峰, 再按原始分辨率 df 做抛物线插值, 避免主频被钉在 bin 上.
  fmin 以下不计入: 抖振信号的建立漂移会淹没主峰."""
  y = np.asarray(y, float) - np.mean(y)
  n = len(y)
  df = 1.0 / (n * dt)
  Nf = 1 << int(np.ceil(np.log2(n * 16)))
  mag = np.abs(np.fft.rfft(y * np.hanning(n), Nf))
  fr = np.fft.rfftfreq(Nf, dt)
  k0 = max(1, int(np.searchsorted(fr, fmin)))
  k = k0 + int(np.argmax(mag[k0:]))
  f0 = fr[k]

  def at(f):
    return mag[int(np.argmin(np.abs(fr - f)))]

  x0, x1, x2 = np.log(at(f0 - df)), np.log(at(f0)), np.log(at(f0 + df))
  den = x0 - 2 * x1 + x2
  if not np.isfinite(den) or abs(den) < 1e-12:
    return f0, df
  f = f0 + 0.5 * (x0 - x2) / den * df
  return (f if np.isfinite(f) and f > 0 else f0), df


cl, cd = u['Cl'] / LREF, u['Cd'] / LREF
f0 = df = float('nan')
if P['phase'] == 'U':
  t = u['time']
  half = len(t) // 2
  f0, df = peak_freq(cl[half:], t[1] - t[0])
  lines += [
    f'物理步数 {len(t)}, t_end = {t[-1]:.4f} s, 平均内迭代 {u["inner"].mean():.1f} 次, '
    f'内迭代残差下降中位数 {np.median(u["res"] / u["res0"]):.2e}',
    f'后半段 Cl: 均值 {cl[half:].mean():.4f}, 范围 [{cl[half:].min():.4f}, {cl[half:].max():.4f}]',
    f'后半段 Cd: 均值 {cd[half:].mean():.4f}, 范围 [{cd[half:].min():.4f}, {cd[half:].max():.4f}]',
    f'Cl 主频 {f0:.2f} Hz (分辨率 {df:.2f} Hz), St = fL/U = {f0 * LREF / U:.4f}, '
    f'k = pi f L/U = {math.pi * f0 * LREF / U:.4f}',
  ]
  tlab, xarr = t, None
else:
  lines += [
    f'定常步数 {len(u)}, 末端残差 {u["res"][-1]:.3e}',
    f'Cl = {cl[-1]:.5f}  (求解器输出 {u["Cl"][-1]:.5f} x 1/{LREF:.6f})',
    f'Cd = {cd[-1]:.5f}  (求解器输出 {u["Cd"][-1]:.5f} x 1/{LREF:.6f})',
  ]
  tlab = u['step']

# ============================ 壁面 Cp: 面中心 + 分环 ============================
# wall_cp.dat 由求解器按 solver.urans.wall_interval 输出; 未开则为空, 此时跳过壁面相关输出
blk, cpn, cp_sides = [], {}, []
wp = os.path.join(out, 'wall_cp.dat')
shock_lines = []
st = sx = rms = None
cp_by_step = {}
if os.path.exists(wp):
  blk = open(wp).read().split('# step ')
  if len(blk) < 2 or not blk[-1].strip().split('\n')[1:]:
    print('警告: wall_cp.dat 为空(求解器未跑到正常结束), 跳过壁面 Cp / 激波统计')
    blk = []
  cur_cp = np.array([r.split() for r in blk[-1].strip().split('\n')[1:]], dtype=float)[:, 2]
  assert len(cur_cp) == len(wallf), (len(cur_cp), len(wallf))
  acc = defaultdict(list)
  for i, (p, q) in enumerate(wallf):
    acc[p].append(cur_cp[i])
    acc[q].append(cur_cp[i])
  cpn = {n: float(np.mean(v)) for n, v in acc.items()}
else:
  print(f'无 {os.path.basename(wp)}(配置未开 wall_interval), 跳过壁面 Cp / 激波统计')


def block_cp(b):
  """wall_cp.dat 的一个 block -> (step, t, 每面 Cp)"""
  rw = b.strip().split('\n')
  return int(rw[0].split()[0]), float(rw[0].split()[2]), \
      np.array([r.split() for r in rw[1:]], dtype=float)[:, 2]


def side_profiles(cpf):
  """每面 Cp -> [(name, x/c, Cp)], 节点值取相邻两面平均"""
  a_ = defaultdict(list)
  for i, (p, q) in enumerate(wallf):
    a_[p].append(cpf[i])
    a_[q].append(cpf[i])
  res = []
  for nm, ids in sides:
    x = (pts[ids - 1, 0] - X0) / CHORD
    cp = np.array([float(np.mean(a_[n])) for n in ids])
    o = np.argsort(x, kind='stable')
    res.append((nm, x[o], cp[o]))
  return res


def shock_x(xb, cpb):
  """dCp/dx 最大处; 钝尾缘小面会让 x 重复, 先去重"""
  kp = np.concatenate([[True], np.diff(xb) > 1e-12])
  xh, ch = xb[kp], cpb[kp]
  g = np.where((xh > 0.05) & (xh < 0.85), np.gradient(ch, xh), -np.inf)
  return xh[np.argmax(g)]


if P['kind'] == 'surface' and cpn:
  cp_sides = side_profiles(cur_cp)
  blocks = blk[1:]
  for b in blocks:
    stp, tb, cp_b = block_cp(b)
    cp_by_step[stp] = (tb, side_profiles(cp_b))
  if P['phase'] == 'U' and len(blocks) > 1:
    st = np.array([block_cp(b)[1] for b in blocks])
    sx = np.array([shock_x(pr[0][1], pr[0][2]) for pr in
                   [cp_by_step[block_cp(b)[0]][1] for b in blocks]])
    sh = sx[len(sx) // 2:]
    fs, _ = peak_freq(sh, st[1] - st[0])
    shock_lines += [
      f'后半段上表面激波位置 x/c: 均值 {sh.mean():.3f}, 范围 [{sh.min():.3f}, {sh.max():.3f}]',
      f'激波位置主频 {fs:.2f} Hz (与 Cl 一致)',
    ]
    # 上表面 Cp 脉动 RMS: 抖振强度指标
    grid = np.linspace(0.0, 1.0, 201)
    tail = [b for b in blocks][len(blocks) // 2:]
    series = [np.interp(grid, cp_by_step[block_cp(b)[0]][1][0][1],
                        cp_by_step[block_cp(b)[0]][1][0][2]) for b in tail]
    rms = np.std(np.array(series), axis=0)
    band = np.where((grid > 0.05) & (grid < 0.95), rms, -1.0)
    shock_lines.append(f'上表面 Cp 脉动 RMS: 最大 {rms[np.argmax(band)]:.3f} '
                       f'@ x/c = {grid[np.argmax(band)]:.3f}')
  else:
    _, xs, cs = cp_sides[0]
    shock_lines.append(f'上表面激波位置 x/c = {shock_x(xs, cs):.3f}')

# ============================ 云图后端 ============================
p_inf = far['p']
if P['structured']:
  S = int(open(os.path.join(base, cfg['io']['structured'])).readline().split()[0])
  files = sorted(glob.glob(os.path.join(out, 'step_*.dat')))
  times = dict(zip(u['step'], zip(u['time'], u['Cl']))) if P['phase'] == 'U' else {}

  def load(fn):
    f = np.loadtxt(fn, skiprows=2)
    d = {k: f[:, i].reshape(-1, S) for i, k in
         enumerate(['x', 'y', 'rho', 'u', 'v', 'T', 'p', 'Ma', 'miubl'])}
    d['Cp'] = (d['p'] - p_inf) / q_inf

    def ds(v):
      return 0.5 * (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1))

    def dn(v):
      return np.gradient(v, axis=0)

    x, y, rr = d['x'], d['y'], d['rho']
    jac = ds(x) * dn(y) - dn(x) * ds(y)
    gx = (ds(rr) * dn(y) - dn(rr) * ds(y)) / jac
    gy = (dn(rr) * ds(x) - ds(rr) * dn(x)) / jac
    grad = np.hypot(gx, gy)
    d['schlieren'] = np.exp(-15 * grad / np.percentile(grad, 99.5))
    vx = (ds(d['v']) * dn(y) - dn(d['v']) * ds(y)) / jac
    uy = (dn(d['u']) * ds(x) - ds(d['u']) * dn(x)) / jac
    d['vorticity'] = (vx - uy) * LREF / U
    return {k: np.hstack([v, v[:, :1]]) for k, v in d.items()}

  def draw(ax, z, levels, cmap, extend):
    c = ax.contourf(zx, zy, z, levels=levels, cmap=cmap, extend=extend, antialiased=False)
    if cmap != 'gray':
      ax.contour(zx, zy, z, levels=levels[1::2], colors='k', linewidths=0.45, alpha=0.7)
    return c

  def sonic(ax, z):
    """翼型 Ma 云图: Ma=1 音速线加粗"""
    ax.contour(zx, zy, z, levels=[1.0], colors='k', linewidths=1.8)

  def outline(ax):
    ax.fill(ww[0], ww[1], color='w', edgecolor='k', linewidth=0.8, zorder=5)
else:
  def load(fn):
    f = np.loadtxt(fn, skiprows=2)
    return {k: f[:, i] for i, k in
            enumerate(['x', 'y', 'rho', 'u', 'v', 'T', 'p', 'Ma', 'miubl'])}

  files = sorted(glob.glob(os.path.join(out, 'step_*.dat')))
  times = dict(zip(u['step'], zip(u['time'], u['Cl']))) if P['phase'] == 'U' else {}

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

  def nodal_avg(cellval):
    s, n = np.zeros(N), np.zeros(N)
    flat = np.array([v - 1 for q in quads for v in q])
    np.add.at(s, flat, np.repeat(cellval, 4))
    np.add.at(n, flat, 1.0)
    return np.divide(s, n, out=np.zeros(N), where=n > 0)

  def nodal_grad(f):
    """三角网线性梯度累加到节点 (矢量化的 2x2 求解)."""
    t0, t1, t2 = tri[:, 0], tri[:, 1], tri[:, 2]
    x1, y1, x2, y2, x3, y3 = (pts[t0, 0], pts[t0, 1], pts[t1, 0], pts[t1, 1],
                              pts[t2, 0], pts[t2, 1])
    det = (x2 - x1) * (y3 - y1) - (x3 - x1) * (y2 - y1)
    f1, f2, f3 = f[t0], f[t1], f[t2]
    gx = ((f2 - f1) * (y3 - y1) - (f3 - f1) * (y2 - y1)) / det
    gy = ((x2 - x1) * (f3 - f1) - (x3 - x1) * (f2 - f1)) / det
    sx, sy, w = np.zeros(N), np.zeros(N), np.zeros(N)
    for t in (t0, t1, t2):
      np.add.at(sx, t, gx)
      np.add.at(sy, t, gy)
      np.add.at(w, t, 1.0)
    return sx / w, sy / w

  def load(fn):
    raw = np.loadtxt(fn, skiprows=2)
    d = {k: nodal_avg(raw[:, i]) for i, k in
         enumerate(['x', 'y', 'rho', 'u', 'v', 'T', 'p', 'Ma', 'miubl'])}
    d['x'] = pts[:, 0]
    d['y'] = pts[:, 1]
    d['Cp'] = (d['p'] - p_inf) / q_inf
    gx, gy = nodal_grad(d['rho'])
    grad = np.hypot(gx, gy)
    d['schlieren'] = np.exp(-15 * grad / np.percentile(grad, 99.5))
    ux, uy = nodal_grad(d['u'])
    vx, vy = nodal_grad(d['v'])
    d['vorticity'] = (vx - uy) * LREF / U
    return d

  def draw(ax, z, levels, cmap, extend):
    # 与结构化(cyc/naca)后端保持同一套等值线风格: 21 级 + 隔级黑色等值线
    xy = pts / LREF
    c = ax.tricontourf(xy[:, 0], xy[:, 1], tri, z, levels=levels, cmap=cmap,
                       extend=extend, antialiased=False)
    if cmap != 'gray':
      ax.tricontour(xy[:, 0], xy[:, 1], tri, z, levels=levels[1::2], colors='k',
                    linewidths=0.45, alpha=0.7)
    return c

  def sonic(ax, z):
    """翼型 Ma 云图: Ma=1 音速线加粗"""
    xy = pts / LREF
    ax.tricontour(xy[:, 0], xy[:, 1], tri, z, levels=[1.0], colors='k', linewidths=1.8)

  def outline(ax):
    poly = pts[[n - 1 for n in order]] / LREF
    ax.fill(poly[:, 0], poly[:, 1], color='w', edgecolor='k', lw=0.7, zorder=5)


def label(fn):
  st = int(os.path.basename(fn)[5:11])
  if st in times:
    return f'step {st}  t={times[st][0]:.4f}s  Cl={times[st][1] / LREF:.4f}'
  return f'step {st}'


def last_valid(files, back=4):
  """最后一个全有限值的场文件(发散后最后几帧可能含 NaN)"""
  for fn in reversed(files[-back:]):
    try:
      if np.isfinite(np.loadtxt(fn, skiprows=2)).all():
        return fn
    except Exception:
      pass
  return None


def body(ax):
  outline(ax)
  ax.set_xlim(P['view'][0], P['view'][1])
  ax.set_ylim(P['view'][2], P['view'][3])
  ax.set_aspect('equal')


# ============================ 历程图 ============================
nrow = 2 + (1 if st is not None else 0)
fig, ax = plt.subplots(nrow, 1, figsize=(9, 2.6 * nrow), sharex=True, squeeze=False)
ax = ax[:, 0]
ax[0].plot(tlab, cl, lw=1)
ax[0].set_ylabel('$C_l$')
ax[1].plot(tlab, cd, lw=1, color='C1')
ax[1].set_ylabel('$C_d$')
if nrow == 3:
  ax[2].plot(st, sx, lw=1, color='C2')
  ax[2].set_ylabel('shock $x/c$')
ax[-1].set_xlabel('t [s]' if P['phase'] == 'U' else 'steady step')
fig.tight_layout()
fig.savefig(os.path.join(out, 'history.png'), dpi=130)
plt.close(fig)

# ============================ 残差曲线 ============================
# history.csv: phase,step,time,inner,res0,res,Cl,Cd
# 定常段 res0 = res = 该步残差; URANS 段 res0/res = 内迭代初/终残差
s_all = h[h['phase'] == 'S']
rows = []
if len(s_all) > 1:
  rows.append(('steady', P['phase'] == 'S', s_all['step'], s_all['res'], None, 'steady step',
               f'steady residual {s_all["res"][0]:.2e} -> {s_all["res"][-1]:.2e}'))
if len(u) > 1 and P['phase'] == 'U':
  rows.append(('urans', False, u['time'], u['res'], u['res0'], 't [s]',
               f'inner residual, median drop {np.median(u["res"] / u["res0"]):.2e}'))
nrow = len(rows)
fig, axs = plt.subplots(nrow, 1, figsize=(9, 2.4 * nrow), squeeze=False)
for k, (_, _, xs, rr, r0, xlab, ttl) in enumerate(rows):
  a_ = axs[k, 0]
  if r0 is not None:
    a_.semilogy(xs, r0, lw=0.8, color='C0', alpha=0.45, label='res0 (inner start)')
  a_.semilogy(xs, rr, lw=1.0, color='C3', label='res (after inner iters)')
  a_.set_ylabel('residual')
  a_.set_xlabel(xlab)
  a_.grid(alpha=0.3, which='both')
  a_.set_title(ttl, fontsize=10)
  a_.legend(fontsize=8)
  if r0 is not None:
    a_.set_ylim(max(rr.min(), 1e-12) * 0.5, max(r0.max(), rr.max()) * 2)
fig.tight_layout()
fig.savefig(os.path.join(out, 'residual.png'), dpi=130)
plt.close(fig)
print('写出 residual.png')

# ============================ Cp 分布 ============================
if not cpn:
  print('无壁面 Cp 数据, 跳过 cp_dist.png')
else:
  fig, z = plt.subplots(figsize=(8.4, 4.4))
  if P['kind'] == 'angle':
    ids = sides[0][1]
    z.plot(ang, np.array([cpn[n] for n in ids]), '.-', ms=3, lw=0.8)
    z.set_xlabel(r'$\theta$ [deg]')
  else:
    for name, xs, cs in cp_sides:
      z.plot(xs, cs, '-', lw=1.0, label=name)
    z.set_xlabel('$x/c$')
    z.legend()
  z.invert_yaxis()
  z.set_ylabel('$C_p$')
  if P['kind'] == 'surface':
    # 与 cp_rms 同一套: 主刻度 0.1 一个数字, 之间 5 个小格(0.02); 不画网格
    z.xaxis.set_major_locator(MultipleLocator(0.1))
    z.xaxis.set_minor_locator(MultipleLocator(0.02))
    z.tick_params(which='both', direction='in', top=True, right=True)
    z.tick_params(which='major', length=6)
    z.tick_params(which='minor', length=3)
  z.set_title(f'{a.case}  Ma={far["Ma"]} alpha={far["alpha"]}deg  (face-centre $C_p$)')
  fig.tight_layout()
  fig.savefig(os.path.join(out, 'cp_dist.png'), dpi=140)
  plt.close(fig)
  print('写出 cp_dist.png')

if rms is not None:
  fig, z = plt.subplots(figsize=(8.4, 4.0))
  z.plot(grid, rms, lw=1.2)
  z.set_xlabel('$x/c$')
  z.set_ylabel(r'$C_p$ RMS (upper surface)')
  # 主刻度 0.1 一个数字, 之间 5 个小格(0.02); 不画网格
  z.xaxis.set_major_locator(MultipleLocator(0.1))
  z.xaxis.set_minor_locator(MultipleLocator(0.02))
  z.tick_params(which='both', direction='in', top=True, right=True)
  z.tick_params(which='major', length=6)
  z.tick_params(which='minor', length=3)
  z.set_title(f'{a.case}  upper-surface $C_p$ fluctuation  Ma={far["Ma"]} '
              f'alpha={far["alpha"]}deg')
  fig.tight_layout()
  fig.savefig(os.path.join(out, 'cp_rms.png'), dpi=140)
  plt.close(fig)
  np.savetxt(os.path.join(out, 'cp_rms.dat'), np.column_stack([grid, rms]),
             header='x/c  Cp_RMS(upper)', fmt='%.6e')
  print('写出 cp_rms.png, cp_rms.dat')

# ============================ 壁面探针: Cp(t) 与 PSD ============================
# 上表面 x/c 探针的 Cp 时间序列 + 单段 Hann 窗 PSD; 抖振主频/St 直接从探针得到
probe_lines = []
if cp_by_step and P['phase'] == 'U' and len(cp_by_step) > 8:
  stp = np.array(sorted(cp_by_step))
  tp = np.array([cp_by_step[s][0] for s in stp])
  sp = np.array([np.interp(a.probe, cp_by_step[s][1][0][1], cp_by_step[s][1][0][2])
                 for s in stp])
  keep = np.isfinite(sp) & np.isfinite(tp)
  stp, tp, sp = stp[keep], tp[keep], sp[keep]
  if len(tp) > 8 and tp[-1] > tp[0]:
    dtp = float(np.median(np.diff(tp)))
    df0 = 1.0 / (len(sp) * dtp)
    # 低频漂移(建立瞬态)会淹没主峰, 所以主峰搜索从 4 个 bin 以上开始
    fpk, dfpk = peak_freq(sp, dtp, fmin=4.0 * df0)
    rms_p = float(np.std(sp))
    win = np.hanning(len(sp))
    Yf = np.fft.rfft((sp - sp.mean()) * win)
    fr = np.fft.rfftfreq(len(sp), dtp)
    S = 2.0 * np.abs(Yf) ** 2 * dtp / np.sum(win ** 2)
    lofrac = float(S[fr <= 20.0].sum() / S.sum()) if S.sum() > 0 else 0.0
    ipk = int(np.argmin(np.abs(fr - fpk)))
    fig2, ax2 = plt.subplots(2, 1, figsize=(9.0, 5.6))
    ax2[0].plot(tp, sp, lw=1.0, color='C0')
    ax2[0].set_xlabel('t [s]')
    ax2[0].set_ylabel(f'$C_p$ @ $x/c={a.probe:g}$')
    ax2[0].set_title(f'{a.case}  Ma={far["Ma"]} alpha={far["alpha"]}deg  '
                     f'wall probe $x/c={a.probe:g}$  ($C_p$ RMS={rms_p:.3f})', fontsize=10)
    ax2[1].loglog(fr[fr > 0], S[fr > 0], lw=1.0, color='C3')
    ax2[1].axvline(fpk, color='k', lw=0.8, ls='--', alpha=0.7)
    ax2[1].set_xlim(fr[fr > 0].min(), min(1000.0, fr[-1]))
    ax2[1].set_xlabel('f [Hz]')
    ax2[1].set_ylabel('PSD [$C_p^2$/Hz]')
    ax2[1].annotate(f'peak {fpk:.1f} Hz\nSt={fpk * LREF / U:.4f}',
                    xy=(fpk, S[ipk]), xytext=(0.72, 0.75), textcoords='axes fraction',
                    fontsize=9, arrowprops=dict(arrowstyle='->', lw=0.8))
    for axx in ax2:
      axx.tick_params(which='both', direction='in', top=True, right=True)
    fig2.tight_layout()
    fig2.savefig(os.path.join(out, 'cp_probe_pt_psd.png'), dpi=140)
    plt.close(fig2)
    # 单独的时域图
    fig4, ax4 = plt.subplots(figsize=(8.4, 3.6))
    ax4.plot(tp, sp, lw=1.1, color='C0')
    ax4.set_xlabel('t [s]')
    ax4.set_ylabel(f'$C_p$ @ $x/c={a.probe:g}$')
    ax4.tick_params(which='both', direction='in', top=True, right=True)
    ax4.set_title(f'{a.case}  wall $C_p$ history @ $x/c={a.probe:g}$   Ma={far["Ma"]}  '
                  f'alpha={far["alpha"]}deg   steps {int(stp[0])}~{int(stp[-1])}   '
                  f'RMS={rms_p:.3f}  ({sp.min():.3f} ~ {sp.max():.3f})', fontsize=10)
    fig4.tight_layout()
    fig4.savefig(os.path.join(out, 'cp_probe_t.png'), dpi=150)
    plt.close(fig4)
    # 独立 PSD 图: 对数-对数坐标, 标出 ±Δf 不确定带与 St 副轴
    nper = (tp[-1] - tp[0]) * fpk
    fig3, ax3 = plt.subplots(figsize=(8.4, 4.8))
    pos = fr > 0
    ax3.loglog(fr[pos], S[pos], lw=1.1, color='C3', zorder=3)
    ax3.axvspan(max(fr[pos].min(), fpk - dfpk), fpk + dfpk, color='0.75', alpha=0.5,
                zorder=1, label=f'±Δf = {dfpk:.1f} Hz')
    ax3.axvline(fpk, color='k', lw=0.9, ls='--', alpha=0.8, zorder=4)
    ax3.set_xlim(fr[pos].min(), min(1000.0, fr[-1]))
    ax3.set_ylim(S.max() * 1e-5, S.max() * 40)
    # Welch 分段平均: 低方差对照(细线=单段高分辨率, 粗线=平均后的趋势)
    nseg = 4
    seg_n = len(sp) // nseg
    if seg_n >= 32:
        ww = np.hanning(seg_n)
        fw = np.fft.rfftfreq(seg_n, dtp)
        Sw = np.zeros_like(fw)
        for i in range(nseg):
            seg = sp[i * seg_n:(i + 1) * seg_n]
            Yw = np.fft.rfft((seg - seg.mean()) * ww)
            Sw += 2.0 * np.abs(Yw) ** 2 * dtp / np.sum(ww ** 2)
        Sw /= nseg
        ax3.loglog(fw[1:], Sw[1:], lw=2.2, color='k', alpha=0.6, zorder=5,
                   label=f'Welch {nseg} segs avg (df={fw[1]:.1f} Hz)')
    ax3.set_xlabel('f [Hz]')
    ax3.set_ylabel('PSD [$C_p^2$/Hz]')
    ax3.tick_params(which='both', direction='in', right=True)
    ax3.tick_params(which='major', length=6)
    ax3.tick_params(which='minor', length=3)
    ax3.legend(fontsize=9, loc='upper right')
    sec = ax3.secondary_xaxis('top',
                              functions=(lambda f: f * LREF / U, lambda st: st * U / LREF))
    sec.set_xlabel('St = fL/U')
    sec.tick_params(direction='in')
    ax3.annotate(f'{fpk:.1f} Hz\nSt = {fpk * LREF / U:.4f}\nPSD = {S[ipk]:.2e}',
                 xy=(fpk, S[ipk]), xytext=(fpk * 2.4, S[ipk] * 0.10), fontsize=9,
                 arrowprops=dict(arrowstyle='->', lw=0.8, shrinkA=2, shrinkB=2))
    ax3.set_title(f'{a.case}  wall $C_p$ PSD @ $x/c={a.probe:g}$   Ma={far["Ma"]}  '
                  f'alpha={far["alpha"]}deg\n'
                  f'steps {int(stp[0])}~{int(stp[-1])}   N={len(tp)}   fs={1 / dtp:.0f} Hz   '
                  f'df={dfpk:.2f} Hz   ({nper:.1f} cycles in window)', fontsize=10)
    fig3.tight_layout()
    fig3.savefig(os.path.join(out, 'psd.png'), dpi=150)
    plt.close(fig3)
    np.savetxt(os.path.join(out, f'cp_probe_x{a.probe:g}.dat'),
               np.column_stack([tp, sp]), header=f't[s]  Cp(x/c={a.probe:g}, upper)',
               fmt='%.6e')
    np.savetxt(os.path.join(out, f'psd_x{a.probe:g}.dat'), np.column_stack([fr, S]),
               header=f'f[Hz]  PSD(Cp, x/c={a.probe:g})', fmt='%.6e')
    probe_lines += [
      f'壁面压力探针 x/c={a.probe:g} (上表面): Cp RMS = {rms_p:.4f}, '
      f'峰谷 {sp.min():.3f} ~ {sp.max():.3f}',
      f'探针 PSD 主频 {fpk:.2f} Hz (分辨率 {dfpk:.2f} Hz), St = {fpk * LREF / U:.4f}'
      f'{"  <- 样本不足, 仅供参考" if len(tp) < 50 else ""}',
    ]
    if lofrac > 0.3:
      probe_lines.append(
        f'  [!] 该点 ≤20 Hz 低频能量占 {lofrac * 100:.0f}%, 抖振峰被建立漂移压制, '
        f'主频/St 仅供参考')
    print(f'写出 psd.png / cp_probe_t.png / cp_probe_pt_psd.png / '
          f'cp_probe_x{a.probe:g}.dat / psd_x{a.probe:g}.dat')
  else:
    print(f'探针数据不足({len(tp)} 帧), 跳过 cp_probe_pt_psd.png')
elif P['kind'] == 'surface' and P['phase'] != 'U':
  print('定常算例无 Cp 时间序列, 跳过探针 PSD')

# ============================ 云图面板 ============================
lastf = last_valid(files)
if lastf is None:
  print('所有场文件均含非有限值, 跳过云图面板')
  files = []
sa = load(lastf) if lastf else None
if sa is not None:
  if P['structured']:
    zx, zy, ww = sa['x'], sa['y'], (sa['x'][0], sa['y'][0])
  else:
    zx = zy = None
    ww = None
specs = [
  ('u', 'RdBu_r', 0.0 if a.u_pos else -1.2 * U, None if a.u_pos else 1.2 * U, r'$u$ [m/s]'),
  ('v', 'RdBu_r', -0.6 * U, 0.6 * U, r'$v$ [m/s]'),
  ('rho', 'jet', None, None, r'$\rho$ [kg/m$^3$]'),
  ('miubl', 'jet', 0.0, None, r'$\tilde{\nu}$ [m$^2$/s]'),
  ('T', 'jet', None, None, r'$T$ [K]'),
  ('Ma', 'jet', None, None, 'Ma'),
  ('Cp', 'jet', -1.6, 1.1, '$C_p$'),
  ('schlieren', 'gray', 0.0, 1.0, 'schlieren'),
  ('vorticity', 'jet', -4.0, 4.0, r'$\omega L/U$'),
]
want = None if a.panels == 'auto' else ([] if a.panels == 'none' else a.panels.split(','))
for key, cmap, lo, hi, bar in (specs if sa is not None else []):
  if want is not None and key not in want:
    continue
  arr = sa[key]
  if not np.isfinite(arr).any():
    print('跳过', key, '(全为非有限值)')
    continue
  if lo is None:
    lo = float(np.floor(np.nanmin(arr) * 20) / 20)
  if hi is None:
    hi = float(np.ceil(np.nanmax(arr) * 20) / 20)
  if not (np.isfinite(lo) and np.isfinite(hi)) or hi <= lo:
    print('跳过', key, '(色标范围无效)')
    continue
  lv = np.linspace(lo, hi, 21)
  try:
    fig, axp = plt.subplots(figsize=(8.2, 7.0))
    c = draw(axp, arr, lv, cmap, 'neither' if cmap == 'gray' else 'both')
    if key == 'Ma' and P['airfoil'] and lo < 1.0 < hi:
      sonic(axp, arr)
    body(axp)
    axp.set_xlabel(P['xlabel'], fontsize=14)
    axp.set_ylabel(P['ylabel'], fontsize=14)
    axp.tick_params(direction='in', top=True, right=True, labelsize=11)
    axp.set_title(label(lastf), fontsize=12)
    cb = fig.colorbar(c, ax=axp, orientation='horizontal', fraction=0.05, pad=0.12)
    cb.ax.tick_params(labelsize=10)
    cb.ax.set_title(bar, fontsize=11, pad=8)
    fig.tight_layout()
    name = ('T' if key == 'T' else key.lower()) + '.png'
    fig.savefig(os.path.join(out, name), dpi=150)
    plt.close(fig)
    print('写出', name)
  except Exception as e:                # 单张出问题不影响其余
    print('跳过', key, ':', str(e)[:80])

# 最后 6 帧纹影 + 壁面 Cp (仅 URANS; 发散后可能只剩前几帧, 只画全有限的)
seq = [f for f in files[-6:] if last_valid([f])] if P['phase'] == 'U' else []
if seq:
    fig, axs = plt.subplots(len(seq), 2, figsize=(14, 3.2 * len(seq)), squeeze=False,
                            gridspec_kw={'width_ratios': [1.6, 1]})
    for row, fn in zip(axs, seq):
      d = load(fn)
      if P['structured']:
        zx, zy, ww = d['x'], d['y'], (d['x'][0], d['y'][0])
        row[0].pcolormesh(d['x'], d['y'], d['schlieren'], cmap='gray', shading='gouraud',
                          vmin=0, vmax=1)
      else:
        xy = pts / LREF
        row[0].tripcolor(xy[:, 0], xy[:, 1], tri, d['schlieren'], cmap='gray',
                         shading='gouraud', vmin=0, vmax=1)
      outline(row[0])
      row[0].set_xlim(P['view'][0], P['view'][1])
      row[0].set_ylim(P['view'][2], P['view'][3])
      row[0].set_aspect('equal')
      row[0].set_title(label(fn), fontsize=9)
      prof = cp_by_step.get(int(os.path.basename(fn)[5:11]), (None, None))[1]
      if prof:
        for nm, xv, cv in prof:
          row[1].plot(xv, -cv, label=nm)
      else:
        xr, cr = d['x'][0], d['Cp'][0]
        masks = ((np.ones(len(xr), bool), 'wall'),) if P['kind'] == 'angle' else (
            (xr > 0, 'upper'), (xr <= 0, 'lower'))
        for mask, nm in masks:
          o = np.argsort(xr[mask])
          row[1].plot(xr[mask][o], -cr[mask][o], label=nm)
      row[1].grid(alpha=0.3)
      row[1].set_ylabel('-Cp')
      if P['kind'] != 'angle':
        row[1].legend(fontsize=8)
    axs[-1][1].set_xlabel('$x/c$' if P['kind'] != 'angle' else '$x/D$')
    fig.tight_layout()
    fig.savefig(os.path.join(out, 'schlieren_sequence.png'), dpi=110)
    plt.close(fig)
    print('写出 schlieren_sequence.png')

# ============================ summary / README ============================
lines = lines + shock_lines + probe_lines
if div_step is not None:
  dt_s = cfg['solver']['urans']['dt'] if P['phase'] == 'U' else 1.0
  lines.insert(0, f'[!] 求解器在 step {div_step} 发散: 以下统计只覆盖前 {len(u)} 步'
                  f'({len(u) * dt_s:.4f} s), 不可作为定量抖振指标')
text = '\n'.join(lines)
open(os.path.join(out, 'summary.txt'), 'w', encoding='utf-8').write(text + '\n')
print(text)

if a.readme:
  sv, uv = cfg['solver'], cfg['solver'].get('urans', {})

  def sutherland(T, t0=273.15, ts=110.4, mu0=1.716e-05):
    return mu0 * (T / t0) ** 1.5 * (t0 + ts) / (T + ts)

  rho_inf = far['p'] / (R * far['T'])
  re_inf = rho_inf * U * CHORD / sutherland(far['T'])
  lab = P['note']
  if a.case == 'oat15a':
    lab = 'OAT15A 定常 RANS' if P['phase'] == 'S' else 'OAT15A 抖振 URANS'
  num = []
  if P['phase'] == 'U':
    num.append(f'物理时间步 dt = {uv["dt"]:g} s（每抖振周期约 {1 / f0 / uv["dt"]:.0f} 步）')
    num.append(f'物理步数 {sv["max_steps"]}，t_end = {uv["dt"] * sv["max_steps"]:.4f} s')
    num.append(f'内迭代最多 {uv["inner"]} 次（tol {uv["inner_tol"]:g}），伪时间 CFL = '
               f'{sv["cfl"]:g}，sweeps = {uv["sweeps"]}')
    num.append(f'定常预热 {uv["steady_iters"]} 步（CFL 2 → {uv["steady_cfl"]:g}）')
  else:
    num.append(f'定常迭代 {uv.get("steady_iters")} 步（CFL 2 → {uv.get("steady_cfl"):g}），'
               f'RANS 段步数 {sv["max_steps"]}（CFL = {sv["cfl"]:g}）')
  num.append(f'输出间隔：场 {sv["dump_interval"]} 步 / 历史 {sv["convergence_interval"]} 步'
             f' / 壁面 Cp {uv.get("wall_interval")} 步')
  rd = [
    f'# {a.case}  {lab}  alpha={far["alpha"]}°' + ('  （求解器发散）' if div_step else ''), '',
    f'- **算例**：{lab}，网格 {CEL} 单元 / {NF} 面（壁面 {len(wallf)} 面）',
    f'- **工况**：Ma = {far["Ma"]}，Re = ρUc/μ = {re_inf:.3e}，alpha = {far["alpha"]}°，'
    f'T = {far["T"]} K，p∞ = {far["p"]:.4f} Pa，U∞ = {U:.4f} m/s',
    f'- **数值设置**：' + '；'.join(num),
    f'- **网格**：`{os.path.join(os.path.dirname(a.config), cfg["io"]["mesh"])}`',
    '- **求解器版本**：PurineCFD-R2 2.3.1，二进制以 `-DPURINE_JAC_FULL` 编译'
    '（隐式求解用完整 4×4 块 Jacobi）',
    f'- **配置与运行**：`{a.config}`（其中 io.log / io.field 为 PLACEHOLDER，'
    f'运行时由脚本改写为指向本目录后执行 `OMP_NUM_THREADS=32 <二进制> <临时配置>`）',
    f'- **壁钟耗时**：{a.wall:.0f} s（进程总时间，含网格读取与初始化）'
    if a.wall == a.wall else '- **壁钟耗时**：-',
    '- **主要结果**：', '', '```', text, '```', '',
  ]
  if P['phase'] == 'U' and f0 == f0:
    t_u = u['time']
    t_lo, t_hi = t_u[len(t_u) // 2], t_u[-1]
    rd.append(f'- **统计窗口**：物理步的后 50%（t = {t_lo:.4f} – {t_hi:.4f} s，'
              f'约 {(t_hi - t_lo) * f0:.1f} 个抖振周期，T = {1 / f0 * 1000:.2f} ms）；'
              f'主频由该窗口 FFT 得到，分辨率 {df:.1f} Hz（补零 + 抛物线插值细化）')
  if div_step is not None:
    rd.append(f'- **警告**：本算例求解器在 step {div_step} 发散，本目录的图与统计只代表发散前的'
              f'前 {len(u)} 步，**不可作为定量结果引用**')
  if LREF != 1.0:
    rd.append(f'- **注意**：力系数参考长度硬编码为 1，输出需乘 1/c = {1 / LREF:.4f}'
              f'（弦长 {CHORD:.6f} m，上表已换算）')
  rd.append(f'- **后处理**：`python tests/post.py {a.config} --case {a.case} --out <本目录>`'
            f'（配置里的 io.field 为 PLACEHOLDER，必须用 --out 指定本目录），输出 '
            f'`history.png`（力系数）`residual.png`（残差）`cp_dist.png`（壁面 Cp）'
            f'`cp_rms.png` 与 `cp_rms.dat`（脉动 RMS）`schlieren_sequence.png`（纹影序列）'
            f'与 9 张云图面板；`cp_probe_t.png`（该点 Cp(t)）与 `psd.png`（对数-对数 PSD，'
            f'含 ±Δf 不确定带与 St 副轴）、`cp_probe_pt_psd.png`、`psd_x{a.probe:g}.dat` '
            f'为 x/c={a.probe:g} 壁面探针（用 `--probe` 改位置）')
  open(os.path.join(out, 'README.md'), 'w', encoding='utf-8').write('\n'.join(rd) + '\n')
  print('写出 README.md')
