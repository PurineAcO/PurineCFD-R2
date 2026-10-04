"""URANS 后处理: 力系数历程与主频、翼型上表面激波位置、马赫数/Cp/数值纹影云图.

用法:
  python tests/post_urans.py config.json
  python tests/post_urans.py cases/naca0012/urans_m07a5.json --airfoil
需要 numpy 与 matplotlib. 云图要求结构化网格(io.structured), 单元按 (n-1)*S_MAX+s 排列.
"""

import argparse
import glob
import json
import math
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

GAMMA, R = 1.4, 287.05

parser = argparse.ArgumentParser()
parser.add_argument('config')
parser.add_argument(
  '--airfoil', action='store_true', help='翼型: 统计上表面激波位置(前缘在 x=-0.5)'
)
args = parser.parse_args()

base = os.path.dirname(os.path.abspath(args.config))
cfg = json.load(open(args.config))
out = os.path.normpath(os.path.join(base, cfg['io']['field']))
far = cfg['farfield']
U = far['Ma'] * math.sqrt(GAMMA * R * far['T'])

# ---------------- 力系数历程 ----------------
h = np.genfromtxt(
  os.path.join(out, 'history.csv'), delimiter=',', names=True, dtype=None, encoding=None
)
u = h[h['phase'] == 'U']
t, cl, cd = u['time'], u['Cl'], u['Cd']
half = len(t) // 2
s = cl[half:] - cl[half:].mean()
freq = np.fft.rfftfreq(len(s), t[1] - t[0])
amp = np.abs(np.fft.rfft(s * np.hanning(len(s))))
f0 = freq[np.argmax(amp[1:]) + 1]
lines = [
  f'物理步数 {len(t)}, t_end = {t[-1]:.4f} s, 平均内迭代 {u["inner"].mean():.1f} 次, '
  f'内迭代残差下降中位数 {np.median(u["res"] / u["res0"]):.2e}',
  f'后半段 Cl: 均值 {cl[half:].mean():.4f}, 范围 [{cl[half:].min():.4f}, {cl[half:].max():.4f}]',
  f'后半段 Cd: 均值 {cd[half:].mean():.4f}, 范围 [{cd[half:].min():.4f}, {cd[half:].max():.4f}]',
  f'Cl 主频 {f0:.3f} Hz, St = fL/U = {f0 / U:.4f}, k = pi f L/U = {math.pi * f0 / U:.4f}',
]

# ---------------- 激波位置: 上表面 dCp/dx 最大处 ----------------
shock_t, shock_x = [], []
if args.airfoil:
  for block in open(os.path.join(out, 'wall_cp.dat')).read().split('# step ')[1:]:
    rows = block.strip().split('\n')
    a = np.array([r.split() for r in rows[1:]], dtype=float)
    up = a[a[:, 1] > 0]
    up = up[np.argsort(up[:, 0])]
    x, cp = up[:, 0], up[:, 2]
    g = np.where((x > -0.45) & (x < 0.35), np.gradient(cp, x), -np.inf)
    shock_t.append(float(rows[0].split()[2]))
    shock_x.append(x[np.argmax(g)] + 0.5)
  sx = np.array(shock_x)[len(shock_x) // 2 :]
  lines.append(
    f'后半段上表面激波位置 x/c: 均值 {sx.mean():.3f}, 范围 [{sx.min():.3f}, {sx.max():.3f}]'
  )

text = '\n'.join(lines)
print(text)
open(os.path.join(out, 'summary.txt'), 'w', encoding='utf-8').write(text + '\n')

n = 3 if args.airfoil else 2
fig, ax = plt.subplots(n, 1, figsize=(9, 2.6 * n), sharex=True)
ax[0].plot(t, cl, lw=1)
ax[0].set_ylabel('Cl')
ax[1].plot(t, cd, lw=1, color='C1')
ax[1].set_ylabel('Cd')
if args.airfoil:
  ax[2].plot(shock_t, shock_x, lw=1, color='C2')
  ax[2].set_ylabel('shock x/c (upper)')
ax[-1].set_xlabel('t [s]')
fig.tight_layout()
fig.savefig(os.path.join(out, 'history.png'), dpi=130)
plt.close(fig)

# ---------------- 云图 ----------------
if 'structured' not in cfg['io']:
  raise SystemExit
S = int(open(os.path.join(base, cfg['io']['structured'])).readline().split()[0])
files = sorted(glob.glob(os.path.join(out, 'step_*.dat')))
times = dict(zip(u['step'], zip(u['time'], u['Cl'])))
p_inf, q_inf = far['p'], 0.5 * GAMMA * far['p'] * far['Ma'] ** 2


def load(fn):
  f = np.loadtxt(fn, skiprows=2)
  # 行为径向层(第0行贴壁), 列为周向
  a = {k: f[:, i].reshape(-1, S) for i, k in enumerate(['x', 'y', 'rho', 'u', 'v', 'T', 'p', 'Ma', 'miubl'])}
  a['Cp'] = (a['p'] - p_inf) / q_inf

  # 周向(列)用周期中心差分, 径向(行)用 np.gradient
  def ds(v):
    return 0.5 * (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1))

  def dn(v):
    return np.gradient(v, axis=0)

  # 曲线坐标链式法则求 |∇ρ|
  r, x, y = a['rho'], a['x'], a['y']
  jac = ds(x) * dn(y) - dn(x) * ds(y)
  gx = (ds(r) * dn(y) - dn(r) * ds(y)) / jac
  gy = (dn(r) * ds(x) - ds(r) * dn(x)) / jac
  grad = np.hypot(gx, gy)
  a['schlieren'] = np.exp(-15 * grad / np.percentile(grad, 99.5))
  vx = (ds(a['v']) * dn(y) - dn(a['v']) * ds(y)) / jac
  uy = (dn(a['u']) * ds(x) - ds(a['u']) * dn(x)) / jac
  a['vorticity'] = (vx - uy) / (far['Ma'] * np.sqrt(GAMMA * 287.05 * far['T']))
  # 末尾补第一列使周向闭合
  return {k: np.hstack([v, v[:, :1]]) for k, v in a.items()}


def label(fn):
  step = int(os.path.basename(fn)[5:11])
  if step in times:
    return f'step {step}  t={times[step][0]:.4f}s  Cl={times[step][1]:.4f}'
  return f'step {step}'


def body(ax, a, view):
  ax.fill(a['x'][0], a['y'][0], color='w', edgecolor='k', linewidth=0.8, zorder=5)
  ax.set_xlim(view[0], view[1])
  ax.set_ylim(view[2], view[3])
  ax.set_aspect('equal')


view = (-0.7, 0.8, -0.35, 0.55) if args.airfoil else (-5.0, 15.0, -10.0, 10.0)
cp_star = (
  2
  / (GAMMA * far['Ma'] ** 2)
  * (((2 + (GAMMA - 1) * far['Ma'] ** 2) / (GAMMA + 1)) ** (GAMMA / (GAMMA - 1)) - 1)
)
a = load(files[-1])


def panel(z, name, levels, cmap, extend, bar):
  fig, ax = plt.subplots(figsize=(8.2, 7.0))
  c = ax.contourf(a['x'], a['y'], z, levels=levels, cmap=cmap, extend=extend, antialiased=False)
  if cmap != 'gray':
    ax.contour(a['x'], a['y'], z, levels=levels[1::2], colors='k', linewidths=0.45, alpha=0.7)
  body(ax, a, view)
  ax.set_xlabel('$x/c$' if args.airfoil else '$x/D$', fontsize=14)
  ax.set_ylabel('$y/c$' if args.airfoil else '$y/D$', fontsize=14)
  ax.tick_params(direction='in', top=True, right=True, labelsize=11)
  ax.set_title(label(files[-1]), fontsize=12)
  cb = fig.colorbar(c, ax=ax, orientation='horizontal', fraction=0.05, pad=0.12)
  cb.ax.tick_params(labelsize=10)
  cb.ax.set_title(bar, fontsize=11, pad=8)
  fig.tight_layout()
  fig.savefig(os.path.join(out, name), dpi=150)
  plt.close(fig)
  print('写出', os.path.join(out, name))


U_inf = far['Ma'] * np.sqrt(GAMMA * 287.05 * far['T'])
specs = [
    ('u', 'RdBu_r', -1.2 * U_inf, 1.2 * U_inf, r'$u$ [m/s]'),
    ('v', 'RdBu_r', -0.6 * U_inf, 0.6 * U_inf, r'$v$ [m/s]'),
    ('rho', 'jet', None, None, r'$\rho$ [kg/m$^3$]'),
    ('miubl', 'jet', 0.0, None, r'$\tilde\nu$ [m$^2$/s]'),
    ('T', 'jet', None, None, r'$T$ [K]'),
    ('Ma', 'jet', None, None, 'Ma'),
    ('Cp', 'jet', -1.6, 1.1, 'Cp'),
    ('schlieren', 'gray', 0.0, 1.0, 'schlieren'),
    ('vorticity', 'jet', -4.0, 4.0, r'$\omega D/U$'),
]
for key, cmap, lo, hi, bar in specs:
    z = a[key]
    if lo is None:
        lo = float(np.floor(z.min() * 20) / 20)
    if hi is None:
        hi = float(np.ceil(z.max() * 20) / 20)
    panel(z, ('T' if key == 'T' else key.lower()) + '.png', np.linspace(lo, hi, 21), cmap,
          'neither' if cmap == 'gray' else 'both', bar)

# 最后 6 帧的纹影与壁面 -Cp
seq = files[-6:]
fig, axs = plt.subplots(
  len(seq), 2, figsize=(14, 3.2 * len(seq)), squeeze=False, gridspec_kw={'width_ratios': [1.6, 1]}
)
for row, fn in zip(axs, seq):
  a = load(fn)
  row[0].pcolormesh(a['x'], a['y'], a['schlieren'], cmap='gray', shading='gouraud', vmin=0, vmax=1)
  body(row[0], a, view)
  row[0].set_title(label(fn), fontsize=9)
  xw, yw, cw = a['x'][0], a['y'][0], a['Cp'][0]
  for side, name in [(yw > 0, 'upper'), (yw <= 0, 'lower')]:
    o = np.argsort(xw[side])
    row[1].plot(xw[side][o], -cw[side][o], label=name)
  row[1].grid(alpha=0.3)
  row[1].set_ylabel('-Cp (1st cell)')
axs[0][1].legend()
axs[-1][1].set_xlabel('x')
fig.tight_layout()
fig.savefig(os.path.join(out, 'schlieren_sequence.png'), dpi=110)
