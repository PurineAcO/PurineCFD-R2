"""NACA0012 后处理: 门类与风格对齐 run/naca0012_ma076_re1e7_c5_a32.

输出 (都在 debugres/ 下): mach.png  cp.png  flow.png  profile.png  convergence.png
用法: python3 cases/naca0012/plot_all.py [debugres]
"""

import json
import math
import sys
from pathlib import Path

import matplotlib
import matplotlib.colors as mcolors
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.colors import BoundaryNorm

BASE = Path(__file__).resolve().parent
OUT = Path(sys.argv[1]) if len(sys.argv) > 1 else BASE / 'debugres'
MESH = BASE / 'naca0012_sa.txt'
NCIR, NRAD = 326, 90
RE = 3.3e6
GAMMA, RGAS, CV = 1.4, 287.05, 717.645
MU0, T0, TS = 1.716e-5, 273.15, 110.4
CFG = json.loads((BASE / 'config.json').read_text())['farfield']
T_INF, MA, ALPHA = CFG['T'], CFG['Ma'], math.radians(CFG['alpha'])
A_INF = math.sqrt(GAMMA * RGAS * T_INF)
U_INF = MA * A_INF
MU_INF = MU0 * (T_INF / T0) ** 1.5 * (T0 + TS) / (T_INF + TS)
RHO_INF = RE * MU_INF / (U_INF * 1.0)
P_INF = RHO_INF * RGAS * T_INF
Q_INF = 0.5 * RHO_INF * U_INF**2
X_REF = 0.25


def sutherland(T):
  return MU0 * (T / T0) ** 1.5 * (T0 + TS) / (T + TS)


def read_nodes(path):
  lines = path.read_text().splitlines()
  start = lines.index('(node)') + 1
  end = lines.index('(end)', start)
  points = [line.split() for line in lines[start:end]]
  return np.array([[float(p[1]), float(p[2])] for p in points]).reshape(NRAD, NCIR, 2)


def latest_field():
  def key(p):
    tag = p.stem[len('field_') :]
    return (0, int(tag)) if tag.isdigit() else (1, 0)

  files = sorted(OUT.glob('field_*.dat'), key=key)
  if not files:
    raise SystemExit(f'{OUT} 下没有 field_*.dat')
  return files[-1]


def read_field(path):
  lines = [line for line in path.read_text().splitlines() if line.strip()]
  return np.array([line.split() for line in lines[2:]], dtype=float)


def read_history(path):
  steps, values = [], []
  if not path.exists():
    return np.array([]), np.array([])
  for line in path.read_text().splitlines():
    parts = line.split(',')
    if len(parts) < 2 or not parts[0].strip().isdigit():
      continue
    steps.append(float(parts[0]))
    values.append(float(parts[1]))
  return np.array(steps), np.array(values)


def wall_frame(nodes, inner):
  wall = nodes[0]
  edge = np.roll(wall, -1, axis=0) - wall
  length = np.linalg.norm(edge, axis=1)
  tangent = edge / length[:, None]
  normal = np.column_stack((tangent[:, 1], -tangent[:, 0]))
  mid = 0.5 * (wall + np.roll(wall, -1, axis=0))
  flip = np.where(np.sum((inner - mid) * normal, axis=1) < 0.0, -1.0, 1.0)
  return mid, length, tangent, normal * flip[:, None]


def wall_distance(centers, nodes):
  wall = nodes[0]
  edge = np.roll(wall, -1, axis=0) - wall
  length = np.linalg.norm(edge, axis=1)
  delta = centers[:, None, :] - wall[None, :, :]
  t = np.clip(np.sum(delta * edge[None], axis=2) / (length[None] ** 2), 0.0, 1.0)
  closest = wall[None] + t[:, :, None] * edge[None]
  return np.linalg.norm(centers[:, None, :] - closest, axis=2).min(axis=1)


def rainbow(name, levels, low=0.18, high=0.55):
  base = plt.get_cmap(name)
  return mcolors.ListedColormap(base(np.linspace(low, high, levels)), name=f'{name}-band')


def frame(ax, window):
  ax.set_aspect('equal')
  ax.set_xlim(window[0], window[1])
  ax.set_ylim(window[2], window[3])
  ax.tick_params(which='both', direction='in', top=True, right=True)
  ax.tick_params(which='major', length=7, width=1.0, labelsize=13)
  ax.tick_params(which='minor', length=4, width=0.8)
  ax.xaxis.set_minor_locator(ticker.AutoMinorLocator(5))
  ax.yaxis.set_minor_locator(ticker.AutoMinorLocator(5))
  for spine in ax.spines.values():
    spine.set_linewidth(1.2)


def draw(ax, x, y, values, levels, cmap):
  plt.rcParams['contour.negative_linestyle'] = 'solid'
  norm = BoundaryNorm(levels, cmap.N)
  clipped = np.clip(values, levels[0], levels[-1])
  filled = ax.contourf(x, y, clipped, levels=levels, cmap=cmap, norm=norm, extend='neither')
  ax.contour(x, y, clipped, levels=levels, colors='k', linewidths=0.6, linestyles='solid')
  bar = plt.colorbar(
    filled, ax=ax, ticks=levels, fraction=0.055, pad=0.03, drawedges=False, format='%g'
  )
  bar.ax.tick_params(length=0, labelsize=12)
  bar.outline.set_linewidth(1.0)
  return bar


def main():
  path = latest_field()
  field = read_field(path)
  # 按 dump 里自带的 (s,n) 重排成 n 主序, 不假设网格文件的单元排列方式
  field = field[np.lexsort((field[:, 10], field[:, 11]))]
  nodes = read_nodes(MESH)
  chord = nodes[0, :, 0].max() - nodes[0, :, 0].min()
  x_le = nodes[0, :, 0].min()
  x, y, rho, u, v, T, p = (field[:, k] for k in range(7))
  ma = np.hypot(u, v) / np.sqrt(GAMMA * RGAS * T)
  cp = (p - P_INF) / Q_INF
  xc = (nodes[0, :, 0] - x_le) / chord
  upper = nodes[0, :, 1] > 0.0

  mid, ds, tangent, normal = wall_frame(nodes, field[:NCIR, :2])
  flow = field[:NCIR, 3:5]
  flip = np.where(np.sum(tangent * flow, axis=1) < 0.0, -1.0, 1.0)
  tangent = tangent * flip[:, None]
  distance = wall_distance(field[:NCIR, :2], nodes)
  mu = sutherland(field[:NCIR, 5])
  tangential = np.abs(field[:NCIR, 3] * tangent[:, 0] + field[:NCIR, 4] * tangent[:, 1])
  tau = mu * tangential / distance
  utau = np.sqrt(tau / field[:NCIR, 2])
  yplus = field[:NCIR, 2] * utau * distance / mu
  cf = 2.0 * tau / Q_INF
  cp_wall = 0.5 * (cp[:NCIR] + np.roll(cp[:NCIR], 1))

  dcp = -cp_wall[:, None] * normal * (ds / chord)[:, None]
  dcf = 0.5 * cf[:, None] * tangent * (ds / chord)[:, None]
  dc = dcp + dcf
  ca_p, cn_p = dcp.sum(axis=0)
  ca_f, cn_f = dcf.sum(axis=0)
  cl_p = cn_p * math.cos(ALPHA) - ca_p * math.sin(ALPHA)
  cd_p = ca_p * math.cos(ALPHA) + cn_p * math.sin(ALPHA)
  cl_f = cn_f * math.cos(ALPHA) - ca_f * math.sin(ALPHA)
  cd_f = ca_f * math.cos(ALPHA) + cn_f * math.sin(ALPHA)
  rx = (mid[:, 0] - (x_le + X_REF * chord)) / chord
  ry = mid[:, 1] / chord
  cm = float(np.sum(rx * dc[:, 1] - ry * dc[:, 0]))

  print(f'文件        {path.name}')
  print(f'工况        Ma={MA:.2f}  Re={RE:.3g}  alpha={math.degrees(ALPHA):.2f} deg')
  print(f'网格        {NRAD - 1} x {NCIR} = {(NRAD - 1) * NCIR} 单元, 弦长 {chord:.6f}')
  print()
  print(f'Ma 范围     {ma.min():.6f} ~ {ma.max():.6f}   超声速单元 {int((ma > 1).sum())}')
  hot = int(ma.argmax())
  print(f'Ma 最大位置 x/c={(x[hot] - x_le) / chord:.4f}, y={y[hot]:.4f}')
  grid = ma.reshape(NRAD - 1, NCIR)
  sup = grid > 1.0
  js = np.where(sup.any(axis=0))[0]
  rows = np.where(sup.any(axis=1))[0]
  if len(js) and len(rows):
    print(
      f'超声速区 x/c = {xc[js].min():.4f} ~ {xc[js].max():.4f}, 径向层 i = {rows.min()} ~ {rows.max()}'
    )
  drop = np.diff(cp_wall)
  print(f'上表面 Cp 最陡下降(x/c) = {xc[np.argmin(drop)]:.4f}')
  print(f'密度范围 {rho.min():.6f} ~ {rho.max():.6f}, 温度 {T.min():.3f} ~ {T.max():.3f} K')
  print(f'nu_tilde {field[:, 8].min():.6e} ~ {field[:, 8].max():.6e}')
  print()
  print('--- 壁面 ---')
  print(f'首层离壁 min {distance.min():.3e}  max {distance.max():.3e}')
  print(f'u_tau min {utau.min():.3f}  max {utau.max():.3f}  mean {utau.mean():.3f} m/s')
  print(
    f'y+    min {yplus.min():.3f}  max {yplus.max():.3f}  mean {yplus.mean():.3f}'
    f'  中位 {np.median(yplus):.3f}'
  )
  print(f'Cf    min {cf.min():.5f}  max {cf.max():.5f}')
  print()
  print(f'--- 力系数 (力矩参考点 x/c={X_REF:.2f}) ---')
  print(f'  Cl     = {cl_p + cl_f:+.5f}    (压差 {cl_p:+.5f}, 摩擦 {cl_f:+.5f})')
  print(f'  Cd     = {cd_p + cd_f:+.5f}    (压差 {cd_p:+.5f}, 摩擦 {cd_f:+.5f})')
  print(f'  Cm_c/4 = {cm:+.5f}')
  print(f'  L/D    = {(cl_p + cl_f) / (cd_p + cd_f):.2f}')

  # (0,0)(0,1) 是分块网格, 这里 reshape 成环形网格
  # 网格里 LE 在 x=x_le, 统一归一化为 x/c
  def norm_x(v):
    return (v - x_le) / chord

  wrap_x = norm_x(np.hstack((nodes[:, :, 0], nodes[:, :1, 0])))
  wrap_y = np.hstack((nodes[:, :, 1], nodes[:, :1, 1])) / chord
  cell_x = norm_x(x.reshape(NRAD - 1, NCIR))
  cell_y = y.reshape(NRAD - 1, NCIR) / chord
  wrap_cx = np.hstack((cell_x, cell_x[:, :1]))
  wrap_cy = np.hstack((cell_y, cell_y[:, :1]))
  outline_x = norm_x(np.append(nodes[0, :, 0], nodes[0, 0, 0]))
  outline_y = np.append(nodes[0, :, 1], nodes[0, 0, 1]) / chord

  def wrapped(values):
    g = values.reshape(NRAD - 1, NCIR)
    return np.hstack((g, g[:, :1]))

  # ---- mach.png ----
  figure, ax = plt.subplots(figsize=(8.4, 7.8))
  levels = np.round(np.arange(0.0, 1.41, 0.1), 2)
  bar = draw(ax, wrap_cx, wrap_cy, wrapped(ma), levels, plt.get_cmap('jet_r', len(levels) - 1))
  bar.ax.set_title('Mach', fontsize=14, pad=8)
  ax.set_xlabel('$x/c$', fontsize=16)
  ax.set_ylabel('$y/c$', fontsize=16)
  frame(ax, (-0.5, 1.5, -1.0, 1.0))
  plt.tight_layout()
  plt.savefig(OUT / 'mach.png', dpi=140)
  plt.close()

  # ---- cp.png ----
  figure, ax = plt.subplots(figsize=(8.4, 7.8))
  levels = np.round(np.arange(-1.4, 0.01, 0.1), 2)
  bar = draw(ax, wrap_cx, wrap_cy, wrapped(cp), levels, rainbow('jet', len(levels) - 1))
  bar.ax.set_title('Cp', fontsize=14, pad=8)
  ax.set_xlabel('$x/c$', fontsize=16)
  ax.set_ylabel('$y/c$', fontsize=16)
  frame(ax, (0.35, 0.55, 0.05, 0.40))
  plt.tight_layout()
  plt.savefig(OUT / 'cp.png', dpi=140)
  plt.close()

  # ---- flow.png ----
  figure, axes = plt.subplots(2, 2, figsize=(13.5, 10))
  for ax, window in zip(axes[0], [(-1.5, 2.5, -2.0, 2.0), (-0.3, 1.3, -0.8, 0.8)]):
    mesh = ax.pcolormesh(wrap_x, wrap_y, grid, cmap='turbo', shading='flat')
    ax.contour(cell_x, cell_y, grid, levels=[1.0], colors='k', linewidths=0.8)
    ax.plot(outline_x, outline_y, 'k-', lw=1.2)
    ax.set_aspect('equal')
    ax.set_xlim(window[0], window[1])
    ax.set_ylim(window[2], window[3])
    figure.colorbar(mesh, ax=ax, shrink=0.9)
  axes[0, 0].set_title(f'Mach number (max={ma.max():.4f}), black line: Ma=1')
  axes[0, 1].set_title('Mach number near the airfoil')

  axes[1, 0].plot(xc[upper], -cp_wall[upper], 'o-', ms=2.5, lw=1.0, label='upper')
  axes[1, 0].plot(xc[~upper], -cp_wall[~upper], 's-', ms=2.5, lw=1.0, label='lower')
  axes[1, 0].axhline(0.0, color='k', lw=0.6)
  axes[1, 0].set_xlabel('x/c')
  axes[1, 0].set_ylabel('-Cp')
  axes[1, 0].set_title('pressure distribution')
  axes[1, 0].invert_yaxis()
  axes[1, 0].grid(True, alpha=0.3)
  axes[1, 0].legend()

  axes[1, 1].plot(xc[upper], cf[upper] * 1e3, 'o-', ms=2.5, lw=1.0, label='upper')
  axes[1, 1].plot(xc[~upper], cf[~upper] * 1e3, 's-', ms=2.5, lw=1.0, label='lower')
  axes[1, 1].axhline(0.0, color='k', lw=0.6)
  axes[1, 1].set_xlabel('x/c')
  axes[1, 1].set_ylabel('Cf x 1e3')
  axes[1, 1].set_title('skin friction')
  axes[1, 1].grid(True, alpha=0.3)
  axes[1, 1].legend()
  plt.tight_layout()
  plt.savefig(OUT / 'flow.png', dpi=140)
  plt.close()

  # ---- profile.png ----
  figure, axes = plt.subplots(1, 2, figsize=(13, 5))
  axes[0].plot(xc, yplus, 'o-', ms=2.5, lw=1.0)
  axes[0].axhline(1.0, color='crimson', ls='--', lw=0.8, label='y+ = 1')
  axes[0].axhline(5.0, color='orange', ls=':', lw=0.8, label='y+ = 5')
  axes[0].set_xlabel('x/c')
  axes[0].set_ylabel('y+')
  axes[0].set_yscale('log')
  axes[0].set_title('y+ at the first cell')
  axes[0].grid(True, which='both', alpha=0.3)
  axes[0].legend()

  mask = np.abs(xc - 0.30) + np.where(nodes[0, :, 1] > 0.0, 0.0, 10.0)
  j = int(mask.argmin())
  for row in (0, 1, 2, 3, 5, 8, 12, 16, 22, 30):
    k = row * NCIR
    dr = wall_distance(field[k : k + NCIR, :2], nodes)[j]
    ut = field[k + j, 3] * tangent[j, 0] + field[k + j, 4] * tangent[j, 1]
    yp = field[k + j, 2] * abs(utau[j]) * dr / sutherland(field[k + j, 5])
    axes[1].plot(yp, abs(ut) / abs(utau[j]), 'o', ms=5, label=f'row {row}')
  law = np.logspace(0, 2.7, 100)
  axes[1].plot(law, law, 'k--', lw=0.8, label='u+ = y+')
  axes[1].set_xscale('log')
  axes[1].set_yscale('log')
  axes[1].set_xlabel('y+')
  axes[1].set_ylabel('u+')
  axes[1].set_xlim(0.1, 300)
  axes[1].set_title('velocity profile at x/c = 0.30 (upper)')
  axes[1].grid(True, which='both', alpha=0.3)
  axes[1].legend(fontsize=8, ncol=2)
  plt.tight_layout()
  plt.savefig(OUT / 'profile.png', dpi=140)
  plt.close()

  # ---- convergence.png ----
  steps, values = read_history(OUT / 'residual_run.csv')
  if len(steps):
    order = np.argsort(steps)
    plt.figure(figsize=(7, 4.5))
    plt.semilogy(steps[order], values[order], lw=1.2)
    plt.axhline(1e-6, color='crimson', ls='--', lw=0.8, label='convergence limit 1e-6')
    plt.xlabel('step')
    plt.ylabel(r'$L_2(\Delta U)$')
    plt.title(f'NACA0012  Ma={MA}  Re={RE:.2g}  alpha={math.degrees(ALPHA):.1f}deg')
    plt.grid(True, which='both', alpha=0.3)
    plt.legend()
    plt.tight_layout()
    plt.savefig(OUT / 'convergence.png', dpi=140)
    plt.close()
  print(f'\n图已写出 {OUT}/mach.png  cp.png  flow.png  profile.png  convergence.png')


if __name__ == '__main__':
  main()
