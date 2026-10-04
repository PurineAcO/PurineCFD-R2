#!/usr/bin/env python3
"""只读已保存的流场文件出图 + 统计, 不重跑求解器.

用法:
    python3 tests/plot_field.py debugres                  # 该目录下所有 field_*.dat
    python3 tests/plot_field.py debugres/field_final.dat  # 只画一帧
"""

import glob
import math
import os
import sys

import numpy as np


def load_field(path):
  var = None
  rows = []
  with open(path, encoding='utf-8', errors='replace') as fp:
    for line in fp:
      if line.startswith('VARIABLES'):
        var = [v.strip().strip('"') for v in line.split('=', 1)[1].split(',')]
      elif line[0].isdigit() or line[0] == '-':
        rows.append([float(v) for v in line.split()])
  return var, np.atleast_2d(np.array(rows))


def build_tri(rows):
  """有结构化 (s,n) 就用结构化三角网, 否则退回 Delaunay(适合 H 型多块网格)."""
  import matplotlib.tri as mtri

  x, y = rows[:, 0], rows[:, 1]
  s_idx, n_idx = rows[:, 10].astype(int), rows[:, 11].astype(int)
  smax, nmax = int(s_idx.max()), int(n_idx.max())
  if smax <= 0 or nmax <= 0:
    return mtri.Triangulation(x, y), 0, 0
  lut = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(s_idx, n_idx))}
  tris = []
  for b in range(1, nmax):
    for a in range(1, smax + 1):
      a2 = a % smax + 1
      q = [lut.get((a, b)), lut.get((a2, b)), lut.get((a2, b + 1)), lut.get((a, b + 1))]
      if any(v is None for v in q):
        continue
      tris.append([q[0], q[1], q[2]])
      tris.append([q[0], q[2], q[3]])
  return mtri.Triangulation(x, y, np.array(tris)), smax, nmax


def report(rows, path):
  u, ma, cp, nu = rows[:, 3], rows[:, 7], rows[:, 9], rows[:, 8]
  wall = rows[rows[:, 11] == 1]
  tag = os.path.basename(path)[len('field_') : -len('.dat')]
  print(
    f'  {tag:>7}  Ma[{ma.min():.4f},{ma.max():.4f}]  Cp[{cp.min():+.3f},{cp.max():+.3f}]  '
    f'u_min={u.min():+8.4f}  回流格={int((u < 0).sum()):5d}/{len(u)}  '
    f'壁面回流={int((wall[:, 3] < 0).sum()):3d}/{len(wall)}  max_nu={nu.max():.4f}'
  )


def plot(dat_path, ma_rng, cp_rng, view=(-5.0, 15.0, -10.0, 10.0), draw_circle=True):
  import matplotlib

  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  from matplotlib.patches import Circle

  var, rows = load_field(dat_path)
  tri, smax, nmax = build_tri(rows)
  ma, cp = rows[:, 7], rows[:, 9]
  title = open(dat_path, encoding='utf-8', errors='replace').readline()
  fig, ax = plt.subplots(1, 2, figsize=(13.6, 5.6))
  for a, val, lab, cm, rng in (
    (ax[0], ma, 'Mach number', 'turbo', ma_rng),
    (ax[1], cp, 'Pressure coefficient $C_p$', 'RdYlBu_r', cp_rng),
  ):
    lo, hi = rng
    lv = np.linspace(lo, hi, 31)
    cf = a.tricontourf(tri, val, levels=lv, cmap=cm, extend='both')
    a.tricontour(tri, val, levels=lv[1::3], colors='k', linewidths=0.25, alpha=0.5)
    a.set_xlim(view[0], view[1])
    a.set_ylim(view[2], view[3])
    a.set_aspect('equal')
    if draw_circle:
      a.add_patch(
        Circle((0.0, 0.0), 0.5, facecolor='white', edgecolor='k', linewidth=0.7, zorder=6)
      )
    a.set_xlabel('x/D', fontsize=14, fontweight='bold')
    a.set_ylabel('y/D', fontsize=14, fontweight='bold')
    a.set_title(lab, fontsize=14)
    a.tick_params(labelsize=11, direction='in', top=True, right=True)
    cb = fig.colorbar(cf, ax=a, fraction=0.046, pad=0.03)
    cb.set_label(f'fixed scale [{lo:.3g}, {hi:.3g}]', fontsize=9)
  fig.suptitle(f'{os.path.basename(dat_path)}   {title.strip()}', fontsize=12)
  fig.tight_layout()
  out = dat_path[:-4] + '.png'
  fig.savefig(out, dpi=130)
  plt.close(fig)
  print(f'           -> {out}')


def frame_key(path):
  """按迭代号排; final 永远排最后."""
  tag = os.path.basename(path)[len('field_') : -len('.dat')]
  try:
    return (0, int(tag))
  except ValueError:
    return (1, 0)


def frames_range(data, col, step, q=0.5):
  """所有帧该列的全局范围(用分位数剔掉瞬态尖峰), 向外取整到 step 的倍数."""
  lo = min(float(np.nanpercentile(d[:, col], q)) for d in data)
  hi = max(float(np.nanpercentile(d[:, col], 100.0 - q)) for d in data)
  return (math.floor(lo / step) * step, math.ceil(hi / step) * step)


def main():
  if len(sys.argv) < 2:
    print(__doc__)
    return 1
  target = sys.argv[1]
  if os.path.isdir(target):
    files = sorted(glob.glob(os.path.join(target, 'field_*.dat')), key=frame_key)
  else:
    files = [target]
  if not files:
    print(f'没有找到流场文件: {target}')
    return 1
  data = [load_field(f)[1] for f in files]
  view = tuple(float(v) for v in sys.argv[2:6]) if len(sys.argv) >= 6 else (-5.0, 15.0, -10.0, 10.0)
  draw_circle = len(sys.argv) < 6
  # 标尺用收敛帧(field_final)定, 早期瞬态帧超出部分直接裁掉, 这样跨 step 才是同一把尺
  ref = next((d for f, d in zip(files, data) if os.path.basename(f) == 'field_final.dat'), None)
  if ref is None:
    ref = data[-1]
  ma_rng = frames_range([ref], 7, 0.05)
  cp_rng = frames_range([ref], 9, 0.25)
  print(
    f'共 {len(files)} 帧, 统一色标(按收敛帧): Ma {ma_rng[0]:.3g}~{ma_rng[1]:.3g}   '
    f'Cp {cp_rng[0]:.3g}~{cp_rng[1]:.3g}'
  )
  for f, rows in zip(files, data):
    report(rows, f)
    try:
      plot(f, ma_rng, cp_rng, view, draw_circle)
    except ImportError:
      print('           (未安装 matplotlib, 只输出统计)')
  return 0


if __name__ == '__main__':
  sys.exit(main())
