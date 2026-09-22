import argparse
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.patches import Polygon

BASE = Path(__file__).resolve().parent

parser = argparse.ArgumentParser(description='画 O 型网格图')
parser.add_argument('mesh', nargs='?', type=Path, default=BASE / 'oat15a_sa.txt')
parser.add_argument('--prefix', default=None, help='输出文件名前缀, 默认取网格名')
args = parser.parse_args()
prefix = args.prefix or args.mesh.stem


def read_mesh(path):
  lines = path.read_text().splitlines()
  header = lines[0].split()
  nnode, ncell = int(header[0]), int(header[2])
  start = lines.index('(node)') + 1
  end = lines.index('(end)', start)
  nodes = np.array([[float(t) for t in line.split()[1:3]] for line in lines[start:end]])
  ncir = nnode - ncell
  return nodes.reshape(nnode // ncir, ncir, 2)


def draw(ax, grid, window, linewidth=0.35, wall_lw=0.9):
  nrad, ncir, _ = grid.shape
  segments = [np.vstack([grid[i], grid[i, :1]]) for i in range(nrad)]
  segments += [grid[:, j] for j in range(ncir)]
  ax.add_collection(LineCollection(segments, linewidths=linewidth, colors='0.3'))
  ax.add_patch(
    Polygon(
      grid[0],
      closed=True,
      facecolor='white',
      edgecolor='0.1',
      linewidth=wall_lw,
      zorder=3,
    )
  )
  ax.set_aspect('equal')
  ax.set_xlim(window[0], window[1])
  ax.set_ylim(window[2], window[3])
  ax.tick_params(which='both', direction='in', top=True, right=True, labelsize=9)
  ax.tick_params(which='major', length=5)
  ax.tick_params(which='minor', length=3)
  ax.xaxis.set_minor_locator(plt.matplotlib.ticker.AutoMinorLocator(5))
  ax.yaxis.set_minor_locator(plt.matplotlib.ticker.AutoMinorLocator(5))
  for spine in ax.spines.values():
    spine.set_linewidth(1.0)


grid = read_mesh(args.mesh)
nrad, ncir, _ = grid.shape
print(f'{args.mesh}: {nrad} x {ncir} = {(nrad - 1) * ncir} 单元')

figure, axes = plt.subplots(1, 2, figsize=(13.5, 6.2))
draw(axes[0], grid, (-10.6, 10.6, -10.6, 10.6), linewidth=0.25)
axes[0].set_title('whole domain (far radius 10)', fontsize=12)
draw(axes[1], grid, (-1.25, 1.65, -1.0, 1.0), linewidth=0.35)
axes[1].set_title('around the airfoil', fontsize=12)
figure.tight_layout()
figure.savefig(BASE / f'{prefix}_overview.png', dpi=150)
plt.close(figure)

figure, axes = plt.subplots(1, 3, figsize=(16.5, 5.4))
draw(axes[0], grid, (-0.58, -0.40, -0.10, 0.10), linewidth=0.8)
axes[0].set_title('leading edge', fontsize=12)
draw(axes[1], grid, (-0.14, 0.16, -0.01, 0.22), linewidth=0.5)
axes[1].set_title(r'upper-surface shock region ($\Delta x/c=0.0015$)', fontsize=11)
draw(axes[2], grid, (0.36, 0.56, -0.10, 0.10), linewidth=0.7)
axes[2].set_title('trailing edge', fontsize=12)
figure.tight_layout()
figure.savefig(BASE / f'{prefix}_detail.png', dpi=150)
plt.close(figure)

# 壁面第一层高度与弦向间距
chord = grid[0, :, 0].max() - grid[0, :, 0].min()
height = np.hypot(*(grid[1] - grid[0]).T)
width = np.hypot(*(np.roll(grid[0], -1, axis=0) - grid[0]).T) / chord
xc = (grid[0, :, 0] - grid[0, :, 0].min()) / chord
print(f'  壁面首层高度 {height.min():.3e} ~ {height.max():.3e} (弦长单位)')
print(f'  周向间距 {width.min():.5f} ~ {width.max():.5f}')
upper = grid[0, :, 1] > 0.0
fine = width < 0.0016
print(f'  Δx<0.0016 的 x/c 区间 {xc[fine].min():.4f} ~ {xc[fine].max():.4f}, '
      f'{int(fine.sum())} 个点 (上表面 {int((fine & upper).sum())}, 下表面 {int((fine & ~upper).sum())})')
print(f'  上表面激波加密段 Δx: {width[fine & upper].min():.5f} ~ {width[fine & upper].max():.5f}')
print(f'  输出 {prefix}_overview.png, {prefix}_detail.png')
