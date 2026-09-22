import argparse
import math
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser(
  description='Convert a PLOT3D single-block 2D grid around an airfoil to the solver mesh format'
)
parser.add_argument(
  'source', nargs='?', type=Path, default=Path(__file__).resolve().parent / 'naca0012.txt'
)
parser.add_argument(
  'target', nargs='?', type=Path, default=Path(__file__).resolve().parent / 'naca0012_cfd.txt'
)
parser.add_argument('--first-height', type=float, default=2.0e-5, help='壁面首层高度, 弦长 1.0')
parser.add_argument('--growth', type=float, default=1.2, help='径向相邻层高度比')
parser.add_argument('--far-radius', type=float, default=30.0, help='远场半径, 弦长单位')
parser.add_argument('--max-cell', type=float, default=0.02, help='径向后段单元厚度上限')
parser.add_argument('--max-cell-height', type=float, default=0.45, help='厚度上限生效到多高')
parser.add_argument('--shock-xc', type=float, default=0.46, help='激波弦向位置')
parser.add_argument('--shock-dx', type=float, default=0.0015, help='激波处弦向目标间距')
parser.add_argument('--shock-width', type=float, default=0.13, help='激波加密窗口半宽')
args = parser.parse_args()
if args.first_height <= 0.0 or args.growth <= 1.0:
  raise SystemExit('--first-height 需为正, --growth 需大于 1')
if args.max_cell <= args.first_height or args.max_cell_height <= 0.0:
  raise SystemExit('--max-cell 需大于首层高度, --max-cell-height 需为正')
if args.shock_dx <= 0.0 or args.shock_width <= 0.0:
  raise SystemExit('--shock-dx 与 --shock-width 需为正')

tokens = args.source.read_text().split()
if int(tokens[0]) != 1:
  raise SystemExit('本脚本只处理单块 PLOT3D 网格')
ni, nj, nk = (int(v) for v in tokens[1:4])
if nk != 1:
  raise SystemExit(f'本脚本只处理二维网格, 实际 nk={nk}')
values = [float(v) for v in tokens[4:]]
if len(values) == 3 * ni * nj:
  values = values[: 2 * ni * nj]  # 丢掉全零的 Z 数组
if len(values) != 2 * ni * nj:
  raise SystemExit(f'数据量 {len(values)} 与 ni*nj*2 = {2 * ni * nj} 不符')

# PLOT3D 按 i 最快存放, 转成 [j][i] 后 i 是径向(体面->远场), j 是周向
grid = [[values[j * ni + i] for i in range(ni)] for j in range(nj)]
radial = [[(grid[j][i], values[ni * nj + j * ni + i]) for j in range(nj)] for i in range(ni)]
for i, line in enumerate(radial):
  if line[0] != line[-1]:
    raise SystemExit(f'第 {i} 条径向线的首尾点不重合, 无法按闭合环处理')

ncir = nj - 1


def stations(points):
  values = [0.0]
  for i in range(len(points) - 1):
    values.append(values[-1] + math.dist(points[i], points[i + 1]))
  return values


def sample(points, station, fraction):
  target = fraction * station[-1]
  for i in range(len(station) - 1):
    width = station[i + 1] - station[i]
    if target <= station[i + 1]:
      weight = 0.0 if width <= 0.0 else (target - station[i]) / width
      return (
        points[i][0] + weight * (points[i + 1][0] - points[i][0]),
        points[i][1] + weight * (points[i + 1][1] - points[i][1]),
      )
  return points[-1]


def span(first, growth, cells):
  return first * cells if growth == 1.0 else first * (growth**cells - 1) / (growth - 1)


def solve_growth(first, cells, total):
  lo, hi = 1.0, 2.0
  while span(first, hi, cells) < total:
    hi *= 1.5
  for _ in range(200):
    mid = 0.5 * (lo + hi)
    if span(first, mid, cells) < total:
      lo = mid
    else:
      hi = mid
  return 0.5 * (lo + hi)


def radial_offsets(total, first, growth, cap, cap_height):
  values = [0.0]
  step = first
  while step < cap:
    values.append(values[-1] + step)
    step *= growth
  while values[-1] < cap_height and values[-1] + cap < total:
    values.append(values[-1] + cap)
  remaining = total - values[-1]
  cells = max(1, math.ceil(math.log(remaining * (growth - 1) / cap + 1, growth)))
  ratio = solve_growth(cap, cells, remaining)
  step = cap
  for _ in range(cells):
    values.append(values[-1] + step)
    step *= ratio
  return values


# 径向: 贴壁几何增长 -> 厚度封顶 -> 再增长铺满整条径向线
lines_of_points = [[radial[i][j] for i in range(ni)] for j in range(ncir)]
stations_list = [stations(points) for points in lines_of_points]
spans = [station[-1] for station in stations_list]
far_ring = radial[ni - 1]
center = (sum(p[0] for p in far_ring) / ncir, sum(p[1] for p in far_ring) / ncir)
reach = [math.dist(points[-1], center) for points in lines_of_points]
outward = [
  ((points[-1][0] - center[0]) / r, (points[-1][1] - center[1]) / r)
  for points, r in zip(lines_of_points, reach)
]
# 远场环按比例外移, 超出原网格的部分沿径向直线外推
far_radius = sum(reach) / ncir
scale = args.far_radius / far_radius
extended = [span + (scale - 1.0) * r for span, r in zip(spans, reach)]
mean_span = sum(extended) / ncir
offsets = radial_offsets(
  mean_span, args.first_height, args.growth, args.max_cell, args.max_cell_height
)
ncell = len(offsets) - 1
fractions = [offset / offsets[-1] for offset in offsets]

rings = []
for points, station, span, direction in zip(lines_of_points, stations_list, extended, outward):
  ring = []
  for fraction in fractions:
    target = fraction * span
    if target <= station[-1]:
      ring.append(sample(points, station, target / station[-1]))
    else:
      length = target - station[-1]
      ring.append((points[-1][0] + direction[0] * length, points[-1][1] + direction[1] * length))
  rings.append(ring)
rings.append(rings[0])
radial = [[rings[j][i] for j in range(ncir + 1)] for i in range(ncell + 1)]
nrad = ncell + 1

# 周向: 上表面激波附近按目标间距加密, 其余保持原分布
wall = [radial[0][j] for j in range(ncir)]
chord = max(p[0] for p in wall) - min(p[0] for p in wall)
x_le = min(p[0] for p in wall)
spacing = [math.dist(wall[j], wall[(j + 1) % ncir]) for j in range(ncir)]
xc = [(p[0] - x_le) / chord for p in wall]
shock = min((j for j in range(ncir) if wall[j][1] > 0.0), key=lambda j: abs(xc[j] - args.shock_xc))
refine = spacing[shock] / args.shock_dx


def density_profile(u, plateau):
  t = abs(u)
  if t >= 1.0:
    return 0.0
  if t <= plateau:
    return 1.0
  return ((1.0 - t) / (1.0 - plateau)) ** 2


density = []
for j in range(ncir):
  u = (xc[j] - args.shock_xc) / args.shock_width
  bump = density_profile(u, 0.5) if wall[j][1] > 0.0 else 0.0
  density.append(1.0 + (refine - 1.0) * bump)
coordinate = [0.0]
for value in density:
  coordinate.append(coordinate[-1] + value)
ncir_new = max(1, round(coordinate[-1]))


def old_index(m):
  target = m
  for k in range(len(coordinate) - 1):
    if target <= coordinate[k + 1]:
      width = coordinate[k + 1] - coordinate[k]
      weight = 0.0 if width <= 0.0 else (target - coordinate[k]) / width
      return k + weight
  return ncir


def ring_point(points, position):
  base = int(math.floor(position)) % ncir
  weight = position - math.floor(position)
  head, tail = points[base], points[(base + 1) % ncir]
  return (head[0] + weight * (tail[0] - head[0]), head[1] + weight * (tail[1] - head[1]))


rings = [
  [ring_point([radial[i][j] for j in range(ncir)], old_index(m)) for m in range(ncir_new)]
  for i in range(nrad)
]
radial = [ring + [ring[0]] for ring in rings]
ncir = ncir_new


def node(i, j):
  return i * ncir + j % ncir + 1


def cell(i, j):
  return i * ncir + j % ncir + 1


def ring_face(i, j):
  return i * ncir + j % ncir + 1


def radial_face(i, j):
  return nrad * ncir + i * ncir + j % ncir + 1


nodes = [radial[i][j] for i in range(nrad) for j in range(ncir)]
cells = [
  (ring_face(i, j), ring_face(i + 1, j), radial_face(i, j), radial_face(i, j + 1))
  for i in range(nrad - 1)
  for j in range(ncir)
]

wall = [(ring_face(0, j), node(0, j), node(0, j + 1), cell(0, j), 0) for j in range(ncir)]
far = [
  (ring_face(nrad - 1, j), node(nrad - 1, j), node(nrad - 1, j + 1), cell(nrad - 2, j), 0)
  for j in range(ncir)
]
interior = [
  (ring_face(i, j), node(i, j), node(i, j + 1), cell(i - 1, j), cell(i, j))
  for i in range(1, nrad - 1)
  for j in range(ncir)
] + [
  (radial_face(i, j), node(i, j), node(i + 1, j), cell(i, j), cell(i, j - 1))
  for i in range(nrad - 1)
  for j in range(ncir)
]

faces = wall + far + interior
face_count = nrad * ncir + (nrad - 1) * ncir
assert len(faces) == face_count, (len(faces), face_count)
assert len(cells) == (nrad - 1) * ncir
assert sorted(face[0] for face in faces) == list(range(1, face_count + 1)), (
  '面编号不是 1..face_count'
)
for face in faces:
  assert 1 <= face[1] <= len(nodes) and 1 <= face[2] <= len(nodes)
for row in cells:
  assert len(set(row)) == 4, '单元引用了重复的面'
refs = Counter(face_id for row in cells for face_id in row)
for face in wall + far:
  assert refs[face[0]] == 1, f'边界面 {face[0]} 被 {refs[face[0]]} 个单元引用'
for face in interior:
  assert refs[face[0]] == 2, f'内部面 {face[0]} 被 {refs[face[0]]} 个单元引用'

lines = [
  f'{len(nodes)} {face_count} {len(cells)} 3',
  'interior = INTER',
  'airfoil = WALL',
  'farfield = FAR',
  '(node)',
]
lines += [f'{i} {x:.16e} {y:.16e}' for i, (x, y) in enumerate(nodes, 1)]
lines += ['(end)', '(edge)']
for name, group in [('interior', interior), ('airfoil', wall), ('farfield', far)]:
  lines.append(name)
  lines += [' '.join(map(str, face)) for face in group]
  lines.append('(end)')
lines += ['(end)', '(cell)']
lines += [' '.join(map(str, (i, *row))) for i, row in enumerate(cells, 1)]

args.target.write_text('\n'.join(lines) + '\n')
print(f'{args.target}: nodes={len(nodes)} faces={face_count} cells={len(cells)}')
print(f'  WALL={len(wall)} FAR={len(far)} INTER={len(interior)}')
print(
  f'  径向 {ncell} 层, 首层 {min(fractions[1] * length for length in extended):.3e} ~'
  f' {max(fractions[1] * length for length in extended):.3e}, 厚度上限 {args.max_cell:.4f}'
)
print(f'  周向 {ncir} 层 (原 {len(spacing)}), 激波 x/c={args.shock_xc:.3f} 处加密 {refine:.2f} 倍')
print(
  f'  弦长 {chord:.6f}, 远场半径 {args.far_radius:.1f} (加密前 {far_radius:.3f} 的'
  f' {scale:.3f} 倍), 周向单元 {ncir}'
)
