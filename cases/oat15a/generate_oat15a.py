import argparse
import math
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent

parser = argparse.ArgumentParser(
  description='由 OAT15A 表面轮廓生成求解器 O 型网格 (径向分布与弦向加密规则同 NACA0012)'
)
parser.add_argument('source', nargs='?', type=Path, default=HERE / 'oat15a.txt')
parser.add_argument('target', nargs='?', type=Path, default=HERE / 'oat15a_sa.txt')
parser.add_argument('--first-height', type=float, default=5.0e-6, help='壁面首层高度, 弦长 1.0')
parser.add_argument('--growth', type=float, default=1.2, help='径向相邻层高度比')
parser.add_argument('--far-radius', type=float, default=10.0, help='远场圆半径, 弦长单位')
parser.add_argument('--max-cell', type=float, default=0.02, help='径向后段单元厚度上限')
parser.add_argument('--max-cell-height', type=float, default=0.45, help='厚度上限生效到多高')
parser.add_argument('--shock-xc', type=float, default=0.5, help='激波弦向位置')
parser.add_argument('--shock-dx', type=float, default=0.0015, help='激波处弦向目标间距')
parser.add_argument('--shock-width', type=float, default=0.13, help='激波加密窗口半宽')
parser.add_argument('--center-x', type=float, default=0.0, help='远场圆心 x (弦中点坐标)')
parser.add_argument('--center-y', type=float, default=0.0, help='远场圆心 y (弦中点坐标)')
parser.add_argument('--bulge', type=float, default=0.08, help='法向控制点允许造成的最大横向偏离(弦长)')
parser.add_argument('--dense', type=int, default=400, help='每条径向线的密集采样点数')
args = parser.parse_args()
if args.first_height <= 0.0 or args.growth <= 1.0:
  raise SystemExit('--first-height 需为正, --growth 需大于 1')
if args.max_cell <= args.first_height or args.max_cell_height <= 0.0:
  raise SystemExit('--max-cell 需大于首层高度, --max-cell-height 需为正')
if args.shock_dx <= 0.0 or args.shock_width <= 0.0:
  raise SystemExit('--shock-dx 与 --shock-width 需为正')

# --- 读表面轮廓 (每行: 块号 x y) ---
raw = []
for line in args.source.read_text().splitlines():
  token = line.split()
  if len(token) >= 3:
    raw.append((float(token[1]), float(token[2])))
if len(raw) < 20:
  raise SystemExit(f'表面点太少: {len(raw)}')
if math.dist(raw[0], raw[-1]) < 1e-9:
  raw.pop()

# --- 归一化到弦长 1, 原点移到弦中点 ---
x_min = min(p[0] for p in raw)
x_max = max(p[0] for p in raw)
chord_raw = x_max - x_min
surf = [((p[0] - x_min) / chord_raw - 0.5, p[1] / chord_raw) for p in raw]

# --- 方向统一成逆时针, 与 NACA0012 的 j 顺序一致 ---
area2 = sum(
  surf[j][0] * surf[(j + 1) % len(surf)][1] - surf[(j + 1) % len(surf)][0] * surf[j][1]
  for j in range(len(surf))
)
reversed_input = area2 < 0.0
if reversed_input:
  surf.reverse()
n = len(surf)


def unit(v):
  length = math.hypot(v[0], v[1])
  return (v[0] / length, v[1] / length) if length > 0.0 else (0.0, 0.0)


seg = [(surf[(j + 1) % n][0] - surf[j][0], surf[(j + 1) % n][1] - surf[j][1]) for j in range(n)]
seg_len = [math.hypot(s[0], s[1]) for s in seg]
tang = [unit(s) for s in seg]
# 角平分线外法向 (逆时针轮廓: 外法向 = 切向顺时针转 90°)
nrm = []
for j in range(n):
  a = tang[(j - 1) % n]
  b = tang[j]
  t = unit((a[0] + b[0], a[1] + b[1]))
  nrm.append((t[1], -t[0]))

stat = [0.0]
for length in seg_len:
  stat.append(stat[-1] + length)
mid = [
  ((surf[j][0] + surf[(j + 1) % n][0]) * 0.5, (surf[j][1] + surf[(j + 1) % n][1]) * 0.5)
  for j in range(n)
]
xc_mid = [p[0] + 0.5 for p in mid]

# --- 周向: 上表面激波附近按目标间距加密, 其余保持原分布 ---
upper = [j for j in range(n) if mid[j][1] > 0.0]
if not upper:
  raise SystemExit('轮廓里找不到上表面点')
shock = min(upper, key=lambda j: abs(xc_mid[j] - args.shock_xc))
refine = seg_len[shock] / args.shock_dx


def density_profile(u, plateau):
  t = abs(u)
  if t >= 1.0:
    return 0.0
  if t <= plateau:
    return 1.0
  return ((1.0 - t) / (1.0 - plateau)) ** 2


density = []
for j in range(n):
  u = (xc_mid[j] - args.shock_xc) / args.shock_width
  bump = density_profile(u, 0.5) if mid[j][1] > 0.0 else 0.0
  density.append(1.0 + (refine - 1.0) * bump)
coordinate = [0.0]
for value in density:
  coordinate.append(coordinate[-1] + value)
ncir = max(1, round(coordinate[-1]))


def old_position(m):
  for k in range(n):
    if m <= coordinate[k + 1]:
      width = coordinate[k + 1] - coordinate[k]
      weight = 0.0 if width <= 0.0 else (m - coordinate[k]) / width
      return k + weight
  return n


def interpolate(points, position):
  base = int(math.floor(position)) % n
  weight = position - math.floor(position)
  head, tail = points[base], points[(base + 1) % n]
  return (head[0] + weight * (tail[0] - head[0]), head[1] + weight * (tail[1] - head[1]))


wall = [interpolate(surf, old_position(m)) for m in range(ncir)]
wall_nrm = [unit(interpolate(nrm, old_position(m))) for m in range(ncir)]

center = (args.center_x, args.center_y)
radius = args.far_radius

# --- 远场角度: 按壁面弧长参数均匀铺满整圈 ---
# NACA0012 原始网格的远场点分布最接近"壁面弧长"(与该规律偏差 26.8°, 按几何角则差 59.5°).
# 弧长参数天然单调、光滑且精确绕整圈, 不会像"几何角平滑"那样在后缘凹陷处长出平台.
theta = []
last = None
for point in wall:
  value = math.atan2(point[1] - center[1], point[0] - center[0])
  if last is not None:
    while value - last > math.pi:
      value -= 2.0 * math.pi
    while value - last < -math.pi:
      value += 2.0 * math.pi
  theta.append(value)
  last = value
wall_seg = [math.dist(wall[j], wall[(j + 1) % ncir]) for j in range(ncir)]
cumulative = [0.0]
for length in wall_seg:
  cumulative.append(cumulative[-1] + length)
total_arc = cumulative[-1]
phi = [theta[0] + 2.0 * math.pi * cumulative[j] / total_arc for j in range(ncir)]
twist = max(abs(math.degrees(p - t)) for p, t in zip(phi, theta))

far = [(center[0] + radius * math.cos(a), center[1] + radius * math.sin(a)) for a in phi]


def span_sum(first, growth, cells):
  return first * cells if growth == 1.0 else first * (growth**cells - 1) / (growth - 1)


def solve_growth(first, cells, total):
  lo, hi = 1.0, 2.0
  while span_sum(first, hi, cells) < total:
    hi *= 1.5
  for _ in range(200):
    middle = 0.5 * (lo + hi)
    if span_sum(first, middle, cells) < total:
      lo = middle
    else:
      hi = middle
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


# --- 径向线: 三次 Bezier, 起点沿壁面法向, 终点沿远场径向, 再按弧长重采样 ---
dense_count = max(50, args.dense)
curves = []
line_length = []
for j in range(ncir):
  point = wall[j]
  end = far[j]
  normal = wall_nrm[j]
  reach = math.dist(point, end)
  radial_dir = unit((end[0] - center[0], end[1] - center[1]))
  straight = unit((end[0] - point[0], end[1] - point[1]))
  # 法向与直达方向夹角越大, 法向控制点就要越短; 否则后缘上下两侧(法向几乎相反)
  # 的曲线会先冲出去再拐回来, 在中层拉出巨大的扇形畸变.
  # 缩短控制点不影响壁面正交性, 因为 B'(0) = 3*a1*n 始终平行于法向.
  cross = abs(normal[0] * straight[1] - normal[1] * straight[0])
  a1 = min(reach / 3.0, args.bulge / max(cross, 1e-6))
  ctrl1 = (point[0] + a1 * normal[0], point[1] + a1 * normal[1])
  ctrl2 = (end[0] - reach / 3.0 * radial_dir[0], end[1] - reach / 3.0 * radial_dir[1])
  samples = []
  for k in range(dense_count + 1):
    t = k / dense_count
    s = 1.0 - t
    weights = (s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t)
    samples.append(
      (
        weights[0] * point[0] + weights[1] * ctrl1[0] + weights[2] * ctrl2[0] + weights[3] * end[0],
        weights[0] * point[1] + weights[1] * ctrl1[1] + weights[2] * ctrl2[1] + weights[3] * end[1],
      )
    )
  arc = [0.0]
  for k in range(1, dense_count + 1):
    arc.append(arc[-1] + math.dist(samples[k], samples[k - 1]))
  curves.append((samples, arc))
  line_length.append(arc[-1])

mean_span = sum(line_length) / ncir
offsets = radial_offsets(
  mean_span, args.first_height, args.growth, args.max_cell, args.max_cell_height
)
ncell = len(offsets) - 1
nrad = ncell + 1
fractions = [offset / offsets[-1] for offset in offsets]

rings = []
for j in range(ncir):
  samples, arc = curves[j]
  total = arc[-1]
  ring = []
  pointer = 0
  for fraction in fractions:
    target = fraction * total
    while pointer < dense_count and arc[pointer + 1] < target:
      pointer += 1
    low = arc[pointer]
    high = arc[pointer + 1] if pointer + 1 <= dense_count else low
    weight = 0.0 if high <= low else (target - low) / (high - low)
    head = samples[pointer]
    tail = samples[pointer + 1] if pointer + 1 <= dense_count else samples[pointer]
    ring.append((head[0] + weight * (tail[0] - head[0]), head[1] + weight * (tail[1] - head[1])))
  rings.append(ring)

radial = [[rings[j][i] for j in range(ncir)] for i in range(nrad)]

# --- 质量自检 ---
signs = []
min_area = float('inf')
worst = None
for i in range(nrad - 1):
  for j in range(ncir):
    quad = (
      radial[i][j],
      radial[i][(j + 1) % ncir],
      radial[i + 1][(j + 1) % ncir],
      radial[i + 1][j],
    )
    value = 0.0
    for k in range(4):
      head = quad[k]
      tail = quad[(k + 1) % 4]
      value += head[0] * tail[1] - tail[0] * head[1]
    value *= 0.5
    signs.append(value)
    if value > min_area:
      pass
    if -value < -min_area or min_area == float('inf'):
      min_area = value
    if worst is None or value > worst[0]:
      worst = (value, i, j)

positive = sum(1 for value in signs if value > 0.0)
direction = [unit((radial[1][j][0] - wall[j][0], radial[1][j][1] - wall[j][1])) for j in range(ncir)]
orth = [
  math.degrees(
    math.acos(max(-1.0, min(1.0, abs(wall_nrm[j][0] * direction[j][0] + wall_nrm[j][1] * direction[j][1]))))
  )
  for j in range(ncir)
]
orth = [min(a, 180.0 - a) for a in orth]

# --- 输出 ---
def node(i, j):
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
wall_group = [(ring_face(0, j), node(0, j), node(0, j + 1), cells[j][0], 0) for j in range(ncir)]
far_group = [
  (
    ring_face(nrad - 1, j),
    node(nrad - 1, j),
    node(nrad - 1, j + 1),
    cells[(nrad - 2) * ncir + j][0],
    0,
  )
  for j in range(ncir)
]
interior = [
  (
    ring_face(i, j),
    node(i, j),
    node(i, j + 1),
    cells[(i - 1) * ncir + j][0],
    cells[i * ncir + j][0],
  )
  for i in range(1, nrad - 1)
  for j in range(ncir)
] + [
  (
    radial_face(i, j),
    node(i, j),
    node(i + 1, j),
    cells[i * ncir + j][0],
    cells[i * ncir + (j - 1) % ncir][0],
  )
  for i in range(nrad - 1)
  for j in range(ncir)
]
faces = wall_group + far_group + interior
face_count = nrad * ncir + (nrad - 1) * ncir
assert len(faces) == face_count, (len(faces), face_count)
assert len(cells) == (nrad - 1) * ncir
assert sorted(face[0] for face in faces) == list(range(1, face_count + 1))
refs = Counter(face_id for row in cells for face_id in row)
for face in wall_group + far_group:
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
for name, group in [('interior', interior), ('airfoil', wall_group), ('farfield', far_group)]:
  lines.append(name)
  lines += [' '.join(map(str, face)) for face in group]
  lines.append('(end)')
lines += ['(end)', '(cell)']
lines += [' '.join(map(str, (i, *row))) for i, row in enumerate(cells, 1)]

args.target.write_text('\n'.join(lines) + '\n')

wall_x = [wall[j][0] + 0.5 for j in range(ncir)]
wall_y = [wall[j][1] for j in range(ncir)]
print(f'{args.target}: nodes={len(nodes)} faces={face_count} cells={len(cells)}')
print(f'  WALL={len(wall_group)} FAR={len(far_group)} INTER={len(interior)}')
print(f'  原始轮廓 {n} 点 (弦长 {chord_raw:.4f}) -> 周向 {ncir} 点'
      f'{", 输入顺序已反转" if reversed_input else ""}')
print(f'  激波 x/c={args.shock_xc:.3f} 处加密 {refine:.2f} 倍, 目标间距 {args.shock_dx:.5f}')
print(f'  径向 {ncell} 层, 首层 {offsets[1] / offsets[-1] * min(line_length):.3e}'
      f' ~ {offsets[1] / offsets[-1] * max(line_length):.3e}, 远场半径 {radius:.3f}')
print(f'  壁面 y/c {min(wall_y):.5f} ~ {max(wall_y):.5f}, x/c {min(wall_x):.5f} ~ {max(wall_x):.5f}')
print(f'  拐向为正的单元 {positive} / {len(signs)}, 单元面积范围 {min(signs):.3e} ~ {max(signs):.3e}')
print(f'  质量最差单元 (i={worst[1]}, j={worst[2]}) 面积 {worst[0]:.3e}')
print(f'  壁面正交偏差: 中位 {sorted(orth)[ncir // 2]:.2f}° 最大 {max(orth):.2f}°')
print(f'  远场角跨度 {math.degrees(phi[-1] - phi[0]):.3f}°, 相对壁面几何角最大扭转 {twist:.3f}°')
