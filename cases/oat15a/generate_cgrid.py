#!/usr/bin/env python3
"""OAT15A C 型网格: 翼面 + 尾迹切口.

内环(逆时针):
    LE -> 下表面 -> 后缘下角 -> 底面下半 -> 中点(1)
       -> 尾迹下侧(顺流, 节点 0..nw-1)
       -> 尾迹上侧(逆流: 末端副本 -> 节点 nw-2 -> ... -> 节点 0 -> 中点(2))
       -> 底面上半 -> 后缘上角 -> 上表面 -> LE

尾迹两侧节点坐标重合(仅末端: 下侧下移 -slit, 上侧副本上移 +slit), 尾迹线上的面只
产生一次, 由上下两侧单元共用 —— 尾迹是内部面, 流动自由穿过.

与 O 型网格的关键区别: O 型把"从底面长出的 5e-6 首层"用于尾迹的顺流方向(于是尾迹
单元长宽比可达 1e3, 底面两个直角处还是异形单元); C 型里同一套首层横跨尾迹 —— 正好
分辨后缘剪切层.

用法:
    python3 generate_cgrid.py oat15a.txt oat15a_c.txt
"""
import argparse
import math
from collections import Counter
from pathlib import Path

HERE = Path(__file__).resolve().parent

ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
ap.add_argument('source', nargs='?', type=Path, default=HERE / 'oat15a.txt')
ap.add_argument('target', nargs='?', type=Path, default=HERE / 'oat15a_c.txt')
ap.add_argument('--first-height', type=float, default=5.0e-6, help='翼面首层高度(弦长单位)')
ap.add_argument('--wake-first', type=float, default=1.5e-3, help='尾迹/底面径向首层(弦长单位)')
ap.add_argument('--growth', type=float, default=1.2, help='径向层高比')
ap.add_argument('--ncell', type=int, default=89, help='径向单元层数')
ap.add_argument('--max-cell', type=float, default=0.02, help='径向后段单元厚度上限')
ap.add_argument('--max-cell-height', type=float, default=0.45, help='厚度上限生效到多高')
ap.add_argument('--far-radius', type=float, default=10.0, help='远场圆半径')
ap.add_argument('--center-x', type=float, default=0.0, help='远场圆心 x')
ap.add_argument('--center-y', type=float, default=0.0, help='远场圆心 y')
ap.add_argument('--bulge', type=float, default=1.0, help='径向线法向控制点长度')
ap.add_argument('--shock-xc', type=float, default=0.5, help='激波弦向位置(x/c)')
ap.add_argument('--shock-dx', type=float, default=0.0015, help='激波处弦向目标间距')
ap.add_argument('--shock-width', type=float, default=0.13, help='激波加密窗口半宽')
ap.add_argument('--te-dx', type=float, default=0.0002, help='后缘底面节点间距')
ap.add_argument('--te-width', type=float, default=0.10, help='后缘加密窗口半宽')
ap.add_argument('--te-plateau', type=float, default=0.02, help='后缘满加密平台宽度(占窗口比例)')
ap.add_argument('--ramp', type=float, default=0.002, help='后缘角处首层过渡弧长')
ap.add_argument('--wake-dx', type=float, default=0.0006, help='尾迹线首段间距')
ap.add_argument('--wake-cell', type=float, default=0.04, help='尾迹线段间距上限')
ap.add_argument('--wake-growth', type=float, default=1.18, help='尾迹线间距增长率')
ap.add_argument('--wake-x', type=float, default=6.0, help='尾迹切口终点 x/c')
ap.add_argument('--wake-tail', type=float, default=0.005, help='尾迹末段长度(缝帽宽度)')
ap.add_argument('--slit', type=float, default=1.0e-4, help='尾迹末端上下侧错开量')
ap.add_argument('--psi-w', type=float, default=40.0, help='尾迹切口在远场占据的半角(度)')
ap.add_argument('--base-span', type=float, default=1.0, help='底面在远场占据的角度(度)')
ap.add_argument('--nrm-smooth', type=float, default=0.05, help='壁面法向按弧长平滑窗口')
ap.add_argument('--dense', type=int, default=400, help='每条径向线的密集采样点数')
args = ap.parse_args()

# ---------------- 表面 ----------------
raw = []
for line in args.source.read_text().splitlines():
  token = line.split()
  if len(token) >= 3:
    raw.append((float(token[1]), float(token[2])))
if len(raw) < 20:
  raise SystemExit(f'表面点太少: {len(raw)}')
if math.dist(raw[0], raw[-1]) < 1e-9:
  raw.pop()

x_min = min(p[0] for p in raw)
x_max = max(p[0] for p in raw)
chord_raw = x_max - x_min
surf = [((p[0] - x_min) / chord_raw - 0.5, p[1] / chord_raw) for p in raw]
n = len(surf)

area2 = sum(
  surf[j][0] * surf[(j + 1) % n][1] - surf[(j + 1) % n][0] * surf[j][1] for j in range(n)
)
if area2 < 0.0:
  surf.reverse()

te_up = max(range(n), key=lambda i: (surf[i][0], surf[i][1]))
te_lo = max(range(n), key=lambda i: (surf[i][0], -surf[i][1]))
le = min(range(n), key=lambda i: surf[i][0])


def walk(a, b):
  out, i = [], a
  while i != b:
    out.append(i)
    i = (i + 1) % n
    if len(out) > n:
      raise SystemExit('轮廓顺序异常')
  return out


lower_idx = walk(le, te_lo)       # LE -> 后缘下角 (含 LE, 不含下角)
upper_idx = walk(te_up, le)       # 后缘上角 -> LE (含上角, 不含 LE)

seg = [(surf[(j + 1) % n][0] - surf[j][0], surf[(j + 1) % n][1] - surf[j][1]) for j in range(n)]
seg_len = [math.hypot(s[0], s[1]) for s in seg]
mid = [
  ((surf[j][0] + surf[(j + 1) % n][0]) * 0.5, (surf[j][1] + surf[(j + 1) % n][1]) * 0.5)
  for j in range(n)
]
xc_mid = [p[0] + 0.5 for p in mid]
upper_seg = [j for j in range(n) if mid[j][1] > 0.0]
shock = min(upper_seg, key=lambda j: abs(xc_mid[j] - args.shock_xc))
refine = seg_len[shock] / args.shock_dx
te_all = max(range(n), key=lambda j: xc_mid[j])
refine_te = seg_len[te_all] / args.te_dx


def density_profile(u, plateau):
  t = abs(u)
  if t >= 1.0:
    return 0.0
  if t <= plateau:
    return 1.0
  return ((1.0 - t) / (1.0 - plateau)) ** 2


def density_at(point):
  u = (point[0] + 0.5 - args.shock_xc) / args.shock_width
  bump = density_profile(u, 0.5) if point[1] > 0.0 else 0.0
  v = (point[0] + 0.5 - 1.0) / args.te_width
  return 1.0 + (refine - 1.0) * bump + (refine_te - 1.0) * density_profile(v, args.te_plateau)


def refine_chain(index_list):
  """按目标间距重采样表面点链, 返回点列表(含起点, 不含终点)."""
  pts = [surf[i] for i in index_list]
  m = len(pts) - 1
  coord = [0.0]
  for k in range(m):
    coord.append(
      coord[-1]
      + density_at(((pts[k][0] + pts[k + 1][0]) * 0.5, (pts[k][1] + pts[k + 1][1]) * 0.5))
    )
  total = max(2, round(coord[-1]))
  out = []
  for mm in range(total):
    k = 0
    while k < m - 1 and mm > coord[k + 1]:
      k += 1
    width = coord[k + 1] - coord[k]
    w = 0.0 if width <= 0.0 else (mm - coord[k]) / width
    out.append((pts[k][0] + w * (pts[k + 1][0] - pts[k][0]), pts[k][1] + w * (pts[k + 1][1] - pts[k][1])))
  return out


# ---------------- 内环 ----------------
x_mid = (surf[te_lo][0] + surf[te_up][0]) * 0.5
y_mid = (surf[te_lo][1] + surf[te_up][1]) * 0.5

loop = []
kind = []

lower_pts = refine_chain(lower_idx)
n_lower = len(lower_pts)
loop += lower_pts
kind += ['wall'] * n_lower

half = abs(surf[te_lo][1] - y_mid)
nb = max(2, round(half / args.te_dx))
for k in range(nb + 1):                       # 后缘下角 -> 中点(1), 含两端
  w = k / nb
  loop.append(
    (surf[te_lo][0] + w * (x_mid - surf[te_lo][0]), surf[te_lo][1] + w * (y_mid - surf[te_lo][1]))
  )
  kind.append('wall')

wake_pts = []                                  # 切口下侧 slit 节点(不含中点)
wx = x_mid
step = args.wake_dx
while True:
  nxt = wx + step
  if args.wake_x - nxt < args.wake_tail:
    break
  wake_pts.append((nxt, y_mid))
  wx = nxt
  step = min(step * args.wake_growth, args.wake_cell)
wake_pts.append((args.wake_x, y_mid))
nw = len(wake_pts)
if nw < 4:
  raise SystemExit(f'尾迹段数太少: {nw}')

# 切口末端用与切口相切的半圆帽(半径 slit): 环是光滑简单曲线, 法向连续旋转,
# 末端单元保持凸. 若末端直接直角收口, 相邻两个位置的径向线会反向张开, 拉出自交单元.
eps = args.slit
c0 = len(loop)                                 # 下侧 slit 起点(中点(1) 之后)
for p in wake_pts[:-1]:
  loop.append(p)
  kind.append('wake')
loop.append(wake_pts[-1])                      # 半圆起点(= slit 末端)
kind.append('cap')                             # 出发的段已经是圆帽首段
c_slit_end = len(loop) - 1

cap_count = 8
for k in range(1, cap_count + 1):              # 半圆: 下 -> 右 -> 上
  a = -math.pi / 2.0 + math.pi * k / cap_count
  loop.append((args.wake_x + eps * math.cos(a), y_mid + eps * (1.0 + math.sin(a))))
  kind.append('cap')
cap_end = len(loop) - 1

d0 = len(loop)                                 # 上侧 slit(回程), 与下侧镜像
for p in reversed(wake_pts[:-1]):
  loop.append((p[0], y_mid + 2.0 * eps))
  kind.append('wake2')
loop.append((x_mid, y_mid))                    # 中点(2)
kind.append('wall')

for k in range(1, nb + 1):                     # 中点(2) -> 后缘上角, 含上角
  w = k / nb
  loop.append((x_mid + w * (surf[te_up][0] - x_mid), y_mid + w * (surf[te_up][1] - y_mid)))
  kind.append('wall')
te_up_pos = len(loop) - 1

upper_pts = refine_chain(upper_idx[1:])        # 后缘上角之后 -> LE
loop += upper_pts
kind += ['wall'] * len(upper_pts)

kind[c0 - 1] = 'wake'          # 中点(1) 出发的段已经是切口首段(内部面)
kind[cap_end] = 'wake2'        # 半圆终点出发的段是上侧切口首段
ncir = len(loop)
nrad = args.ncell + 1
assert len(kind) == ncir
assert kind[c_slit_end] == 'cap' and kind[cap_end] == 'wake2'


# ---------------- 弧长 / 法向 ----------------
def unit(v):
  length = math.hypot(v[0], v[1])
  return (v[0] / length, v[1] / length) if length > 0.0 else (0.0, 0.0)


wall_seg = [math.dist(loop[j], loop[(j + 1) % ncir]) for j in range(ncir)]
cumulative = [0.0]
for length in wall_seg:
  cumulative.append(cumulative[-1] + length)
total_arc = cumulative[-1]

tang = []
for j in range(ncir):
  a = unit((loop[j][0] - loop[(j - 1) % ncir][0], loop[j][1] - loop[(j - 1) % ncir][1]))
  b = unit((loop[(j + 1) % ncir][0] - loop[j][0], loop[(j + 1) % ncir][1] - loop[j][1]))
  tang.append(unit((a[0] + b[0], a[1] + b[1])))
nrm = [(t[1], -t[0]) for t in tang]

if args.nrm_smooth > 0.0:
  window = args.nrm_smooth
  smoothed = []
  for j in range(ncir):
    sx = sy = 0.0
    for d in range(ncir):
      delta = abs(cumulative[d] - cumulative[j])
      delta = min(delta, total_arc - delta)
      if delta <= window:
        weight = 1.0 - delta / window
        sx += weight * nrm[d][0]
        sy += weight * nrm[d][1]
    smoothed.append(unit((sx, sy)) if (sx, sy) != (0.0, 0.0) else nrm[j])
  nrm = smoothed

# ---------------- 远场角度 ----------------
psi = math.radians(args.psi_w)
span = math.radians(args.base_span)
eps = psi / (4.0 * nw)
theta = [0.0] * ncir

# 用整环弧长上的"分段线性断点"映射: 每段端点只在一个断点上, 相邻两个不同的环位置
# 一定拿到不同角度 —— 否则远场点会重合, 外边界出现零长度边.
breaks = [
  (0.0, -math.pi),
  (cumulative[n_lower], -(psi + span)),
  (cumulative[c0 - 1], -psi),
  (cumulative[c_slit_end], -eps),
  (cumulative[cap_end], eps),
  (cumulative[d0 + nw - 1], psi),
  (cumulative[te_up_pos], psi + span),
  (total_arc, math.pi),
]
for k in range(len(breaks) - 1):
  assert breaks[k][0] < breaks[k + 1][0], f'断点弧长未递增: {breaks[k]}'

for j in range(ncir):
  s = cumulative[j]
  k = 0
  while k < len(breaks) - 2 and s > breaks[k + 1][0]:
    k += 1
  s0, a0 = breaks[k]
  s1, a1 = breaks[k + 1]
  w = 0.0 if s1 <= s0 else (s - s0) / (s1 - s0)
  theta[j] = a0 + w * (a1 - a0)

radius = args.far_radius
center = (args.center_x, args.center_y)
far = [
  (center[0] + radius * math.cos(theta[j]), center[1] + radius * math.sin(theta[j]))
  for j in range(ncir)
]

# ---------------- 径向线 ----------------
dense_count = max(50, args.dense)
curves = []
line_length = []
for j in range(ncir):
  point = loop[j]
  end = far[j]
  normal = nrm[j]
  reach = math.dist(point, end)
  radial_dir = unit((end[0] - center[0], end[1] - center[1]))
  straight = unit((end[0] - point[0], end[1] - point[1]))
  cross = abs(normal[0] * straight[1] - normal[1] * straight[0])
  a1 = min(reach / 3.0, args.bulge / max(cross, 1e-6))
  ctrl1 = (point[0] + a1 * normal[0], point[1] + a1 * normal[1])
  ctrl2 = (end[0] - reach / 3.0 * radial_dir[0], end[1] - reach / 3.0 * radial_dir[1])
  samples = []
  for k in range(dense_count + 1):
    t = k / dense_count
    s = 1.0 - t
    w = (s * s * s, 3.0 * s * s * t, 3.0 * s * t * t, t * t * t)
    samples.append(
      (
        w[0] * point[0] + w[1] * ctrl1[0] + w[2] * ctrl2[0] + w[3] * end[0],
        w[0] * point[1] + w[1] * ctrl1[1] + w[2] * ctrl2[1] + w[3] * end[1],
      )
    )
  arc = [0.0]
  for k in range(1, dense_count + 1):
    arc.append(arc[-1] + math.dist(samples[k], samples[k - 1]))
  curves.append((samples, arc))
  line_length.append(arc[-1])

mean_span = sum(line_length) / ncir


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


def resample(values, count):
  out = []
  m = len(values) - 1
  for k in range(count + 1):
    t = k / count * m
    i = min(int(t), m - 1)
    w = t - i
    out.append(values[i] * (1.0 - w) + values[i + 1] * w)
  return out


offsets_a = radial_offsets(mean_span, args.first_height, args.growth, args.max_cell, args.max_cell_height)
offsets_b = radial_offsets(mean_span, args.wake_first, args.growth, args.max_cell, args.max_cell_height)
frac_a = resample([v / offsets_a[-1] for v in offsets_a], args.ncell)
frac_b = resample([v / offsets_b[-1] for v in offsets_b], args.ncell)

# 底面/切口/半圆帽用"尾迹首层"(底面上没有边界层, 首层不该被压到 5e-6), 翼面用
# "壁面首层". 过渡放在底面两半上线性完成: 两端分别与翼面/切口连续, 避免首层高度
# 突变把后缘角的单元拉成退化斜四边形.
def blend(j):
  if j < n_lower or j > te_up_pos:
    return 0.0
  if j <= c0 - 1:
    return (j - n_lower) / max(1, c0 - 1 - n_lower)
  if j <= d0 + nw - 2:
    return 1.0
  return max(0.0, (te_up_pos - j) / max(1, te_up_pos - (d0 + nw - 1)))


weight = [blend(j) for j in range(ncir)]
fractions = [[(1.0 - w) * frac_a[k] + w * frac_b[k] for k in range(args.ncell + 1)] for w in weight]

rings = []
for j in range(ncir):
  samples, arc = curves[j]
  total = arc[-1]
  ring = []
  pointer = 0
  for fraction in fractions[j]:
    target = fraction * total
    while pointer < dense_count and arc[pointer + 1] < target:
      pointer += 1
    low = arc[pointer]
    high = arc[pointer + 1] if pointer + 1 <= dense_count else low
    w = 0.0 if high <= low else (target - low) / (high - low)
    head = samples[pointer]
    tail = samples[pointer + 1] if pointer + 1 <= dense_count else samples[pointer]
    ring.append((head[0] + w * (tail[0] - head[0]), head[1] + w * (tail[1] - head[1])))
  rings.append(ring)

radial = [[rings[j][i] for j in range(ncir)] for i in range(nrad)]

# ---------------- 拓扑 ----------------
# 下侧位置 c0+k (节点 k) 与上侧位置 cap_end+(nw-1-k) 是同一节点; 下侧边 (c0-1+k)
# 与上侧边 (cap_end+(nw-1-k)) 是同一条几何边, 只在下侧产生面, 上下两侧单元共用.
share_face = {cap_end + (nw - 1 - k): c0 - 1 + k for k in range(nw)}
assert len(share_face) == nw


# 环 0(内边界)的面: 尾迹上侧那一批被映射到共享 id, 所以自己的编号不能预留 ——
# 这里只对"实际存在"的内边界面做紧凑编号, 否则面 id 会出现空洞.
canon = {}
for j in range(ncir):
  if kind[j] in ('wall', 'wake', 'cap'):
    canon[j] = len(canon)
n0 = len(canon)
face_count = n0 + 2 * (nrad - 1) * ncir


def node(i, j):
  return i * ncir + j + 1


def ring_face(i, j):
  if i == 0:
    if j in share_face:
      j = share_face[j]
    return canon[j] + 1
  return n0 + (i - 1) * ncir + (j % ncir) + 1


def radial_face(i, j):
  return n0 + (nrad - 1) * ncir + i * ncir + (j % ncir) + 1


nodes = [radial[i][j] for i in range(nrad) for j in range(ncir)]
cells = [
  (ring_face(i, j), ring_face(i + 1, j), radial_face(i, j), radial_face(i, j + 1))
  for i in range(nrad - 1)
  for j in range(ncir)
]

wall_index = [j for j in range(ncir) if kind[j] in ('wall', 'cap')]
wall_group = [(ring_face(0, j), node(0, j), node(0, (j + 1) % ncir), j + 1, 0) for j in wall_index]
far_group = [
  (
    ring_face(nrad - 1, j),
    node(nrad - 1, j),
    node(nrad - 1, (j + 1) % ncir),
    (nrad - 2) * ncir + j + 1,
    0,
  )
  for j in range(ncir)
]
interior = [
  (
    ring_face(i, j),
    node(i, j),
    node(i, (j + 1) % ncir),
    (i - 1) * ncir + j + 1,
    i * ncir + j + 1,
  )
  for i in range(1, nrad - 1)
  for j in range(ncir)
] + [
  (
    radial_face(i, j),
    node(i, j),
    node(i + 1, j),
    i * ncir + j + 1,
    i * ncir + (j - 1) % ncir + 1,
  )
  for i in range(nrad - 1)
  for j in range(ncir)
]
wake_link = [
  (
    ring_face(0, c0 - 1 + k),
    node(0, c0 - 1 + k),
    node(0, c0 + k),
    c0 + k,
    cap_end + (nw - 1 - k) + 1,
  )
  for k in range(nw)
]

faces = interior + wall_group + far_group + wake_link

assert len(cells) == (nrad - 1) * ncir, (len(cells), (nrad - 1) * ncir)
assert len(set(f[0] for f in faces)) == face_count, '面 id 不完整/重复'
assert sorted(f[0] for f in faces) == list(range(1, face_count + 1)), '面 id 不连续'
refs = Counter(face_id for row in cells for face_id in row)
for face in wall_group + far_group:
  assert refs[face[0]] == 1, f'边界面 {face[0]} 被引用 {refs[face[0]]} 次'
for face in interior + wake_link:
  assert refs[face[0]] == 2, f'内部面 {face[0]} 被引用 {refs[face[0]]} 次'

# ---------------- 质量自检 ----------------
signs = []
for i in range(nrad - 1):
  for j in range(ncir):
    quad = (radial[i][j], radial[i][(j + 1) % ncir], radial[i + 1][(j + 1) % ncir], radial[i + 1][j])
    value = 0.0
    for k in range(4):
      head, tail = quad[k], quad[(k + 1) % 4]
      value += head[0] * tail[1] - head[1] * tail[0]
    signs.append(0.5 * value)
negative = sum(1 for v in signs if v <= 0.0)

aspects = []
for j in range(ncir):
  for i in range(nrad - 1):
    a = math.dist(radial[i][j], radial[i][(j + 1) % ncir])
    b = math.dist(radial[i + 1][j], radial[i + 1][(j + 1) % ncir])
    c = math.dist(radial[i][j], radial[i + 1][j])
    d = math.dist(radial[i][(j + 1) % ncir], radial[i + 1][(j + 1) % ncir])
    aspects.append(max(a, b) / max(min(a, b, c, d), 1e-30))
aspects.sort()

direction = [unit((radial[1][j][0] - loop[j][0], radial[1][j][1] - loop[j][1])) for j in range(ncir)]
orth = []
for j in range(ncir):
  dot = abs(nrm[j][0] * direction[j][0] + nrm[j][1] * direction[j][1])
  a = math.degrees(math.acos(max(-1.0, min(1.0, dot))))
  orth.append(min(a, 180.0 - a))

lines = [
  f'{len(nodes)} {face_count} {len(cells)} 3',
  'interior = INTER',
  'airfoil = WALL',
  'farfield = FAR',
  '(node)',
]
lines += [f'{i} {x:.16e} {y:.16e}' for i, (x, y) in enumerate(nodes, 1)]
lines += ['(end)', '(edge)']
for name, group in [('interior', interior + wake_link), ('airfoil', wall_group), ('farfield', far_group)]:
  lines.append(name)
  lines += [' '.join(map(str, face)) for face in group]
  lines.append('(end)')
lines += ['(end)', '(cell)']
lines += [' '.join(map(str, (i, *row))) for i, row in enumerate(cells, 1)]
args.target.write_text('\n'.join(lines) + '\n')

print(f'{args.target}: nodes={len(nodes)} faces={face_count} cells={len(cells)}')
print(f'  内环 {ncir} 点: 下表面 {n_lower}, 底面 {nb}x2, 尾迹 {nw}x2, 上表面 {len(upper_pts)}')
print(f'  WALL={len(wall_group)} FAR={len(far_group)} 尾迹内部面={len(wake_link)}')
print(f'  激波 x/c={args.shock_xc:.3f} 加密 {refine:.2f} 倍 (目标 {args.shock_dx:.5f})')
print(f'  后缘 x/c=1.0 加密 {refine_te:.2f} 倍 (目标 {args.te_dx:.5f}, 窗口 {args.te_width:.3f})')
print(f'  尾迹切口 x/c {x_mid - 0.5:.2f} -> {args.wake_x:.2f}, {nw} 段/侧, 远场半角 {args.psi_w:.1f}°')
print(f'  径向 {args.ncell} 层, 线长 {min(line_length):.3f} ~ {max(line_length):.3f} (均值 {mean_span:.3f})')
print(f'  翼面首层 {frac_a[1] * mean_span:.3e}, 尾迹首层 {frac_b[1] * mean_span:.3e}')
print(f'  单元: {len(signs)} 个, 非正拐向 {negative} 个, 面积 {min(signs):.3e} ~ {max(signs):.3e}')
print(f'  长宽比: 中位 {aspects[len(aspects) // 2]:.1f}, 99% {aspects[int(len(aspects) * 0.99)]:.1f}, 最大 {aspects[-1]:.1f}')
print(f'  壁面正交偏差: 中位 {sorted(orth)[ncir // 2]:.2f}° 最大 {max(orth):.2f}°')
