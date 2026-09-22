import argparse
import math
from collections import Counter
from pathlib import Path

parser = argparse.ArgumentParser(
  description='Convert a ring point list (mesh2.txt) to the solver mesh format'
)
parser.add_argument(
  'source', nargs='?', type=Path, default=Path(__file__).resolve().parent / 'mesh2.txt'
)
parser.add_argument(
  'target', nargs='?', type=Path, default=Path(__file__).resolve().parent / 'mesh2_cfd.txt'
)
args = parser.parse_args()

rows = args.source.read_text().split('\n')
na, nr = (int(v) for v in rows[0].split())
points = [tuple(float(v) for v in row.split()) for row in rows[1:] if row.strip()]
if len(points) != na * nr:
  raise SystemExit(f'点数为 {len(points)}, 与头部声明的 {na}*{nr} 不符')

# 每 na 个点是一圈, 首尾点重合, 因此一圈只有 na-1 个不同节点
rings = [points[j * na : (j + 1) * na] for j in range(nr)]
for j, ring in enumerate(rings):
  if ring[0] != ring[-1]:
    raise SystemExit(f'第 {j} 圈首尾点不重合, 无法按闭合环处理')
nt = na - 1


def node(j, k):
  return j * nt + k % nt + 1


def cell(j, k):
  return j * nt + k % nt + 1


def ring_face(j, k):
  return j * nt + k % nt + 1


def radial_face(j, k):
  return nr * nt + j * nt + k % nt + 1


nodes = [rings[j][k] for j in range(nr) for k in range(nt)]
cells = [
  (ring_face(j, k), ring_face(j + 1, k), radial_face(j, k), radial_face(j, k + 1))
  for j in range(nr - 1)
  for k in range(nt)
]

wall = [(ring_face(0, k), node(0, k), node(0, k + 1), cell(0, k), 0) for k in range(nt)]
far = [
  (ring_face(nr - 1, k), node(nr - 1, k), node(nr - 1, k + 1), cell(nr - 2, k), 0)
  for k in range(nt)
]
interior = [
  (ring_face(j, k), node(j, k), node(j, k + 1), cell(j - 1, k), cell(j, k))
  for j in range(1, nr - 1)
  for k in range(nt)
] + [
  (radial_face(j, k), node(j, k), node(j + 1, k), cell(j, k), cell(j, k - 1))
  for j in range(nr - 1)
  for k in range(nt)
]

faces = wall + far + interior
face_count = nr * nt + (nr - 1) * nt
assert len(faces) == face_count, (len(faces), face_count)
assert len(cells) == (nr - 1) * nt
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
  'cylinder = WALL',
  'farfield = FAR',
  '(node)',
]
lines += [f'{i} {x:.16e} {y:.16e}' for i, (x, y) in enumerate(nodes, 1)]
lines += ['(end)', '(edge)']
for name, group in [('interior', interior), ('cylinder', wall), ('farfield', far)]:
  lines.append(name)
  lines += [' '.join(map(str, face)) for face in group]
  lines.append('(end)')
lines += ['(end)', '(cell)']
lines += [' '.join(map(str, (i, *row))) for i, row in enumerate(cells, 1)]

args.target.write_text('\n'.join(lines) + '\n')
centre = (25.0, 0.0)
print(f'{args.target}: nodes={len(nodes)} faces={face_count} cells={len(cells)}')
print(f'  WALL={len(wall)} FAR={len(far)} INTER={len(interior)}')
print(
  f'  圆柱半径 {sum(math.dist(p, centre) for p in rings[0][:-1]) / nt:.6f}, '
  f'远场半径 {sum(math.dist(p, centre) for p in rings[-1][:-1]) / nt:.6f}, 周向单元 {nt}'
)
