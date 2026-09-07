import argparse
import json
import math
from pathlib import Path

parser = argparse.ArgumentParser(description='Generate the Ma=0.2, Re=45 cylinder mesh')
parser.add_argument('output', nargs='?', type=Path, default=Path(__file__).resolve().parent)
out = parser.parse_args().output.resolve()
out.mkdir(parents=True, exist_ok=True)
nt, nr = 128, 96
radius, outer_radius, first_height = 0.5, 30.0, 0.006
lo, hi = 1.0, 1.2
for _ in range(80):
  q = (lo + hi) / 2
  if first_height * (q**nr - 1) / (q - 1) > outer_radius - radius:
    hi = q
  else:
    lo = q
q = (lo + hi) / 2
radii = [radius + first_height * (q**j - 1) / (q - 1) for j in range(nr + 1)]
nodes = [
  (r * math.cos(2 * math.pi * i / nt), r * math.sin(2 * math.pi * i / nt))
  for r in radii
  for i in range(nt)
]


def node(j, i):
  return j * nt + i % nt + 1


faces, cells, lookup = [], [], {}
for j in range(nr):
  for i in range(nt):
    vertices = [node(j, i), node(j + 1, i), node(j + 1, i + 1), node(j, i + 1)]
    cell_faces = []
    for k in range(4):
      a, b = vertices[k], vertices[(k + 1) % 4]
      key = tuple(sorted((a, b)))
      if key not in lookup:
        lookup[key] = len(faces) + 1
        faces.append([a, b, len(cells) + 1, 0])
      else:
        faces[lookup[key] - 1][3] = len(cells) + 1
      cell_faces.append(lookup[key])
    cells.append(cell_faces)
groups = {'interior': [], 'cylinder': [], 'farfield': []}
for i, face in enumerate(faces, 1):
  a, b, c, d = face
  group = 'interior' if d else ('cylinder' if a <= nt and b <= nt else 'farfield')
  groups[group].append((i, *face))
with (out / 'mesh.txt').open('w') as f:
  f.write(
    f'{len(nodes)} {len(faces)} {len(cells)} 3\ninterior = INTER\ncylinder = WALL\nfarfield = FAR\n(node)\n'
  )
  for i, (x, y) in enumerate(nodes, 1):
    f.write(f'{i} {x:.16e} {y:.16e}\n')
  f.write('(end)\n(edge)\n')
  for group, rows in groups.items():
    f.write(group + '\n')
    for row in rows:
      f.write(' '.join(map(str, row)) + '\n')
    f.write('(end)\n')
  f.write('(end)\n(cell)\n')
  for i, row in enumerate(cells, 1):
    f.write(' '.join(map(str, (i, *row))) + '\n')
T, Ma, Re, D, R = 300.0, 0.2, 45.0, 1.0, 287.05
mu = 1.716e-5 * (T / 273.15) ** 1.5 * (273.15 + 110.4) / (T + 110.4)
U = Ma * math.sqrt(1.4 * R * T)
rho = Re * mu / (U * D)
p = rho * R * T
(out / 'parameters.json').write_text(
  json.dumps(
    dict(
      Ma=Ma,
      Re=Re,
      T_K=T,
      D_m=D,
      U_m_s=U,
      rho_kg_m3=rho,
      p_Pa=p,
      mu_Pa_s=mu,
      nt=nt,
      nr=nr,
      outer_radius_D=outer_radius,
      first_height_D=first_height,
      growth_ratio=q,
    ),
    indent=2,
  )
  + '\n'
)
print((out / 'parameters.json').read_text())
