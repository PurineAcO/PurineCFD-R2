"""重新生成 tests/reference_200.json。

物理模型或求解器改动后, 回归基准需要重跑。用法:

    python3 tests/make_reference.py [工作目录]
"""

import json
import os
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PROBES = [0, 16, 32, 48, 64, 96, 127, 128, 512, 1024, 2048, 4096, 8192, 12287]
COLUMNS = ['x', 'y', 'rho', 'u', 'v', 'T', 'p', 'Ma', 'nu_tilde']
DESCRIPTION = '200 iterations of the steady SA-RANS solver with the C5 compressibility term, 128 x 96 cylinder mesh, CFL=1'


def main():
  work = Path(sys.argv[1]).resolve() if len(sys.argv) > 1 else ROOT / 'build' / 'reference-work'
  work.mkdir(parents=True, exist_ok=True)
  subprocess.run(
    [sys.executable, str(ROOT / 'cases/cylinder/generate_case.py'), str(work)], check=True
  )
  config = json.loads((ROOT / 'config.json').read_text())
  config['io'] = {
    'mesh': str(work / 'mesh.txt'),
    'log': 'run.log',
    'field': 'field',
  }
  config['solver']['max_steps'] = 200
  (work / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
  environment = dict(os.environ)
  environment['OMP_NUM_THREADS'] = '1'
  subprocess.run(
    [str(ROOT / 'build' / 'purinecfd'), str(work / 'config.json')],
    cwd=work,
    env=environment,
    check=True,
  )
  rows = [
    [float(value) for value in line.split()]
    for line in (work / 'field' / 'step_000200.dat').read_text().splitlines()[2:]
  ]
  reference = {
    'description': DESCRIPTION,
    'columns': COLUMNS,
    'cells': len(rows),
    'probes': {str(index): rows[index] for index in PROBES},
    'minimum': [min(row[k] for row in rows) for k in range(9)],
    'maximum': [max(row[k] for row in rows) for k in range(9)],
  }
  (ROOT / 'tests/reference_200.json').write_text(json.dumps(reference, indent=2) + '\n')
  print(f'{len(rows)} 单元, 已写入 tests/reference_200.json')


if __name__ == '__main__':
  main()
