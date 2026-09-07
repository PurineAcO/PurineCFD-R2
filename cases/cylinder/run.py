"""Generate and run the bundled cylinder case with an explicit physical model."""

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def main():
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('--model', choices=['laminar', 'sa'], required=True)
  parser.add_argument('--re', type=float, default=47)
  parser.add_argument('--steps', type=int, default=80000)
  parser.add_argument('--output', type=Path, required=True)
  args = parser.parse_args()
  if not math.isfinite(args.re) or args.re <= 0 or args.steps < 2:
    parser.error('Re must be finite and positive; steps must be at least 2')
  binary = ROOT / 'build/purinecfd'
  if not binary.is_file():
    parser.error('Build first: cmake -S . -B build && cmake --build build -j 4')
  target = args.output.resolve()
  if target.exists():
    parser.error(f'Output already exists; choose a new directory: {target}')
  subprocess.run(
    [sys.executable, str(Path(__file__).with_name('generate_case.py')), str(target)], check=True
  )
  parameters = json.loads((target / 'parameters.json').read_text())
  parameters['Re'] = args.re
  parameters['model'] = args.model
  parameters['converged'] = False
  parameters['rho_kg_m3'] = (
    args.re * parameters['mu_Pa_s'] / (parameters['U_m_s'] * parameters['D_m'])
  )
  parameters['p_Pa'] = parameters['rho_kg_m3'] * 287.05 * parameters['T_K']
  (target / 'parameters.json').write_text(json.dumps(parameters, indent=2) + '\n')
  config = json.loads((ROOT / 'config.json').read_text())
  config['io'] = {'mesh': 'mesh.txt', 'log': 'run.log', 'field': 'field'}
  config['solver'].update(model=args.model, max_steps=args.steps)
  config['farfield'] = {
    'Ma': parameters['Ma'],
    'T': parameters['T_K'],
    'p': parameters['p_Pa'],
  }
  (target / 'config.json').write_text(json.dumps(config, indent=2) + '\n')
  print(f'Running; progress is written to {target / "run.log"}', flush=True)
  subprocess.run([str(binary), str(target / 'config.json')], check=True)
  latest = max((target / 'field').glob('step_*.dat'), key=lambda f: int(f.stem.split('_')[1]))
  shutil.copy2(latest, target / 'baseflow.dat')
  parameters['converged'] = 'Converged at step ' in (target / 'run.log').read_text()
  (target / 'parameters.json').write_text(json.dumps(parameters, indent=2) + '\n')
  print(f'Finished: {target}; converged={parameters["converged"]}')
  if not parameters['converged']:
    print('This field is for installation checks only. Use more steps before stability analysis.')


if __name__ == '__main__':
  main()
