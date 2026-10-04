#!/usr/bin/env python3
"""近壁剖面 vs Thwaites 层流边界层解析预测.

用法: python3 tests/bl_profile.py debugres/_prev/field_final.dat
"""

import math
import os
import sys

import numpy as np

U = 69.44379
NU = 1.1573
R = 0.5
S_MAX, N_MAX = 192, 128
BAND = 1.0  # 积分带宽 y/D, 足够包含厚边界层


def load(path):
  rows = []
  with open(path, encoding='utf-8', errors='replace') as fp:
    for line in fp:
      if line.startswith(('VARIABLES', 'TITLE')):
        continue
      p = line.split()
      if len(p) >= 12:
        try:
          rows.append([float(v) for v in p])
        except ValueError:
          pass
  return np.array(rows)


def radial(rows, s):
  sel = rows[rows[:, 10] == s]
  sel = sel[np.argsort(sel[:, 11])]
  return sel


def thwaites(npts=400):
  """Thwaites 方法: 绕圆柱前缘的层流动量厚度 theta(phi)."""
  phi = np.linspace(1e-4, math.pi, npts)
  x = R * phi
  ue = 2.0 * U * np.sin(phi)
  due_dx = 2.0 * U * np.cos(phi) / R
  integ = np.concatenate([[0.0], np.cumsum(0.5 * (ue[1:] ** 5 + ue[:-1] ** 5) * np.diff(x))])
  theta = np.sqrt(np.maximum(0.45 * NU * integ / np.maximum(ue, 1e-12) ** 6, 0.0))
  lam = theta**2 / NU * due_dx
  h = 2.61 - 3.75 * lam + 5.24 * lam**2  # Thwaites 标准拟合
  h = np.where(lam < -0.09, 2.99, h)
  return phi, theta, lam, h, ue


def main():
  path = sys.argv[1] if len(sys.argv) > 1 else 'debugres/_prev/field_final.dat'
  rows = load(path)
  phi_t, th_t, lam_t, h_t, ue_t = thwaites()
  sep = lam_t < -0.09
  iphi_sep = math.degrees(phi_t[np.argmax(sep)]) if sep.any() else None
  print(f'file {path}   cells {len(rows)}   r_max {np.hypot(rows[:, 0], rows[:, 1]).max():.1f} m')
  print('Thwaites analytic (Re=60 cylinder):')
  for d in (30, 60, 90, 103):
    i = int(np.argmin(abs(np.degrees(phi_t) - d)))
    print(
      f'   phi={d:3d}deg  theta_m={th_t[i] * 1e3:7.2f} mm  H={h_t[i]:5.2f} '
      f' delta*={h_t[i] * th_t[i] * 1e3:7.2f} mm  lambda={lam_t[i]:+.4f}'
    )
  print(f'   predicted separation phi = {iphi_sep:.1f}deg')

  rec = []
  for s in range(1, S_MAX + 1):
    prof = rows[rows[:, 10] == s]
    if len(prof) < 20:
      continue
    prof = prof[np.argsort(prof[:, 11])]
    ang = math.atan2(prof[0, 1], prof[0, 0])
    phi = abs(math.pi - abs(ang))
    ue = 2.0 * U * math.sin(phi)
    if ue < 0.15 * U:
      continue
    y = np.hypot(prof[:, 0], prof[:, 1]) - R
    u = prof[:, 3]
    keep = y <= BAND
    yy, uu = y[keep], u[keep]
    ds = float(np.sum((1.0 - uu / ue) * np.gradient(yy)))
    rec.append((math.degrees(phi), ue / U, ds, prof[0, 9], uu[0] / ue))

  rec.sort()
  arr = np.array(rec)
  print(
    f'\n{"phi":>6} {"Ue/U":>6} {"d*_CFD[mm]":>11} {"d*_Thw[mm]":>11} {"ratio":>7} {"Cp_wall":>9} {"u1/Ue":>7}'
  )
  for d in (15, 30, 45, 60, 75, 90, 105, 120, 135, 150):
    i = int(np.argmin(abs(arr[:, 0] - d)))
    ph = math.radians(arr[i, 0])
    tt = float(np.interp(ph, phi_t, h_t * th_t))
    print(
      f'{arr[i, 0]:6.1f} {arr[i, 1]:6.3f} {arr[i, 2] * 1e3:11.3f} {tt * 1e3:11.3f} '
      f'{arr[i, 2] / tt if tt > 0 else 0:7.2f} {arr[i, 3]:+9.4f} {arr[i, 4]:7.3f}'
    )

  try:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(1, 3, figsize=(18.5, 5.4))
    ax[0].plot(arr[:, 0], arr[:, 2] * 1e3, 'o-', ms=4, label='CFD')
    ax[0].plot(np.degrees(phi_t), h_t * th_t * 1e3, '-', lw=2, label='Thwaites')
    if iphi_sep:
      ax[0].axvline(iphi_sep, color='r', ls=':', label=f'Thwaites sep {iphi_sep:.0f}deg')
    ax[0].set_xlabel(r'$\varphi$ [deg]  (0 = front stagnation)', fontsize=12)
    ax[0].set_ylabel(r'$\delta^*$ [mm]  (ref = potential $U_e$)', fontsize=12)
    ax[0].set_title('boundary layer thickness', fontsize=13)
    ax[0].grid(alpha=0.4)
    ax[0].legend(fontsize=10)

    ax[1].plot(arr[:, 0], arr[:, 3], 'o-', ms=4, color='tab:red', label='CFD')
    ax[1].plot(
      np.degrees(phi_t),
      1 - (2 * np.sin(phi_t)) ** 2,
      '-',
      color='tab:blue',
      label=r'potential $1-(2\sin\varphi)^2$',
    )
    ax[1].axhline(0.0, color='k', lw=0.6)
    if iphi_sep:
      ax[1].axvline(iphi_sep, color='r', ls=':')
    ax[1].set_ylim(-3.2, 1.4)
    ax[1].set_xlabel(r'$\varphi$ [deg]', fontsize=12)
    ax[1].set_ylabel(r'$C_p$ on wall', fontsize=12)
    ax[1].set_title('wall pressure: no separation plateau', fontsize=13)
    ax[1].grid(alpha=0.4)
    ax[1].legend(fontsize=10)

    for d, c in ((45, 'tab:blue'), (75, 'tab:green'), (90, 'tab:red')):
      i = int(np.argmin(abs(arr[:, 0] - d)))
      phi = math.radians(arr[i, 0])
      s = int(round((math.pi - phi) / (2 * math.pi) * S_MAX - 0.5)) % S_MAX + 1
      prof = rows[rows[:, 10] == s]
      prof = prof[np.argsort(prof[:, 11])]
      y = np.hypot(prof[:, 0], prof[:, 1]) - R
      ue = 2.0 * U * math.sin(phi)
      ax[2].plot(prof[:, 3] / ue, y, '-', color=c, lw=2, label=f'CFD {arr[i, 0]:.0f}deg')
      tt = float(np.interp(phi, phi_t, h_t * th_t))
      ax[2].axhline(tt, color=c, ls=':', lw=1.2)
    ax[2].axvline(1.0, color='k', lw=0.8)
    ax[2].set_xlim(0, 1.05)
    ax[2].set_ylim(0, 1.6)
    ax[2].set_xlabel(r'$u/U_e^{pot}$', fontsize=12)
    ax[2].set_ylabel('y/D', fontsize=12)
    ax[2].set_title(r'profiles (dotted = Thwaites $\delta^*$)', fontsize=13)
    ax[2].grid(alpha=0.4)
    ax[2].legend(fontsize=10)

    fig.tight_layout()
    out_png = os.path.join(os.path.dirname(path) or '.', 'bl_profile.png')
    fig.savefig(out_png, dpi=135)
    print(f'\nfigure -> {out_png}')
  except ImportError:
    print('(matplotlib missing)')


if __name__ == '__main__':
  main()
