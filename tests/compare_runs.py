"""新旧两次 URANS 算例对比: 历程、内迭代收敛、流场云图差值.

用法:
  python tests/compare_runs.py --a run/3-cyc-... --b run/4-cyc-... \
      --config config.json --out run/7-cmp-... [--airfoil] [--labels 标量耗散 精确Jacobi]

两个目录须含 history.csv 与 step_*.dat. --config 提供网格拓扑与来流参数.
"""

import argparse
import glob
import json
import math
import os

import matplotlib

matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

GAMMA, R = 1.4, 287.05

parser = argparse.ArgumentParser()
parser.add_argument('--a', required=True)
parser.add_argument('--b', required=True)
parser.add_argument('--config', required=True)
parser.add_argument('--out', required=True)
parser.add_argument('--airfoil', action='store_true')
parser.add_argument('--labels', nargs=2, default=['A', 'B'])
args = parser.parse_args()

cfg = json.load(open(args.config))
base = os.path.dirname(os.path.abspath(args.config))
far = cfg['farfield']
U = far['Ma'] * math.sqrt(GAMMA * R * far['T'])
os.makedirs(args.out, exist_ok=True)


def hist(d):
    h = np.genfromtxt(
        os.path.join(d, 'history.csv'), delimiter=',', names=True, dtype=None, encoding=None
    )
    return h[h['phase'] == 'U']


def spectral(h):
    t, cl = h['time'], h['Cl']
    half = len(t) // 2
    s = cl[half:] - cl[half:].mean()
    freq = np.fft.rfftfreq(len(s), t[1] - t[0])
    amp = np.abs(np.fft.rfft(s * np.hanning(len(s))))
    return freq, amp, cl[half:], t, peak(s, t[1] - t[0])


def peak(y, dt):
  """补零 + 按原始分辨率做抛物线插值, 避免峰值被钉在 bin 上"""
  y = np.asarray(y, float) - np.mean(y)
  n = len(y)
  df = 1.0 / (n * dt)
  N = 1 << int(np.ceil(np.log2(n * 16)))
  mag = np.abs(np.fft.rfft(y * np.hanning(n), N))
  fr = np.fft.rfftfreq(N, dt)
  k = int(np.argmax(mag[1:])) + 1
  f0 = fr[k]

  def at(f):
    return mag[int(np.argmin(np.abs(fr - f)))]

  a, b, c = np.log(at(f0 - df)), np.log(at(f0)), np.log(at(f0 + df))
  return f0 + 0.5 * (a - c) / (a - 2 * b + c) * df


HA, HB = hist(args.a), hist(args.b)
fA, aA, cA, tA, pA = spectral(HA)
fB, aB, cB, tB, pB = spectral(HB)
tb = tA[: min(len(tA), len(tB))]


def shock(d, airfoil):
    if not airfoil:
        return None
    ts, xs = [], []
    for block in open(os.path.join(d, 'wall_cp.dat')).read().split('# step ')[1:]:
        rows = block.strip().split('\n')
        a = np.array([r.split() for r in rows[1:]], dtype=float)
        up = a[a[:, 1] > 0]
        up = up[np.argsort(up[:, 0])]
        x, cp = up[:, 0], up[:, 2]
        g = np.where((x > -0.45) & (x < 0.35), np.gradient(cp, x), -np.inf)
        ts.append(float(rows[0].split()[2]))
        xs.append(x[np.argmax(g)] + 0.5)
    return np.array(ts), np.array(xs)


SA, SB = shock(args.a, args.airfoil), shock(args.b, args.airfoil)

# ---------------- 历程对比 ----------------
n = 5 if args.airfoil else 4
fig, ax = plt.subplots(n, 1, figsize=(9.5, 2.4 * n), sharex=False)
ax[0].plot(HA['time'], HA['Cl'], lw=1, label=args.labels[0])
ax[0].plot(HB['time'], HB['Cl'], lw=1, label=args.labels[1])
ax[0].set_ylabel('Cl')
ax[0].legend(fontsize=9, ncol=2)
ax[1].plot(HA['time'], HA['Cd'], lw=1)
ax[1].plot(HB['time'], HB['Cd'], lw=1)
ax[1].set_ylabel('Cd')
ax[2].plot(fA, aA, lw=1, label=f"St={pA / U:.4f}")
ax[2].plot(fB, aB, lw=1, label=f"St={pB / U:.4f}")
ax[2].set_xlim(0, min(60.0, fA.max()))
ax[2].set_xlabel('f [Hz]')
ax[2].set_ylabel(r'$|Cl|$ spectrum')
ax[2].legend(fontsize=8)
ax[3].semilogy(tb, HA['res'][: len(tb)] / HA['res0'][: len(tb)], lw=1)
ax[3].semilogy(tb, HB['res'][: len(tb)] / HB['res0'][: len(tb)], lw=1)
ax[3].set_ylabel(r'$res/res_0$')
ax[3].set_xlim(0.0, float(max(HA['time'][-1], HB['time'][-1])))
ax[3].axhline(1e-3, color='k', lw=0.6, ls='--')
if args.airfoil and SA is not None:
    ax[4].plot(SA[0], SA[1], lw=1)
    ax[4].plot(SB[0], SB[1], lw=1)
    ax[4].set_ylabel('shock x/c')
    ax[4].set_xlim(0.0, float(max(HA['time'][-1], HB['time'][-1])))
ax[-1].set_xlabel('t [s]')
fig.tight_layout()
fig.savefig(os.path.join(args.out, 'history_cmp.png'), dpi=130)
plt.close(fig)

# ---------------- 汇总 ----------------
lines = [f'{"指标":<22}{args.labels[0]:>16}{args.labels[1]:>16}']
rows = [
    ('物理步数', len(HA), len(HB)),
    ('平均内迭代次数', HA['inner'].mean(), HB['inner'].mean()),
    ('内迭代残差下降中位数', np.median(HA['res'] / HA['res0']), np.median(HB['res'] / HB['res0'])),
    ('Cl 均值', cA.mean(), cB.mean()),
    ('Cl 幅值', cA.max() - cA.min(), cB.max() - cB.min()),
    ('Cd 均值', HA['Cd'][len(HA) // 2 :].mean(), HB['Cd'][len(HB) // 2 :].mean()),
    ('St', pA / U, pB / U),
]
if SA is not None:
    rows.append(('激波 x/c 均值', SA[1][len(SA[1]) // 2 :].mean(), SB[1][len(SB[1]) // 2 :].mean()))
for name, va, vb in rows:
    if isinstance(va, (int, np.integer)):
        lines.append(f'{name:<22}{va:>16d}{vb:>16d}')
    else:
        lines.append(f'{name:<22}{va:>16.4f}{vb:>16.4f}')
open(os.path.join(args.out, 'summary_cmp.txt'), 'w', encoding='utf-8').write('\n'.join(lines) + '\n')
print('\n'.join(lines))

# ---------------- 流场云图 (末帧 A | B | A-B) ----------------
S = int(open(os.path.join(base, cfg['io']['structured'])).readline().split()[0])


def load(d):
    fn = sorted(glob.glob(os.path.join(d, 'step_*.dat')))[-1]
    f = np.loadtxt(fn, skiprows=2)
    keys = ['x', 'y', 'rho', 'u', 'v', 'T', 'p', 'Ma', 'miubl']
    a = {k: f[:, i].reshape(-1, S) for i, k in enumerate(keys)}

    def ds(v):
        return 0.5 * (np.roll(v, -1, axis=1) - np.roll(v, 1, axis=1))

    def dn(v):
        return np.gradient(v, axis=0)

    rho, x, y = a['rho'], a['x'], a['y']
    jac = ds(x) * dn(y) - dn(x) * ds(y)
    gx = (ds(rho) * dn(y) - dn(rho) * ds(y)) / jac
    gy = (dn(rho) * ds(x) - ds(rho) * dn(x)) / jac
    grad = np.hypot(gx, gy)
    a['schlieren'] = np.exp(-15 * grad / np.percentile(grad, 99.5))
    vx = (ds(a['v']) * dn(y) - dn(a['v']) * ds(y)) / jac
    uy = (dn(a['u']) * ds(x) - ds(a['u']) * dn(x)) / jac
    a['vorticity'] = (vx - uy) / U
    return {k: np.hstack([v, v[:, :1]]) for k, v in a.items()}, fn


AA, fnA = load(args.a)
AB, fnB = load(args.b)
view = (-0.2, 1.4, -0.5, 0.6) if args.airfoil else (-1.0, 5.0, -3.0, 3.0)
fields = [('Ma', 'jet', None, 'Ma'), ('vorticity', 'RdBu_r', 4.0, r'$\omega D/U$')]
fig, axs = plt.subplots(len(fields), 3, figsize=(15, 4.0 * len(fields)))
for row, (key, cmap, lim, bar) in zip(axs, fields):
    za, zb = AA[key], AB[key]
    if key == 'Ma':
        lo = float(np.floor(min(za.min(), zb.min()) * 20) / 20)
        hi = float(np.ceil(max(za.max(), zb.max()) * 20) / 20)
        levels = np.linspace(lo, hi, 25)
    else:
        if lim is None:
            lim = float(np.ceil(max(np.abs(za).max(), np.abs(zb).max()) * 20) / 20)
        levels = np.linspace(-lim, lim, 25)
    dv = zb - za
    dlim = np.percentile(np.abs(dv), 99.5) or 1.0
    for k, (z, ttl, cm, lv) in enumerate(
        [(za, args.labels[0], cmap, levels), (zb, args.labels[1], cmap, levels),
         (dv, 'B - A', 'coolwarm', np.linspace(-dlim, dlim, 25))]
    ):
        c = row[k].contourf(AA['x'], AA['y'], z, levels=lv, cmap=cm, extend='both')
        row[k].fill(AA['x'][0], AA['y'][0], color='w', edgecolor='k', lw=0.8, zorder=5)
        row[k].set_xlim(view[0], view[1])
        row[k].set_ylim(view[2], view[3])
        row[k].set_aspect('equal')
        row[k].set_title(f'{bar}  {ttl}', fontsize=10)
        fig.colorbar(c, ax=row[k], orientation='horizontal', fraction=0.05, pad=0.1)
    row[0].set_ylabel('y/c' if args.airfoil else 'y/D')
fig.tight_layout()
fig.savefig(os.path.join(args.out, 'field_cmp.png'), dpi=150)
plt.close(fig)
print('写出', os.path.join(args.out, 'history_cmp.png'), os.path.join(args.out, 'field_cmp.png'))
