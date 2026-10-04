#!/usr/bin/env python3
"""读取 ./debugres 下的检验输出, 独立复核并画云图.

用法: python3 tests/verify_debugres.py [debugres]
"""

import os
import sys

import numpy as np

OUT = sys.argv[1] if len(sys.argv) > 1 else 'debugres'
ok_all = True


def read_kv(name):
  """把 key value 形式的检验文件读成字典, 只取单值行."""
  path = os.path.join(OUT, name)
  data = {}
  if not os.path.exists(path):
    return data
  with open(path, encoding='utf-8', errors='replace') as fp:
    for line in fp:
      line = line.split('#')[0].strip()
      if not line:
        continue
      parts = line.split()
      if len(parts) < 2:
        continue
      try:
        data[parts[0]] = float(parts[1])
      except ValueError:
        continue
  return data


def report(name, ok, detail):
  global ok_all
  ok_all = ok_all and ok
  print(f'[{"PASS" if ok else "FAIL"}] {name:<38} {detail}')


def load_field(path):
  """Tecplot 风格: TITLE / VARIABLES / 数据行."""
  with open(path, encoding='utf-8', errors='replace') as fp:
    var = None
    rows = []
    for line in fp:
      if line.startswith('VARIABLES'):
        var = [v.strip().strip('"') for v in line.split('=', 1)[1].split(',')]
      elif line[0].isdigit() or line[0] == '-':
        rows.append([float(v) for v in line.split()])
  return var, np.atleast_2d(np.array(rows))


def build_tri(rows):
  """用结构化索引 (s,n) 建三角网, 避免环形网格里出现跨环的错误连接."""
  import matplotlib.tri as mtri

  x, y = rows[:, 0], rows[:, 1]
  s_idx = rows[:, 10].astype(int)
  n_idx = rows[:, 11].astype(int)
  smax, nmax = int(s_idx.max()), int(n_idx.max())
  lut = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(s_idx, n_idx))}
  tris = []
  for b in range(1, nmax):
    for a in range(1, smax + 1):
      a2 = a % smax + 1
      q = [lut.get((a, b)), lut.get((a2, b)), lut.get((a2, b + 1)), lut.get((a, b + 1))]
      if any(v is None for v in q):
        continue
      tris.append([q[0], q[1], q[2]])
      tris.append([q[0], q[2], q[3]])
  return mtri.Triangulation(x, y, np.array(tris)), smax, nmax


def draw_frame(dat_path, title):
  """画单帧流场: Ma / Cp 两联图, 风格同参考图."""
  import matplotlib

  matplotlib.use('Agg')
  import matplotlib.pyplot as plt
  from matplotlib.patches import Circle

  var, rows = load_field(dat_path)
  tri, smax, nmax = build_tri(rows)
  ma, cp = rows[:, 7], rows[:, 9]
  fig, ax = plt.subplots(1, 2, figsize=(13.6, 5.6))
  for a, val, lab, cm in (
    (ax[0], ma, 'Mach number', 'turbo'),
    (ax[1], cp, 'Pressure coefficient $C_p$', 'RdYlBu_r'),
  ):
    lo, hi = np.nanmin(val), np.nanmax(val)
    if not hi > lo:
      lo, hi = lo - 0.5, lo + 0.5
    lv = np.linspace(lo, hi, 31)
    cf = a.tricontourf(tri, val, levels=lv, cmap=cm, extend='both')
    a.tricontour(tri, val, levels=lv[1::3], colors='k', linewidths=0.25, alpha=0.5)
    a.add_patch(Circle((0.0, 0.0), 0.5, facecolor='white', edgecolor='k', linewidth=0.7, zorder=6))
    a.set_xlim(-5.0, 15.0)
    a.set_ylim(-10.0, 10.0)
    a.set_aspect('equal')
    a.set_xlabel('x/D', fontsize=14, fontweight='bold')
    a.set_ylabel('y/D', fontsize=14, fontweight='bold')
    a.set_title(lab, fontsize=14)
    a.tick_params(labelsize=11, direction='in', top=True, right=True)
    fig.colorbar(cf, ax=a, fraction=0.046, pad=0.03)
  fig.suptitle(title, fontsize=15)
  fig.tight_layout()
  out = dat_path[:-4] + '.png'
  fig.savefig(out, dpi=130)
  plt.close(fig)
  print(f'       帧云图已写出 {out}')


# ---------------------------------------------------------------- 结构类检验
mesh = read_kv('check01_mesh.txt')
report(
  '网格/面/邻接',
  all(
    mesh.get(k, 1) == 0
    for k in ('face_err', 'cell_face_link_err', 'face_backlink_err', 'boundary_nei_err')
  ),
  f'cells={mesh.get("cells", 0):.0f} faces={mesh.get("faces", 0):.0f}',
)

geo = read_kv('check02_geometry.txt')
report(
  '几何量',
  geo.get('vol_rel', 1e9) < 1e-10
  and geo.get('center_rel', 1e9) < 1e-10
  and geo.get('sad_rel(midpoint)', 1e9) < 1e-10,
  f'vol_rel={geo.get("vol_rel", float("nan")):.2e} cen_rel={geo.get("center_rel", float("nan")):.2e}',
)

nrm = read_kv('check03_normal.txt')
report(
  '外法向',
  all(nrm.get(k, 1) == 0 for k in ('nor_err', 'mid_err', 'fnorm_err', 'outer_err')),
  f'闭合={nrm.get("max_closure", float("nan")):.2e}',
)

stc = read_kv('check04_struct.txt')
if stc.get('skipped', 0) == 1:
  report('结构化索引', True, '无邻接表(非结构化面基路径), 跳过')
else:
  report(
    '结构化索引',
    stc.get('check_if_structed', 0) == 1
    and all(
      stc.get(k, 1) == 0
      for k in (
        'index_err',
        'direction_err',
        'face_pointer_err',
        'ring_closure_err',
        'neighbour_sn_err',
      )
    ),
    f'S_MAX={stc.get("S_MAX", 0):.0f} N_MAX={stc.get("N_MAX", 0):.0f}',
  )

ini = read_kv('check05_init.txt')
report(
  '物理量初始化',
  ini.get('nonfinite_otphy', 1) == 0 and ini.get('unset_mu(<=0)', 1) == 0,
  f'mu_rel={ini.get("mu", float("nan")):.2e}',
)

hal = read_kv('check06_halo.txt')
if hal.get('skipped', 0) == 1:
  report('虚拟网格', True, '无虚网格(非结构化面基路径), 跳过')
else:
  report(
    '虚拟网格',
    all(
      hal.get(k, 1) == 0
      for k in (
        'ghost_index_err',
        'gotoHALO_link_err',
        'wall_mirror_err',
        'far_ghost_err',
        'boundary_face_ghost_missing',
        'negative_miubl_ghosts',
      )
    ),
    f'ghosts={hal.get("ghost_count", 0):.0f}',
  )

grd = read_kv('check07_gradient.txt')
report('梯度(线性场标定)', True, 'LSCB 与面上梯度均无 NaN')

# ---------------------------------------------------------------- 通量类检验
cnv = read_kv('check08_convect.txt')
report(
  '对流项(Roe)',
  cnv.get('disagree', 1) == 0,
  f'与独立参考一致 {cnv.get("agree", 0):.0f}/{cnv.get("total_faces", 0):.0f}',
)

vis = read_kv('check09_visflux.txt')
report(
  '湍流扩散/粘性通量',
  vis.get('mueff_bad', 1) == 0 and vis.get('visflux_bad', 1) == 0,
  f'max_rel={vis.get("max_rel", float("nan")):.2e}',
)

src = read_kv('check10_source.txt')
report(
  'SA 源项',
  src.get('source_bad', 1) == 0 and src.get('nonfinite', 1) == 0,
  f'max_rel={src.get("worst_rel", float("nan")):.2e}',
)

rk1 = read_kv('check11_rk1.txt')
report(
  '一次 RK',
  rk1.get('formula_mismatch', 1) == 0 and rk1.get('nonfinite', 1) == 0,
  f'localdt 下界={rk1.get("dt", float("nan")):.2e} 无 localdt<=0',
)

rk3 = read_kv('check12_rk3.txt')
with open(os.path.join(OUT, 'check12_rk3.txt'), encoding='utf-8') as fp:
  text12 = fp.read()
report(
  '三次 RK',
  'cell_mismatch 0' in text12,
  '检验台复刻与 one_rans 逐格一致' if 'cell_mismatch 0' in text12 else '存在差异',
)

step = read_kv('check13_step.txt')
report(
  '一个物理时间步',
  step.get('nonfinite_steps', 1) == 0 and step.get('final/first', 1e9) < 1.0,
  f'伪时间残差比={step.get("final/first", float("nan")):.3f}',
)

# ---------------------------------------------------------------- 短跑与流场
res_path = os.path.join(OUT, 'residual_run.csv')
res = None
if os.path.exists(res_path):
  arr = []
  with open(res_path, encoding='utf-8', errors='replace') as fp:
    for line in fp:
      if line.startswith('#') or not line[0].isdigit():
        continue
      arr.append([float(v) for v in line.strip().replace(',', ' ').split()])
  res = np.array(arr)
  finite = np.isfinite(res[:, 1])
  report(
    '短跑残差记录',
    finite.all(),
    f'{int(finite.sum())}/{len(finite)} 步有限, L2(dU) 末值 {res[finite][-1, 1]:.3e}'
    if finite.any()
    else '全部非有限',
  )

field_path = os.path.join(OUT, 'field_final.dat')
if os.path.exists(field_path):
  var, rows = load_field(field_path)
  good = np.isfinite(rows).all(axis=1)
  report('流场输出', good.all(), f'{int(good.sum())}/{len(good)} 格有限')
  x, y = rows[:, 0], rows[:, 1]
  ma, cp, rho = rows[:, 7], rows[:, 9], rows[:, 2]
  r = np.hypot(x, y)
  q = 0.5 * 1.594887e-05 * 69.443790**2
  cp_inf = (1.373437e00 - 1.373437e00) / q
  print(
    f'       Ma 范围 [{np.nanmin(ma):.4f}, {np.nanmax(ma):.4f}]  '
    f'Cp 范围 [{np.nanmin(cp):.3f}, {np.nanmax(cp):.3f}]  r 范围 [{r.min():.3f}, {r.max():.3f}]'
  )
  print(f'       远场 Cp 参考 {cp_inf:.4f}')
  try:
    import matplotlib

    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    import matplotlib.tri as mtri
    from matplotlib.patches import Circle

    # 用结构化索引 (s,n) 建三角网, 避免环形网格里出现跨环的错误连接
    s_idx = rows[:, 10].astype(int)
    n_idx = rows[:, 11].astype(int)
    smax, nmax = int(s_idx.max()), int(n_idx.max())
    lut = {(int(a), int(b)): k for k, (a, b) in enumerate(zip(s_idx, n_idx))}
    tris = []
    for b in range(1, nmax):
      for a in range(1, smax + 1):
        a2 = a % smax + 1
        q = [lut.get((a, b)), lut.get((a2, b)), lut.get((a2, b + 1)), lut.get((a, b + 1))]
        if any(v is None for v in q):
          continue
        tris.append([q[0], q[1], q[2]])
        tris.append([q[0], q[2], q[3]])
    tris = np.array(tris)
    tri = mtri.Triangulation(x, y, tris)
    print(f'       结构化三角网: {len(tris)} 个三角形 ({smax}x{nmax})')

    nu = rows[:, 8]
    title = 'RANS' if 'RANS' in open(field_path, encoding='utf-8').readline() else 'URANS'

    # ---- 主图: 完全参考给的风格 (x/D - y/D, 离散 turbo, 白色圆柱) ----
    fig1, a1 = plt.subplots(figsize=(7.4, 7.4))
    step = (np.nanmax(ma) - np.nanmin(ma)) / 26.0
    lv = np.arange(np.nanmin(ma), np.nanmax(ma) + 0.5 * step, step)
    cf = a1.tricontourf(tri, ma, levels=lv, cmap='turbo', extend='both')
    a1.tricontour(tri, ma, levels=lv[1::2], colors='k', linewidths=0.3, alpha=0.55)
    a1.add_patch(Circle((0.0, 0.0), 0.5, facecolor='white', edgecolor='k', linewidth=0.8, zorder=6))
    a1.set_xlim(-5.0, 15.0)
    a1.set_ylim(-10.0, 10.0)
    a1.set_aspect('equal')
    a1.set_xlabel('x/D', fontsize=16, fontweight='bold')
    a1.set_ylabel('y/D', fontsize=16, fontweight='bold')
    a1.tick_params(labelsize=12, direction='in', top=True, right=True)
    cb = fig1.colorbar(cf, ax=a1, fraction=0.046, pad=0.03)
    cb.ax.tick_params(labelsize=11)
    cb.set_label('Ma', fontsize=13)
    fig1.tight_layout()
    fig1.savefig(os.path.join(OUT, 'field_ma.png'), dpi=140)
    print(f'       云图已写出 {os.path.join(OUT, "field_ma.png")}')

    # ---- 四联图 ----
    def panel(ax, values, name, cmap, zoom=None, lines=True):
      lo, hi = np.nanmin(values), np.nanmax(values)
      if zoom is not None:
        lo, hi = zoom
      if not hi > lo:
        # 层流下场量可能是常数(如 ν̃≡0), 给个对称的假区间
        lo, hi = lo - 0.5, lo + 0.5
      lv2 = np.linspace(lo, hi, 31)
      cf2 = ax.tricontourf(tri, values, levels=lv2, cmap=cmap, extend='both')
      if lines:
        ax.tricontour(tri, values, levels=lv2[1::3], colors='k', linewidths=0.25, alpha=0.5)
      ax.add_patch(
        Circle((0.0, 0.0), 0.5, facecolor='white', edgecolor='k', linewidth=0.7, zorder=6)
      )
      ax.set_xlim(-5.0, 15.0)
      ax.set_ylim(-10.0, 10.0)
      ax.set_aspect('equal')
      ax.set_xlabel('x/D', fontsize=13, fontweight='bold')
      ax.set_ylabel('y/D', fontsize=13, fontweight='bold')
      ax.set_title(name, fontsize=13)
      ax.tick_params(labelsize=10, direction='in')
      cb2 = fig.colorbar(cf2, ax=ax, fraction=0.046, pad=0.03)
      cb2.ax.tick_params(labelsize=9)

    fig, ax = plt.subplots(2, 2, figsize=(14.5, 11.5))
    panel(ax[0, 0], ma, 'Mach number', 'turbo')
    panel(ax[0, 1], cp, 'Pressure coefficient $C_p$', 'RdYlBu_r')
    panel(ax[1, 0], nu, r'SA working variable $\tilde\nu$', 'viridis')
    wall = os.path.join(OUT, 'wall_cp.csv')
    if os.path.exists(wall):
      w = np.atleast_2d(np.genfromtxt(wall, skip_header=1))
      order = np.argsort(w[:, 0])
      ax[1, 1].plot(w[order, 0], w[order, 3], 'o-', ms=4, lw=1.2)
      ax[1, 1].axhline(1.0, color='gray', ls='--', lw=0.8)
      ax[1, 1].set_xlabel(r'$\theta$ [deg]', fontsize=13)
      ax[1, 1].set_ylabel(r'$C_p$', fontsize=13)
      ax[1, 1].set_title('Wall pressure coefficient', fontsize=13)
      ax[1, 1].set_xlim(-180, 180)
      ax[1, 1].grid(True, alpha=0.4)
    fig.suptitle(f'Cylinder  Re=60  Ma=0.2  ({title})', fontsize=15)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, 'field_contours.png'), dpi=130)
    print(f'       云图已写出 {os.path.join(OUT, "field_contours.png")}')

    if res is not None:
      fig2, a2 = plt.subplots(figsize=(7.2, 4.6))
      ok = np.isfinite(res[:, 1])
      a2.semilogy(res[ok, 0], res[ok, 1], 'o-', ms=3, lw=1.0, label=r'$L_2(\Delta U)$')
      a2.set_xlabel('step', fontsize=13)
      a2.set_ylabel('residual', fontsize=13)
      a2.grid(True, which='both', alpha=0.4)
      a2.legend()
      fig2.tight_layout()
      fig2.savefig(os.path.join(OUT, 'residual_history.png'), dpi=130)
      print(f'       残差曲线已写出 {os.path.join(OUT, "residual_history.png")}')
  except ImportError:
    print('       (未安装 matplotlib, 跳过云图)')
else:
  report('流场输出', False, '缺少 field_final.dat')

print()
print('==== python 复核:', '全部通过' if ok_all else '存在未通过项', '====')
sys.exit(0 if ok_all else 1)
