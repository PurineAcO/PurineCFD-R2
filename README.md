# PurineCFD-R2

二维有限体积 URANS 求解器，当前版本2.3.3

- **2.3.3** 增强了求解器的鲁棒性，即使用包括边界在内完全解析的Jacobi进行迭代，还增加了自适应CFL稳定，抑制了远场发散。改版还调整了求解器源码文件的组织。
- **2.3.1** 起，双时间步的隐式求解使用**完整的 4×4 块 Jacobi**（[dualtime_full.hpp](src/dualtime_full.hpp)）：对角块与相邻块由对流 Jacobi $A_n$ 、Roe 耗散矩阵 $\hat{A}_n$ 、黏性梯度响应 $T_f$ 、冻结梯度 $B_f$ 组装，不再用单个标量谱半径近似。编译期由 `-DPURINE_JAC_FULL` 选择（对应 CMake 选项 `PURINE_JAC_FULL`，默认 ON）；关掉则退回 2.3.0 的标量版本 [dualtime.hpp](src/dualtime.hpp)。

## 构建与运行

### 依赖

| 用途 | 依赖 | 声明处 |
| --- | --- | --- |
| 编译求解器 | C++17 编译器 + OpenMP、CMake ≥ 3.16 | 系统包（Linux 用 GCC；Windows 推荐 MSYS2 的 MinGW-w64 GCC） |
| JSON 解析 | 无 | nlohmann/json 3.11.3 已随源码在 [third_party/](third_party/) |
| 后处理 | Python ≥ 3.10 + numpy、matplotlib、scipy | [pyproject.toml](pyproject.toml) 的 `post` 组 |
| 涡量动画 | pillow | `anim` 组 |
| 图形界面 | Python ≥ 3.10 + PySide6 | `gui` 组 |
| 代码检查（可选） | ruff | `dev` 组 |

四个组都在 `[tool.uv] default-groups` 里，所以 `uv sync` 一条命令就能装齐全部（见「Python 环境」）。
只想编译求解器、不跑后处理与 GUI 的话，Python 部分整个可以跳过。

### 编译

提供了快速编译脚本：

```sh
bash fastbuild.sh
```

CMake 选项：

| 选项 | 默认 | 作用 |
| --- | --- | --- |
| `PURINE_JAC_FULL` | ON | 定义 `PURINE_JAC_FULL`，启用完整 4×4 块 Jacobi |
| `PURINE_NATIVE` | OFF | 加 `-march=native` |
| `PURINE_SANITIZE` | OFF | 加 ASan/UBSan（另有一个 `debug` 目标编 [tests/debug.cpp](tests/debug.cpp) 检验台） |

### 使用CLI运行

根目录的 [config.json](config.json) 已指向圆柱算例，输入以下命令运行：

```sh
python3 cases/cylinder/generate_case.py
OMP_NUM_THREADS=32 ./build/purinecfd config.json
```

输出写在 `io.field` 指向的目录（这里是 `run/field/`）：`history.csv` 每步一行，`step_*.dat` 每
`dump_interval` 步一个，配了 `urans.wall_interval` 时还有 `wall_cp.dat`；屏幕输出进 `run/run.log`。

注意 `io.field` 是输出目录，重复跑同一个 `io.field` 会覆盖上一次的结果

### 使用GUI运行

当你确认准备好Python环境后，执行以下命令：

```sh
./gui.sh
```

GUI的使用详见 [gui/README.md](gui/README.md) 与下文「图形界面」一节。

### Python 环境

在仓库根目录依次执行：

```sh

curl -LsSf https://astral.sh/uv/install.sh | sh
# Windows / PowerShell 用: powershell -ExecutionPolicy ByPass -c "irm https://astral.sh/uv/install.ps1 | iex"
exec $SHELL
uv --version
cd /path/to/PurineCFD-R2
uv sync
uv run python -c "import numpy, matplotlib, scipy, PIL, PySide6; print('ok')"
```

使用求解器后处理时，需激活Python环境。

WSL 用户需要 WSLg（`echo $DISPLAY` 有输出即可）。WSL 发行版默认不带中文字体，
[gui/main.py](gui/main.py) 会自动从 `/mnt/c/Windows/Fonts/` 寻找；找不到时会提示
`sudo apt install fonts-noto-cjk`。

### 从中间流场启动

`io` 里加一个可选的 `restart`，指向任意一次运行写出的 `step_*.dat` 即可接着那个流场算：

```json
"io": {
  ...
  "restart": "../../run/25-oat-a37-20261016-1632/step_001500.dat"
}
```

中断与发散两种情况都会把当前这一物理步的流场落进 `io.field`，文件名为 `step_%06d.dat`，当受到坏单元时，退出码为1，正常打断时退出码为0

**目前该技术尚未成熟，舍入误差控制的不好，不建议使用**

## 算例

三个算例的完整记录都在 `run/` 下。以下显示了在仓库开发者机上运行结果

### 圆柱绕流（[config.json](config.json)）


- 工况： Ma=0.2， $\mathrm{Re}=\rho UD/\mu=1000$ ，T=300 K，α=0°，来流沿 +x 
- 物性： Sutherland 黏度；由 Re 反算 p∞≈22.89061 Pa，U∞≈69.44379 m/s 
- 网格： D=1 m，192×128 单元，远场半径 30D，第一层高度 0.0015D 
- 边界： 无滑移绝热壁面（ $\tilde\nu=0$ ），亚声速特征远场（ $\tilde\nu_\infty=3\nu_\infty$ ） 
- 数值： Δt=1.5e-3 s，800 物理步，内迭代最多 20，伪时间 CFL=100 
- 结果： $\mathrm{St}=fD/U=0.197$（13.67 Hz），Cd 均值 **1.144**（1.089–1.201），Cl 范围 **−0.687–0.691**；壁钟 164 s 

### NACA0012（[urans_m07a5.json](cases/naca0012/urans_m07a5.json)）


- 工况： Ma=0.7，Re=1e7，α=5°，T=300 K 
- 网格： 弦长 1 m，O 型网格 326×89，远场约 9c 
- 数值： 先 1500 步定常（steady_cfl=50）给初场，再 3000 步 URANS，Δt=2e-4 s，内迭代 40 
- 结果： 激波自持振荡（抖振）：上表面激波 x/c **0.207–0.314**（均值 0.274），Cl **0.439–0.674**（均值 0.584），Cd 均值 0.0354，主频 **12.92 Hz**（ $k=\pi fc/U=0.167$ ），上表面 Cp RMS 峰值 **0.40**（位于 x/c=0.285 激波处）；壁钟 2106 s 

两个算例每个物理步都建议让内迭代降到 `inner_tol`（连续方程残差下降约两个量级）；NACA0012 的内迭代加到 80 次后振荡频率与幅值基本不变。结果尚未做网格独立性验证。

### OAT15A 抖振（[urans_m073a35_re33e6.json](cases/oat15a/urans_m073a35_re33e6.json)）

网格用 `cases/oat15a/icm_mesh.txt`（62591 单元、474 个壁面面、不带结构化邻接表）。配制为「1500 步定常预热
（CFL 2→50）+ 2000 物理步」，壁钟约 2 小时（96 线程，每物理步 3.5–4.0 s）。

| 工况 | dt | 伪时间 CFL | 后半段 Cl | 后半段 Cd | 主频 / St | 激波 x/c |
| --- | --- | --- | --- | --- | --- | --- |
| α=3.5°，Ma=0.73，Re=3.3e6 | 1e-4 s | 100 | 均值 **0.896**（0.731–1.027）| 均值 **0.0439** | **71.7 Hz** / 0.0651（$\Delta f$=10 Hz）| 均值 0.455 |
| α=3.9°，Ma=0.73，Re=3.3e6 | 5e-5 s | 200 | 均值 **0.854**（0.627–1.067）| 均值 **0.0505** | **74.5 Hz** / 0.0676（$\Delta f$=20 Hz）| 均值 0.395 |

（力系数输出参考长度硬编码为 1，标准弦长都是 0.23 m，实际 $C_l$ 需再除以 0.23。）

这两个工况正是 2.3.3 之前blow-up两个算例；修好后跑满全程且**从未触发过 CFL 回退**。

## 以开发者身份运行

### 关于json

配置文件分 `io`、`solver`、`farfield` 三部分。`io` 下 `mesh`、`log`、`field` 必填，`structured`（结构化邻接表）与 `restart`（起始流场）可选：给出结构化邻接表路径时启用虚拟单元（PLOT3D 单块，拓扑由 `S_MAX`/`N_MAX` 隐含），缺省表示网格不带结构化信息。

`farfield.alpha` 是以度为单位、从 +x 方向逆时针为正的来流迎角，取值需满足 $|\alpha|<90$。`farfield.p` 是来流静压，要按雷诺数给工况就反算。

线程数由 `OMP_NUM_THREADS` 决定（**优先于** `solver.threads`），未设置时取 `omp_get_max_threads()`（进程可用核数，SMT 按硬件线程计数），二者都受 `OMP_THREAD_LIMIT` 限制。本机实测最佳线程数：圆柱与 NACA0012 都是 32（8–128 扫描）。

未知字段、重复键和无效值都会报错；**所有文件路径都相对配置文件所在目录解析**。

线程数**可能**改变数值结果（机制见「数值约定」里关于分块的那段），所以对比两次运行前先确认线程数一致；
每次运行的 `run.log` 第一行就是 `OpenMP threads: N`。

### 代码检查

```sh
uv run ruff check cases tests gui            # Python 静态检查
uv run ruff format --check cases tests gui    # Python 格式（缩进 2、单引号，见 pyproject.toml）
clang-format --dry-run --Werror src/*.hpp src/*.cpp
```

构建启用严格编译警告（`-Wall -Wextra -Wpedantic -Werror`，未使用参数/变量仍为警告）；
`-DPURINE_SANITIZE=ON` 启用 ASan/UBSan，`-DPURINE_NATIVE=ON` 加 `-march=native`。

### 检验台

```sh
bash tests/fdb.sh cases/oat15a/urans_m073a31_re33e6.json 10 40
```

[fdb.sh](tests/fdb.sh) 会编出 `debug` 目标（[tests/debug.cpp](tests/debug.cpp)）并跑 OAT15A 3.1° 标定，
输出 `debugres/summary.txt`（逐项检验与标定结论）与 `debugres/calib.csv`（CFL × 内迭代收敛曲线）。
这个检验台与主程序走同一套双时间步路径：逐项检查网格/几何/梯度/通量/SA 源项，再用远场边界的数值微分
Jacobi 对拍。

**一般建议普通使用时不运行检验台，仅在开发时使用**，我们不保证提供的`debug.cpp`是最新的

## 代码结构

`src/` 下每个模块一个 `.hpp`，数据定义与实现合并，`main.cpp` 是唯一的翻译单元；除 `main.cpp` 外全部为头文件，按包含关系分层：

| 层 | 文件 | 职责 |
| --- | --- | --- |
| 数据与控制量 | [config.hpp](src/config.hpp)、[classconfig.hpp](src/classconfig.hpp) | 基本类型（`vec2`/`vec4`/`mat5`）、单元与面的数据结构、全局控制量 |
| 网格与几何 | [readmesh.hpp](src/readmesh.hpp)、[geometry.hpp](src/geometry.hpp)、[HALO.hpp](src/HALO.hpp) | 读入与连接拓扑、面积/中心/带长面法向；结构化网格的虚拟单元与远场插值 |
| 物性与初边值 | [physic.hpp](src/physic.hpp)、[initialize.hpp](src/initialize.hpp)、[udf.hpp](src/udf.hpp)、[boundary.hpp](src/boundary.hpp) | Sutherland 黏度等物性、来流初始化、无滑移绝热壁面与亚声速特征远场 |
| 空间离散 | [interpolate.hpp](src/interpolate.hpp)、[grad.hpp](src/grad.hpp)、[convect.hpp](src/convect.hpp)、[dissipation.hpp](src/dissipation.hpp)、[SA.hpp](src/SA.hpp) | MUSCL 重构与面插值、最小二乘单元梯度与面梯度、Roe 通量差分裂（JST 人工耗散可选）、SA 湍流通量与源项 |
| 时间推进 | [dualtime.hpp](src/dualtime.hpp)、[dualtime_full.hpp](src/dualtime_full.hpp) | 双时间步：BDF1 外迭代与伪时间内迭代；块 Jacobi 组装与对称 Gauss-Seidel 求解 |
| 驱动与输出 | [main.cpp](src/main.cpp)、[io.hpp](src/io.hpp) | 配置解析、主循环、日志与流场输出 |

`timarch.hpp`（显式时间推进，2.2.x版本）、`residual.hpp`、`parallel.hpp` 目前不被引用，为保留模块。

## 数值约定

**空间离散**：有限体积、单元中心格式。对流项用 Roe 通量差分裂，面左右状态由 MUSCL 重构给出；梯度用单元最小二乘，黏性通量由面梯度计算；SA 方程对流项一阶迎风、扩散项两点格式，破坏项线性化后进对角。`dissipation.hpp` 的 JST 人工耗散为可选（`cc::scheme`）。

**时间推进**：双时间步。物理时间一阶后向欧拉（BDF1）

$$\frac{V}{\Delta t}\left(Q^{n+1}-Q^{n}\right)+R\left(Q^{n+1}\right)=0$$

伪时间隐式推进，每个内迭代解

$$\left[\frac{V}{\Delta\tau}+\frac{V}{\Delta t}+\frac{\partial R}{\partial Q}\right]\Delta Q=-\left[R\left(Q^{k}\right)+\frac{V\left(Q^{k}-Q^{n}\right)}{\Delta t}\right],\qquad \frac{V}{\Delta\tau}=\frac{1}{\mathrm{cfl}}\sum_{f}\lambda_{f}$$

$\Delta t\to 0$ 时退化为定常隐式迭代。线性系统用按单元编号分块的对称 Gauss-Seidel 扫描（`block_sgs`）求解：块内取本遍最新值、块间取上一遍的值。
**注意：分块数等于 OpenMP 线程数**（[dualtime.hpp](src/dualtime.hpp) 里 `lo = cell_num*tid/nt`、`hi = cell_num*(tid+1)/nt`），
所以改线程数改的是**迭代算子本身**，而不只是浮点归约顺序——残差下降、收敛路径、乃至混沌的 URANS 轨迹都会随之变化。
线程越多块越小，块间耦合越弱（越接近纯 Jacobi），内迭代越难收敛；线程越少块越大，隐式耦合越强、通常更稳但更慢。
要复现某个结果，必须用同样的线程数。2.3.1 起 $\partial R/\partial Q$ 取**完整的 4×4 块 Jacobi**，对角块与相邻块分别为

$$D_{c}=\frac{V}{\Delta\tau}+\frac{1}{\mathrm{cfl}}\sum_{f}\lambda_{f}+\frac{1}{2}\sum_{f}s_{f}A_{n,f}(c)+\frac{1}{2}\sum_{f}\left|\bar{A}_{n,f}\right|+\sum_{f}T_{f}-\frac{1}{2}\sum_{f}s_{f}B_{f}$$

$$L_{cj}=\frac{1}{2}s_{f}A_{n,f}(j)-\frac{1}{2}\left|\bar{A}_{n,f}\right|-T_{f}-\frac{1}{2}s_{f}B_{f}$$

这里 $\lambda_{f}$ 是面上对流与黏性谱半径之和，$T_{f}$ 是黏性通量对法向梯度的响应（内部面取两点薄层近似，壁面取一阶单侧、距离用 `tur.sad`）， $B_{f}$ 是冻结梯度下的黏性通量 Jacobi。

<!-- 
**容易被忽略的几点**：

- 守恒量是 $Q=[\rho,\rho u,\rho v,\rho E]$ ， $E=C_vT+(u^2+v^2)/2$ ，`vec4` 的成员名为 `.c/.x/.y/.e`；而 `phy.e` 是**单位质量**总能量，不是 $\rho E$ 。
- 面的 `nor` 为 $\boldsymbol{n}\,\Delta s$ （带面长）；`fnorm` 是 0/1 标志， $s=2\,\mathrm{fnorm}-1$ 把固定面法向翻成单元外法向。组装 Jacobi 时 $A_{n}$ 、 $B_{f}$ 随法向反号而变号，而 $|\bar{A}_{n}|$ 与 $T_{f}$ 不变号。
- 面变量（`phy`、`phgrad`、`tur`）是重构后的面值。内部面梯度取两侧单元最小二乘梯度的加权平均；**壁面**用一阶单侧两点差分，因此壁面梯度对单元自身状态的响应是单侧的。
- `config.hpp` 中 $C_v$ 与 $R/(\gamma-1)$ 并不严格相等（717.645 与 717.625），用 $\frac{\partial Q}{\partial W}$ 做变量变换时会残留约 1e-4 的相对误差。
- `other_cell()` 对非内部面返回 `nullptr`；结构化网格的虚拟单元存在独立的 `GhostList` 中，其 `index > cell_num`，不能用来索引单元数组。
- 每个面只由 `nei[0]` 所属的单元写一次，这是并行下避免同一面被多线程同时写的前提。 -->

收敛量取相邻两次检查之间 $\rho$ 、 $u$ 、 $v$ 、 $E$ 、 $\tilde\nu$ 的最大归一化状态更新量 $\max|\Delta\varphi|/\varphi_{\mathrm{ref}}$ ，相对首次检查降到 1e-4 且绝对值小于 1e-6 才停止定常迭代。

## 出错定位

**每个物理步**做完都扫一遍全场（`bad_cells()`）：密度、总能、T 或 ν̃ 出现非有限值，或密度/T 非正，
就打印出错单元的编号与坐标：

```text
  bad cell #59 (0.0257452,0.0510509)
step    860 发散: 12 个坏单元, 回退上一流场, CFL x1.0000 -> x0.5000
```

然后**回退到上一个有效流场，并把伪时间 CFL 减半重试这一步**；CFL 一路减到 1e-3 仍不收敛才放弃，
放弃前会先把回退后的流场写盘（见「续算与中断保存」）再返回非零码。所以：

- 日志里出现 `发散: ... 回退` 不等于这次运行失败，只是这一步自动降了 CFL
- 只要没有 `CFL 已降到下限仍不收敛, 放弃`，整趟就是正常跑完的

配置本身的错误在启动阶段打印一行直接报错并退出，求解器不使用异常。

## URANS（双时间步）

`./build/purinecfd <配置>` 运行 URANS，与 Fluent 密度基隐式求解器的做法一致；标量版算法在 [dualtime.hpp](src/dualtime.hpp)，2.3.1 的完整 Jacobi 版在 [dualtime_full.hpp](src/dualtime_full.hpp)：

- 外迭代：物理时间用一阶后向欧拉（BDF1），方程见“数值约定”。
- 内迭代：伪时间隐式推进；2.3.1 起 $\partial R/\partial Q$ 取完整的 4×4 块 Jacobi（见“数值约定”），按单元编号分块做对称 Gauss-Seidel 扫描。
- SA 方程在流动方程之后单独求解，同样 BDF1 + 隐式伪时间；对流一阶迎风，破坏项进对角。
- 残差 $R$ 复用 Roe + MUSCL + 最小二乘梯度 + SA 黏性通量。

配置在 `solver` 下增加 `urans` 一节，`end_step`（或 `max_steps`）为物理时间步数，`cfl` 为内迭代伪时间 CFL，`convergence_interval` 为屏幕输出间隔：

| 键 | 含义 | 默认 |
| --- | --- | --- |
| `dt` | 物理时间步长（s），必填 | — |
| `inner` | 每个物理步最多内迭代次数，必填 | — |
| `inner_tol` | 连续方程残差相对首次内迭代降到该值即结束内迭代 | 1e-3 |
| `sweeps` | 每次内迭代的对称 GS 扫描次数 | 4 |
| `steady_iters` | URANS 前的定常隐式迭代次数（Δt→∞），用来给初场 | 0 |
| `steady_cfl` | 定常迭代的 CFL 终值（前 200 步从 2 线性增大） | 50 |
| `wall_interval` | 壁面 Cp 输出间隔，0 不输出 | 0 |
| `seed` | 初场反对称涡扰动幅值（相对 U∞），对称问题起振用 | 0 |

## 后处理

求解器在 `io.field` 目录下写出 `history.csv`（每步 Cl、Cd、内迭代次数与残差）、`wall_cp.dat`、
`step_*.dat`（列依次为 x、y、ρ、u、v、T、p、Ma、ν̃）。两个脚本都只读「配置 + run 目录」：

```sh
# 曲线 / 云图 / PSD / README，一次出 15 张图
uv run --with numpy --with matplotlib --with scipy python tests/post.py \
    cases/oat15a/urans_m073a35_re33e6.json --case oat15a --out run/<目标目录> --readme --wall 7022

# 涡量动画（GIF）
uv run --with numpy --with matplotlib --with pillow python tests/vort_anim.py \
    cases/oat15a/urans_m073a35_re33e6.json --out run/<目标目录> --case oat15a --field vorticity
```

`tests/post.py` 的要点：

- `--case` 必填（`cylinder` / `naca0012` / `oat15a`），决定网格结构化与否、视图范围、力系数参考长度
- `--out` 指定 run 目录；配置里的 `io.field` 是占位符时**必须**用它
- 统计量：力系数均值与范围、主频与 St（Hann 窗 + 补零 FFT 定位峰，再按原始分辨率 $1/T$ 抛物线插值）、
  壁面探针 PSD、上表面激波位置序列、Cp 脉动 RMS
- `--panels u,v,rho,...` 挑云图面板；`--no-history` / `--no-residual` / `--no-cp-dist` / `--no-cp-rms` /
  `--no-probe` / `--no-sequence` 关掉对应输出（默认全部开启）
- 求解器非正常退出也能跑：会插一句「求解器在 step N 发散」，统计只覆盖已有的步数
- 结构化算例（`cylinder` / `naca0012`）的配置若缺 `io.structured`，会打印提示并退回非结构后端渲染

`tests/vort_anim.py` 把 `step_*.dat` 序列做成动画；`--field` 可选
`vorticity` / `schlieren` / `Ma` / `Cp` / `u` / `v`，`--fps` / `--stride` / `--level` / `--size`
分别控制帧率、抽帧、涡量色标幅度与图幅。

算例配置与最新结果见上面的“算例”一节。

## 图形界面

```sh
./gui.sh          # 等价于 uv run --with pyside6 python gui/main.py
```

窗口是六个选项卡，按操作顺序排：

| 选项卡 | 内容 |
| --- | --- |
| ① 算例 | 算例（cylinder / naca0012 / oat15a）、配置模板（首项是内置默认模板）、网格文件、结构化邻接表、起始流场（可留空，填了就接着那个场算）、输出目录 |
| ② 来流 | 马赫数、雷诺数、静温 T、参考长度 L、攻角 α；下方实时显示换算出的 μ / a / U∞ / ρ∞ / p∞ / q∞ |
| ③ 求解 | 求解参数**表单**，字段由所选模板的 `solver` 段动态生成，不暴露 JSON |
| ④ 后处理 | 曲线/统计与 9 个云图面板的复选框、探针位置、u 色标 |
| ⑤ 运行 | 可执行文件、生成配置 JSON、自动测线程数、开始运行 / 仅跑后处理 / 停止 |
| ⑥ 控制台 | 求解器日志 + stderr + 后处理输出 |

<!-- 设计要点：

- **② 只需要填雷诺数**：按 Sutherland 公式反算 `p∞`，与 [cases/cylinder/generate_case.py](cases/cylinder/generate_case.py) 同一套公式
- **没有模板也能跑**：① 页模板下拉首项是「（内置默认）」，`cases/<算例>/` 下没有可用 json 时用内置 `solver` 模板兜底；模板里的 `io` 路径会被重写成真实 run 目录
- **不改模板**：每次运行在 `cases/<算例>/_run_<编号>.json` 生成一份临时配置，模板文件字节不变
- **能从中间流场续算**：① 页「起始流场」填一个 `step_*.dat` 即写入 `io.restart`；求解结束时控制台会直接把最新场文件路径列出来提醒
- **停止不会丢数据**：⑤ 页「停止」先发 SIGTERM，求解器跑完当前物理步、写下该步流场再退出（30 秒无响应才强杀）
- **不改动现有命令行流程**：GUI 只做「生成配置 + 起进程 + 收日志」，求解器与 `tests/post.py` 的调用方式完全不变
- 输出目录沿用 `run/<编号>-<标签>-<日期>-<时分>/`，编号自动取 `run/` 下最大值加一
- 控制台靠**轮询 `run.log`**（求解器内部用 `freopen` 把 stdout 重定向进了日志文件），`stderr` 走管道直连
- 求解结束自动接后处理，非正常退出也接
- **自动测线程数**：在 `run/_bench/` 里用各档线程数各跑 25 步定常预热 + 1 个物理步，比壁钟与定常残差，
  给出「最快」与「残差下降最多」两档（线程数会改变迭代算子，二者未必一致）；只走临时目录，不碰正式 run -->
<!-- 
依赖只有 PySide6（[pyproject.toml](pyproject.toml) 的 `gui` 组）。 -->

细节见 [gui/README.md](gui/README.md)。
