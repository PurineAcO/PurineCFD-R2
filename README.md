# PurineCFD-R2

代码版本 2.3.1

2.3.1 起，双时间步的隐式求解使用**完整的 4×4 块 Jacobi**（[dualtime_full.hpp](src/dualtime_full.hpp)）：对角块与相邻块由对流 Jacobi `A_n`、Roe 耗散矩阵 `|Ā_n|`、黏性梯度响应 `T_f`、冻结梯度 `B_f` 组装，不再用单个标量谱半径近似；2.3.0 的标量版本原样保留在 [dualtime.hpp](src/dualtime.hpp)。

## 构建与运行

依赖 C++17 编译器与 OpenMP（Linux 用 GCC；Windows 推荐 MSYS2 的 MinGW-w64 GCC）、CMake ≥3.16、Python ≥3.10。nlohmann/json 3.11.3 随源码提供。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPURINE_NATIVE=ON -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build build -j 8
python3 cases/cylinder/generate_case.py
./build/purinecfd config.json
```

配置见 [config.json](config.json)，分为 `io`、`solver`、`farfield` 三部分。`io` 下 `mesh`、`log`、`field` 必填，`structured` 可选：给出结构化邻接表路径时启用虚拟单元（PLOT3D 单块，拓扑由 `S_MAX`/`N_MAX` 隐含），缺省表示网格不带结构化信息。未知字段、重复键和无效值都会报错；文件路径相对于配置文件解析。`farfield.alpha` 是以度为单位、从 +x 方向逆时针为正的来流迎角，取值需满足 $|\alpha|<90$。

线程数由 `OMP_NUM_THREADS` 决定，未设置时取 `omp_get_max_threads()`（进程可用核数，SMT 按硬件线程计数），也可用 `solver.threads` 指定，二者都受 `OMP_THREAD_LIMIT` 限制。本机实测最佳线程数：圆柱与 NACA0012 都是 32（8~128 扫描）。

代码检查与检验台：

```sh
uv run ruff check cases tests
uv run ruff format --check cases tests
clang-format --dry-run --Werror src/*.hpp src/*.cpp
bash tests/fdb.sh 400 60   # 构建 debug 检验台并跑圆柱, 结果写到 debugres/
```

构建启用严格编译警告（`-Wall -Wextra -Wpedantic -Werror`，未使用参数/变量仍为警告）；`-DPURINE_SANITIZE=ON` 启用 ASan/UBSan，`-DPURINE_NATIVE=ON` 加 `-march=native`。

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

`timarch.hpp`（显式时间推进）、`residual.hpp`、`parallel.hpp` 目前不被主程序引用：第一个由 `tests/debug.cpp` 检验台使用，后两个为保留模块。

## 数值约定

**空间离散**：有限体积、单元中心格式。对流项用 Roe 通量差分裂，面左右状态由 MUSCL 重构给出；梯度用单元最小二乘，黏性通量由面梯度计算；SA 方程对流项一阶迎风、扩散项两点格式，破坏项线性化后进对角。`dissipation.hpp` 的 JST 人工耗散为可选（`cc::scheme`）。

**时间推进**：双时间步。物理时间一阶后向欧拉（BDF1）

$$\frac{V}{\Delta t}\left(Q^{n+1}-Q^{n}\right)+R\left(Q^{n+1}\right)=0$$

伪时间隐式推进，每个内迭代解

$$\left[\frac{V}{\Delta\tau}+\frac{V}{\Delta t}+\frac{\partial R}{\partial Q}\right]\Delta Q=-\left[R\left(Q^{k}\right)+\frac{V\left(Q^{k}-Q^{n}\right)}{\Delta t}\right],\qquad \frac{V}{\Delta\tau}=\frac{1}{\mathrm{cfl}}\sum_{f}\lambda_{f}$$

$\Delta t\to 0$ 时退化为定常隐式迭代。线性系统用按单元编号分块的对称 Gauss-Seidel 扫描（`block_sgs`）求解：块内取本遍最新值、块间取上一遍的值，结果与线程调度无关。2.3.1 起 $\partial R/\partial Q$ 取**完整的 4×4 块 Jacobi**，对角块与相邻块分别为

$$D_{c}=\frac{V}{\Delta\tau}+\frac{1}{\mathrm{cfl}}\sum_{f}\lambda_{f}+\frac{1}{2}\sum_{f}s_{f}A_{n,f}(c)+\frac{1}{2}\sum_{f}\left|\bar{A}_{n,f}\right|+\sum_{f}T_{f}-\frac{1}{2}\sum_{f}s_{f}B_{f}$$

$$L_{cj}=\frac{1}{2}s_{f}A_{n,f}(j)-\frac{1}{2}\left|\bar{A}_{n,f}\right|-T_{f}-\frac{1}{2}s_{f}B_{f}$$

这里 $\lambda_{f}$ 是面上对流与黏性谱半径之和，$T_{f}$ 是黏性通量对法向梯度的响应（内部面取两点薄层近似，壁面取一阶单侧、距离用 `tur.sad`），$B_{f}$ 是冻结梯度下的黏性通量 Jacobi。

**容易被忽略的几点**：

- 守恒量是 $Q=[\rho,\rho u,\rho v,\rho E]$，$E=C_vT+(u^2+v^2)/2$，`vec4` 的成员名为 `.c/.x/.y/.e`；而 `phy.e` 是**单位质量**总能量，不是 $\rho E$。
- 面的 `nor` $=\boldsymbol{n}\,\Delta s$ 带面长；`fnorm` 是 0/1 标志，$s=2\,\mathrm{fnorm}-1$ 把固定面法向翻成单元外法向。组装 Jacobi 时 $A_{n}$、$B_{f}$ 随法向反号而变号，而 $|\bar{A}_{n}|$ 与 $T_{f}$ 不变号。
- 面变量（`phy`、`phgrad`、`tur`）是重构后的面值。内部面梯度取两侧单元最小二乘梯度的加权平均；**壁面**用一阶单侧两点差分，因此壁面梯度对单元自身状态的响应是单侧的。
- `config.hpp` 中 $C_v$ 与 $R/(\gamma-1)$ 并不严格相等（717.645 与 717.625），用 $\partial Q/\partial W$ 做变量变换时会残留约 $10^{-4}$ 的相对误差。
- `other_cell()` 对非内部面返回 `nullptr`；结构化网格的虚拟单元存在独立的 `GhostList` 中，其 `index > cell_num`，不能用来索引单元数组。
- 每个面只由 `nei[0]` 所属的单元写一次，这是并行下避免同一面被多线程同时写的前提。

收敛量取相邻两次检查之间 $\rho$、$u$、$v$、$E$、$\tilde\nu$ 的最大归一化状态更新量 $\max|\Delta\varphi|/\varphi_{\mathrm{ref}}$，相对首次检查降到 $10^{-4}$ 且绝对值小于 $10^{-6}$ 才停止定常迭代。

## 出错定位

每个 `convergence_interval` 步做一次全场有效性检查（串行扫描）：密度、总能或 ν̃ 出现 NaN，或密度非正，就打印首个出错单元的编号与坐标，并以非零码退出。

```text
Error: invalid state in cell #59 (0.0257452,0.0510509)
```

求解器不使用异常：出错时打印一行 `Error: ...` 并返回 `false`，由 `main` 统一返回非零退出码。

## 算例

两个算例的完整记录都在 `run/` 下（`README.md`、`summary.txt`、`history.csv`、逐帧 `step_*.dat` 与各类云图）。下面数字为 2.3.1 的结果，32 线程。

### 圆柱绕流（[config.json](config.json)）

| 项目 | 设置 |
| --- | --- |
| 工况 | Ma=0.2，$\mathrm{Re}=\rho UD/\mu=1000$，T=300 K，α=0°，来流沿 +x |
| 物性 | Sutherland 黏度；由 Re 反算 p∞≈22.89061 Pa，U∞≈69.44379 m/s |
| 网格 | D=1 m，192×128 单元，远场半径 30D，第一层高度 0.0015D |
| 边界 | 无滑移绝热壁面（$\tilde\nu=0$），亚声速特征远场（$\tilde\nu_\infty=3\nu_\infty$） |
| 数值 | Δt=1.5e-3 s，800 物理步，内迭代最多 20，伪时间 CFL=100 |
| 结果 | $\mathrm{St}=fD/U=$**0.197**（13.67 Hz），Cd 均值 **1.144**（1.089~1.201），Cl 范围 **−0.687~0.691**；壁钟 164 s |

### NACA0012（[urans_m07a5.json](cases/naca0012/urans_m07a5.json)）

| 项目 | 设置 |
| --- | --- |
| 工况 | Ma=0.7，Re=1e7，α=5°，T=300 K |
| 网格 | 弦长 1 m，O 型网格 326×89，远场约 9c |
| 数值 | 先 1500 步定常（steady_cfl=50）给初场，再 3000 步 URANS，Δt=2e-4 s，内迭代 40 |
| 结果 | 激波自持振荡（抖振）：上表面激波 x/c **0.207~0.314**（均值 0.274），Cl **0.439~0.674**（均值 0.584），Cd 均值 0.0354，主频 **12.92 Hz**（$k=\pi fc/U=$**0.167**），上表面 Cp RMS 峰值 **0.40**（位于 x/c=0.285 激波处）；壁钟 2106 s |

两个算例每个物理步都建议让内迭代降到 `inner_tol`（连续方程残差下降约两个量级）；NACA0012 的内迭代加到 80 次后振荡频率与幅值基本不变。结果尚未做网格独立性验证。

## URANS（双时间步）

`./build/purinecfd <配置>` 运行 URANS，与 Fluent 密度基隐式求解器的做法一致；标量版算法在 [dualtime.hpp](src/dualtime.hpp)，2.3.1 的完整 Jacobi 版在 [dualtime_full.hpp](src/dualtime_full.hpp)：

- 外迭代：物理时间用一阶后向欧拉（BDF1），方程见“数值约定”。
- 内迭代：伪时间隐式推进；2.3.1 起 $\partial R/\partial Q$ 取完整的 4×4 块 Jacobi（见“数值约定”），按单元编号分块做对称 Gauss-Seidel 扫描。
- SA 方程在流动方程之后单独求解，同样 BDF1 + 隐式伪时间；对流一阶迎风，破坏项进对角。
- 残差 $R$ 复用 Roe + MUSCL + 最小二乘梯度 + SA 黏性通量。

配置在 `solver` 下增加 `urans` 一节，`max_steps` 为物理时间步数，`cfl` 为内迭代伪时间 CFL，`convergence_interval` 为屏幕输出间隔：

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

`io.field` 目录下输出 `history.csv`（每步 Cl、Cd、内迭代次数与残差）、`wall_cp.dat`、`step_*.dat`（列依次为 x、y、ρ、u、v、T、p、Ma、ν̃）。后处理：

```sh
uv run --with matplotlib --with numpy python tests/post_urans.py cases/naca0012/urans_m07a5.json --airfoil
```

给出主频、激波位置统计、力系数历程和马赫数 / Cp / 数值纹影云图。主频用 Hann 窗 + 补零的 FFT，并按原始分辨率做抛物线插值（同时给出分辨率）；翼型还会用上表面激波位置序列独立做一次 FFT 交叉验证。

算例配置与最新结果见上面的“算例”一节。
