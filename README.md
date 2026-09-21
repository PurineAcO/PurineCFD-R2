# PurineCFD-R2

二维定常 SA-RANS 有限体积求解器，使用凸四边形网格、JST 人工耗散、五级 RK 伪时间推进和 OpenMP。

## 构建与运行

依赖支持 C++17 与 OpenMP 的编译器（GCC ≥9、Clang、MSVC）、CMake ≥3.16、Python ≥3.10。nlohmann/json 3.11.3 随源码提供。除 OpenMP 外无第三方依赖；目录创建与物理核心探测在 POSIX 和 Windows 下分别走各自的分支，其余代码平台无关。GCC 8 的 `<filesystem>` 实现不完整，`config.cpp` 需要 GCC ≥9。

POSIX 平台（GCC/Clang）：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPURINE_NATIVE=ON -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build build -j 8
python3 cases/cylinder/generate_case.py
./build/purinecfd config.json
```

Windows（Clang + libomp；`CMAKE_PREFIX_PATH` 指向 LLVM 安装目录，使其找到 `libomp`）：

```powershell
cmake -S . -B build -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release `
  -DCMAKE_CXX_COMPILER="<LLVM>/bin/clang++.exe" -DCMAKE_PREFIX_PATH="<LLVM>"
cmake --build build -j 8
python cases/cylinder/generate_case.py
.\build\purinecfd.exe config.json
```

配置见 [config.json](config.json)，分为 `io`、`solver`、`farfield` 三部分。字段全部必填；未知字段、重复键和无效值会报错。文件路径相对于配置文件解析。`farfield.alpha` 是以度为单位、从 +x 方向逆时针为正的来流迎角，取值需满足 `|alpha| < 90`。

默认按进程可用的物理核心数分配线程，SMT 不重复计数；拓扑读取失败时退回 `omp_get_num_procs()`。`OMP_NUM_THREADS` 可显式设置线程数，自动设置遵守 `OMP_THREAD_LIMIT`。

## 按什么顺序读代码

`include/` 放数据定义和函数声明，`src/` 放实现。建议先读数据结构，再沿迭代主线阅读公式。

| 顺序 | 文件 | 关注的问题 |
| --- | --- | --- |
| 1 | [config.h](include/config.h)、[classconfig.h](include/classconfig.h) | `vec2`、`physics`、`turbulence`、`dissipation` 各存什么？一个单元、一个面分别保存什么？ |
| 2 | [main.cpp](src/main.cpp) | 程序如何启动？`solve()` 里一个伪时间步和一个 RK 阶段做哪些事？ |
| 3 | [readmesh.cpp](src/readmesh.cpp)、[geometry.cpp](src/geometry.cpp) | 如何读入并连接网格，计算面积、中心和面法向？ |
| 4 | [physic.cpp](src/physic.cpp)、[initialize.cpp](src/initialize.cpp)、[boundary.cpp](src/boundary.cpp) | 如何计算空气物性，初始化来流，施加壁面及远场条件？ |
| 5 | [interpolate.cpp](src/interpolate.cpp)、[grad.cpp](src/grad.cpp)、[convect.cpp](src/convect.cpp) | 如何由面值求梯度，并汇总穿过单元边界的通量？ |
| 6 | [dissipation.cpp](src/dissipation.cpp)、[SA.cpp](src/SA.cpp) | 人工耗散和湍流模型分别加入哪些项？ |
| 7 | [timarch.cpp](src/timarch.cpp)、[residual.cpp](src/residual.cpp)、[io.cpp](src/io.cpp) | 时间步怎样确定？何时停止？结果怎样输出？ |

配置解析和线程设置分别在 `config.cpp`、`parallel.cpp`；前者用 `std::filesystem` 归一化配置中的相对路径，后者按平台统计物理核心（POSIX 下读 `sched_getaffinity` 与 `/sys` 拓扑，Windows 下用 `GetProcessAffinityMask` 与 `GetLogicalProcessorInformationEx`），SMT 不重复计数。边界条件参数在 `udf.h`。

## 数值约定

- 单元守恒量为 `Q = [ρ, ρu, ρv, ρE]`，`E = Cv*T + (u²+v²)/2`。代码中的 `phy.e` 是单位质量总能量。守恒量存在 `conser`，本步 RK 的基准存在 `conserformer`。
- 面的 `nor = n*Δs` 包含面长 Δs。单元用 `fnorm` 把固定面法向转换成自己的外法向。
- Green–Gauss 梯度为 `∇φ = Σf(φf*n_f*Δs_f)/V`，这里的 `V` 是二维单元面积（`vol`），按单位厚度计。
- 每阶段按 `Qᵏ = Qⁿ − αk*Δt*R(Qᵏ⁻¹)` 更新，`Qⁿ` 始终是本步开始时的状态。系数依次为 `1/4、1/6、3/8、1/2、1`。
- `tur.miubl` 是 SA 工作变量 ν̃，湍流运动黏度为 `νt = ν̃*fv1`；`tur.sad` 是到最近壁面中点的距离。`SA.cpp` 中的源项按产生、破坏、梯度平方三项展开。

`main.cpp` 的 `solve()` 里，每个 RK 阶段依次：恢复原始量（`reform`/`form_physic`）→ 壁面与远场边界 → 面插值 → 单元梯度与 JST 激波检测 → 面梯度与通量 → 单元汇总与湍流方程 → 统一推进守恒量。OpenMP 循环末尾的同步保证相邻单元读取的是同一阶段数据；面插值只由 `nei[0]` 所属的单元负责，保证每个面只被一个线程写一次。壁面距离 `tur.sad` 在初始化时算一次。

收敛量取相邻两次检查之间 ρ、u、v、E、ν̃ 的最大归一化状态更新量 `max|Δφ|/φ_ref`。相对于首次检查降低至 `1e-4`，且绝对值小于 `1e-6`，才停止迭代。

## 出错定位

求解过程中出现 NaN 或非正的密度、温度时，会报告出错时执行的函数、步数、网格编号和该网格的物理量：

```text
Error: Invalid flow state at step 37, function rk_stage, cell #1 (0.502703,0.012341)
 rho=-7.713843e-03 u=5.473518e+02 v=2.984908e+00 T=1.946697e+02 p=1.030077e+00 miubl=0.000000e+00
```

标记在并行循环内完成，只有真出错时才串行扫描定位，报告内容与线程数无关。

## 代码风格

数据结构与函数命名沿用仓库早期版本：`cc::` 命名空间、`vec2`、`phy`/`tur`/`diss`、`conser`/`conserformer`、`#define allcell`、中文注释、4 空格缩进。不使用异常：出错时打印一行 `Error: ...` 并返回 `false`，由 `main` 统一返回非零退出码。

## 圆柱算例与检查

| 项目 | 设置 |
| --- | --- |
| 来流 | Ma=0.2，Re=45，T=300 K，沿 +x 方向 |
| 尺寸和网格 | D=1 m，128×96 单元，远场半径 30D，第一层高度 0.006D |
| 物性 | Sutherland 黏度；由 Re=ρUD/μ 反算 p∞≈1.03007749 Pa，U∞≈69.44379 m/s |
| 边界 | 无滑移绝热壁面，ν̃=0；亚声速特征远场，ν̃∞=3ν∞ |
| 迭代 | CFL=1，SA 包含 ft2 项；参考结果在 32,800 步收敛，Cd≈1.49 |

该结果尚未验证网格独立性。日志默认写入 `run/run.log`；`Converged at step` 表示达到收敛判据，步数用尽会单独提示。错误返回非零退出码。
流场为 `run/field/step_XXXXXX.dat`，列依次为 x、y、ρ、u、v、T、p、Ma、ν̃。

```sh
uv run ruff check cases tests
uv run ruff format --check cases tests
clang-format --dry-run --Werror include/*.h src/*.cpp
PURINECFD_BIN="$PWD/build/purinecfd" uv run pytest -q
```

`PURINECFD_BIN` 省略时按平台取 `build/purinecfd`（POSIX）或 `build/purinecfd.exe`（Windows）。两个断言自动线程数与 CPU 亲和性关系的用例依赖 Linux 的亲和接口，在其它平台会被跳过。

构建启用严格编译警告。`-DPURINE_SANITIZE=ON` 启用 ASan/UBSan。测试覆盖数值回归、线程一致性、输入校验和错误报告。
