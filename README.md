# PurineCFD-R2

二维定常 SA-RANS 有限体积求解器，使用凸四边形网格、JST 人工耗散、五级 RK 伪时间推进和 OpenMP。

## 构建与运行

依赖 Linux、GCC（C++17/OpenMP）、CMake ≥3.16、Python ≥3.10。nlohmann/json 3.11.3 随源码提供。

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DPURINE_NATIVE=ON -DCMAKE_INTERPROCEDURAL_OPTIMIZATION=ON
cmake --build build -j 8
python3 cases/cylinder/generate_case.py
./build/purinecfd config.json
```

配置见 [config.json](config.json)，分为 `io`、`solver`、`farfield` 三部分。字段全部必填；未知字段、重复键和无效值会报错。文件路径相对于配置文件解析。

默认按进程可用的物理核心数分配线程，SMT 不重复计数；CPU 拓扑读取失败会报错。`OMP_NUM_THREADS` 可显式设置线程数，自动设置遵守 `OMP_THREAD_LIMIT`。

## 按什么顺序读代码

`include/` 放数据定义和函数声明，`src/` 放实现。建议先读数据结构，再沿迭代主线阅读公式。

| 顺序 | 文件 | 关注的问题 |
| --- | --- | --- |
| 1 | [flow.h](include/flow.h)、[mesh.h](include/mesh.h) | 一个单元、一个面分别保存什么？变量的物理含义和单位是什么？ |
| 2 | [main.cpp](src/main.cpp)、[solver.cpp](src/solver.cpp) | 程序如何启动？一个伪时间步和一个 RK 阶段做哪些事？ |
| 3 | [mesh_reader.cpp](src/mesh_reader.cpp)、[geometry.cpp](src/geometry.cpp) | 如何读入并连接网格，计算面积、中心和面法向？ |
| 4 | [physics.cpp](src/physics.cpp)、[initialization.cpp](src/initialization.cpp)、[boundary.cpp](src/boundary.cpp) | 如何计算空气物性，初始化来流，施加壁面及远场条件？ |
| 5 | [gradients.cpp](src/gradients.cpp)、[fluxes.cpp](src/fluxes.cpp) | 如何由面值求梯度，并汇总穿过单元边界的通量？ |
| 6 | [jst.cpp](src/jst.cpp)、[spalart_allmaras.cpp](src/spalart_allmaras.cpp) | 人工耗散和湍流模型分别加入哪些项？ |
| 7 | [timestep.cpp](src/timestep.cpp)、[convergence.cpp](src/convergence.cpp)、[io.cpp](src/io.cpp) | 时间步怎样确定？何时停止？结果怎样输出？ |

配置解析和线程设置分别在 `config.cpp`、`parallel.cpp`。

## 数值约定

- 单元守恒量为 `Q = [ρ, ρu, ρv, ρE]`，`E = Cv*T + (u²+v²)/2`。代码中的 `flow.e` 是单位质量总能量。
- 面的 `area_normal = n*Δs` 包含面长 Δs。单元用 `normal_points_outward` 将固定面法向转换成自己的外法向。
- Green–Gauss 梯度为 `∇φ = Σf(φf*n_f*Δs_f)/V`，这里的 `V` 是二维单元面积，按单位厚度计。
- 每阶段按 `Qᵏ = Qⁿ − αk*Δt*R(Qᵏ⁻¹)` 更新，`Qⁿ` 始终是本步开始时的状态。系数依次为 `1/4、1/6、3/8、1/2、1`。
- `turbulence.nu_tilde` 是 SA 工作变量 ν̃，湍流运动黏度为 `νt = ν̃*fv1`。`spalart_allmaras.cpp` 中的源项按产生、破坏、梯度平方三项展开。

`solver.cpp` 中，每阶段依次恢复状态与边界、计算单元梯度、计算面通量、汇总右端项、统一更新状态。OpenMP 循环末尾的同步保证相邻单元读取的是同一阶段数据。固定几何量在初始化时缓存，面通量每阶段更新一次。

收敛量取相邻两次检查之间 ρ、u、v、E、ν̃ 的最大归一化状态更新量 `max|Δφ|/φ_ref`。相对于首次检查降低至 `1e-4`，且绝对值小于 `1e-6`，才停止迭代。

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

构建启用严格编译警告。`-DPURINE_SANITIZE=ON` 启用 ASan/UBSan。测试覆盖数值回归、线程一致性、输入校验和错误报告。
