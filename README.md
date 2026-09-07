# PurineCFD-R2

二维定常可压缩流动求解器，支持层流和 Spalart–Allmaras（SA-RANS），使用 C++17、有限体积法、JST 人工耗散和 OpenMP。

**第一次使用请从[安装与圆柱算例教程](docs/quickstart.md)开始。**教程覆盖 Windows/WSL、工具安装、编译、运行、收敛检查和结果文件，命令可逐行复制。

- [数值方法、代码结构与配置](docs/numerics.md)
- [后续全局稳定性分析：Crouch-Py](https://github.com/PurineAcO/crouch-hopf)

已有 Linux 开发环境时，在仓库目录执行：

```sh
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
OMP_NUM_THREADS=4 python3 cases/cylinder/run.py --model laminar --re 47 --output run/cylinder-laminar-re47
```

`--model sa` 切换到 SA。底层配置的 `solver.model` 必填，不接受旧模型开关。圆柱运行入口会自动生成网格、设置 Re、保存模型与收敛标记，并拒绝覆盖已有结果。

附带算例为 Ma=0.2、T=300 K 的圆柱；NACA 0012 需要另行准备和验证网格及边界条件。当前程序求定常基流，稳定性分析由另一个程序完成。
