# 从安装到第一个圆柱算例

这里的“终端”是输入命令的窗口。每个代码框中的命令按顺序复制，一行执行完再执行下一行；不要复制代码框外的文字。第一次建议使用默认圆柱，不修改网格或物性。

## 1. 准备 Ubuntu

本程序支持 Linux，使用 GCC、CMake 和 OpenMP。Windows 用户先安装 WSL：在开始菜单搜索 **PowerShell**，选择“以管理员身份运行”，输入：

```powershell
wsl --install -d Ubuntu
```

按提示重启，从开始菜单打开 **Ubuntu**，首次启动时设置用户名和密码。以后本文所有命令都在 Ubuntu 窗口执行，不能放进 PowerShell。Ubuntu 用户从下一步开始。

详细安装说明见 [Microsoft WSL 文档](https://learn.microsoft.com/en-us/windows/wsl/install)。本文命令面向 Ubuntu 22.04/24.04；计算流程在 Ubuntu 22.04、Python 3.10、GCC 11 上检查。

## 2. 安装工具并下载程序

```sh
sudo apt update
sudo apt install -y git build-essential cmake python3 curl
mkdir -p ~/cfd
cd ~/cfd
git clone https://github.com/PurineAcO/PurineCFD-R2.git
cd ~/cfd/PurineCFD-R2
```

`sudo` 要求输入刚设置的 Ubuntu 密码；输入时屏幕不显示字符，这是正常行为。`~/cfd` 表示你的 Ubuntu 用户目录里的 cfd 文件夹。第二次使用无需重复安装或 clone，直接执行最后一行 `cd` 即可。

若使用尚未合并的 draft PR，请在 GitHub 的 PR 页面查看提交分支，按 PR 描述中的试用命令下载；上面的命令下载已合并的 main。

## 3. 编译

编译是把 C++ 源文件转换为可运行程序，只需在首次使用或源码改变后执行。

```sh
cd ~/cfd/PurineCFD-R2
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build -j 4
```

看到 `Built target purinecfd` 表示成功。可执行文件位于 `build/purinecfd`。不需要 CUDA 或显卡。

## 4. 先做 200 步安装检查

```sh
cd ~/cfd/PurineCFD-R2
OMP_NUM_THREADS=4 python3 cases/cylinder/run.py --model laminar --re 47 --steps 200 --output run/install-check
```

结束时会打印 `converged=False`。200 步仅检查安装和输出，不够得到定常解。结果目录里应有 `run.log`、`parameters.json` 和 `baseflow.dat`。

命令中的 `OMP_NUM_THREADS=4` 表示使用四个 CPU 线程；机器较忙时可改成 2。该设置只影响本次命令。

## 5. 计算收敛的层流基流

```sh
cd ~/cfd/PurineCFD-R2
OMP_NUM_THREADS=4 python3 cases/cylinder/run.py --model laminar --re 47 --output run/cylinder-laminar-re47
```

运行期间窗口可能没有新增文字，进度写在日志中。另开一个 Ubuntu 窗口查看：

```sh
tail -f ~/cfd/PurineCFD-R2/run/cylinder-laminar-re47/run.log
```

在查看日志的窗口按 `Ctrl+C` 只会停止查看，不会停止另一个窗口的计算。不要关闭正在计算的窗口。

运行脚本会生成网格、根据 Re 计算来流压力、写入配置、启动计算，并将最后一个场文件复制成 `baseflow.dat`。该算例固定 Ma=0.2、T=300 K、D=1 m，网格为 128×96、远场半径 30D。达到判据后日志出现 `Converged at step`，脚本打印 `converged=True`。默认最多运行 80,000 步，速度取决于机器。

## 6. 切换 SA 或 Re

```sh
OMP_NUM_THREADS=4 python3 cases/cylinder/run.py --model sa --re 47 --output run/cylinder-sa-re47
```

`--model laminar` 是层流，`--model sa` 是 SA-RANS。改变 `--re` 后必须使用新的输出目录，例如 `--re 49 --output run/cylinder-laminar-re49`。脚本拒绝覆盖已有目录。

这是圆柱入口，不会生成 NACA 0012 网格。翼型算例需要另外准备经过核对的网格、来流和边界条件。底层求解器可用 `./build/purinecfd config.json` 读取配置，`solver.model` 必填；现有来流沿 +x，未提供攻角配置入口。

## 7. 找到结果并交给稳定性程序

| 文件 | 用途 |
| --- | --- |
| `run.log` | 进度与收敛信息 |
| `config.json` | 本次运行的设置 |
| `mesh.txt` | 本次网格 |
| `parameters.json` | Re、模型、量纲尺度与收敛标记 |
| `baseflow.dat` | 最后一个流场，文本列为 x、y、ρ、u、v、T、p、Ma、ν̃ |
| `field/step_*.dat` | 运行过程保存的场文件 |

将整个 `run/cylinder-laminar-re47` 文件夹交给稳定性程序，不要只复制 `baseflow.dat`。[继续进行稳定性分析](https://github.com/PurineAcO/crouch-hopf/blob/main/docs/quickstart.md)。

Windows 用户可在 Ubuntu 执行 `explorer.exe .` 打开当前目录，或在资源管理器访问 `\\wsl.localhost\Ubuntu\home` 找到文件。

## 常见问题

- `command not found`：确认是在 Ubuntu 中操作，并完成第 2 步。
- `Build first`：完成第 3 步，检查是否存在 `build/purinecfd`。
- `Output already exists`：换一个输出目录，避免覆盖上一组结果。
- `converged=False`：检查日志；若只是步数耗尽，用新的目录增加 `--steps`。有 `Error` 或无效流场提示时不要继续稳定性分析。
- 停止计算：在计算窗口按 `Ctrl+C`。被中断的运行不能当作收敛结果。

求助时请发送完整命令、`config.json`、`parameters.json`、`run.log` 和终端错误文字。
