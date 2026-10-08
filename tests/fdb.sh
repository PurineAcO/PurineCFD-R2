# fdb means fast-debug-build.sh
# please run this bash in dir PurineCFD-R2
# 用法: bash tests/fdb.sh [配置] [步数] [inner] [CFL] [线程数] [中间帧间隔]
# 缺省跑 OAT15A 3.1° 的隐式求解检验台, 结果写到 debugres/
# 例: bash tests/fdb.sh cases/oat15a/urans_m073a31_re33e6.json 10 40
set -e

cd "$(dirname "$0")/.."

CFG="${1:-cases/oat15a/urans_m073a31_re33e6.json}"

# 首次运行需要 configure, 之后这行直接跳过
[ -f debug/CMakeCache.txt ] || cmake -S . -B debug

cmake --build debug --target debug

# 不删: 上一轮的结果挪到 debugres/_prev, 随时可以翻回去看
mkdir -p debugres/_prev
mv -f debugres/*.txt debugres/*.csv debugres/_prev/ 2>/dev/null || true
./debug/debug "$CFG" "${@:2}"

echo "==================== 结果 ===================="
echo "debugres/summary.txt  逐项检验与标定结论"
echo "debugres/calib.csv    CFL × 内迭代收敛曲线 (cfl,step,k,res,ratio)"
