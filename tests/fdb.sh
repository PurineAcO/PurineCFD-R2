# fdb means fast-debug-build.sh
# please run this bash in dir PurineCFD-R2
# 用法: bash tests/fdb.sh [步数] [子迭代] [dt] [CFL] [rans] [laminar|trace] [线程数] [中间帧间隔]
#       bash tests/fdb.sh plot     # 不重跑, 只用 debugres 里已保存的流场出图
set -e

cd "$(dirname "$0")/.."

plot_debugres() {
    if python3 -c "import matplotlib" 2>/dev/null; then
        python3 tests/plot_field.py debugres
    elif command -v uv >/dev/null 2>&1; then
        uv run --offline --with matplotlib python tests/plot_field.py debugres
    else
        python3 tests/plot_field.py debugres
    fi
}

if [ "${1:-}" = "plot" ]; then
    plot_debugres
    exit 0
fi

# 首次运行需要 configure, 之后这行直接跳过
[ -f debug/CMakeCache.txt ] || cmake -S . -B debug

cmake --build debug --target debug

# 不删: 上一轮的帧挪到 debugres/_prev, 随时可以翻回去看
mkdir -p debugres/_prev
mv -f debugres/field_*.dat debugres/field_*.png debugres/forces_*.csv \
      debugres/wall_cp*.csv debugres/_prev/ 2>/dev/null || true
./debug/debug "$@"

echo "==================== 流场出图(只读已保存的文件) ===================="
plot_debugres

echo "==================== python 复核 ===================="
if python3 -c "import matplotlib" 2>/dev/null; then
    python3 tests/verify_debugres.py debugres
elif command -v uv >/dev/null 2>&1; then
    uv run --offline --with matplotlib python tests/verify_debugres.py debugres
else
    python3 tests/verify_debugres.py debugres
fi
