# fdb means fast-debug-build.sh
# please run this bash in dir PurineCFD-R2


set -e
cd "$(dirname "$0")/.."
[ -f debug/CMakeCache.txt ] || cmake -S . -B debug
cmake --build debug --target debug
echo ====================================================
./debug/debug