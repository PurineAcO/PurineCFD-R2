#!/usr/bin/env bash
set -e
cd "$(dirname "$(readlink -f "$0")")"
if [ -z "${DISPLAY:-}${WAYLAND_DISPLAY:-}" ]; then
  echo "无图形环境: DISPLAY 与 WAYLAND_DISPLAY 均未设置" >&2
  exit 1
fi
exec uv run --with pyside6 python gui/main.py
