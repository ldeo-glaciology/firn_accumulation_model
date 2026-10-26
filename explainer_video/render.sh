#!/bin/bash
# usage: ./render.sh <scene_file.py> <SceneClass> [ql|qm|final]
set -e
cd "$(dirname "$0")"
Q=${3:-ql}
export PYTHONPATH="$(pwd):$PYTHONPATH"
if [ "$Q" = "final" ]; then
  FLAGS="-r 1920,1080 --frame_rate 30"
else
  FLAGS="-$Q"
fi
/home/user/venv/bin/manim $FLAGS --media_dir media scenes/$1 $2
