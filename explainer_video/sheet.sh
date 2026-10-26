#!/bin/bash
# contact sheet of frames every N seconds: ./sheet.sh video.mp4 out.png [every_s]
EVERY=${3:-4}
ffmpeg -v error -y -i "$1" -vf "fps=1/$EVERY,scale=427:-1,tile=4x5:padding=4" -frames:v 1 "$2"
