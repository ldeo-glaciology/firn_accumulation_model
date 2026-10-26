#!/bin/bash
# usage: ./queue.sh QUALITY file:Class file:Class ...
Q=$1; shift
for item in "$@"; do
  f=${item%%:*}; c=${item##*:}
  ./render.sh $f $c $Q > test/log_${f%.py}.txt 2>&1 && echo "OK $f" >> test/queue_status.txt || echo "FAIL $f" >> test/queue_status.txt
done
