#!/bin/bash
# usage: time_run.sh <fa/anc-base-path-without-ext> <workdir> <resultlog>
base=$1; work=$2; log=$3
mkdir -p "$work"; cd "$work"
cp "${base}.fa" "${base}.anc" .
name=$(basename "$base")
cp /home/claude/work/dialign-mmap/dialign-mmap/src/dialign2-2 ./d
export DIALIGN2_DIR=/home/claude/work/dialign_package/dialign2_dir
t0=$(date +%s.%N)
./d -anc -n "${name}.fa" > run.out 2> run.err &
pid=$!
peak=0
while kill -0 $pid 2>/dev/null; do
  rss=$(awk '/VmRSS/{print $2}' /proc/$pid/status 2>/dev/null)
  [ -n "$rss" ] && [ "$rss" -gt "$peak" ] && peak=$rss
  sleep 0.2
done
wait $pid; rc=$?
t1=$(date +%s.%N)
dt=$(python3 -c "print(f'{$t1-$t0:.2f}')")
echo "$name rc=$rc wall_s=$dt peak_rss_kb=$peak" | tee -a "$log"
