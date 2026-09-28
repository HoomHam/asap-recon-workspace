#!/bin/zsh
# Scratch batch (s14): soft-bin stability, lower-SNR sweep, registration proxies. One stream per subject.
KEY=$1; SIG=$2
PY=/opt/homebrew/Caskroom/miniforge/base/bin/python
run() { $PY tune_gas.py --key $KEY --sigma $SIG "$@" 2>&1 | grep --line-buffered -E "^[0-9:]{8} (prod|gplb|MS|k[0-9]|done|non-prod|noise x)|Traceback|Error" | sed "s/^/[$KEY $*] /"; }
run --set reg
for bd in 0.5 1 4; do run --set stab --bindist0sq $bd; done
run --set stab --nbins 8 --bins 0,2,4,6
run --set stab --nbins 24 --bins 0,6,12,18
for nm in 1 1.5 2 3 4; do run --set snr --noise-mult $nm; done
echo "[$KEY] ALL DONE"
