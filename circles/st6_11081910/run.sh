#!/bin/zsh
# st6 (mathematics-52): on the 267 open circles, (a) refined descent with 900 s per descent (stage_sub.py, pmax 1e6, four-factor last),
# (b) point search with t height <= 3*10^6 (circ_ptsearch.gp). 4 + 4 shards, 1 thread each, nice 10. DONE when all end.
cd ~/52_trial/circle_search
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
P=../venv/bin/python; D=st6_11081910
for s in 0 1 2 3; do nice -n 10 $P stage_sub.py $D/open267.lst 4 $s $D/de_sh$s.txt 1000000 1 900 > $D/de_sh$s.log 2>&1 & done
for s in 0 1 2 3; do nice -n 10 $P rungp.py $D/ps_sh$s.gp > $D/ps_sh$s.log 2>&1 & done
wait
touch $D/DONE
