#!/bin/zsh
# st5 (mathematics-52): redo of st4 with stage_sub.py (each descent in a child process, killed after 120 s).
cd ~/52_trial/circle_search
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
P=../venv/bin/python; D=st5_11081121
for s in 0 1; do nice -n 10 $P stage_sub.py $D/re406.lst 2 $s $D/re406_sh$s.txt 400 0 120 > $D/re406_sh$s.log 2>&1 & done
for s in 0 1 2 3 4 5; do nice -n 10 $P stage_sub.py $D/open897.lst 6 $s $D/st3_sh$s.txt 1000000 1 120 > $D/st3_sh$s.log 2>&1 & done
wait
touch $D/DONE
