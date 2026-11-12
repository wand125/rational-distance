# mw5 (mathematics-52): other pairs for 569/720 and 769/360 (pair-1 curves have rank 2 that mwrank cannot prove)
cd ~/52_trial/circle_search
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 MW_DUMP=mw5_11101225/dump_$1.txt
PY=$(cat mw5_11101225/PY)
shift
for rj in "$@"; do
  set -- ${=rj//,/ }
  nice -n 10 zsh mw5_11101225/to.sh 1500 $PY mw_one.py $1 $2 3000 2 > mw5_11101225/one_${1//\//_}_$2.txt 2>> mw5_11101225/err.txt
  l=$(tail -n 1 mw5_11101225/one_${1//\//_}_$2.txt 2>/dev/null)
  [[ -n $l ]] && echo $l >> mw5_11101225/all.txt || echo "$1 $2 timeout" >> mw5_11101225/all.txt
done
