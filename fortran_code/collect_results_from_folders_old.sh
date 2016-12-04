# collect results

rm -f ../results/results_time.dat
rm -f ../results/results_min.dat

touch ../results/results_time.dat
touch ../results/results_min.dat

declare -i length n imp
#declare -a SobolList
SobolList=(100 200 300 400 500 600 700 800 900 1000 1250 1500 1750 2000)
length=${#SobolList[@]}
length=$((length-1))
imp=50   # number of implementations

for j in $(seq 0 $length); do
 n=${SobolList[$j]}
 for i in $(seq 1 $imp); do
   cat ../results/bobyqa_sobol$n/startpoint_$i/fortran_code/monteCarloGOPAmin.dat >> ../results/results_min.dat
   cat ../results/bobyqa_sobol$n/startpoint_$i/fortran_code/monteCarloGOPA.dat >> ../results/results_time.dat
 done
done
