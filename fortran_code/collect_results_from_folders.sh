# collect results

#rm -f ../results/results_time.dat
#rm -f ../results/results_min.dat

#touch ../results/results_time.dat
#touch ../results/results_min.dat

declare -i length n
#declare -a SobolList
#SobolList=(250 500 750 1000 1500 2000)
SobolList=(4000 5000)

length=${#SobolList[@]}
length=$((length-1))

for j in $(seq 0 $length); do
 n=${SobolList[$j]}
	touch ../results/min_bobyqa_sob$n.dat
	touch ../results/time_bobyqa_sob$n.dat
#	touch ../results/min_amoeba_sob$n.dat
#	touch ../results/time_amoeba_sob$n.dat
 for i in $(seq 1 25); do
   cat ../results/bobyqa_sobol$n/startpoint_$i/fortran_code/monteCarloGOPAmin.dat >> ../results/min_bobyqa_sob$n.dat
   cat ../results/bobyqa_sobol$n/startpoint_$i/fortran_code/monteCarloGOPA.dat >> ../results/time_bobyqa_sob$n.dat
#   cat ../results/amoeba_sobol$n/startpoint_$i/fortran_code/monteCarloGOPAmin.dat >> ../results/min_amoeba_sob$n.dat
#   cat ../results/amoeba_sobol$n/startpoint_$i/fortran_code/monteCarloGOPA.dat >> ../results/time_amoeba_sob$n.dat
 done
done
