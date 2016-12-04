#generate random starting points using bounds from config file
#gfortran GOPA_init_point.f90 -o initpoint
#./initpoint


###########################################################################################################
echo "housekeeping..."

# Remove any locked file
rm -f *lock 
# Remove old libraries
rm -f *.mod
# Remove Montecarlo dat file
rm -f monteCarloGOPAmin.dat
rm -f monteCarloGOPA.dat
rm -f startingPoints.dat


###########################################################################################################
echo "compiling..."

#ifort -O3 -g -heap-arrays nrtype.f90 myParams.f90 stateControl.f90 genericParams.f90 utilities.f90 simplex.f90 objective.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch
#ifort -m64 -g -debug all -heap-arrays nrtype.f90 stateControl.f90 genericParams.f90 simplex.f90 global.f90 objective_griewank.f90 utilities.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch
gfortran -O3               nrtype.f90 myParams.f90 stateControl.f90 genericParams.f90 utilities.f90 simplex.f90 objective.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch

if [ $? -ne 0 ]; then
  echo "Errors compiling GlobalSearch"
  exit
fi


###########################################################################################################
echo "run Monte Carlo..."

#need to allow modifying txt file, otherwise cannot rewrite
chmod u+x config.txt

#declare variables
declare -i n m dim lineout linein imp
#declare -a SobolList
dim=7
SobolList=(100 200 300 400 500 600 700 800 900 1000 1250 1500 1750 2000)
length=${#SobolList[@]}
length=$((length-1)) # length of SobolList (minus 1 because indexation starts at 0)
imp=50   # number of implementations
queue="shared"
timelimit="300:00"

# create folder results
mkdir ../results

# create one folder by algo and number of Sobol points
#for j in 250 500 750 1000 1500 2000
for j in $(seq 0 $length); do
 n=${SobolList[$j]}
 mkdir ../results/bobyqa_sobol$n
 mkdir ../results/amoeba_sobol$n
done

# create one folder by starting point
for j in $(seq 0 $length); do
 n=${SobolList[$j]}
 for i in $(seq 1 $imp); do
  mkdir ../results/bobyqa_sobol$n/startpoint_$i
  mkdir ../results/amoeba_sobol$n/startpoint_$i
 done
done


# copy everything in each folder and change config file appropriately
for i in $(seq 1 $imp); do
	
	# overwrite starting point on config file 
	#(using XX random starting points but always the same point for each step in performance profile)
	for k in $(seq 1 $dim); do
	  lineout=37+$dim+$k
	  linein=($i-1)*$dim+$k
	  sed -i "${lineout}s/.*/`sed -n "${linein}p" init_points.dat`/g" config.txt
	done

	# varying number of Sobol points
	for j in $(seq 0 $length); do
 		n=${SobolList[$j]}   # number of sobols generated
	 	m=$n/10	# number of sobols kept	

	 #################################### bobyqa ###########################################
	 # copy files in new folder
	 cd ..
	 cp -r fortran_code results/bobyqa_sobol$n/startpoint_$i 			#	 cp -r ./* ../results/bobyqa_$j
	 cp -r data results/bobyqa_sobol$n/startpoint_$i
	 cp -r SWEout results/bobyqa_sobol$n/startpoint_$i

	 # go to new folder
	 cd results/bobyqa_sobol$n/startpoint_$i/fortran_code

      # overwrite number of Sobol points generated / kept on config file
	echo "bobyqa dim7 sobol is $n, kept is $m and iter is" $i
	sed -i "28s/.*/7, 297, -1, ${n[@]}, ${m[@]}, 0, -1/g" config.txt
	 	 
	# go back to original folder
	cd ../../../../fortran_code



	 #################################### amoeba ###########################################
	 # copy files in new folder
	 cd ..
	 cp -r fortran_code results/amoeba_sobol$n/startpoint_$i
	 cp -r data results/amoeba_sobol$n/startpoint_$i
	 cp -r SWEout results/amoeba_sobol$n/startpoint_$i

	 # go to new folder
	 cd results/amoeba_sobol$n/startpoint_$i/fortran_code

      # overwrite number of Sobol points generated / kept on config file
	echo "amoeba dim7 sobol is $n, kept is $m and iter is" $i
	sed -i "28s/.*/7, 297, -1, ${n[@]}, ${m[@]}, 0, -1/g" config.txt
	 	 
	# go back to original folder
	cd ../../../../fortran_code


     done
done

cd ../results

# run optimimization in each folder
for i in $(seq 1 $imp); do
	for j in $(seq 0 $length); do
 		n=${SobolList[$j]}
	
	 # go to new folder
	 cd bobyqa_sobol$n/startpoint_$i/fortran_code

	# run minimization with Bobyqa algo on Yale cluster
	bsub -q $queue -W $timelimit ./GlobalSearch 0 config.txt b &
	 
	# go back to folder results
	cd ../../../

	# go to new folder
	cd amoeba_sobol$n/startpoint_$i/fortran_code

	# run minimization with Amoeba algo on Yale cluster
	bsub -q $queue -W $timelimit ./GlobalSearch 0 config.txt a &
	 
	# go back to folder results
	cd ../../../

     done

done


# to print the results in a file, use the following after on the previous line:
#> simulations_output_$(date +%Y%m%d).txt

echo "done. check simulations_output file"