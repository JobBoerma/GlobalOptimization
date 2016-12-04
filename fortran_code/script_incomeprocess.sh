#generate random starting points using bounds from config file
#gfortran GOPA_init_point.f90 -o initpoint
#./initpoint

# if there is any locked file
rm -f *lock 
# Remove old libraries
rm -f *.mod
# Remove Montecarlo dat file
rm -f monteCarloGOPAmin.dat
rm -f monteCarloGOPA.dat
rm -f startingPoints.dat


###########################################################################################################
echo "Compiling"
#ifort -O3 -g -heap-arrays nrtype.f90 myParams.f90 stateControl.f90 genericParams.f90 utilities.f90 simplex.f90 objective.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch
#ifort -m64 -g -debug all -heap-arrays nrtype.f90 stateControl.f90 genericParams.f90 simplex.f90 global.f90 objective_griewank.f90 utilities.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch
gfortran -O3               nrtype.f90 myParams.f90 stateControl.f90 genericParams.f90 utilities.f90 simplex.f90 objective.f90 minimize.f90 GlobalSearch.f90 -o GlobalSearch

