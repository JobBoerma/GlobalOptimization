MODULE myParams
    !====================================================================
    USE nrtype
    IMPLICIT NONE

    ! Pre (post=0) or post government (post=1)
    integer(i4b), parameter :: post = 0
    character(len=*), dimension(2), parameter :: casename = (/ 'pre', 'pos'/)

    ! Choose country (1- GERMANY, 2- US, 3-SWEDEN)
    integer(i4b), parameter :: country = 3
    character(len=*), dimension(3), parameter :: countryname = (/ 'GER', 'USA', 'SWE'/)

! -------------------------------------------------

    ! First and last year in the original data sample (in levels)
    integer(i4b), dimension(3), parameter :: YearF = (/1984, 1976, 1979/), &
                                             YearL = (/2011, 2010, 2010/)

    ! Number of years for simulation (all years from beg to end, at annual freq)
    integer(i4b), parameter :: Tlev_s = YearL(country)-YearF(country)+1

    ! Number of years in levels in the data (same as simulation, except for the US)
    integer(i4b), dimension(3), parameter :: TdatC = (/Tlev_s,28,Tlev_s/)
    integer(i4b), parameter :: Tlev_d = TdatC(country)

    ! Short-run horizon
    integer(i4b), dimension(3), parameter :: L_srvec = (/1, 2, 1/)
    integer(i4b), parameter :: L_sr = L_srvec(country), L_mr=L_sr+2, L_lr=L_sr+4 

    ! Number of years to compute changes, at each horizon.
    ! Again, Tdat and Tsim should coincide except for the US
    integer(i4b), parameter :: Tdat_sr = Tlev_d-L_sr, Tdat_mr = Tlev_d-L_mr, Tdat_lr = Tlev_d-L_lr
    integer(i4b), parameter :: Tsim_sr = Tlev_s-L_sr, Tsim_mr = Tlev_s-L_mr, Tsim_lr = Tlev_s-L_lr

    ! Directory Paths (no need to change, as long as all the subdirectories are created)
    character(len=*), parameter :: rootdir = "../" 
    character(len=*), parameter :: datadir = rootdir//"data/"//countryname(country)//"/"//casename(post+1)
    character(len=*), parameter :: savedir = rootdir//countryname(country)//"out/"//casename(post+1) 

    real(dp), dimension(7):: storeTheta

    ! Simulation parameters
    integer(i4b), parameter  :: nsim = 10               ! Number of simulations
    integer(i4b), parameter  :: nhhsim = 10000!10000     ! Number of HHs per simulation
    integer(i4b), parameter  :: ntargets = 3 ! N of parameters
    integer(i4b), parameter  :: age1= 25, ageT=60,nage= ageT-age1+1 ! N of parameters
    integer(i4b), parameter :: NmomCS1 = Tsim_sr+Tsim_mr+Tsim_lr, NmomLC = 36
    integer(i4b), parameter :: NmomCS = NmomCS1*ntargets, NmomCS_all = NmomCS1*4
    integer(i4b), parameter :: Nmom = NmomCS + nage

    ! GLOBALS
    ! Read-in from matlab
    REAL(DP), DIMENSION(Tlev_s) :: xt            ! Data mean and gdp growth
    REAL(DP), DIMENSION(Tsim_sr) :: mt1            ! Data mean and gdp growth
    REAL(DP) :: VAR25            ! Data mean and gdp growth
    integer :: idsr(Tsim_sr), idmr(Tsim_mr), idlr(Tsim_lr)

    ! dimension is # moments, to be read in in cofig.txt
    REAL(DP), DIMENSION(NmomCS_all)   :: Cx_csA, gama_csA            ! Data moments - vector
    REAL(DP), DIMENSION(NmomCS)   :: Cx_cs, gama_cs            ! Data moments - vector
    REAL(DP), DIMENSION(nage)   :: Cx_lc, gama_lc            ! Data moments - vector
    REAL(DP), DIMENSION(Nmom)       :: Cx, gama            ! Data moments - targetted
    
    ! random vars
    REAL(DP), DIMENSION(nhhsim*5,Tlev_s,nsim) :: randva ! matrix of all random vars to be used
    


END MODULE myParams
