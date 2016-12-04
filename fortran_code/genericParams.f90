!=========================================================================
! 2. GENERICPARAMS: Define variables we need for the global optimization.
    !Note, you should not put function specific parameters here. Use
    !a separate module.
!=========================================================================
MODULE genericParams
USE nrtype
USE stateControl
IMPLICIT NONE

    INTEGER(I4B) :: p_default_alg = 0 !default algorithm. 0-bobyqa, 1-amoeba, 2-dfpmin
    INTEGER(I4B) :: p_nx, p_nmom, p_iprint=2, p_maxeval, p_ninterppt
    INTEGER(I4B):: p_qr_ndraw, p_maxpoints
    INTEGER(I4B) :: p_searchType, p_lotteryPoints
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p_wspace
    REAL(DP), DIMENSION(:), ALLOCATABLE :: p_init
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: p_range, p_bound
    INTEGER(I4B), PARAMETER  :: p_nsmplx=100
        ! initialize the seed for random numbers
    INTEGER(I4B), PARAMETER :: p_SEED=314159265
    REAL(DP), parameter :: p_tolf_amoeba= 1.0d-8 ! Amoeba variables
    INTEGER(I4B), parameter :: p_itmax_amoeba=1000 ! Amoeba variables

    INTEGER(I4B) :: fe_counter = 0
    REAL(DP) :: cpuTime = 0

CONTAINS
    SUBROUTINE initialize(seqNu, update, config)
        !Initialize the instance. This includes parsing the configuration file,
        !setting parameter values, etc.
        INTEGER(I4B), INTENT(IN) :: seqNu
        LOGICAL, INTENT(IN), OPTIONAL :: update
        CHARACTER(LEN=25), INTENT(IN), OPTIONAL :: config
        LOGICAL :: isUpdate
        INTEGER(I4B) :: i, n, fileDesc, openStat
        INTEGER, DIMENSION(:), ALLOCATABLE :: seed

        !if config is present,
        IF (PRESENT(config)) THEN
            isUpdate = .FALSE.

            IF (PRESENT(update)) THEN
                isUpdate = update
            END IF

            IF (isUpdate .EQV. .FALSE.) THEN
                !if not updating, then reset all files
                call reset()
                call myopen(unit=fileDesc, file='internalConfig.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
                call myopen(unit=fileDesc, file='sobol.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
                call myopen(unit=fileDesc, file='lastSobol.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
                call myopen(unit=fileDesc, file='lastParam.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
                call myopen(unit=fileDesc, file='sobolFnVal.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
                call myopen(unit=fileDesc, file='x_starts.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
                call myopen(unit=fileDesc, file='searchResults.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
                call myopen(unit=fileDesc, file='searchStart.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
            ELSE
                !if updating, then don't reset state stuff, lastSobol, previously
                !calculated function values, or previously used parameters
                call myopen(unit=fileDesc, file='sobol.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
                call myopen(unit=fileDesc, file='x_starts.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
                call myopen(unit=fileDesc, file='lastParam.dat',status='replace', IOSTAT=openStat, ACTION='write')
                call myclose(fileDesc)
            END IF
            !set initTerm to this one
            call setInit(seqNu)

            ! set state to 0
            IF(.not. setState(0, seqNu)) THEN
                print *,"<genericParams:initialize> Error. Unable to set state to 0 when initializing."
                stop 0
            END IF

            !initialize
            IF (isUpdate) THEN
                ! if we are updating, then store old values first
                call parseConfig(.FALSE., 'internalConfig.dat')
            END IF
            call parseConfig(isUpdate, config)

            !print to internal config file so we have values for warm starts.
            call myopen(UNIT=fileDesc, FILE='internalConfig.dat', STATUS='replace', IOSTAT=openStat, ACTION='write')
            write(fileDesc,1001) p_nx, p_nmom, p_maxeval, p_qr_ndraw, p_maxpoints-1, p_searchType, p_lotteryPoints
            DO i=1,p_nx
                write(fileDesc,*) p_range(i,1),p_range(i,2)
            END DO
            DO i=1,p_nx
                write(fileDesc,*) p_init(i)
            END DO
            DO i=1,p_nx
                write(fileDesc,*) p_bound(i,1),p_bound(i,2)
            END DO
            call myclose(fileDesc)
        ELSE
            !Otherwise, read from internal config file
            call parseConfig(.FALSE., 'internalConfig.dat')
        END IF
1001 format (7(i6,x))

        IF(allocated(p_wspace)) THEN
            deallocate(p_wspace)
        END IF
        allocate(p_wspace((p_ninterppt+5)*(p_ninterppt+p_nx)+3*p_nx*(p_nx+5)/2))

        call random_seed(size = n)
        allocate(seed(n))
        seed(1)=123456
        call random_seed(put = seed)
        deallocate(seed)
    END SUBROUTINE initialize

    SUBROUTINE parseConfig(isUpdate, configFile)
        !parse the configuration file
        LOGICAL, INTENT(IN) :: isUpdate
        CHARACTER(LEN=*), INTENT(IN) :: configFile
        INTEGER(I4B) :: openStat, i
        INTEGER(I4B) :: nx, nmom, maxeval, ninterppt
        INTEGER(I4B) :: qr_ndraw,maxpoints
        INTEGER(I4B) :: searchType, lotteryPoints, fileDesc
        CHARACTER(LEN=100) :: line
        LOGICAL :: parsedLine

        call myopen(UNIT=fileDesc, FILE=configFile, STATUS='old', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            print *,"<genericParams.parseConfig()> Error: unable to open file:",configFile
            stop 0
        END IF

        !Parse the basic parameters in the first line
        parsedLine = .FALSE.
        DO WHILE (parsedLine .eqv. .FALSE.)
            read(fileDesc,'(A100)',END=10) line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.
            read(line,*) nx, nmom, maxeval, qr_ndraw, maxpoints, searchType, lotteryPoints

            IF (maxeval < 0) THEN
                maxeval=40*(nx+1)
            END IF

            ninterppt=2*nx+1

            IF (qr_ndraw < 0) THEN
                qr_ndraw =500
            END IF

            IF (maxpoints < 0 ) THEN
                maxpoints=min(200,qr_ndraw)
            ELSE IF(maxpoints > qr_ndraw) THEN
                print *, "<genericParams:parseConfig> Error in config file. maxpoints > # sobol points"
                stop 0
            END IF

            IF (lotteryPoints < 0) THEN
                lotteryPoints = maxpoints
            END IF

        END DO
        IF (isUpdate) THEN
            !if updating, don't allow changes in nx or nmom
            IF (nx /= p_nx) THEN
                print *, "<genericParams:parseConfig> Cannot update the number of parameters.",nx,p_nx
                stop 0
            END IF

            IF (nmom /= p_nmom) THEN
                print *, "<genericParams:parseConfig> Cannot update the number of moments."
                stop 0
            END IF
        ELSE
            IF (allocated(p_range) .eqv. .FALSE.) THEN
                allocate(p_range(nx,2))
            END IF
            IF (allocated(p_init) .eqv. .FALSE.) THEN
                allocate(p_init(nx))
            END IF
            IF (allocated(p_bound) .eqv. .FALSE.) THEN
                allocate(p_bound(nx,2))
            END IF
        END IF

        p_nx=nx
        p_nmom=nmom
        p_maxeval=maxeval
        p_qr_ndraw=qr_ndraw
        p_maxpoints = maxpoints+1 !we add one, because we want to store an extra line
                                  !in the output file. This is the easiest way to do it
        p_ninterppt = ninterppt
        p_searchType = searchType
        p_lotteryPoints = lotteryPoints

        ! now parse the range for the sobol points
        parsedLine = .FALSE.
        DO WHILE ( (parsedLine .eqv. .FALSE.) .AND. (.not. isUpdate) )
            read(fileDesc,'(A100)',END=11) line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.

            read(line,*) p_range(1,1), p_range(1,2)

            DO i=2,p_nx
                read(fileDesc,*,END=11) p_range(i,1), p_range(i,2)
            END DO
        END DO

        ! now parse the initial guess for the parameters
        parsedLine = .FALSE.
        DO WHILE ( (parsedLine .eqv. .FALSE.) .AND. (.not. isUpdate) )
            read(fileDesc,'(A100)',END=12) line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.

            read(line,*) p_init(1)

            DO i=2,p_nx
                read(fileDesc,*,END=12) p_init(i)
            END DO
        END DO

        ! now parse the bounds for the parameters
        parsedLine = .FALSE.
        DO WHILE ( (parsedLine .eqv. .FALSE.) .AND. (.not. isUpdate) )
            read(fileDesc,'(A100)',END=13) line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.
            read(line,*) p_bound(1,1), p_bound(1,2)

            DO i=2,p_nx
                read(fileDesc,*,END=13) p_bound(i,1), p_bound(i,2)
            END DO
        END DO

        call myclose(fileDesc)
        return

10      print *,"<genericParams:parseConfig> Error. ",configFile," Unable to read first line of parameters"
11      print *,"<genericParams:parseConfig> Error. ",configFile," Unable to read range of values for all parameters"
12      print *,"<genericParams:parseConfig> Error. ",configFile," Unable to read initial guesses for all parameters"
13      print *,"<genericParams:parseConfig> Error. ",configFile," Unable to read bounds for all parameters"
!100     format(i10)

    END SUBROUTINE parseConfig
END MODULE genericParams
