PROGRAM GlobalSearch
    ! The main program. It is written as a pseudo state machine, and
    ! multiple instances can be run at once. The main driver (which is
    ! the initially the cold start) will set states for all other
    ! instances of the program.
    !
    ! The states are as follows:
    !     -1: exit state - informs all processes to terminate
    !      1: initial state. Initializes all variables, parameters, etc.
    !           Prev States: None
    !           Next States: 2 (main driver), 3 (others)
    !      2: This state generates the sobol points. Only the "main driver"
    !         program can be in this state.
    !           Prev States: 1
    !           Next States: 3
    !      3: This state begins solving the objective function at the sobol
    !         points.
    !           Prev States: 1, 2
    !           Next States: 8
    !      4: This state tries to minimize the objective function starting at
    !         the smallest calculated sobol value, and iterating through the
    !         number of points specified in the config file.
    !           Prev States: 6
    !           Next States: 5
    !      5: This state checks if any minimization points have been missed.
    !         If so, make note of the missing points and minimize from them.
    !           Prev State: 4
    !           Next States: 4 (if missing points), 7(otherwise)
    !      6: This state sorts the sobol points by minimum objective function
    !         value
    !           Prev States: 8
    !           Next States: 4
    !      7: This state takes the minimum objective value found, and restarts
    !         the minimization function at this point
    !           Prev States: 5
    !           Next States: -1
    !      8: This state checks if any sobol points have not been calculated.
    !         If this is true, make note of the missing points and calculate
    !         them. Otherwise, continue through the process
    !           Prev States: 3
    !           Next States: 3 (if missing points), 6 (otherwise)
    USE nrtype
    USE genericParams
    USE myParams
    USE stateControl
    USE utilities
    USE minimize
    USE OBJECTIVE, only : objFun, dfovec, initial0 !,GetData

    IMPLICIT NONE

    ! command line options
    LOGICAL ::isWarm           ! if this is a warm start or not
    LOGICAL :: updatePoints     ! if this instance invokation is to update the points
    INTEGER(I4B) :: alg         ! the default algorithm this instance will use for
                                ! minimizing the objective function

    ! temporary variables
    LOGICAL :: complete         ! if the number of points to evaluate (either sobol or
                                ! minimization) has been completed
    CHARACTER(LEN=1000) :: errorString  ! variable used for writing error messages
    INTEGER(I4B) :: temp, temp2, i, driver, newDriver, currentState, seqNo
    CHARACTER(LEN=25) :: arg1, config
    CHARACTER(LEN=1) :: arg2
    REAL(DP), DIMENSION(:), ALLOCATABLE  :: qrdraw
    REAL(DP), DIMENSION(:,:), ALLOCATABLE :: trial

    ! variables for file manipulation
    INTEGER(I4B) :: fileDesc, iostat, ierr

    ! For timing
    REAL(SP):: start, finish

    ! added by AA TK
    REAL(DP), DIMENSION(:) , ALLOCATABLE	:: bestParam
    REAL(DP)			:: bestValue

    call cpu_time(start)

    alg = p_default_alg
    print*,"calling initial0"
    call initial0
    call parseCommandLine()

    ! Random sequence
    !allocate(Cx(p_nmom),gama(p_nmom))
    !call GetData

    IF (isWarm .eqv. .FALSE.) THEN
        IF (updatePoints .eqv. .FALSE.) THEN
            ! If this is a true cold start, then reset the sequence
            seqNo = 1
        ELSE
            !if we are updating sobol points, we need to make sure that
            !the current set of sobol points are already generated
            call waitState(3)
            seqNo = getSequenceNumber()
        END IF

        !If cold start or update, then initialize with config file (mandatory parameter)
        call initialize(seqNo, updatePoints, config)

        ! we have a bit of a timing issue. Need to initialize with sequence number 1,
        ! but then we don't create the seq.dat file. So this takes care of it. Ugly,
        ! but don't know a better way just yet.
        IF(seqNo == 1) THEN
            seqNo = getSequenceNumber()
        END IF

        !now let's begin the state machine
        !Note: setState actually changes the state, and returns if it
        !      was successful. So, this is not a "pure" function
        IF (setState(1, seqNo) .eqv. .FALSE.) THEN
            write(errorString,*) "<main>: Error. Cannot set state to 1 in cold start."
            call writeToLog(errorString)
            stop 0
        END IF
    ELSE
        ! we are doing a warm start. Let's make sure that all the cold state initialization
        ! is complete before trying to get a new sequence number.
        call waitState(3)
        seqNo = getSequenceNumber()
        call initialize(seqNo)
    END IF

    !keep track of which sequence is the "driver" state when we initialized
    !need to know this because a future update could require reallocation
    !of trial (sobol points) or other variables
    driver = getInit()
    allocate(qrdraw(p_nx), trial(p_nx, p_qr_ndraw))

    ! addded AA TK
    allocate(bestParam(p_nx))

    IF (setState(2, seqNo)) THEN
        !If this process is able to set the state, then it must be the "driver"
        !So generate sobol points, etc.
        call insobl(p_nx, p_qr_ndraw)
        DO i = 1,p_qr_ndraw
            CALL I4_SOBOL(p_nx,qrdraw)
            trial(:,i) = p_range(:,1)+qrdraw*(p_range(:,2)-p_range(:,1))
        END DO
        call myopen(unit=fileDesc, file='sobol.dat',status='replace',IOSTAT=iostat, ACTION='write') ! save initial trial points
        write(fileDesc,*) trial
        call myclose(fileDesc)
    ELSE
        ! not driving. Let's wait until points are generated so we can start
        ! evaluating the functions at those points.
        call waitState(3)
    END IF

    ! now, let's go through the state machine
    currentState = 0
    DO WHILE(currentState .NE. p_exitState)

        ! If the driver has been updated, that means we've updated
        ! points so we need to reinitialize
        newDriver = getInit()
        IF (driver /= newDriver) THEN
            call initialize(seqNo)
            driver = newDriver
            deallocate(trial)
            allocate(trial(p_nx, p_qr_ndraw))
        END IF

        ! Get the current state.
        currentState = getState()

        ! perform the action as per the current state. The details of each state
        ! are given at the top of this file.
        SELECT CASE (currentState)
            CASE (p_exitState)
                cycle
            CASE (1)
                !we must have updated the number of sobol points. So,
                !let's wait until initialization is done
                call waitState(3)
            CASE (2)
                IF (setState(3,seqNo)) THEN
                    !we are the driver, and we have now generated all the sobol
                    !points. Lets move onto evaluating them.
                    cycle
                ELSE
                    !Otherwise, let's wait until we can start solving.
                    call waitState(3)
                END IF
            CASE (3)
                !Now everyone should be solving the function at the
                !sobol points.
                call solveAtSobolPoints(seqNo, complete)
                IF (complete) THEN
                    ! If we have completed solving the objective function
                    ! at each sobol point, then check  if everything is
                    ! done.
                    IF(setState(8,seqNo)) THEN
                        cycle
                    ELSE
                        call statePause(currentState)
                    END IF
                END IF
            CASE (4)
                ! all instances search for a minimum over the given point
                call search(seqNo, alg, complete)
                IF (complete) THEN
                    ! If we have minimized at every point, look for any
                    ! missing points
                    IF (setState(5,seqNo)) THEN
                        cycle
                    ELSE
                        call statePause(currentState)
                    END IF
                END IF
            CASE (5)
                IF (isInit(seqNo)) THEN
                    ! have only one instance (the driver) search for the missing
                    ! points.
                    call findMissingSearch(complete)
                    IF (complete) THEN
                        ! if everything was completed, then continue on to the
                        ! last minimization.
                        IF(setState(7,seqNo)) THEN
                            cycle
                        END IF
                    ELSE
                        ! otherwise, go back and minimize at the points that
                        ! were not solved.
                        IF(setState(4,seqNo)) THEN
                            cycle
                        END IF
                    END IF
                ELSE
                    call statePause(currentState)
                END IF
            CASE (6)
                ! sort the sobol points - the algorithm will use the best n points
                ! as specified in the configuration file.
                IF(isInit(seqNo)) THEN
                    call setupSobol()
                ELSE
                    call statePause(currentState)
                END IF
                IF(setState(4,seqNo)) THEN
                    cycle
                END IF
            CASE (7)
                ! perform the last optimization
                IF(setState(p_exitState,seqNo)) THEN
                    call lastSearch(seqNo)
                    !call cpu_time(finish)
                    !cpuTime = finish-start
                    !open(unit = 111, file = "monteCarloGOPA.dat", position = "append", STATUS='unknown')
                    !we write out algorithm indicator (writes 1 or 0 where 1=amoeba and 0=bobyqa, then we write function evaluation counter, cputime,
                    !number of sobols that are kept to start local searches, and number of sobols overall that are drawn.
                    !write(111,*) alg, fe_counter, cpuTime, p_maxpoints, p_qr_ndraw, nsim, nhhsim, p_nx
                    !close(111)
                ELSE
                    call statePause(currentState)
                END IF
            CASE (8)
                IF (isInit(seqNo)) THEN
                    ! find any missing points at which the algorithm hasn't solved
                    ! and store them.
                    call findMissingSobol(complete)
                    IF (complete) THEN
                        ! if nothing is missing, continue to the next state
                        ! which performs the minimization algorithm.
                        IF(setState(6,seqNo)) THEN
                            cycle
                        END IF
                    ELSE
                        ! otherwise, return to solving at each of the sobol points
                        IF(setState(3,seqNo)) THEN
                            cycle
                        END IF
                    END IF
                ELSE
                    call statePause(currentState)
                END IF
            CASE DEFAULT
                write(errorString, *) "<main> : Error, unknown state: ", currentState
                call writeToLog(errorString)
                stop 0
        END SELECT
        currentState = getState()
    END DO
    write(errorString, *) seqNo," has no more tasks."
    call writeToLog(errorString)

    call cpu_time(finish)
    cpuTime = finish-start
    open(unit = 111, file = "monteCarloGOPA.dat", position = "append", STATUS='unknown')
    !we write out algorithm indicator (writes 1 or 0 where 1=amoeba and 0=bobyqa, then we write function evaluation counter, cputime,
    !number of sobols that are kept to start local searches, and number of sobols overall that are drawn.
    write(111,*) alg, fe_counter, cpuTime, p_maxpoints, p_qr_ndraw, nsim, nhhsim, p_nx, p_init(1)
    close(111)
contains
    SUBROUTINE solveAtSobolPoints(seqNo, complete)
        !This routine solves the given function at the sobol point
        !indicated in lastSobol.dat If we have solved all points,
        !then complete is returned as true.
        implicit none
        INTEGER(I4B), INTENT(IN) :: seqNo
        LOGICAL, INTENT(OUT) :: complete
        INTEGER(I4B) :: whichPoint, openStat
        REAL(DP) :: fval

        whichPoint = getSobolPoint()

        IF (whichPoint < 0) THEN
            complete = .TRUE.
        ELSE
            complete = .FALSE.

            write(*, *) seqNo," solving sobol point ",whichPoint
            write(errorString, *) seqNo," solving sobol point ",whichPoint
            call writeToLog(errorString)

            call myopen(unit=fileDesc, file='sobol.dat',status='old', IOSTAT=openStat, ACTION='read') ! save initial trial points
            read(fileDesc,*) trial
            call myclose(fileDesc)

            fval=objFun(trial(:,whichPoint))
            call myopen(unit=fileDesc, file='sobolFnVal.dat', STATUS='old', IOSTAT=openStat, ACTION='write', position='append')
            write(fileDesc,7000) whichPoint, fval, trial(:,whichPoint)
            call myclose(fileDesc)
        END IF

7000    format(i6, 200f40.20)
    END SUBROUTINE solveAtSobolPoints

    SUBROUTINE search(seqNo, algor, complete)
        !This routine searches for a minimum at the next
        !point, as given by lastParam.dat, using the algorithm
        !specified in algor. If we have gone through all the
        !parameter guesses, complete is set to TRUE. Otherwise,
        !it is set to .FALSE.
        implicit none
        INTEGER(I4B), INTENT(IN) :: seqNo, algor
        LOGICAL, INTENT(OUT) :: complete
        INTEGER(I4B) :: i,k, whichPoint, lotPoint
        REAL(DP), DIMENSION(:,:), ALLOCATABLE :: x_starts
        REAL(DP), DIMENSION(p_nx) :: evalParam

        allocate(x_starts(p_nx+1,p_maxpoints))

        !get which point we want to evaluate
        whichPoint = getParamPoint()

        IF (whichPoint < 0) THEN
            complete = .TRUE.
        ELSE
            complete = .FALSE.

            !read the starting points. We need to do this every time, because we
            !could have updated the sobol points in which case this would have
            !changed
            call myopen(UNIT=fileDesc, FILE='x_starts.dat', STATUS='old', IOSTAT=iostat, ACTION='read')
            DO i=1,p_maxpoints-1
                read(fileDesc,*) (x_starts(k,i), k=1,p_nx+1)
            END DO
            call myclose(fileDesc)

            !get parameter value based on previous solutions. Note there are multiple search types,
            !   0: the best point
            !   1: a lottery using n points (where n is given in the config file)
            evalParam = getModifiedParam(whichPoint, x_starts(2:,whichPoint), p_searchType, lotPoint)
            write(*, *) seqNo," searching using sobol point ",whichPoint," and best point ",lotPoint
            write(errorString, *) seqNo," searching using sobol point ",whichPoint," and best point ",lotPoint
            call writeToLog(errorString)

            !let's store the best sobol point, just so we focus our searches
            !do this after getModifiedParam so that on first execution we can use initial
            !guess
            IF (whichPoint == 1) THEN
                i=0
                call myopen(UNIT=fileDesc, FILE='searchResults.dat', STATUS='unknown', IOSTAT=ioStat,&
                    ACTION='write',position='append')
                write(fileDesc,7452) seqNo, i, x_starts(1,1), x_starts(2:,1)
                call myclose(fileDesc)

                call myopen(UNIT=fileDesc, FILE='searchStart.dat', STATUS='unknown', IOSTAT=ioStat, &
                    ACTION='write',position='append')
                write(fileDesc,7452) seqNo, i, x_starts(2:,1)
                call myclose(fileDesc)
            END IF

            !We now have the point at which to solve. So, finish the search
            call completeSearch(seqNo, algor, whichPoint, evalParam)
        END IF

        deallocate(x_starts)
        return

7452    format(2i10, 200f40.20)

    END SUBROUTINE search

    SUBROUTINE completeSearch (seqNo, algor, whichPoint, evalParam)
        !This subroutine completes the search given a specific point at which to evaluate.
        INTEGER(I4B), INTENT(IN) :: seqNo, algor, whichPoint
        INTEGER(I4B) :: openStat
        REAL(DP), DIMENSION(p_nx), INTENT(INOUT) :: evalParam
        REAL(DP) :: itratio, rhobeg, rhoend, fn_val
        REAL(DP), DIMENSION(p_nx) :: evalBack
        REAL(DP), DIMENSION(p_nx+3) :: currentBest

        evalBack = evalParam

        SELECT CASE (algor)
        !minimize according to the algorithm
            CASE (0)
                !search
                itratio=REAL(whichPoint)/REAL(p_maxpoints)
                IF (itratio>0.5) THEN
                    rhobeg  = (minval(p_bound(:,2)-p_bound(:,1))/2.5_DP)/(4*itratio)
                    rhoend  = 1.0D-3/(4*itratio)
                ELSE
                    rhobeg = minval(p_bound(:,2)-p_bound(:,1))/2.50_DP
                    rhoend  = 1.0D-3
                END IF
                call bobyqa_h(p_nx,p_ninterppt,evalParam,p_bound(:,1),p_bound(:,2),rhobeg,rhoend,p_iprint,p_maxeval,p_wspace,p_nmom)
            CASE (1)
                call runAmoeba(evalParam,p_tolf_amoeba,p_itmax_amoeba)
            CASE DEFAULT
                write(errorString,*) "<main:search> Error, unknown search method: ",algor
                call writeToLog(errorString)
                stop 0
        END SELECT
        fn_val = objFun(evalParam)

        !output result
        currentBest = getBestLine(openStat)
        i = INT(currentBest(2))
        IF (fn_val < currentBest(3)) THEN
            i = i+1
        END IF

        !save the results
        call myopen(UNIT=fileDesc, FILE='searchResults.dat', STATUS='unknown', IOSTAT=openStat, ACTION='write',position='append')
        write(fileDesc,7453) seqNo, i, fn_val, evalParam
        call myclose(fileDesc)

        call myopen(UNIT=fileDesc, FILE='searchStart.dat', STATUS='unknown', IOSTAT=openStat, ACTION='write',position='append')
        write(fileDesc,7453) seqNo, whichPoint, evalBack
        call myclose(fileDesc)

7453    format(2i10, 200f40.20)

    END SUBROUTINE completeSearch

    SUBROUTINE runAmoeba(amoeba_pt,tolf,itmax)
        !This subroutine executes the amoeba search. It takes the point passed in
        !as a parameter, and generates a simplex centered on it. The points of the
        !simplex are proportional to the number of sobol points generated and the
        !dimension of the parameter space. For example, if we have 2 parameters and
        !100 sobol points, then the simplex distance is approximately 1/10 of the
        !range of each parameter.
        use simplex, only : simplex_coordinates2
        IMPLICIT NONE
        REAL(DP), intent(in) :: tolf
        INTEGER(I4B), intent(in) :: itmax
        INTEGER(I4B) :: ia
        REAL(DP), DIMENSION(p_nx), intent(inout) :: amoeba_pt
        REAL(DP), DIMENSION(p_nx,p_nx+1) ::  x
        REAL(DP), DIMENSION(p_nx+1,p_nx) :: xT
        REAL(DP), DIMENSION(p_nx+1) :: fnVals
        REAL(DP), DIMENSION(p_nx) :: temp

        call simplex_coordinates2 ( p_nx, x )

        do ia=1,p_nx+1
            temp = amoeba_pt+x(1:p_nx,ia)*(p_range(:,2)-p_range(:,1))/(DBLE(p_maxpoints)**(1.0_dp/p_nx))
            temp=max(min(temp,p_bound(:,2)),p_bound(:,1))
            fnVals(ia)=objFun(temp)
            xT(ia,:)=temp
        end do

        ia=itmax
        call amoeba(xT,fnVals,tolf,objFun,ia)
        ia=minloc(fnVals,1)
        amoeba_pt = xT(ia,:)
    END SUBROUTINE runAmoeba


    SUBROUTINE lastSearch(seqNo)
        !This routine solves the given function for the last time,
        !taking as the initial point the best point we have found
        !so far. It always uses the same algorithm, bobyqa.
        implicit none
        INTEGER(I4B), INTENT(IN) :: seqNo
        INTEGER(I4B) :: i, openStat
        REAL(DP), DIMENSION(p_nx) :: evalParam
        REAL(DP) :: rhobeg, rhoend, fn_val

	! added by AA and TK to keep result of last search only if better than best point so far
!	REAL(DP), DIMENSION(p_nx) 	:: bestParam
!	REAL(DP)			:: bestValue

        evalParam = getBestPoint(openStat)

	! added by AA and TK to keep result of last search only if better than best point so far
    	bestParam = evalParam
	    bestValue = objFun(bestParam)

        call myopen(UNIT=fileDesc, FILE=savedir//'/THETASOL.dat', STATUS='replace', IOSTAT=openStat, ACTION='write')
        write(fileDesc,*) evalParam
        call myclose(fileDesc)

        i=-1
        call myopen(UNIT=fileDesc, FILE='searchStart.dat', STATUS='unknown', IOSTAT=openStat, ACTION='write',position='append')
        write(fileDesc,7451) seqNo, i, evalParam
        call myclose(fileDesc)

        write(errorString, *) "Executing final search at ", evalParam, " max times: ",p_maxeval
        call writeToLog(errorString)

        rhobeg  = minval(p_bound(:,2)-p_bound(:,1))/10.0_DP
        rhoend  = 1.0D-3/4.0_dp
        call bobyqa_h(p_nx,p_ninterppt,evalParam,p_bound(:,1),p_bound(:,2),rhobeg,rhoend,p_iprint,p_maxeval,p_wspace,p_nmom)

        fn_val = objFun(evalParam)

        call myopen(UNIT=fileDesc, FILE='searchResults.dat', STATUS='unknown', IOSTAT=openStat, ACTION='write',position='append')
        write(fileDesc,7451) i, i, fn_val, evalParam
        call myclose(fileDesc)

        ! addeed this to print in a file the global argmin and fmin (AA and TK)
!    	IF (fn_val > bestValue) THEN
!       	    fn_val = bestValue
!            evalParam = bestParam
!	    ENDIF
    	IF (fn_val < bestValue) THEN
       	    bestValue = fn_val
            bestParam = evalParam
	    ENDIF

        ! addeed this to print in a file the global argmin and fmin (AA and TK): note that the last entry (line) in searchResults might therefore be different from the final result reported below
        open(unit = 111, file = "monteCarloGOPAmin.dat", position = "append", STATUS='unknown')
        	!write(111,8111) fn_val, evalParam
          write(111,8111) bestValue, bestParam
        close(111)
8111    format(200f40.20) ! (aa)


        return

7451    format(2i10, 200f40.20)

    END SUBROUTINE lastSearch

    FUNCTION getModifiedParam(whichPoint, evalPoint, whichMethod, lotPoint) RESULT (y)
        !This routine returns the next parameter at which to evaluate. It combines
        !the next sobol point with a previously evaluated point specified by whichMethod.
        !Method 1: simply get the best point, and scale to that
        !Method 2: each point has a probability of being selected, based
        !          on rank. Pick one and then scale to that
        INTEGER(I4B), INTENT(IN) :: whichPoint
        REAL(DP), DIMENSION(p_nx), INTENT(IN) :: evalPoint
        INTEGER(I4B), INTENT(IN) :: whichMethod
        INTEGER(I4B), INTENT(OUT) :: lotPoint
        REAL(DP), DIMENSION(p_nx) :: y
        REAL(DP), DIMENSION(p_nx) :: basePoint
        REAL(DP) :: w11

        SELECT CASE (whichMethod)
            CASE (0)
                basePoint = getBestPoint(lotPoint)
            CASE (1)
                basePoint = getLotteryPoint(lotPoint)
            CASE DEFAULT
                write(errorString, *) "<genericSearch:getModifiedParam> : Error, unknown selectionMethod: ", whichMethod
                call writeToLog(errorString)
                stop 0
        END SELECT

        w11=min(max(0.10, (sqrt(real(whichPoint)/real(p_maxpoints)))),0.995)
        y= w11*basePoint+(1.0_DP-w11)*evalPoint
    END FUNCTION getModifiedParam

    FUNCTION getLotteryPoint(lotPointSelect) RESULT (y)
        !This routine returns a parameter where the probability of that parameter is
        !proportional to the inverse of its rank amongst all sorted points.
        INTEGER(I4B), INTENT(OUT) :: lotPointSelect
        REAL(DP), DIMENSION(p_nx) :: y
        REAL(DP), DIMENSION(p_maxpoints,p_nx) :: sortedPoints
        INTEGER(I4B) :: numPoints, sumVal, i
        REAL(DP), DIMENSION(p_nx) :: breakArray
        REAL(DP), DIMENSION(p_maxpoints) :: probOfPoint
        REAL(DP) :: sum2

        breakArray = -1.0_dp
        !sort the points and find out how many we have
        call getSortedResults(sortedPoints)
        DO numPoints = 1,p_maxpoints
            IF(ALL(sortedPoints(numPoints,:) .EQ. breakArray))THEN
                exit
            END IF
        END DO
        numPoints = numPoints - 1

        IF(numPoints .eq. 0) THEN
            y = p_init
            lotPointSelect=0
            return
        END IF
        !but, only use the number of points we are allowed to as specified
        !in the configuration file
        numPoints = min(numPoints, p_lotteryPoints)

        sumVal = (numPoints * (numPoints+1)) / 2

        !calculate the probability of each point
        probOfPoint = 0.0_dp
        DO i=1,numPoints
            probOfPoint(i) = DBLE(sumVal)/DBLE(i)
        END DO
        sum2 = sum(probOfPoint)
        probOfPoint = probOfPoint/sum2

        DO i=2,numPoints
            probOfPoint(i) = probOfPoint(i-1)+probOfPoint(i)
        END DO

        !now, let's get a random number and find the appropriate point
        call random_number(sum2)
        DO i=1,numPoints
            IF (sum2 < probOfPoint(i)) THEN
                exit
            END IF
        END DO
        y=sortedPoints(i,:)
        lotPointSelect=i
    END FUNCTION getLotteryPoint

    FUNCTION getBestPoint(lotPoint) RESULT (y)
        !This routine returns the best set of parameters we have solved so far, as
        !stored in searchSaved.dat
        INTEGER(I4B), INTENT(OUT) :: lotPoint
        REAL(DP), DIMENSION(p_nx) :: y
        REAL(DP), DIMENSION(p_nx+3) :: y2

        y2=getBestLine(lotPoint)
        y=y2(4:)
    END FUNCTION getBestPoint

    FUNCTION getBestLine(lotPoint) RESULT(y)
        !This function returns the best line that we have calculateod
        !so far (in terms of the minimum function value)
        INTEGER(I4B), INTENT(OUT) :: lotPoint
        REAL(DP), DIMENSION(p_nx+3) :: y
        REAL(DP), DIMENSION(p_maxpoints) :: fval_in
        INTEGER(I4B), DIMENSION(p_maxpoints) :: seqN, LocMinNo
        REAL(DP), DIMENSION(p_maxpoints,p_nx) :: pest_saved
        INTEGER(I4B) :: count1, ilo, openStat

        call myopen(UNIT=fileDesc, FILE='searchResults.dat', STATUS='unknown', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            call myclose(fileDesc)
            write(errorString, *) "<genericSearch.getBestLine()> Error: unable to open searchResults.dat file."
            call writeToLog(errorString)
            stop 0
        END IF
        read(fileDesc,7460,end=10) seqN(1), LocMinNo(1), fval_in(1), pest_saved(1,:)

        count1=2
        DO WHILE (count1 .le. p_maxpoints)
            ! read until we run out of lines. Possible that file is really large (by mistake), so
            ! stop if we've reached the maximum number of iterations
            read(fileDesc,7460,end=492) seqN(count1), LocMinNo(count1), fval_in(count1), pest_saved(count1,:)
            count1=count1+1
        END DO
492     call myclose(fileDesc)
        count1 = count1-1
        ilo=minloc(fval_in(1:count1),1)
        y(1) = DBLE(seqN(ilo))
        y(2) = DBLE(LocMinNo(ilo))
        y(3) = fval_in(ilo)
        y(4:)=pest_saved(ilo,:) ! the best parameter that I have found so far
        lotPoint = 1
        return

        !if we are here, file is empty, so return initial guess
10      call myclose(fileDesc)
        y(1) = 0.0_dp
        y(2) = 0.0_dp
        y(3) = 100.0_dp
        y(4:) = p_init
        lotPoint = 0
        return

7460    format(2i10, 200f40.20)

    END FUNCTION getBestLine

    SUBROUTINE getSortedResults(y)
        !This function sorts all the results so far, and returns them in an array.
        !We use this to determine the top points for the lottery selection.
        REAL(DP), DIMENSION(p_maxpoints,p_nx),INTENT(OUT) :: y
        REAL(DP), DIMENSION(p_maxpoints) :: fval_in
        INTEGER(I4B), DIMENSION(p_maxpoints) :: seqN, LocMinNo, fval_init_index
        REAL(DP), DIMENSION(p_maxpoints,p_nx) :: pest_saved
        INTEGER(I4B) :: count1, openStat

        call myopen(UNIT=fileDesc, FILE='searchResults.dat', STATUS='old', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            call myclose(fileDesc)
            write(errorString,*)"<genericSearch.getSortedResults()> Error: unable to open searchResults.dat file."
            call writeToLog(errorString)
            stop 0
        END IF
        read(fileDesc,7452,end=10) seqN(1), LocMinNo(1), fval_in(1), pest_saved(1,:)

        count1=2
        DO WHILE (count1 .le. p_maxpoints)
            ! read until we run out of lines. Possible that file is really large (by mistake), so
            ! stop if we've reached the maximum number of iterations
            read(fileDesc,7452,end=492) seqN(count1), LocMinNo(count1), fval_in(count1), pest_saved(count1,:)
            count1=count1+1
        END DO
492     call myclose(fileDesc)
        count1 = count1-1

        !Sort the function values. fval_init_index(1) gives the smallest
        !element in fval, fval_init_index(2) gives the next smallest, etc.
        CALL indexx(count1,fval_in(1:count1),fval_init_index(1:count1))

        ! Get the points over which we want to minimize. Note that x_starts is
        ! now in order of the evaluated function
        DO i=1,count1
            y(i,:) = pest_saved(i,:)
        END DO
        y(i,:) = -1.0_dp
        return

        !if we are here, file is empty, so return -1s
10      call myclose(fileDesc)
        y(1,:)=-1.0_dp
        return
7452    format(2i10, 200f40.20)

    END SUBROUTINE getSortedResults

    SUBROUTINE setupSobol()
        !This routine returns the best set of parameters we have solved so far, as
        !stored in searchSaved.dat
        INTEGER(I4B) :: i,j
        REAL(DP), DIMENSION(p_nx+1,p_maxpoints) :: x_starts
        REAL(DP), DIMENSION(p_qr_ndraw+1) :: fval
        INTEGER(I4B), DIMENSION(p_qr_ndraw+1) :: fval_init_index, whichPoints

        call myopen(UNIT=fileDesc, FILE='sobolFnVal.dat', STATUS='old', IOSTAT=ioStat, ACTION='read')
        DO i=1,p_qr_ndraw
            read(fileDesc,7460,end=45) whichPoints(i),fval(i), trial(:,i)
        END DO
        call myclose(fileDesc)

        !Sort the function values. fval_init_index(1) gives the smallest
        !element in fval, fval_init_index(2) gives the next smallest, etc.
        CALL indexx(p_qr_ndraw,fval(1:p_qr_ndraw),fval_init_index(1:p_qr_ndraw))

        ! Get the points over which we want to minimize. Note that x_starts is
        ! now in order of the evaluated function
        ! p_maxpoints is actually one greater than expected because we keep an extra
        ! point (the best sobol point) in the output file.
        DO i=1,p_maxpoints-1
            x_starts(1,i) = fval(fval_init_index(i))
            x_starts(2:p_nx+1,i) = trial(:, fval_init_index(i))
        END DO

        call myopen(UNIT=fileDesc, FILE='x_starts.dat', STATUS='replace', IOSTAT=ioStat, ACTION='write')
        DO i=1,p_maxpoints-1
            write(fileDesc,7462) (x_starts(j,i), j=1,p_nx+1)
        END DO
        call myclose(fileDesc)

        call myopen(unit=fileDesc, file='lastParam.dat',status='replace', IOSTAT=ioStat, ACTION='write')
        call myclose(fileDesc)
        return

        ! if here, didn't read enough values, which shouldn't happen
45      write(errorString, *) "<main:setupSobol> Error. Couldn't read ",p_qr_ndraw, " sobol evaluations."
        call writeToLog(errorString)
        stop 0

7460    format(i10, 200f40.20)
7462    format(200f40.20)

    END SUBROUTINE setupSobol

    FUNCTION getSobolPoint() RESULT (y)
        !get the next sobol point that the algorithm needs to solve
        INTEGER(I4B) :: y, openStat

        call myopen(UNIT=fileDesc, FILE='lastSobol.dat', STATUS='unknown', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            write(errorString,*) "<genericSearch.getSobolPoint()> Error: ", seqNo," unable to open lastSobol.dat file."
            call writeToLog(errorString)
            stop 0
        END IF
        read(fileDesc, 269, END=10) y
        call myclose(fileDesc)

        IF (y < p_qr_ndraw) THEN
            !update if we need to read more points
            call deleteFirstLine('lastSobol.dat')
            IF(isEmptyFile('lastSobol.dat'))THEN
                call setInit(seqNo)
            END IF
        ELSE IF (y == p_qr_ndraw) THEN
            call setInit(seqNo)
            call deleteFirstLine('lastSobol.dat')
        ELSE
            !if we are here, we changed the number of sobol points to calculate and now have too many. So delete
            !all lines in the file and indicate that we are done
            y = -1
            call myopen(unit=fileDesc, file='lastSobol.dat',status='replace', IOSTAT=openStat, ACTION='write')
            call myclose(fileDesc)
        END IF

        return

        !if we are here, file is empty, so return a -1
10      y = -1
        call myclose(fileDesc)
269     format(i10)
    END FUNCTION getSobolPoint

    FUNCTION getParamPoint() RESULT (y)
        !get the next point that the algorithm needs to minimize
        INTEGER(I4B) :: y, openStat

        call myopen(UNIT=fileDesc, FILE='lastParam.dat', STATUS='unknown', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            write(errorString,*) "<genericSearch.getParamPoint()> Error: unable to open lastParam.dat file."
            call writeToLog(errorString)
            stop 0
        END IF
        read(fileDesc, 269, END=10) y
        call myclose(fileDesc)

        IF (y < p_maxpoints) THEN
            !update if we need to read more points
            call deleteFirstLine('lastParam.dat')
            IF(isEmptyFile('lastParam.dat'))THEN
                call setInit(seqNo)
            END IF
        ELSE IF (y == p_maxpoints) THEN
            call setInit(seqNo)
            call deleteFirstLine('lastParam.dat')
        ELSE
            !if we are here, we changed the number of sobol points to calculate and now have too many. So delete
            !all lines in the file and indicate that we are done
            y = -1
            call myopen(unit=fileDesc, file='lastParam.dat',status='replace', IOSTAT=openStat, ACTION='write')
            call myclose(fileDesc)
        END IF

        return

        !if we are here, file is empty, so return a -1
10      y = -1
        call myclose(fileDesc)

269     format(i10)
    END FUNCTION getParamPoint

    SUBROUTINE findMissingSobol(allFinished)
        !Called by the main instance. Waits until function values have been derived
        !for all sobol points. If some points are missing (often due to a warm start
        !instance having been killed), solve for it.
        LOGICAL, INTENT(OUT) :: allFinished

        INTEGER(I4B) :: openStat, i
        LOGICAL, DIMENSION(p_qr_ndraw) :: solvedPoints
        CHARACTER(LEN=1000) :: errorString

        allFinished = .TRUE.
        solvedPoints = .FALSE.

        !we aren't actually guaranteed that all points are complete. Just that we have tried everything. So
        !find out which ones are missing, and have all processes go back and solve them
        solvedPoints = .FALSE.
        call myopen(UNIT=fileDesc, FILE='sobolFnVal.dat', STATUS='unknown', IOSTAT=openStat, ACTION='read')
        DO
            read(fileDesc,270, END=10) i
            IF (i > p_qr_ndraw) THEN
                !This shouldn't happen. let's note the error and stop.
                write(errorString, *) seqNo, " found point: ",i,"greater than max: ",p_qr_ndraw
                call writeToLog(errorString)
                call setInit(seqNo)
                IF(setState(-1,seqNo))THEN
                    stop 0
                ELSE
                    stop 0
                END IF
            END IF
            solvedPoints(i) = .TRUE.
        END DO
10      call myclose(fileDesc)

        !Add the sobol points that need to be solved.
        call myopen(UNIT=fileDesc, FILE='lastSobol.dat', STATUS='old', IOSTAT=openStat, ACTION='write')
        DO i=1,p_qr_ndraw
            IF(solvedPoints(i))THEN
                cycle
            END IF
            write(fileDesc,270) i
            allFinished = .FALSE.
        END DO
        call myclose(fileDesc)

270     format(i10)

    END SUBROUTINE findMissingSobol

    SUBROUTINE findMissingSearch(allFinished)
        !Called by the main instance. Waits until function values have been derived
        !for all sobol points. If some points are missing (often due to a warm start
        !instance having been killed), solve for it.
        LOGICAL, INTENT(OUT) :: allFinished

        INTEGER(I4B) :: openStat, lastP, i
        LOGICAL, DIMENSION(p_maxpoints) :: solvedPoints

        allFinished = .TRUE.
        solvedPoints = .FALSE.

        call myopen(UNIT=fileDesc, FILE='searchStart.dat', STATUS='unknown', IOSTAT=openStat, ACTION='read')
        DO
            READ (fileDesc,7460, END=10) openStat, lastP
            IF(lastP /= 0) THEN
                IF (lastP > p_maxpoints) THEN
                    !This shouldn't happen. let's note the error and stop.
                    write(errorString, *) seqNo, " found point: ",lastP,"greater than max: ",p_maxpoints
                    call writeToLog(errorString)
                    call setInit(seqNo)
                    IF(setState(-1,seqNo))THEN
                        stop 0
                    ELSE
                        stop 0
                    END IF
                END IF
                solvedPoints(lastP) = .TRUE.
            END IF
        END DO
10      call myclose(fileDesc)

        !Add the points that need to be solved.
        call myopen(UNIT=fileDesc, FILE='lastParam.dat', STATUS='old', IOSTAT=openStat, ACTION='write')
        DO i=1,p_maxpoints-1
            IF(solvedPoints(i))THEN
                cycle
            END IF
            write(fileDesc,7470) i
            allFinished = .FALSE.
        END DO
        call myclose(fileDesc)

7460    format(2i10, 200f40.20)
7470    format(i10)
    END SUBROUTINE findMissingSearch

    SUBROUTINE deleteFirstLine(fileName)
        !delete the first line of the given file
        CHARACTER(LEN=*), INTENT(IN) :: fileName
        CHARACTER(80) Name
        INTEGER(I4B) :: IO, iter, fileDesc

        call myopen(UNIT=fileDesc, FILE=fileName, STATUS='unknown', IOSTAT=IO, ACTION='readwrite')
        OPEN( 2, STATUS = 'SCRATCH' )
        IO = 0
        iter = 0
        DO WHILE (IO == 0)
            READ( fileDesc, *, IOSTAT = IO ) Name
            iter = iter+1
            IF ((IO == 0).and.(iter > 1)) THEN
                WRITE( 2, * ) Name
            END IF
        END DO
        REWIND( 2 )                     ! back to the beginning of SCRATCH
        CLOSE( UNIT=fileDesc, STATUS='DELETE' )   ! delete original
        OPEN(UNIT=fileDesc, FILE=fileName, STATUS='unknown', IOSTAT=IO, ACTION='write')      ! recreate original
        IO = 0
        DO WHILE (IO == 0)
            READ( 2, *, IOSTAT = IO ) Name
            IF (IO == 0) WRITE( fileDesc, * ) Name
        END DO
        call myclose(fileDesc)          ! keep
        CLOSE( 2 )                      ! delete
    END SUBROUTINE deleteFirstLine

    SUBROUTINE parseCommandLine()
        !As the name suggests, parse the command line

        temp=COMMAND_ARGUMENT_COUNT()

        IF (temp == 0) THEN
            print *,"Invoke with:"
            print *,"               ./genericSearch < -1 | <0|2> configfile | 1>  [(a)moeba | (d)fpmin | (b)obyqa]"
            print *,"where"
            print *,"              -1 = stop all instances"
            print *,"               0 = cold start"
            print *,"               1 = warm start"
            print *,"               2 = update points"
            print *,"      configfile = intial parameters. Mandatory if cold start or point update"
            print *,"           a,d,b = search algorithm. Optional, default is bobyqa"
            stop 0
        END IF

        IF (temp .GE. 1) THEN
            call GET_COMMAND_ARGUMENT(1, arg1, ierr)
            read (arg1,*) temp2

            isWarm = .FALSE.
            updatePoints = .FALSE.
            IF (temp2 == -1) THEN
                call setInit(p_exitState)
                IF(setState(p_exitState, -1) .neqv. .TRUE.) THEN
                    print *,"error ending all states."
                    stop 0
                END IF
                write(errorString, *) "ending all starts."
                call writeToLog(errorString)
                stop 1
            ELSE IF (temp2 == 1) THEN
                isWarm = .TRUE.
            ELSE IF (temp2 == 2) THEN
                updatePoints = .TRUE.
            END IF
        END IF

        IF ( (temp == 1) .AND. (isWarm .EQV. .FALSE.) ) THEN
            print *, "<main> Error. Must specify a config file if updating points or starting cold."
            stop 0
        END IF

        IF (temp .GE. 2) THEN
            IF (isWarm) THEN
                CALL GET_COMMAND_ARGUMENT(2, arg2, ierr)
                SELECT CASE (arg2)
                    CASE ('a')
                        alg = 1
                    CASE ('d')
                        alg = 2
                    CASE ('b')
                        alg = 0
                    CASE default
                        alg = p_default_alg
                END SELECT
            ELSE
                call GET_COMMAND_ARGUMENT(2, config, ierr)
                alg = p_default_alg
            END IF
        END IF

        IF (temp .GE. 3) THEN
            IF (isWarm) THEN
                print *,"<main>: Error. Cannot specify a config file for warm starts."
                STOP 0
            ELSE
                call GET_COMMAND_ARGUMENT(3, arg2, ierr)
                SELECT CASE (arg2)
                    CASE ('a')
                        alg = 1
                    CASE ('d')
                        alg = 2
                    CASE ('b')
                        alg = 0
                    CASE default
                        alg = p_default_alg
                END SELECT
            END IF
        END IF

        IF (temp .GE. 4) THEN
            print *, "<main>: Error. Unknown parameters."
            STOP 0
        END IF
    END SUBROUTINE parseCommandLine

    SUBROUTINE recover(state)
        !it seems like the main driver state has died. So, set this state to
        !be the driver program
        INTEGER(I4B), INTENT(IN) :: state
        CHARACTER(LEN=1000) :: errorString

        call setInit(seqNo)
        IF (.NOT. setState(state,seqNo)) THEN
            write(errorString, *) seqNo, " set as initial but unable to set state. This should never happen."
            call writeToLog(errorString)
            stop 0
        END IF
    END SUBROUTINE recover

    SUBROUTINE statePause(state)
        !wait in this state
        INTEGER(I4B), INTENT(IN) :: state
        CHARACTER(LEN=1000) :: errorString

        IF(.not. isAlive())THEN
            write(errorString, *) seqNo, " found dead state: ", state,". Trying to recover"
            call writeToLog(errorString)
            call recover(state)
        ELSE
            write(errorString, *) seqNo, " in state ", state,", waiting for next state."
            call writeToLog(errorString)
            call SLEEP(100)
        END IF
    END SUBROUTINE statePause

    FUNCTION isAlive() RESULT (y)
        !check if the driver process is alive. This is defined as if the state has
        !been changed sometime in the last 2 hours.
        LOGICAL :: y
        INTEGER(I4B) :: currentTime, stateTime, ticks
        INTEGER(I4B), PARAMETER :: maxWait = 7200 !seconds - ie. if no change in state, there is probably a problem
        call SYSTEM_CLOCK(COUNT=currentTime,COUNT_RATE=ticks)
        stateTime = getStateTime()
        y = ((currentTime-stateTime)/ticks).LE.maxWait
    END FUNCTION isAlive

    FUNCTION isEmptyFile(fileName) RESULT (y)
        !check if the given file is empty.
        CHARACTER(len=*), INTENT(IN) :: fileName
        LOGICAL :: y
        INTEGER(I4B) :: fileDesc,openStat
        call myopen(UNIT=fileDesc, FILE=fileName, STATUS='unknown', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            write(errorString,*) "<genericSearch.isEmptyFile()> Error: ", seqNo," unable to open ", fileName
            call writeToLog(errorString)
            stop 0
        END IF
        read(fileDesc, *, END=100)
        call myclose(fileDesc)
        y=.FALSE.
        return

100     y = .TRUE.
        call myclose(fileDesc)
        return
    END FUNCTION isEmptyFile
END PROGRAM GlobalSearch
