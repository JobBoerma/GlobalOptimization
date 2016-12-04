PROGRAM GOPA_init_point

	IMPLICIT NONE

        CHARACTER(LEN=100) 	:: line
	REAL			:: mynumber
	INTEGER			:: p_nx, nmom, maxeval, qr_ndraw, maxpoints, searchType, lotteryPoints
	REAL, ALLOCATABLE	:: p_bound(:,:)
	LOGICAL			:: parsedLine
	INTEGER			:: i, jj
	

	! Read dimension and bounds from configfile
        open(UNIT=222, FILE="config.txt", STATUS='old', ACTION='read')
        
	!Parse the basic parameters in the first line
        parsedLine = .FALSE.
	DO WHILE (parsedLine .eqv. .FALSE.)
            read(222,'(A100)') line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.
            read(line,*) p_nx, nmom, maxeval, qr_ndraw, maxpoints, searchType, lotteryPoints
	END DO

	allocate(p_bound(p_nx,2))

	!Parse the bounds
       	parsedLine = .FALSE.
        DO WHILE (parsedLine .eqv. .FALSE.)
            read(222,'(A100)') line

            IF (line(1:1)=='!') THEN
                cycle
            END IF

            parsedLine = .TRUE.

            read(line,*) p_bound(1,1), p_bound(1,2)

            DO i=2,p_nx
                read(222,*) p_bound(i,1), p_bound(i,2)
            END DO
        END DO

	

	! Generate 1000 initial points scaled to the bounds
	do jj = 1,1000
	  do i = 1,p_nx
		! generate a random number between 0 and 1
        CALL RANDOM_NUMBER(mynumber)        	
		! scale the number to the bounds
		mynumber = mynumber * (p_bound(i,2) - p_bound(i,1)) + p_bound(i,1)
		! write the scaled number in dat file
        	open(unit = 111, file = "init_points.dat", position = "append", STATUS='unknown')
        		write(111,*) mynumber
        	close(111)
            close(222)
	enddo
      enddo
END PROGRAM


        !INTEGER :: i_seed
        !INTEGER, DIMENSION(:), ALLOCATABLE :: a_seed
        !INTEGER, DIMENSION(1:8) :: dt_seed

        ! overwrites initial points with random points
        ! source: web.ph.surrey.ac.uk/fortweb/glossary/random_seed.html
        !CALL RANDOM_SEED(size=i_seed)
        !ALLOCATE(a_seed(1:i_seed))
        !CALL RANDOM_SEED(get=a_seed)
        !CALL DATE_AND_TIME(values=dt_seed)
        !a_seed(i_seed)=dt_seed(8); a_seed(1)=dt_seed(8)*dt_seed(7)*dt_seed(6)
        !CALL RANDOM_SEED(put=a_seed)
        !DEALLOCATE(a_seed)
        

