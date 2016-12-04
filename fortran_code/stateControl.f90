!=================================================================
! 1. STATECONTROL: This module does all the file io stuff that manages 
!                  control of states.
!=================================================================
MODULE stateControl
    USE nrtype
    IMPLICIT NONE

    INTEGER(I4B), parameter :: p_exitState=-1
CONTAINS

    SUBROUTINE myopen(unit, file, status, iostat, action, position)
        !open a file. Locks files if they are opened for any purpose other
        !than reading. Locked files cannot be read/written until the
        !lock is removed.
        INTEGER(I4B), INTENT(OUT) :: unit
        CHARACTER(LEN=*), INTENT(IN) :: file, status, action
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: position
        INTEGER(I4B), INTENT(OUT) :: iostat
        CHARACTER(LEN=50) :: lockFile
        CHARACTER(LEN=200) :: msg
        LOGICAL :: done
        INTEGER(I4B), SAVE :: fileUnit = 20
        INTEGER(I4B) :: waitMax

        lockFile = '.lock'
        waitMax = 10

        done = .FALSE.
        fileUnit = 20+mod(fileUnit+2,10)
        unit = fileUnit

        DO WHILE (done .eqv. .FALSE.)
            open(fileUnit+1,file=trim(trim(file)//lockFile),status='new',err=99)
            !don't lock files if we are only opening for reading
            IF (action .EQ. 'read') THEN
                close(fileUnit+1, STATUS='delete')
            END IF
            IF (present(position)) THEN
                open(UNIT=fileUnit,FILE=file,STATUS=status,IOSTAT=iostat,ACTION=action,POSITION=position)
            ELSE
                open(UNIT=fileUnit,FILE=file,STATUS=status,IOSTAT=iostat,ACTION=action)
            END IF
            done = .TRUE.
            cycle
99          print *, file, " is locked. waiting. "
            waitMax = waitMax - 1
            IF(waitMax < 0) THEN
                msg = trim(file)//'. Going to break lock and hope things work.'
                msg = 'Waited too long for '//trim(msg)
                call writeToLog(msg)
                close(fileUnit+1, STATUS='delete')
            END IF
            call SLEEP(1)
        END DO
    END SUBROUTINE myopen

    SUBROUTINE myclose(unit, stat)
        !closes a file. Removes the lock on any locked file.
        INTEGER(I4B), INTENT(IN) :: unit
        CHARACTER(LEN=*), INTENT(IN), OPTIONAL :: stat

        close(unit+1, STATUS='delete')
        IF (present(stat)) THEN
            close(UNIT=unit, STATUS=stat)
        ELSE
            close(UNIT=unit)
        END IF
    END SUBROUTINE myclose

    SUBROUTINE reset()
        !delete any of the files managing the state
        close(41)
        open(UNIT=41, FILE='initTerm.dat', STATUS='replace');close(41)
        open(UNIT=41, FILE='seq.dat', STATUS='replace');close(41)
        open(UNIT=41, FILE='state.dat', STATUS='replace');close(41)
        open(UNIT=41, FILE='stateTime.dat', STATUS='replace');close(41)
        open(UNIT=41, FILE='logFile.txt', STATUS='replace');close(41)
    END SUBROUTINE reset

    SUBROUTINE setInit(initTerm)
        !define which state is the initial one
        INTEGER(I4B), INTENT(IN) :: initTerm
        INTEGER(I4B) :: fileDesc,iostat

        call myopen(UNIT=fileDesc, FILE='initTerm.dat', STATUS='replace', IOSTAT=iostat, ACTION="write")
        write(fileDesc,*) initTerm
        call myclose(fileDesc)
    END SUBROUTINE setInit

    FUNCTION isInit(seqNo) RESULT (y)
        !Function that returns TRUE if the instance calling it is the "main driver"
        INTEGER(I4B), INTENT(IN) :: seqNo
        LOGICAL :: y

        y = (seqNo == getInit())
    END FUNCTION isInit

    FUNCTION getInit() RESULT (y)
        !Function that returns the instance number of the "main driver"
        INTEGER(I4B) :: y, openStat
        INTEGER(I4B) :: fileDesc
        call myopen(UNIT=fileDesc, FILE='initTerm.dat', STATUS='old', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            print *,"<stateControl.getInit()> Error: unable to open init file."
            stop 0
        END IF
        read(fileDesc,*) y
        call myclose(fileDesc)
    END FUNCTION getInit

    FUNCTION getSequenceNumber() RESULT (y)
        ! This function returns a sequence number. If seqNo file exists,
        ! just increments from previous
        INTEGER(I4B) :: openStat, y, fileDesc

        call myopen(UNIT=fileDesc, FILE='seq.dat', STATUS='unknown', IOSTAT=openstat, ACTION='read')
        read(fileDesc, *, END=10) y
        call myclose(fileDesc)
        y = y+1
        call myopen(UNIT=fileDesc, FILE='seq.dat', STATUS='replace', IOSTAT = openStat, ACTION='write')
        write(fileDesc, *) y
        call myclose(fileDesc)
        return

10      y = 1
        call myclose(fileDesc)
        call myopen(UNIT=fileDesc, FILE='seq.dat', STATUS='replace', IOSTAT = openStat, ACTION='write')
        write(fileDesc, *) y
        call myclose(fileDesc)
    END FUNCTION getSequenceNumber

    FUNCTION setState(s, seqNu) RESULT (y)
        !Sets the current state. This can only be done by the main driver program.
        !if another instance tries, FALSE is returned and the state is not changed.
        INTEGER(I4B), INTENT(IN) :: s, seqNu
        LOGICAL :: y
        INTEGER(I4B) :: fileDesc, iostat
        INTEGER(I4B) :: currentTime

        IF (getInit() == seqNu) THEN
            y = .TRUE.
            call myopen(UNIT=fileDesc, FILE='state.dat', STATUS='replace', IOSTAT=iostat, ACTION='write')
            write(fileDesc,*) s
            call myclose(fileDesc)
            call SYSTEM_CLOCK(currentTime)
            call myopen(UNIT=fileDesc, FILE='stateTime.dat', STATUS='replace', IOSTAT=iostat, ACTION='write')
            write(fileDesc,*) currentTime
            call myclose(fileDesc)
        ELSE
            y = .FALSE.
        END IF
    END FUNCTION setState

    FUNCTION getState() RESULT (y)
        !returns the current state value
        INTEGER(I4B) :: y, openStat, fileDesc

        call myopen(UNIT=fileDesc, FILE='state.dat', STATUS='old', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            print *,"<stateControl.getState()> Error: unable to open state file."
            stop 0
        END IF
        read(fileDesc,*) y
        call myclose(fileDesc)
    END FUNCTION getState

    FUNCTION getStateTime() RESULT (y)
        !returns the current state value
        INTEGER(I4B) :: y, openStat, fileDesc

        call myopen(UNIT=fileDesc, FILE='stateTime.dat', STATUS='old', IOSTAT=openStat, ACTION='read')
        IF (openStat /= 0) THEN
            ! file doesn't exist. This is an error
            print *,"<stateControl.getStateTime()> Error: unable to open stateTime file."
            stop 0
        END IF
        read(fileDesc,*) y
        call myclose(fileDesc)
    END FUNCTION getStateTime

    SUBROUTINE waitState(s)
        !instance waits until the current state is the one requested.
        INTEGER(I4B), INTENT(IN) :: s
        INTEGER(I4B) :: currentstate

        currentstate = getState()
        DO WHILE ( (currentstate < s) .AND. (currentstate .ne. p_exitState) )
            call SLEEP(5)
            currentstate = getState()
        END DO
    END SUBROUTINE waitState

    SUBROUTINE writeToLog(str)
        !writes a string to the log file
        CHARACTER(LEN=*), INTENT(IN) :: str
        INTEGER(I4B) ::  fileDesc, iostat

        call myopen(unit=fileDesc, file='logFile.txt',status='unknown',IOSTAT=iostat, ACTION='write',position="append")
        write(fileDesc,*) trim(str)
        call myclose(fileDesc)
    END SUBROUTINE writeToLog

END MODULE stateControl
