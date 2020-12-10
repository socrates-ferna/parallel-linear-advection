MODULE misc_subroutines
    USE MPI
    USE numerical_schemes
    USE input_functions
    IMPLICIT NONE
!---------------------------------------------------------------------!
!                                                                     !
!   MISCELLANEOUS SUBROUTINES TO MAKE THE MAIN PROGRAM MORE READABLE  !
!                                                                     !
!---------------------------------------------------------------------!
    CONTAINS
    SUBROUTINE SCHEMESELECTION(schstr,schpt,spstenc,timestenc,pres,fut,pst,logicalinit)
        PROCEDURE(UPWIND), POINTER, INTENT(OUT) :: schpt
        LOGICAL :: logicalinit
        INTEGER :: spstenc, timestenc, pres, fut, pst
	CHARACTER(*) :: schstr
        SELECT CASE (schstr) 
            CASE ('upw')
                schpt => UPWIND 
                WRITE(*,*) 'Chosen scheme is upwind, stencilSize=', spstenc,timestenc
            CASE ('cnt')
                schpt => CENTRAL
                WRITE(*,*) 'Chosen scheme is central, stencilSize=', spstenc,timestenc
            CASE ('lax')
                schpt => LAX
                WRITE(*,*) 'Chosen scheme is lax, stencilSize=',spstenc,timestenc
            CASE ('lpf')
                timestenc = 2
                pst = 0
                pres = 1
                fut = 2
                logicalinit = .TRUE.
                schpt => LEAPFROG
                WRITE(*,*) 'Chosen scheme is leapfrog, stencilSize=', spstenc,timestenc
            CASE ('lxw')
                schpt => LAXWENDROFF
                WRITE(*,*) 'Chosen scheme is lax, stencilSize=', spstenc,timestenc
            CASE ('mcc')
                WRITE(*,*) 'Chosen scheme is MacCormack'
            CASE ('tow')
                spstenc = 2
                schpt => THIRDORDERUPWIND
                WRITE(*,*) 'Chosen scheme is Third Order Upwind', spstenc,timestenc
        END SELECT
    END SUBROUTINE SCHEMESELECTION

    SUBROUTINE FUNCTIONSELECTION(fn,fnpoint)
        PROCEDURE(SGN), POINTER :: fnpoint
        CHARACTER(*) :: fn
        SELECT CASE (fn)
	    CASE ('sgn')
                fnpoint => SGN
            CASE ('exp')
                fnpoint => EXPONENTIAL
            CASE ('lin')
                fnpoint => LINEAR !USED DURING DEVELOPMENT TO SEE DIFFERENT VALUES AT EACH GRID POINT
        END SELECT
    END SUBROUTINE FUNCTIONSELECTION

    SUBROUTINE INITIAL_READ(id,vel,sch,func,courant,dim,leftb,rightb,nwrites,finaltime,strid,buf,pos)
        INTEGER, INTENT(OUT) :: pos, strid, nwrites, dim
        INTEGER :: ierr, id, status
	CHARACTER(300) :: msg
        CHARACTER(*) :: buf
        CHARACTER(*), INTENT(OUT) :: sch, func
        REAL(KIND=8) :: vel, courant, leftb, rightb, finaltime

        IF (id == 0) THEN
            OPEN(UNIT=1,FILE='input.ini',STATUS='old', ACTION='READ', IOSTAT=status,IOMSG=msg)
            READ(1,'(2X,F10.3)') vel
            CALL MPI_PACK(vel,1,MPI_DOUBLE_PRECISION,buf,300,pos,MPI_COMM_WORLD,ierr)
        
            READ(1,'(7X,A15)') sch
            CALL MPI_PACK(sch,15,MPI_CHARACTER,buf,300,pos,MPI_COMM_WORLD,ierr)
        
            READ(1,'(11X,A15)') func
            CALL MPI_PACK(func,15,MPI_CHARACTER,buf,300,pos,MPI_COMM_WORLD,ierr)
        
            READ(1,'(4X,F7.4)') courant
            CALL MPI_PACK(courant,1,MPI_DOUBLE_PRECISION,buf,300,pos,MPI_COMM_WORLD,ierr)
        
            READ(1,'(8X,I6)') dim
            CALL MPI_PACK(dim,1,MPI_INTEGER,buf,300,pos,MPI_COMM_WORLD,ierr)
        
            READ(1,'(3X,F10.3)') leftb
            CALL MPI_PACK(leftb,1,MPI_DOUBLE_PRECISION,buf,300,pos,MPI_COMM_WORLD,ierr)
        
            READ(1,'(3X,F10.3)') rightb
            CALL MPI_PACK(rightb,1,MPI_DOUBLE_PRECISION,buf,300,pos,MPI_COMM_WORLD,ierr)
        
            READ(1,'(14X,I3)') nwrites
            CALL MPI_PACK(nwrites,1,MPI_INTEGER,buf,300,pos,MPI_COMM_WORLD,ierr)
        
            READ(1,'(9X,F10.3)') finaltime
            CALL MPI_PACK(finaltime,1,MPI_DOUBLE_PRECISION,buf,300,pos,MPI_COMM_WORLD,ierr)
            READ(1,'(9X,I3)') strid !!! External simulation numbering system to group saved times properly when writing tecplot files
        
            CLOSE(1)
        
            CALL MPI_BCAST(buf,300,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
        
        ELSE
            CALL MPI_BCAST(buf,300,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
            CALL MPI_UNPACK(buf,300,pos,vel,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives u ',vel
            CALL MPI_UNPACK(buf,300,pos,sch,15,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives scheme ',sch
            CALL MPI_UNPACK(buf,300,pos,func,15,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives infunction ',func
            CALL MPI_UNPACK(buf,300,pos,courant,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives CFL ',courant
            CALL MPI_UNPACK(buf,300,pos,dim,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives npoints ',dim
            CALL MPI_UNPACK(buf,300,pos,leftb,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives leftb ',leftb
            CALL MPI_UNPACK(buf,300,pos,rightb,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives rightb ',rightb
            CALL MPI_UNPACK(buf,300,pos,nwrites,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives ncontrolTimes',nwrites
            CALL MPI_UNPACK(buf,200,pos,finaltime,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        END IF
        
    END SUBROUTINE INITIAL_READ
END MODULE misc_subroutines
