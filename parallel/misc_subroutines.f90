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
                !WRITE(*,*) 'Chosen scheme is upwind, stencilSize=', spstenc,timestenc
            CASE ('cnt')
                schpt => CENTRAL
                !WRITE(*,*) 'Chosen scheme is central, stencilSize=', spstenc,timestenc
            CASE ('lax')
                schpt => LAX
                !WRITE(*,*) 'Chosen scheme is lax, stencilSize=',spstenc,timestenc
            CASE ('lpf')
                timestenc = 2
                pst = 0
                pres = 1
                fut = 2
                logicalinit = .TRUE.
                schpt => LEAPFROG
                !WRITE(*,*) 'Chosen scheme is leapfrog, stencilSize=', spstenc,timestenc
            CASE ('lxw')
                schpt => LAXWENDROFF
                !WRITE(*,*) 'Chosen scheme is lax, stencilSize=', spstenc,timestenc
            CASE ('mcc')
                !WRITE(*,*) 'Chosen scheme is MacCormack'
            CASE ('tow')
                spstenc = 2
                schpt => THIRDORDERUPWIND
                !WRITE(*,*) 'Chosen scheme is Third Order Upwind', spstenc,timestenc
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

    SUBROUTINE INITIAL_READ(id,u,scheme,infunction,CFL,npoints,xl,xr,ncontrolTimes,stopTime,strandid,pack_buf,position)
        INTEGER, INTENT(OUT) :: position, strandid, ncontrolTimes, npoints
        INTEGER :: ierr, id, status
	    CHARACTER(300) :: msg
        CHARACTER(*) :: pack_buf
        CHARACTER(15), INTENT(OUT) :: scheme, infunction
        REAL(KIND=8) :: u, CFL, xl, xr, stopTime
        NAMELIST / input_list / u, scheme, infunction, CFL, npoints, xl, xr, ncontrolTimes, stopTime, strandid
        IF (id == 0) THEN

            OPEN(UNIT=1,FILE='input.nml',STATUS='old', ACTION='READ', IOSTAT=status,IOMSG=msg)
            READ(UNIT=1, NML=input_list)
            CLOSE(1)

            CALL MPI_PACK(u,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(scheme,15,MPI_CHARACTER,pack_buf,300,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(infunction,15,MPI_CHARACTER,pack_buf,300,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(CFL,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(npoints,1,MPI_INTEGER,pack_buf,300,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(xl,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(xr,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(ncontrolTimes,1,MPI_INTEGER,pack_buf,300,position,MPI_COMM_WORLD,ierr)
            CALL MPI_PACK(stopTime,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)
        
            CALL MPI_BCAST(pack_buf,300,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
        
        ELSE
            CALL MPI_BCAST(pack_buf,300,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
            CALL MPI_UNPACK(pack_buf,300,position,u,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives u ',u
            CALL MPI_UNPACK(pack_buf,300,position,scheme,15,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives scheme ',scheme
            CALL MPI_UNPACK(pack_buf,300,position,infunction,15,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives infunction ',infunction
            CALL MPI_UNPACK(pack_buf,300,position,CFL,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives CFL ',CFL
            CALL MPI_UNPACK(pack_buf,300,position,npoints,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives npoints ',npoints
            CALL MPI_UNPACK(pack_buf,300,position,xl,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives xl ',xl
            CALL MPI_UNPACK(pack_buf,300,position,xr,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives xr ',xr
            CALL MPI_UNPACK(pack_buf,300,position,ncontrolTimes,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
            !WRITE(*,*) 'Proc',id,'receives ncontrolTimes',ncontrolTimes
            CALL MPI_UNPACK(pack_buf,200,position,stopTime,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
        END IF
        
        
    END SUBROUTINE INITIAL_READ

    SUBROUTINE FLUXDETERMINATION(u,id,pp,np,imsbs,imsbe,imrbs,imrbe,ipsbs,ipsbe,iprbs,iprbe,bstart,bend, &
                                 sst,istart,iend,commsize,increment)
        INTEGER, INTENT(IN) :: id, istart, iend, commsize, sst
        INTEGER, INTENT(OUT) :: pp, np, imsbs, imsbe, imrbs, imrbe, ipsbs, ipsbe, iprbs, iprbe, &
                                bstart, bend, increment
        INTEGER :: ierr
        REAL(KIND=8) :: u

        IF ( u > 0.0 ) THEN
            IF(id == 0) THEN
                pp = MPI_PROC_NULL !THIS WAY there are no IF's inside the time loop to distinguish processors
                np = id +1
                IF(commsize == 1) np= MPI_PROC_NULL
            ELSE IF(id == commsize - 1) THEN
                pp = id - 1
                np = MPI_PROC_NULL
            ELSE
                pp = id - 1 !previous proc
                np = id + 1 !next proc
            END IF
            
            !IMPORTANT: THE I-1 AND I+1 CONCEPTS IN THE FOLLOWING VARIABLES CORRESPONDS TO &
            !THE UPSTREAM AND DOWNSTREAM CELL WITH RESPECT TO THE RECEIVER PROCESSOR
            !THIS WILL ALLOW US TO KEEP THE TIME LOOP COMM STRUCTURE INDEPENDENT OF THE FLOW DIRECTION
            !THE SCHEMES' SUBROUTINES ALSO ADAPT TO THE FLOW DIRECTION IN THIS SENSE BY TAKING ABSOLUTE VALUE OF THE VELOCITY AND REVERSING THE SEARCHED CELL (turns i-1 into i+1)
            !"i-1" is upstream cell
            !"i+1" is downstream cell
            increment = 1
            imsbs = iend - sst + 1    ! "i-1" sent buffer start index
            imsbe = iend              ! "i-1" sent buffer end index
            imrbs = istart - sst      ! "i-1" received buffer start index
            imrbe = istart - 1        ! "i-1" received buffer end index
        
            ipsbs = istart            ! "i+1" sent buffer start index
            ipsbe = istart + sst - 1  ! "i+1" sent buffer end index
            iprbs = iend + 1          ! "i+1" received buffer start index
            iprbe = iend + sst        ! "i+1" received buffer end index
        
            !WRITE(*,100) 'u+',' Processor',id,'sends indices',imsbs,':',imsbe
            !WRITE(*,100) 'u+',' Processor',id,'recvs indices',imrbs,':',imrbe
            
            bstart = istart
            bend = istart + sst -1
        
        ELSE IF ( u < 0.0 ) THEN
            IF(id == 0) THEN
                np = MPI_PROC_NULL
                pp = id +1
                IF(commsize == 1) pp = MPI_PROC_NULL
            ELSE IF(id == commsize - 1) THEN
                np = id - 1
                pp = MPI_PROC_NULL
            ELSE
                pp = id + 1 ! rank-1 sends info to this process
                np = id - 1 ! this process sends info to np
            END IF
            
            !WRITE(*,*) 'Flux goes from right to left'
            increment = -1
            imsbs = istart            ! "i-1" sent buffer start index
            imsbe = istart + sst - 1  ! "i-1" sent buffer end index
            imrbs = iend + 1          ! "i-1" received buffer start index
            imrbe = iend + sst        ! "i-1" received buffer end index
        
            ipsbs = iend - sst + 1    ! "i+1" sent buffer start index
            ipsbe = iend              ! "i+1" sent buffer end index
            iprbs = istart - sst      ! "i+1" received buffer start index
            iprbe = istart - 1        ! "i+1" received buffer end index
        
            !WRITE(*,100) 'u-',' Processor',id,'sends indices',imsbs,':',imsbe
            !WRITE(*,100) 'u-',' Processor',id,'recvs indices',imrbs,':',imrbe
            bstart = iend - sst +1
            bend = iend
        !    intstart = istart                    !!!THESE ARE IF I USE BUFFERED SEND
        !    intend = iend - sst
        ELSE
            WRITE(*,*) 'Flux is neither positive nor negative. It is either zero or an error has ocurred'
            CALL MPI_FINALIZE(ierr)
            STOP
        END IF
        100 FORMAT(A3,A10,I3,A15,I4,A1,I4)
    END SUBROUTINE FLUXDETERMINATION

    SUBROUTINE DOUBLESIZEARRAY(array, aux_array)
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: array
        REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: aux_array
        !INTEGER :: d1,d2
        ALLOCATE(aux_array(2*SIZE(array))); aux_array = 0.0D0
        aux_array(1:SIZE(array)) = array(:)
        DEALLOCATE(array)
        ALLOCATE(array(SIZE(aux_array)))
        array(:) = aux_array(:)
        DEALLOCATE(aux_array)
    END SUBROUTINE DOUBLESIZEARRAY
END MODULE misc_subroutines
