MODULE input_functions
!---------------------------------------------------------------------------------------\
! Purpose:                                                                              |
!   Module containing all the input functions required in the assignment.               |
!   The pointer function_pointer in the main program will point to one of these         |
!   Record of revisions:                                                                |
!       Date        Programmer      Description                                         |
!      ======      ===========     =============                                        |
!   01/11/2020         SFF               -                                              |
!---------------------------------------------------------------------------------------\
IMPLICIT NONE


CONTAINS
    REAL(KIND=8) FUNCTION SGN(xin, u, t)
        REAL(KIND=8), INTENT(IN) :: xin, u, t
        SGN = 0.5*(DSIGN(1.D0, xin - u*t) + 1.0D0)
        !PRINT*, 'sgn called, x:',xin,'u',u,'t',t,'SGN',SGN
    END FUNCTION SGN
    
    REAL(KIND=8) FUNCTION EXPONENTIAL(xin, u, t)
        REAL(KIND=8), INTENT(IN) :: xin, u, t
        EXPONENTIAL = 0.5D0*DEXP(-(xin - u*t)**2.0D0)
    END FUNCTION EXPONENTIAL

    REAL(KIND=8) FUNCTION LINEAR(xin, u, t)
        REAL(KIND=8), INTENT(IN) :: xin, u, t
        LINEAR = xin
    END FUNCTION LINEAR

    REAL(KIND=8) FUNCTION ERRORF(ref,num,type)
        REAL(KIND=8), INTENT(IN) :: ref,num
        INTEGER :: ierr
        CHARACTER(3) :: type
        IF(type == 'abs') THEN
            ERRORF = ref-num
        ELSE IF(type == 'rel') THEN
            ERRORF = (ref-num)/num
        ELSE
            WRITE(*,*) 'Error in type passed to ERROR function. STOP NOW'
            ERRORF = 0.D0
            CALL MPI_FINALIZE(ierr)
            STOP
        END IF
    END FUNCTION ERRORF


    SUBROUTINE ANALYTICAL(flux, start, end,fn,timeindex,x_arr,vel,time,d1,d2,fut,id)
        INTEGER, INTENT(IN) :: start, end, timeindex, fut, d1, d2,id
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2) :: flux
        REAL(KIND=8), INTENT(IN), DIMENSION(d1:d2) :: x_arr
        PROCEDURE(SGN),POINTER :: fn
        REAL(KIND=8), INTENT(IN) :: time, vel
        !WRITE(*,'(I2,A20,2I2)') id, 'Has called analytical', d1,d2
        !WRITE(*,*) 'CurrentTime:', time
        DO i=start,end
            !WRITE(*,'(A2,I0,A2,F8.2,I5)') 'x(',i,')=',x_arr(i), id
            flux(i) = fn(x_arr(i), vel, time)
            !WRITE(*,'(A4,I0,A2,F8.2,I5)') 'phi(',i,')=',flux(i), id
        END DO

    END SUBROUTINE ANALYTICAL

    SUBROUTINE ERROR(flux,start,end,fn,timeindex,x_arr,vel,time,d1,d2,fut,id,& 
                    error_array,local_analytical)
        INTEGER, INTENT(IN) :: start, end, timeindex, fut, d1, d2,id
        INTEGER :: i
        REAL(KIND=8), INTENT(OUT), DIMENSION(d1:d2) :: local_analytical
        REAL(KIND=8), INTENT(OUT), DIMENSION(d1:d2) :: error_array
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2) :: flux
        REAL(KIND=8), INTENT(IN), DIMENSION(d1:d2) :: x_arr
        PROCEDURE(SGN),POINTER :: fn
        REAL(KIND=8), INTENT(IN) :: time, vel
        
        CALL ANALYTICAL(local_analytical,start,end,fn,timeindex, &
                        x_arr,vel,time,d1,d2,fut,id)
        error_array = local_analytical - flux
        !WRITE(*,'(3(A5,F6.2,A9,I2,F6.2,A3,F6.2))') ('flux',flux(i), & 
        !    'local_ana',i,local_analytical(i),'err',error_array(i), i=start,end)
        
    END SUBROUTINE ERROR

    SUBROUTINE NORMS(error_array,L1,L2,LINF)
        REAL(KIND=8), INTENT(IN), DIMENSION(:) :: error_array
        REAL(KIND=8), INTENT(OUT) :: L1,L2,LINF
        
        L1 = SUM(error_array)
        L2 = DSQRT(SUM(error_array**2)) !! Fortran also offers the intrinsic function NORM2()
        LINF = MAXVAL(DABS(error_array))
    END SUBROUTINE NORMS


    !comment
END MODULE input_functions
