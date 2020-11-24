MODULE input_functions
!---------------------------------------------------------------------------------------\
! Purpose:                                                                              |
!   Module containing all the input functions required in the assignment.               |
!   The pointer function_pointer in the main program will point to one of these         |
!   Record of revisions:                                                                |
!       Date        Programmer      Description                                         |
!      ======      ===========     =============                                        |
!   01/11/2020         SFF         Empty structure                                        |
!---------------------------------------------------------------------------------------\
IMPLICIT NONE

CONTAINS
    REAL(KIND=8) FUNCTION SGN(xin, u, t)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: xin, u, t
        SGN = 0.5*(DSIGN(1.D0, xin - u*t) + 1.0D0)
    END FUNCTION SGN
    !PENDING TO CHECK WHETHER POLYMORPHISM HERE IS POSSIBLE (TO PASS THE WHOLE ARRAY AND RECEIVE IT BACK)
    !I THINK POLYMORPH IS NOT POSSIBLE BECAUSE OF THE DATA TYPE RIGIDITY
    !I THINK I HAVE TO INDICATE THAT THE INPUTS ARE ARRAYS TO DO A WHOLE ARRAY OPERATION
    
    REAL(KIND=8) FUNCTION EXPONENTIAL(xin, u, t)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: xin, u, t
        EXPONENTIAL = 0.5D0*DEXP(-(xin - u*t)**2.0D0)
    END FUNCTION EXPONENTIAL

    REAL(KIND=8) FUNCTION LINEAR(xin, u, t)
    IMPLICIT NONE
    REAL(KIND=8), INTENT(IN) :: xin, u, t
        LINEAR = xin
    END FUNCTION LINEAR

    SUBROUTINE ANALYTICAL(flux, start, end,fn,timeindex,x_arr,vel,time,d1,d2,fut,id)
        INTEGER, INTENT(IN) :: start, end, timeindex, fut, d1, d2,id
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: flux
        REAL(KIND=8), INTENT(IN), DIMENSION(d1:d2) :: x_arr
        PROCEDURE(SGN),POINTER :: fn
        REAL(KIND=8), INTENT(IN) :: time, vel
        WRITE(*,'(I2,A20,2I2)') id, 'Has called analytical', d1,d2
        DO i=start,end
            WRITE(*,'(A2,I2,A2,F8.2,I5)') 'x(',i,')=',x_arr(i), id
            flux(i,timeindex) = fn(x_arr(i), vel, time)
        END DO

    END SUBROUTINE ANALYTICAL

    !SUBROUTINE ERROR(flux,)
    !comment
END MODULE input_functions