MODULE numerical_schemes
!---------------------------------------------------------------------------------------\
! Purpose:                                                                              |
!   Module containing all the implemented numerical schemes required in the assignment. |
!   The pointer scheme_pointer in the main program will point to one of these           |
!   Record of revisions:                                                                |
!       Date        Programmer      Description                                         |
!      ======      ===========     =============                                        |
!   04/11/2020         SFF         Upwind scheme                                        |
!---------------------------------------------------------------------------------------\
IMPLICIT NONE
    CONTAINS
    SUBROUTINE UPWIND(f, spacing, tstep, vel, start, end,id,d1,d2,pres,fut,inc) !INCOMPLETE SCHEME

        INTEGER, INTENT(IN) :: start, end, id, d1, d2,pres,fut,inc
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: f
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        ! TE FALTA PASARLE INCREMENT=1 O INCREMENT=-1 para que haga bien la toma de info en i-1
        !IF (id == 0) THEN
        !WRITE(*,'(A5,I4,A30,2I4)') 'Proc',id,'has called upwind for indices', start,end
        !    WRITE(*,111) (i, f(i,0), i=start,end)
        !    111 FORMAT('Values before:', /,&
        !               I2, F8.2)
        !END IF
        DO i = start, end
            f(i,fut) = f(i,pres) - vel*tstep/spacing * (f(i,pres) - f(i-inc,pres))
        END DO
        !comment
    END SUBROUTINE UPWIND
    
    SUBROUTINE CENTRAL(f, spacing, tstep, vel, start, end,id,d1,d2,pres,fut,inc)

        INTEGER, INTENT(IN) :: start, end, id, d1, d2,pres,fut,inc
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: f
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        DO i = start, end
            f(i,fut) = f(i,pres) - vel*tstep/(2*spacing) * (f(i+inc,pres) - f(i-inc,pres))
        END DO

    END SUBROUTINE CENTRAL

    SUBROUTINE LAX(f, spacing, tstep, vel, start, end, id, d1, d2,pres,fut,inc)

        INTEGER, INTENT(IN) :: start, end, id, d1, d2,pres,fut,inc
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: f
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        DO i = start, end
            f(i,fut) =  0.5D0 * (f(i+inc,pres) + f(i-inc,pres)) - vel*tstep/(2*spacing) * (f(i+inc,pres) - f(i-inc,pres))
        END DO

    END SUBROUTINE LAX
    
    SUBROUTINE LEAPFROG(f, spacing, tstep, vel, start, end, id, d1, d2,pres,fut,inc)

        INTEGER, INTENT(IN) :: start, end, id, d1, d2,pres,fut,inc
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: f
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        DO i = start, end
            f(i,fut) = f(i,pres-1) - vel*tstep/spacing * (f(i+inc,pres) - f(i-inc,pres))
        END DO

    END SUBROUTINE LEAPFROG

    SUBROUTINE LAXWENDROFF(f, spacing, tstep, vel, start, end, id, d1, d2,pres,fut,inc)
        INTEGER, INTENT(IN) :: start, end, id, d1, d2,pres,fut,inc
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: f
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        DO i = start, end
            f(i,fut) = f(i,pres) - vel*tstep/(2*spacing) * (f(i+inc,pres) - f(i-inc,pres)) &
            + 0.5D0 * (vel**2) * (tstep**2) * (f(i+inc,pres) - 2*f(i,pres) + f(i-inc,pres)) / (spacing**2)
        END DO

    END SUBROUTINE LAXWENDROFF

    SUBROUTINE THIRDORDERUPWIND(f, spacing, tstep, vel, start, end, id, d1, d2,pres,fut,inc)
        INTEGER, INTENT(IN) :: start, end, id, d1, d2,pres,fut,inc
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: f
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel

        DO i=start,end
            f(i,fut) = f(i,pres) - vel*tstep*((f(i,pres) - f(i-inc,pres))/(2*spacing) - (f(i+inc,pres) &
                       - 3*f(i,pres) + 3*f(i-inc,pres) - f(i-2*inc,pres))/(6*spacing))
        END DO
    
    END SUBROUTINE THIRDORDERUPWIND
    
    !SUBROUTINE MACCORMACK IS IMPLEMENTED INSIDE THE TIMELOOP (GOTO 510 IN MAIN PROGRAM)


END MODULE numerical_schemes
