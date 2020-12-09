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
    SUBROUTINE UPWIND(flux, spacing, tstep, vel, starti, endi,id,d1,d2,pres,fut) !INCOMPLETE SCHEME

        INTEGER, INTENT(IN) :: starti, endi, id, d1, d2,pres,fut
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: flux
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        ! TE FALTA PASARLE INCREMENT=1 O INCREMENT=-1 para que haga bien la toma de info en i-1
        !IF (id == 0) THEN
        !WRITE(*,'(A5,I4,A30,2I4)') 'Proc',id,'has called upwind for indices', starti,endi
        !    WRITE(*,111) (i, flux(i,0), i=starti,endi)
        !    111 FORMAT('Values before:', /,&
        !               I2, F8.2)
        !END IF
        DO i = starti, endi
            flux(i,fut) = flux(i,pres) - vel*tstep/spacing * (flux(i,pres) - flux(i-1,pres))
        END DO
        !comment
    END SUBROUTINE UPWIND
    
    SUBROUTINE CENTRAL(flux, spacing, tstep, vel, starti, endi,id,d1,d2,pres,fut)

        INTEGER, INTENT(IN) :: starti, endi, id, d1, d2,pres,fut
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: flux
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        DO i = starti, endi
            flux(i,fut) = flux(i,pres) - vel*tstep/(2*spacing) * (flux(i+1,pres) - flux(i-1,pres))
        END DO

    END SUBROUTINE CENTRAL

    SUBROUTINE LAX(flux, spacing, tstep, vel, starti, endi, id, d1, d2,pres,fut)

        INTEGER, INTENT(IN) :: starti, endi, id, d1, d2,pres,fut
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: flux
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        DO i = starti, endi
            flux(i,fut) =  0.5D0 * (flux(i+1,pres) - flux(i-1,pres)) - &
            vel*tstep/(2*spacing) * (flux(i+1,pres) - flux(i-1,pres))
        END DO

    END SUBROUTINE LAX
    
    SUBROUTINE LEAPFROG(flux, spacing, tstep, vel, starti, endi, id, d1, d2,pres,fut)

        INTEGER, INTENT(IN) :: starti, endi, id, d1, d2,pres,fut
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: flux
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        DO i = starti, endi
            flux(i,fut) = flux(i,pres-1) - vel*tstep/spacing * (flux(i+1,pres) - flux(i-1,pres))
        END DO

    END SUBROUTINE LEAPFROG

    SUBROUTINE LAXWENDROFF(flux, spacing, tstep, vel, starti, endi, id, d1, d2,pres,fut)
        INTEGER, INTENT(IN) :: starti, endi, id, d1, d2,pres,fut
        INTEGER :: i
        REAL(KIND=8), INTENT(INOUT), DIMENSION(d1:d2,0:fut) :: flux
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        DO i = starti, endi
            flux(i,fut) = flux(i,pres) - vel*tstep/(2*spacing) * (flux(i+1,pres) - flux(i-1,pres)) &
            + 0.5D0 * (vel**2) * (tstep**2) * (flux(i+1,pres) - 2*flux(i,pres) + flux(i-1,pres)) / (spacing**2)
        END DO

    END SUBROUTINE LAXWENDROFF


    !SUBROUTINE MACCORMACK IS IMPLEMENTED INSIDE THE TIMELOOP (GOTO 510 IN MAIN PROGRAM)


END MODULE numerical_schemes