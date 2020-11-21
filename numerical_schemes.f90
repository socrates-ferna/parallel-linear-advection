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
    SUBROUTINE upwind(flux, spacing, tstep, vel, starti, endi,id,d1,d2) !INCOMPLETE SCHEME
        INTEGER, INTENT(IN) :: starti, endi, id, d1, d2
        INTEGER :: i
        REAL, INTENT(INOUT), DIMENSION(d1:d2,0:1) :: flux
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        ! TE FALTA PASARLE INCREMENT=1 O INCREMENT=-1 para que haga bien la toma de info en i-1
        IF (id == 0) THEN
            WRITE(*,'(A40,2I2)') 'Proc 0 has called upwind for indices', starti,endi
            WRITE(*,111) (i, flux(i,0), i=starti,endi)
            111 FORMAT('Values before:', /,&
                       I2, F8.2)
        END IF
        DO i = starti, endi
            flux(i,1) = flux(i,0) - vel*tstep/spacing * (flux(i,0) - flux(i-1,0))
        END DO

    END SUBROUTINE upwind
    
    SUBROUTINE central(flux, spacing, tstep, vel, starti, endi,id, d1, d2) !INCOMPLETE SCHEME
        INTEGER, INTENT(IN) :: starti, endi, id, d1, d2
        INTEGER :: i
        REAL, INTENT(INOUT), DIMENSION(d1:d2,0:1) :: flux
        !REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL(KIND=8), INTENT(IN) :: spacing, tstep, vel
        DO i = starti, endi
            flux(i,1) = flux(i,0) - vel*tstep/(2*spacing) * (flux(i+1,0) - flux(i-1,0))
        END DO

        flux(starti:endi, 0) = flux(starti:endi, 1)

    END SUBROUTINE central

END MODULE numerical_schemes