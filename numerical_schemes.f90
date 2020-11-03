MODULE numerical_schemes
!---------------------------------------------------------------------------------------\
! Purpose:                                                                              |
!   Module containing all the implemented numerical schemes required in the assignment. |
!   The pointer scheme_pointer in the main program will point to one of these           |
!   Record of revisions:                                                                |
!       Date        Programmer      Description                                         |
!      ======      ===========     =============                                        |
!   01/11/2020         SFF         Upwind scheme                                        |
!---------------------------------------------------------------------------------------\
CONTAINS
    SUBROUTINE upwind(flux, grid, spacing, tstep, vel, starti, endi) !INCOMPLETE SCHEME
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: starti, endi
        INTEGER :: i
        REAL, INTENT(INOUT), DIMENSION(0:,0:) :: flux
        REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL, INTENT(IN) :: spacing, tstep, vel

        do i = starti, endi
            flux(i,1) = flux(i,0) - vel*tstep/spacing * (flux(i,0) - flux(i-1,0))
        end do

    END SUBROUTINE upwind
    SUBROUTINE central(flux, grid, spacing, tstep, vel, starti, endi) !INCOMPLETE SCHEME
        IMPLICIT NONE
        INTEGER, INTENT(IN) :: starti, endi
        INTEGER :: i
        REAL, INTENT(INOUT), DIMENSION(0:,0:) :: flux
        REAL, INTENT(IN), DIMENSION(:) :: grid
        REAL, INTENT(IN) :: spacing, tstep, vel

        do i = starti, endi
            flux(i,1) = flux(i,0) - vel*tstep/(2*spacing) * (flux(i+1,0) - flux(i-1,0))
        end do

    END SUBROUTINE central
END MODULE numerical_schemes