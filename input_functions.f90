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
    !comment
END MODULE input_functions