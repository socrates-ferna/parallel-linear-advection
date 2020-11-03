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
    REAL FUNCTION sgn(xin, u, t)
    IMPLICIT NONE
    REAL, INTENT(IN) :: xin, u, t
        sgn = 0.5*(SIGN(1., xin - u*t) + 1.)
    END FUNCTION sgn
    !PENDING TO CHECK WHETHER POLYMORPHISM HERE IS POSSIBLE (TO PASS THE WHOLE ARRAY AND RECEIVE IT BACK)
    !I THINK I HAVE TO INDICATE THAT THE INPUTS ARE ARRAYS TO DO A WHOLE ARRAY OPERATION
    REAL FUNCTION exponential(xin, u, t)
    IMPLICIT NONE
    REAL, INTENT(IN) :: xin, u, t
        exponential = 0.5*EXP(-(xin - u*t)**2)
    END FUNCTION exponential

END MODULE input_functions