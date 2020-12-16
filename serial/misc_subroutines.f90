MODULE misc_subroutines
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
