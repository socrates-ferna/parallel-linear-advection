MODULE shared_declarations
IMPLICIT NONE
SAVE
INTEGER :: startiindex, endindex
INTEGER :: i_counter
REAL, DIMENSION(:,:), ALLOCATABLE :: flujo
REAL, DIMENSION(:), ALLOCATABLE :: malla
REAL :: espaciado, timstep, veloc
!comment
END MODULE shared_declarations