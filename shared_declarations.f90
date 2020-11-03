MODULE shared_declarations
IMPLICIT NONE
SAVE
INTEGER, POINTER :: starti, endi
INTEGER :: i
REAL, DIMENSION(:,:), POINTER :: flux
REAL, DIMENSION(:), POINTER :: grid
REAL, POINTER :: spacing, tstep, vel

END MODULE shared_declarations