PROGRAM parallel_linear_advection
USE MPI
USE numerical_schemes, only: upwind, central
USE input_functions, only: sgn, exponential
!-----------------------------------------------------------------------
! Purpose: 1D code for parallel resolution of the linear advection equation
! using several numerical schemes. The analytical solution is used to verify the numerical results.
!
! Licence: This code is distributed under GNU GPL Licence
! Author: Sócrates Fernández Fernández, s(dot)fernaferna(at)gmail(dot)com
!
! LAST MOD: 03/11/2020
!-----------------------------------------------------------------------
! BLOCK 0: DECLARATIONS AND MPI INIT
!-----------------------------------------------------------
IMPLICIT NONE
INTEGER :: stat(MPI_STATUS_SIZE), nprocs, id, ierr !parallel vars
INTEGER :: i, j, npoints, nperproc, spatialStencil, timeStencil, istart, iend
REAL(KIND=8) :: u, CFL, xl, xr, dx, dt, current_time
REAL, DIMENSION(:), ALLOCATABLE :: controlTimes, x
REAL, DIMENSION(:,:), ALLOCATABLE :: phi
CHARACTER(15) :: scheme, infunction
PROCEDURE(upwind), POINTER :: scheme_pointer
PROCEDURE(sgn), POINTER :: function_pointer
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)

!BLOCK I: INPUTS. WILL BE READ FROM A FILE IN LATER COMMITS
!-----------------------------------------------------------
u = 1.0D0
scheme = 'upwind'
infunction = 'sign'
CFL = 1.0D0
npoints = 15
xl = -40.0D0
xr = 40.0D0 
ALLOCATE(x(0:npoints-1))
ALLOCATE(controlTimes(0:2)); controlTimes = [0.D0, 5.0D0, 10.0D0] !user should be able to define start, end, writeStep
current_time = controlTimes(0)

!BLOCK I.I: AUXILIARY CALCULATIONS DUE TO INPUT
!-----------------------------------------------------------
SELECT CASE (scheme) !remember to add the function in numerical_schemes.f90
    CASE ('upwind')
        spatialStencil = 1 !we will have to define a timeStencil and a spatialStencil in future versions (LeapFrog)
        timeStencil = 1
        scheme_pointer => upwind !MODULE NOT WORKING
        WRITE(*,*) 'Chosen scheme is upwind, stencilSize=', spatialStencil
END SELECT

SELECT CASE (infunction)
    CASE ('sign')
        function_pointer => sgn
    CASE ('exponential')
        function_pointer => exponential
END SELECT


dx = (xr - xl) / REAL(npoints-1) !a segment of npoints has npoints-1 interior intervals
dt = CFL * dx / abs(u) ! CFL=u*dt/dx
nperproc = npoints / nprocs !Number of points assigned to each proc before remainder

!BLOCK II: DOMAIN SPLITTING AND INITIALISATION
!-----------------------------------------------------------
IF (id == 0) THEN 
    istart = 0
    iend = nperproc * (id + 1)
    WRITE(*,*) 'MY RANK IS:', id, 'istart=', istart, 'iend=', iend
    ALLOCATE(phi(istart:iend+spatialStencil, 0:timeStencil))
    !PENDING SOME WRITES OVER HERE REGARDING ARRAY SIZE
ELSE IF (id == nprocs -1 ) THEN 
    istart = nperproc * id + 1
    iend = npoints - 1
    WRITE(*,*) 'MY RANK IS:', id, 'istart=', istart, 'iend=', iend
    ALLOCATE(phi(istart-spatialStencil:iend, 0:timeStencil))
ELSE
    istart = nperproc * id + 1
    iend = nperproc * (id + 1)
    WRITE(*,*) 'MY RANK IS:', id, 'istart=', istart, 'iend=', iend
    ALLOCATE(phi(istart-spatialStencil:iend+spatialStencil, 0:timeStencil))
END IF

DO i= istart, iend
    x(i) = xl + i*dx
    phi(i) = function_pointer(x(i), u, current_time) !function still undefined
    WRITE(*,*) 'i=', i, 'Function value at x=', x(i), 'phi(i)=', phi(i)
END DO

do i=0,nprocs-1
    if ( id == i) then
        open(10, file= 'initialised' // id // '.dat', form='formatted', status='replace')
        do j = istart,iend
            WRITE(10,*) phi(i)
        end do
        close(10)
    end if
end do

if ( id == 0 ) then 
    WRITE(*,*) 'WAIT FOR USER INPUT'
    READ(*,*)
END IF
CALL MPI_FINALIZE(ierr)
END PROGRAM parallel_linear_advection