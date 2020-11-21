PROGRAM parallel_linear_advection
USE MPI
USE numerical_schemes, only: upwind, central
USE input_functions, only: sgn, exponential, linear
USE shared_declarations
!-----------------------------------------------------------------------
! Purpose: 1D code for parallel resolution of the linear advection equation
! using several numerical schemes. The analytical solution is used to verify the numerical results.
!
! Licence: This code is distributed under GNU GPL Licence
! Author: Sócrates Fernández Fernández, s(dot)fernaferna(at)gmail(dot)com
!
! LAST MOD: 04/11/2020
!-----------------------------------------------------------------------
! BLOCK 0: DECLARATIONS AND MPI INIT
!-----------------------------------------------------------
IMPLICIT NONE
INTEGER :: stat(MPI_STATUS_SIZE), nprocs, id, ierr, reqs(2), nbstat !parallel vars
INTEGER :: i, j, npoints, nperproc, spatialStencil, timeStencil, istart, iend, iteration
INTEGER :: senderproc, receiverproc, sentbufferstart, sentbufferend, receivedbufferstart, receivedbufferend
INTEGER :: bstart, bend, intstart,intend ! boundary (receive info form other proc) indices and interior indices
REAL(KIND=8) :: u, CFL, xl, xr, dx, dt, current_time
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: controlTimes, x
REAL, DIMENSION(:,:), ALLOCATABLE :: phi
LOGICAL :: send, receive
CHARACTER(15) :: scheme, infunction
PROCEDURE(upwind), POINTER :: scheme_pointer
PROCEDURE(sgn), POINTER :: function_pointer

CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)

!-----------------------------------------------------------
!BLOCK I: INPUTS. WILL BE READ FROM A FILE IN LATER COMMITS
!-----------------------------------------------------------
u = 1.0D0
scheme = 'upwind'
infunction = 'linear'
CFL = 1.0D0
npoints = 10
xl = -40.0D0
xr = 40.0D0 
ALLOCATE(x(0:npoints-1))
ALLOCATE(controlTimes(0:2)); controlTimes = [0.D0, 5.0D0, 10.0D0] !user should be able to define start, end, writeStep
current_time = controlTimes(0)
iteration = 0

!-----------------------------------------------------------
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
    CASE ('linear')
        function_pointer => linear
END SELECT


dx = (xr - xl) / REAL(npoints-1) !a segment of npoints has npoints-1 interior intervals
dt = CFL * dx / abs(u) ! CFL=u*dt/dx
nperproc = npoints / nprocs !Number of points assigned to each proc before remainder

WRITE(*,005) id,'it time dt', iteration, current_time, dt

!-----------------------------------------------------------
!BLOCK II: DOMAIN SPLITTING AND INITIALISATION
!-----------------------------------------------------------
IF (id == 0) THEN
    IF (u > 0.0) THEN
        send=.TRUE.
        receive=.FALSE.
    ELSE IF (u < 0.0) THEN
        send=.FALSE.
        receive=.TRUE.
    END IF
    istart = 0
    iend = nperproc * (id + 1)
    WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
!    ALLOCATE(flujo(starti:endi, 3)); flujo(:,:) = 0.0
!    WRITE(*,*) flujo (2,1)
    !PENDING SOME WRITES OVER HERE REGARDING ARRAY SIZE
ELSE IF (id == nprocs -1 ) THEN
    IF (u < 0.0) THEN
        send=.TRUE.
        receive=.FALSE.
    ELSE IF (u > 0.0) THEN
        send=.FALSE.
        receive=.TRUE.
    END IF
    istart = nperproc * id + 1
    iend = npoints - 1
    WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
ELSE
    send=.TRUE.
    receive=.TRUE.
    istart = nperproc * id + 1
    iend = nperproc * (id + 1)
    WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
    
END IF

001 FORMAT(3(A7,I3))
002 FORMAT(A2,I2,2(A10,F8.3))
ALLOCATE(phi(istart-spatialStencil:iend+spatialStencil, 0:timeStencil));phi(:,:)=0.0

!THIS IS NOT GOOD PROGRAMMING PRACTICE, DO X INIT FIRST AND THEN PHI INIT
DO i= istart-spatialStencil, iend+spatialStencil !procs on domain limits will have "ghost elements"
    x(i) = xl + i*dx
    phi(i,0) = function_pointer(x(i), u, current_time)
    WRITE(*,002) 'i=', i, 'phi at x=', x(i), 'phi(i)=', phi(i,0)
END DO

!CALL A WRITING SUBROUTINE HERE

!-----------------------------------------------------------
!BLOCK III: FLUX DETERMINATION AND BOUNDARY EXCHANGE SETUP
!-----------------------------------------------------------
IF ( u >= 0.0) THEN
    senderproc = id - 1 ! rank-1 sends info to this process
    receiverproc = id + 1 ! this process sends info to receiverproc
    sentbufferstart = iend - spatialStencil + 1
    sentbufferend = iend
    receivedbufferstart = istart - spatialStencil
    receivedbufferend = istart - 1
    WRITE(*,100) 'u+',' Processor',id,'sends indices',sentbufferstart,':',sentbufferend
    WRITE(*,100) 'u+',' Processor',id,'recvs indices',receivedbufferstart,':',receivedbufferend
    bstart = istart
    bend = istart + spatialStencil -1
    intstart = istart + spatialStencil
    intend = iend
ELSE IF ( u <= 0.0) THEN
    WRITE(*,*) 'Flux goes from right to left'
    senderproc = id +1 ! rank-1 sends info to this process
    receiverproc = id - 1 ! this process sends info to receiverproc
    sentbufferstart = istart
    sentbufferend = istart + spatialStencil - 1
    receivedbufferstart = iend + spatialStencil
    receivedbufferend = iend + 1
    WRITE(*,100) 'u-',' Processor',id,'sends indices',sentbufferstart,':',sentbufferend
    WRITE(*,100) 'u-',' Processor',id,'recvs indices',receivedbufferstart,':',receivedbufferend
    bstart = iend - spatialStencil +1
    bend = iend
    intstart = istart
    intend = iend - spatialStencil
ELSE
    WRITE(*,*) 'Flux is neither positive nor negative. It is either zero or an error has ocurred'
    CALL MPI_FINALIZE(ierr)
    STOP
END IF

100 FORMAT(A3,A10,I3,A15,I4,A1,I4)

!intstart = istart + spatialStencil
!intend = iend - spatialStencil

!-----------------------------------------------------------
!BLOCK IV:
!-----------------------------------------------------------
010 FORMAT(A9,I2,A25)
005 FORMAT(I3,A20,I5,2F10.2)
020 FORMAT(I4, A25, 2I4, F6.2)

WRITE(*,010) 'Processor',id,'waiting for the barrier'
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

DO 
    !dt=min(dt,controlTimes(2)-current_time)
    WRITE(*,005) id,'INIT it time dt', iteration, current_time, dt
    IF (send) THEN
        CALL MPI_ISEND(phi(sentbufferstart:sentbufferend,0),spatialStencil, &
        MPI_REAL,receiverproc,0,MPI_COMM_WORLD,reqs(1),ierr)
        WRITE(*,020) id, 'Sent buffer index is',sentbufferstart,sentbufferend,&
        phi(sentbufferstart:sentbufferend,0)
        IF (id == 0) THEN
            CALL scheme_pointer(phi, dx, dt, u, istart+1, iend,id, &
            istart-spatialStencil, iend+spatialStencil) !compute interior
            CALL MPI_WAIT(reqs(1),MPI_STATUS_IGNORE,ierr)
            WRITE(*,10) 'Proc',id,'has waited for the receival ISEND'
            
        ELSE IF (id == nprocs - 1) THEN
            CALL scheme_pointer(phi, dx, dt, u, istart, iend-1,id, &
            istart-spatialStencil, iend+spatialStencil) !compute interior
            CALL MPI_WAIT(reqs(1),MPI_STATUS_IGNORE,ierr)
        END IF
    END IF
    !PROBLEMA no puedo enviar antes?
    IF (receive) THEN
        WRITE(*,010) 'Processor', id, 'receives'
        
        CALL MPI_IRECV(phi(receivedbufferstart:receivedbufferend,0),spatialStencil, &
        MPI_REAL,senderproc,0,MPI_COMM_WORLD,reqs(2),ierr)
        
        CALL scheme_pointer(phi, dx, dt, u, intstart, intend,id, &
        istart-spatialStencil, iend+spatialStencil) !compute interior
        
        CALL MPI_WAIT(reqs(2),MPI_STATUS_IGNORE,ierr)
        
        WRITE(*,020) id, 'Received buffer index is',receivedbufferstart,&
        receivedbufferend, phi(receivedbufferstart:receivedbufferend,0)
        
        CALL scheme_pointer(phi, dx, dt, u, bstart, bend,id, &
        istart-spatialStencil, iend+spatialStencil) !compute interio
    
!    ELSE IF( id == 0) THEN !TIENES QUE VER COMO OBLIGAS A QUE LA BC SE MANTENGA
!        CALL scheme_pointer(phi, dx, dt, u, istart+1, iend) !compute interior without affecting
!    
!    ELSE IF( id == nprocs -1) THEN
!        CALL scheme_pointer(phi, dx, dt, u, istart, iend-1)
    END IF
    phi(:,0)=phi(:,1)
    iteration = iteration + 1;
    current_time = current_time + dt
    WRITE(*,005) id,'it time dt', iteration, current_time, dt
    IF (current_time >= controlTimes(2))EXIT
END DO

!400 FORMAT(A2,1X,F10.4)

WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
DO i=istart,iend
    WRITE(*,*) 'id',id, i, phi(i,1)
END DO
!aquí va el boundary update

!do i=0,nprocs-1
!    if ( id == i) then
!        open(10, file= 'initialised' // CHAR(id) // '.dat', form='formatted', status='replace')
!        do j = istart,iend
!            WRITE(10,*) phi(i,0)
!        end do
!        close(10)
!    end if
!end do

!if ( id == 0 ) then 
!    WRITE(*,*) 'WAIT FOR USER INPUT'
!    READ(*,*)
!END IF
CALL MPI_FINALIZE(ierr)
END PROGRAM parallel_linear_advection