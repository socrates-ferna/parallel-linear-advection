PROGRAM parallel_linear_advection
USE MPI
USE numerical_schemes, only: UPWIND, CENTRAL, LAX, LEAPFROG, LAXWENDROFF!, ANALYTICAL
USE input_functions, only: SGN, EXPONENTIAL, LINEAR, ANALYTICAL, ERROR
!USE shared_declarations
!-----------------------------------------------------------------------
! Purpose: 1D code for parallel resolution of the linear advection equation
! using several numerical schemes. The analytical solution is used to verify the numerical results.
!
! Licence: This code is distributed under GNU GPL Licence
! Author: Sócrates Fernández Fernández, s(dot)fernaferna(at)gmail(dot)com
!
! LAST MOD: 25/11/2020
!-----------------------------------------------------------------------
! BLOCK 0: DECLARATIONS AND MPI INIT
!-----------------------------------------------------------
IMPLICIT NONE
INTEGER :: stat(MPI_STATUS_SIZE), nprocs, id, ierr, send_req, recv_req, nbstat !parallel vars
INTEGER :: i, j, npoints, nperproc, spatialStencil, timeStencil, istart, iend, iteration, &
             senderproc, receiverproc, sentbufferstart, sentbufferend, receivedbufferstart, &
             receivedbufferend, past, present, future
INTEGER :: bstart, bend, intstart,intend ! boundary (receive info form other proc) indices and interior indices
REAL(KIND=8) :: u, CFL, xl, xr, dx, dt, current_time !, recvbuf=0.0D0, sentbuf=53.0D0
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: controlTimes, x, error_array
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: phi
LOGICAL :: send, receive, special_init
CHARACTER(15) :: scheme, infunction
PROCEDURE(UPWIND), POINTER :: SCHEME_POINTER
PROCEDURE(SGN), POINTER :: FUNCTION_POINTER

!comment
CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)

!-----------------------------------------------------------
!BLOCK I: INPUTS. WILL BE READ FROM A FILE IN LATER COMMITS
!-----------------------------------------------------------
u = 1.0D0
scheme = 'upwind'
infunction = 'sgn'
CFL = 0.90D0
npoints = 20
xl = -40.0D0
xr = 40.0D0 
!ALLOCATE(x(0:npoints-1))
ALLOCATE(controlTimes(0:2)); controlTimes = [0.D0, 5.0D0, 10.0D0] !user should be able to define start, end, writeStep
current_time = controlTimes(0)
iteration = 0

!-----------------------------------------------------------
!BLOCK I.I: AUXILIARY CALCULATIONS DUE TO INPUT
!-----------------------------------------------------------
!past = 0
!present = 1
!future = 2

SELECT CASE (scheme) !remember to add the function in numerical_schemes.f90
    CASE ('upwind')
        spatialStencil = 1 !we will have to define a timeStencil and a spatialStencil in future versions (LeapFrog)
        timeStencil = 1
        present = 0
        future = 1
        special_init = .FALSE.
        SCHEME_POINTER => UPWIND !MODULE NOT WORKING
        !WRITE(*,*) 'Chosen scheme is upwind, stencilSize=', spatialStencil, timeStencil
    CASE ('central')
        spatialStencil = 1
        timeStencil = 1
        present = 0
        future = 1
        special_init = .FALSE.
        SCHEME_POINTER => CENTRAL
        WRITE(*,*) 'Chosen scheme is central, stencilSize=', spatialStencil, timeStencil
    CASE ('lax')
        spatialStencil = 1
        timeStencil = 1
        present = 0
        future = 1
        special_init = .FALSE.
        SCHEME_POINTER => LAX
        WRITE(*,*) 'Chosen scheme is lax, stencilSize=', spatialStencil, timeStencil
    CASE ('leapfrog')
        spatialStencil = 1
        timeStencil = 2
        past = 0
        present = 1
        future = 2
        special_init = .TRUE.
        SCHEME_POINTER => LEAPFROG
        WRITE(*,*) 'Chosen scheme is leapfrog, stencilSize=', spatialStencil, timeStencil
    CASE ('lax-wendroff')
        spatialStencil = 1
        timeStencil = 1
        present = 0
        future = 1
        special_init = .FALSE.
        SCHEME_POINTER => LAXWENDROFF
        WRITE(*,*) 'Chosen scheme is lax, stencilSize=', spatialStencil, timeStencil
END SELECT

SELECT CASE (infunction)
    CASE ('sgn')
        FUNCTION_POINTER => SGN
    CASE ('exponential')
        FUNCTION_POINTER => EXPONENTIAL
    CASE ('linear')
        FUNCTION_POINTER => LINEAR
END SELECT

!ALLOCATE(x(0-spatialStencil:npoints-1+spatialStencil));x(:)=0


dx = (xr - xl) / REAL(npoints-1) !a segment of npoints has npoints-1 interior intervals
dt = CFL * dx / abs(u) ! CFL=u*dt/dx
nperproc = npoints / nprocs !Number of points assigned to each proc before remainder
!WRITE(*,*) 'nperproc:', nperproc
!!WRITE(*,005) id,'it time dt', iteration, current_time, dt

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

001 FORMAT(3(A7,I5))
002 FORMAT(A2,I2,2(A10,F8.3))

ALLOCATE(x(istart-spatialStencil:iend+spatialStencil));x(:)=0
ALLOCATE(phi(istart-spatialStencil:iend+spatialStencil, 0:timeStencil));phi(:,:)=0.0 !missing poner las BCs en todos los tsteps


DO i= istart-spatialStencil, iend+spatialStencil !procs on domain limits will have "ghost elements"
    x(i) = xl + i*dx
    WRITE(*,'(A2,I2,A2,F8.2,I10)') 'x(',i,')=',x(i), id
END DO
    
CALL ANALYTICAL(phi(istart-spatialStencil:istart+spatialStencil,present), &
                istart-spatialStencil,iend+spatialStencil,FUNCTION_POINTER, &
                present, x, u, current_time, & 
                istart-spatialStencil,iend+spatialStencil,future,id)
CALL ANALYTICAL(phi(istart-spatialStencil:istart+spatialStencil,future), &
                istart-spatialStencil,iend+spatialStencil,FUNCTION_POINTER, &
                present, x, u, current_time, & 
                istart-spatialStencil,iend+spatialStencil,future,id)
    !phi(i,present) = FUNCTION_POINTER(x(i), u, current_time)
    !WRITE(*,002) 'i=', i, 'phi at x=', x(i), 'phi(i)=', phi(i,present)

!do i=istart-spatialStencil, iend+spatialStencil
!    WRITE(*,'(A4,I2,A2,F8.2)') 'phi(',i,')=',phi(i,present)
!    WRITE(*,'(A2,I2,A2,F8.2)') 'x(',i,')=',x(i)
!end do

!CALL A WRITING SUBROUTINE HERE

!-----------------------------------------------------------
!BLOCK III: FLUX DETERMINATION AND BOUNDARY EXCHANGE SETUP
!-----------------------------------------------------------
IF ( u >= 0.0 ) THEN
    IF(id == 0) THEN
        senderproc = MPI_PROC_NULL !THIS WAY there are no IF's inside the time loop to distinguish processors
        receiverproc = id +1
        IF(nprocs == 1) receiverproc=MPI_PROC_NULL
    ELSE IF(id == nprocs - 1) THEN
        senderproc = id - 1
        receiverproc = MPI_PROC_NULL
    ELSE
        senderproc = id - 1 ! rank-1 sends info to this process
        receiverproc = id + 1 ! this process sends info to receiverproc
    END IF
    
    sentbufferstart = iend - spatialStencil + 1
    sentbufferend = iend
    receivedbufferstart = istart - spatialStencil
    receivedbufferend = istart - 1
    !WRITE(*,100) 'u+',' Processor',id,'sends indices',sentbufferstart,':',sentbufferend
    !WRITE(*,100) 'u+',' Processor',id,'recvs indices',receivedbufferstart,':',receivedbufferend
    bstart = istart
    bend = istart + spatialStencil -1
!    intstart = istart + spatialStencil     !!!THESE ARE IF I USE BUFFERED SEND
!    intend = iend
ELSE IF ( u <= 0.0 ) THEN
    IF(id == 0) THEN
        receiverproc = MPI_PROC_NULL
        senderproc = id +1
        IF(nprocs == 1) senderproc = MPI_PROC_NULL
    ELSE IF(id == nprocs - 1) THEN
        receiverproc = id - 1
        senderproc = MPI_PROC_NULL
    ELSE
        senderproc = id + 1 ! rank-1 sends info to this process
        receiverproc = id - 1 ! this process sends info to receiverproc
    END IF
    
    WRITE(*,*) 'Flux goes from right to left'
    sentbufferstart = istart
    sentbufferend = istart + spatialStencil - 1
    receivedbufferstart = iend + spatialStencil
    receivedbufferend = iend + 1
    !WRITE(*,100) 'u-',' Processor',id,'sends indices',sentbufferstart,':',sentbufferend
    !WRITE(*,100) 'u-',' Processor',id,'recvs indices',receivedbufferstart,':',receivedbufferend
    bstart = iend - spatialStencil +1
    bend = iend
!    intstart = istart                    !!!THESE ARE IF I USE BUFFERED SEND
!    intend = iend - spatialStencil
ELSE
    WRITE(*,*) 'Flux is neither positive nor negative. It is either zero or an error has ocurred'
    CALL MPI_FINALIZE(ierr)
    STOP
END IF

100 FORMAT(A3,A10,I3,A15,I4,A1,I4)

intstart = istart + spatialStencil
intend = iend - spatialStencil

!-----------------------------------------------------------
!BLOCK IV:
!-----------------------------------------------------------
010 FORMAT(A9,I2,A25)
005 FORMAT(I3,A20,I5,2F10.2)
020 FORMAT(I4, A25, 2I4, F7.2,A10,I3)

!!IF (id == 1) THEN
!    phi(receivedbufferstart:receivedbufferend,present)= 123.88D0
!END IF


WRITE(*,010) 'Processor',id,'waiting for the barrier'
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

DO 
    !dt=min(dt,controlTimes(2)-current_time)
    WRITE(*,005) id,'INIT it time dt', iteration, current_time, dt

    !IF (send .and. receive) THEN !sospecho que hay riesgo de que si se descompensan los procesadores los mensajes terminen siendo machacados incluso si son recibiods
        !sentbuf = phi(sentbufferstart,present)
        CALL MPI_ISEND(phi(sentbufferstart:sentbufferend,present),spatialStencil, &
        MPI_DOUBLE_PRECISION,receiverproc,0,MPI_COMM_WORLD,send_req,ierr)
        WRITE(*,020) id, 'Sent buffer index is',sentbufferstart,sentbufferend, & 
        phi(sentbufferstart:sentbufferend,present),'iteration',iteration

        !WRITE(*,010) 'Processor', id, 'receives'
        
        CALL MPI_IRECV(phi(receivedbufferstart:receivedbufferend,present),spatialStencil, &
        MPI_DOUBLE_PRECISION,senderproc,0,MPI_COMM_WORLD,recv_req,ierr)
        

        CALL SCHEME_POINTER(phi, dx, dt, u, intstart, intend,id, &
        istart-spatialStencil, iend+spatialStencil,present,future) !compute interior
        
        CALL MPI_WAIT(recv_req,MPI_STATUS_IGNORE,ierr)
        WRITE(*,'(A15,F8.2)') 'recvbuf is', phi(receivedbufferstart:receivedbufferend,present)
        !phi(receivedbufferstart:receivedbufferend,present)=recvbuf
        
        !WRITE(*,*) 'Proc',id,'recv_req:',recv_req,'after receiving'
        WRITE(*,020) id, 'Received buffer index is',receivedbufferstart,receivedbufferend, & 
        phi(receivedbufferstart:receivedbufferend,present),'iteration',iteration
        
        CALL SCHEME_POINTER(phi, dx, dt, u, bstart, bend,id, &
        istart-spatialStencil, iend+spatialStencil,present,future)

        CALL MPI_WAIT(send_req,MPI_STATUS_IGNORE,ierr) !ojalá esto baste para que no se machaquen mensajes si el procesador se adelanta  &
        !y manda de nuevo otro bufer antes de que se haya computado el boundary del receptor

        CALL SCHEME_POINTER(phi, dx, dt, u, sentbufferstart, sentbufferend,id, &
        istart-spatialStencil, iend+spatialStencil,present,future)

    !ELSE IF ( send ) THEN
    !    CALL MPI_ISEND(phi(sentbufferstart:sentbufferend,present),spatialStencil, &
    !    MPI_DOUBLE_PRECISION,receiverproc,0,MPI_COMM_WORLD,send_req,ierr)
    !    WRITE(*,020) id, 'Sent buffer index is',sentbufferstart,sentbufferend, & 
    !    phi(sentbufferstart:sentbufferend,present),'iteration',iteration
        
        !!!!!!IF (id == 0) THEN
    !    CALL SCHEME_POINTER(phi, dx, dt, u, intstart, intend,id, &
    !        istart-spatialStencil, iend+spatialStencil,present,future) !compute interior
    !    CALL MPI_WAIT(send_req,MPI_STATUS_IGNORE,ierr)
    !    WRITE(*,'(A5,I4,A22,I3)') 'Proc',id,'has waited receival',iteration
    !    CALL SCHEME_POINTER(phi, dx, dt, u, sentbufferstart,sentbufferend,id, &
    !        istart-spatialStencil,iend+spatialStencil,present,future)
        !!!!!ELSE IF (id == nprocs - 1) THEN
        !!!!    CALL SCHEME_POINTER(phi, dx, dt, u, istart, iend-1,id, &
        !!!!    istart-spatialStencil, iend+spatialStencil,present,future) !compute interior
        !!!!    CALL MPI_WAIT(send_req,MPI_STATUS_IGNORE,ierr)
        !END IF
    !ELSE IF ( receive ) THEN
    !!!!!!    !WRITE(*,010) 'Processor', id, 'receives'
    !    
    !    CALL MPI_IRECV(phi(receivedbufferstart:receivedbufferend,present),spatialStencil, &
    !    MPI_DOUBLE_PRECISION,senderproc,0,MPI_COMM_WORLD,recv_req,ierr)
    !    
    !    CALL SCHEME_POINTER(phi, dx, dt, u, intstart, intend,id, &
    !    istart-spatialStencil, iend+spatialStencil,present,future) !compute interior without computing last point
    !    
    !    CALL MPI_WAIT(recv_req,MPI_STATUS_IGNORE,ierr)
    !    
    !    WRITE(*,020) id, 'Received buffer index is',receivedbufferstart, receivedbufferend, & 
    !    phi(receivedbufferstart:receivedbufferend,present),'iteration',iteration
    !    
    !    CALL SCHEME_POINTER(phi, dx, dt, u, bstart, bend,id, &
    !    istart-spatialStencil, iend+spatialStencil,present,future)
!
    !ELSE
    !    WRITE(*,'(A4,1X,I2,A50)') 'Proc',id,'reached a logical error sending or receiving info'
    !    CALL MPI_FINALIZE(ierr)
    !    STOP
    !END IF

    phi(:,present)=phi(:,future)
    iteration = iteration + 1;
    current_time = current_time + dt
    
    !WRITE(*,005) id,'it time dt', iteration, current_time, dt
    !CALL MPI_WAITALL(2,[send_req, recv_req],MPI_STATUSES_IGNORE,ierr)
    IF (current_time >= 5.0D0 )THEN
        ALLOCATE(error_array(istart:iend));error_array(:)=-10.0D0
        CALL ERROR(phi(istart:iend,present), &
                istart,iend,FUNCTION_POINTER, &
                present, x, u, current_time, & 
                istart,iend,future,id,error_array)
        EXIT
    END IF
END DO

!400 FORMAT(A2,1X,F10.4)

!WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
!DO i=istart,iend
WRITE(*,121) (i, phi(i,present), i=istart,iend)
121 FORMAT(6(I5,F10.2),/)
!END DO

CALL MPI_FINALIZE(ierr)
END PROGRAM parallel_linear_advection


!no puedo modificar el sent bufer hasta que no esté explícitamente terminado. No vale la de llamar al scheme pointer solo una vez
!el request identifica la operación iniciada por la nonblocking call. guarda el tag, el buffer, communicator, source and dest

!!!!!!!!!PENDING REMINDERS!!!!!!!!
!TAKING CARE THAT BC IS ALWAYS IMPOSED
