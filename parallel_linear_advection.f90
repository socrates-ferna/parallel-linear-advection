PROGRAM parallel_linear_advection
USE MPI
USE numerical_schemes, only: UPWIND, CENTRAL, LAX, LEAPFROG, LAXWENDROFF!, ANALYTICAL
USE input_functions, only: SGN, EXPONENTIAL, LINEAR, ANALYTICAL, ERROR, NORMS
!USE shared_declarations
!-----------------------------------------------------------------------
! Purpose: 1D code for parallel resolution of the linear advection equation
! using several numerical schemes. The analytical solution is used to verify the numerical results.
!
! Licence: This code is distributed under GNU GPL Licence
! Author: Sócrates Fernández Fernández, s(dot)fernaferna(at)gmail(dot)com
!
! LAST MOD: 27/11/2020
!-----------------------------------------------------------------------
! BLOCK 0: DECLARATIONS AND MPI INIT
!-----------------------------------------------------------
IMPLICIT NONE
INTEGER :: stat(MPI_STATUS_SIZE), nprocs, id, ierr, send_req, recv_req, nbstat, position !parallel vars
INTEGER :: i, j, npoints, nperproc, spatialStencil, timeStencil, istart, iend, iteration, &
             senderproc, receiverproc, sentbufferstart, sentbufferend, receivedbufferstart, &
             receivedbufferend, past, present, future, sendsize,ncontrolTimes
INTEGER :: bstart, bend, intstart,intend, status ! boundary (receive info form other proc) indices and interior indices
REAL(KIND=8) :: u, CFL, xl, xr, dx, dt, current_time,stopTime, L1, L2, LINF,&
                L1_tot,L2_tot,LINF_overall !, recvbuf=0.0D0, sentbuf=53.0D0
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: controlTimes, x, error_array, receive_arr, x_tot, x_nodes,&  !we need the nodes array for writing the file into tecplot format
                                            analytical_res, recv_err_arr, recv_ana_arr
INTEGER, DIMENSION(:), ALLOCATABLE :: displacements_arr, sendsizes_arr
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: phi
LOGICAL :: send, receive, special_init
CHARACTER(15) :: scheme, infunction
CHARACTER(300) :: msg,pack_buf
PROCEDURE(UPWIND), POINTER :: SCHEME_POINTER
PROCEDURE(SGN), POINTER :: FUNCTION_POINTER
!TYPE :: input_vector
!    REAL(KIND=8) :: velocity
!    CHARACTER(15) :: scheme_string
!    CHARACTER(15) :: infunction_string
!    REAL(KIND=8) :: Courant_number
!    INTEGER :: mesh_size
!    REAL(KIND=8) :: left_boundary
!    REAL(KIND=8) :: right_boundary
!    INTEGER :: n_control_times
!    REAL(KIND=8) :: last_time
!END TYPE input_vector
!TYPE(input_vector) :: readinput


CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)

!-----------------------------------------------------------
!BLOCK I: INPUTS. WILL BE READ FROM A FILE IN LATER COMMITS
!-----------------------------------------------------------
ALLOCATE(displacements_arr(0:nprocs-1))
ALLOCATE(sendsizes_arr(0:nprocs-1))
position=0

IF (id == 0) THEN
    OPEN(UNIT=1,FILE='input.ini',STATUS='old', ACTION='READ', IOSTAT=status,IOMSG=msg)
    READ(1,'(2X,F10.3)') u
    CALL MPI_PACK(u,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(7X,A15)') scheme
    CALL MPI_PACK(scheme,15,MPI_CHARACTER,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(11X,A15)') infunction
    CALL MPI_PACK(infunction,15,MPI_CHARACTER,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(4X,F10.3)') CFL
    CALL MPI_PACK(CFL,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(8X,I6)') npoints
    CALL MPI_PACK(npoints,1,MPI_INTEGER,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(3X,F10.3)') xl
    CALL MPI_PACK(xl,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(3X,F10.3)') xr
    CALL MPI_PACK(xr,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(14X,I3)') ncontrolTimes
    CALL MPI_PACK(ncontrolTimes,1,MPI_INTEGER,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(9X,F10.3)') stopTime
    CALL MPI_PACK(stopTime,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    CLOSE(1)
    WRITE(*,'(2X,F5.3)') u
    WRITE(*,'(7X,A15)') scheme
    WRITE(*,'(11X,A15)') infunction
    WRITE(*,'(4X,F10.3)') CFL
    WRITE(*,'(8X,I6)') npoints
    WRITE(*,'(3X,F10.3)') xl
    WRITE(*,'(3X,F10.3)') xr
    WRITE(*,'(14X,I3)') ncontrolTimes
    WRITE(*,'(9X,F10.3)') stopTime

    CALL MPI_BCAST(pack_buf,300,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
ELSE
    CALL MPI_BCAST(pack_buf,300,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(pack_buf,300,position,u,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    WRITE(*,*) 'Proc',id,'receives u ',u
    CALL MPI_UNPACK(pack_buf,300,position,scheme,15,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
    WRITE(*,*) 'Proc',id,'receives scheme ',scheme
    CALL MPI_UNPACK(pack_buf,300,position,infunction,15,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
    WRITE(*,*) 'Proc',id,'receives infunction ',infunction
    CALL MPI_UNPACK(pack_buf,300,position,CFL,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    WRITE(*,*) 'Proc',id,'receives CFL ',CFL
    CALL MPI_UNPACK(pack_buf,300,position,npoints,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    WRITE(*,*) 'Proc',id,'receives npoints ',npoints
    CALL MPI_UNPACK(pack_buf,300,position,xl,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    WRITE(*,*) 'Proc',id,'receives xl ',xl
    CALL MPI_UNPACK(pack_buf,300,position,xr,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    WRITE(*,*) 'Proc',id,'receives xr ',xr
    CALL MPI_UNPACK(pack_buf,300,position,ncontrolTimes,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    WRITE(*,*) 'Proc',id,'receives ncontrolTimes',ncontrolTimes
    CALL MPI_UNPACK(pack_buf,200,position,stopTime,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
END IF




ALLOCATE(controlTimes(0:ncontrolTimes));controlTimes(0)=0.0
DO i=1,ncontrolTimes
    controlTimes(i)=stopTime/REAL(ncontrolTimes)*i
END DO

current_time = controlTimes(0)
iteration = 0
j = 1

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



dx = (xr - xl) / REAL(npoints-1) !a segment of npoints has npoints-1 interior intervals
dt = CFL * dx / abs(u) ! CFL=u*dt/dx
nperproc = npoints / nprocs !Number of points assigned to each proc before remainder


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
ALLOCATE(phi(istart-spatialStencil:iend+spatialStencil, 0:timeStencil));phi(:,:)=0.0

sendsize = SIZE(phi(istart:iend,present))
CALL MPI_GATHER(sendsize,1,MPI_INTEGER,sendsizes_arr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

IF(id == 0) THEN
    PRINT*,sendsizes_arr
    ALLOCATE(x_nodes(0:npoints))
    ALLOCATE(x_tot(0:npoints-1))
    displacements_arr(0)=0
    DO i=1,nprocs-1
        displacements_arr(i) = sendsizes_arr(i) + displacements_arr(i-1)
    END DO
    DO i=0,npoints
        x_nodes(i) = xl + (i-0.5D0)*dx
    END DO
    DO i=0,npoints-1
        x_tot(i) = xl + i*dx
    END DO
END IF

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
!BLOCK IV: SIMULATION EXECUTION
!-----------------------------------------------------------
010 FORMAT(A9,I2,A25)
005 FORMAT(I3,A20,I5,2F10.2)
020 FORMAT(I4, A25, 2I4, F7.2,A10,I3)


WRITE(*,010) 'Processor',id,'waiting for the barrier'
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

DO 
    dt=min(dt,controlTimes(j)-current_time)
    WRITE(*,005) id,'INIT it time dt', iteration, current_time, dt

    !IF (send .and. receive) THEN 
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
    
    IF (current_time >= controlTimes(j) )THEN
        j = j+1
        WRITE(*,*) 'Control point reached, writing...'
        ALLOCATE(error_array(istart:iend));error_array(:)=-100.0D0
        ALLOCATE(analytical_res(istart:iend));analytical_res(:)= 23.0D0
        CALL ERROR(phi(istart:iend,present), &
                istart,iend,FUNCTION_POINTER, &
                present, x, u, current_time, & 
                istart,iend,future,id,error_array,analytical_res)
        CALL NORMS(error_array,L1,L2,LINF)

        IF (id == 0) THEN
            ALLOCATE(receive_arr(0:npoints-1))
	    ALLOCATE(recv_err_arr(0:npoints-1))
	    ALLOCATE(recv_ana_arr(0:npoints-1))
        END IF
        CALL MPI_GATHERV(analytical_res(istart:iend),sendsize,MPI_DOUBLE_PRECISION, & 
        recv_ana_arr, sendsizes_arr, displacements_arr,MPI_DOUBLE_PRECISION,& 
        0,MPI_COMM_WORLD,ierr)
        CALL MPI_GATHERV(phi(istart:iend,present),sendsize,MPI_DOUBLE_PRECISION, & 
        receive_arr, sendsizes_arr, displacements_arr,MPI_DOUBLE_PRECISION,& 
        0,MPI_COMM_WORLD,ierr)
        CALL MPI_GATHERV(error_array(istart:iend),sendsize,MPI_DOUBLE_PRECISION, & 
        recv_err_arr, sendsizes_arr, displacements_arr,MPI_DOUBLE_PRECISION,& 
        0,MPI_COMM_WORLD,ierr)

        CALL MPI_REDUCE(L1,L1_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(L2,L2_tot,1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
        CALL MPI_REDUCE(LINF,LINF_overall,1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)

        IF (id == 0) THEN

            WRITE(*,*) 'Writing output files'
            WRITE(*,*) 'L1=', L1,'L2=',L2,'LINF=',LINF
            WRITE(*,*) (i, receive_arr(i), i=0,npoints-1)
            
            OPEN(UNIT=100,FILE='result.dat',STATUS='old',ACTION='WRITE', IOSTAT=status,IOMSG=msg) !pending to vary the name
            
            222 FORMAT(A6,A3,A1,A3,A1,F4.2,A1,I3,A1,F4.2,A1) !!!YOU NEED TO GET THE NUMBER OF DIGITS OF NPOINTS INTEGER
            WRITE(100,222) 'TITLE="', infunction, '_', scheme, '_',  CFL, '_',  npoints,'_',current_time, '"'
            WRITE(100,*) 'VARIABLES = "x-nodes", "x-centres", "Analytical", "Numerical", "Error"'
            WRITE(100,*) 'ZONE'
            WRITE(100,222) 'T=   "', infunction, '_', scheme, '_',  CFL, '_',  npoints, '"'
            
            333 FORMAT('I=',I3,', DATAPACKING=BLOCK, VARLOCATION=([2,3,4]=CELLCENTERED)')
            WRITE(100,333) npoints+1
            
            444 FORMAT(101(E13.6))
            WRITE(100,444) (x_nodes(i), i=0,npoints)
            WRITE(100,444) (x_tot(i), i=0,npoints-1)
            WRITE(100,444) (recv_ana_arr(i), i=0,npoints-1)
            WRITE(100,444) (receive_arr(i), i=0,npoints-1)
            WRITE(100,444) (recv_err_arr(i), i=0,npoints-1)

            CLOSE(100)
        END IF
        IF(j >= SIZE(controlTimes)) EXIT
        dt = CFL*dx/u
    END IF
END DO


!DO i=istart,iend
WRITE(*,121) (i, phi(i,present), i=istart,iend)
121 FORMAT(6(I5,F10.2),/)
!END DO

CALL MPI_FINALIZE(ierr)
END PROGRAM parallel_linear_advection

!!PENDING TO DEFINE NEW MPI DATATYPE AND DO THE ISEND/IRECV USING IT. IF YOU DON'T THEN THE SPATIALST LARGER THAN 1 WILL NOT BE SENT
