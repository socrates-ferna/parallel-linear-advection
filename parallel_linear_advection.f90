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
! GitHub: socrates-ferna
! LAST MOD: 7/12/2020
!-----------------------------------------------------------------------
! BLOCK 0: DECLARATIONS AND MPI INIT
!-----------------------------------------------------------
IMPLICIT NONE
INTEGER :: stat(MPI_STATUS_SIZE), commsize, id, ierr, send_req, recv_req, nbstat, position, reqs(4) !parallel vars
INTEGER :: i, j, npoints, nperproc, sst, tst, istart, iend, iteration, &
             pp, np, imsbs, imsbe, imrbs,ipsbs, ipsbe, iprbs, iprbe, &
             imrbe, past, present, future, sendsize,ncontrolTimes, i1
INTEGER :: bstart, bend, intstart,intend, status
REAL(KIND=8) :: u, CFL, xl, xr, dx, dt, currentTime,stopTime
                
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: controlTimes, x,  x_tot, L1_tot,L2_tot,LINF_overall, L1, L2, LINF
                                            
INTEGER, DIMENSION(:), ALLOCATABLE :: displacements, ssizes
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: phi, error_array, analytical_res, receive_arr, r_err, r_ana, saved_results
LOGICAL :: send, receive, special_init
CHARACTER(15) :: scheme, infunction,CFL_str,time_str,npoints_str, aux_str, i1_str
CHARACTER(300) :: msg,pack_buf,filename, format_str
PROCEDURE(UPWIND), POINTER :: SCHEME_POINTER
PROCEDURE(SGN), POINTER :: FUNCTION_POINTER


CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, commsize, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)


001 FORMAT(3(A7,I5))
002 FORMAT(A2,I2,2(A10,F8.3))

!----------------------------------------------------------------!
!BLOCK I: INPUTS. READ FROM A FILE, PACK & BCAST TO ALL PROCS    !
!----------------------------------------------------------------!
ALLOCATE(displacements(0:commsize-1))
ALLOCATE(ssizes(0:commsize-1))
position=0

IF (id == 0) THEN
    OPEN(UNIT=1,FILE='input.ini',STATUS='old', ACTION='READ', IOSTAT=status,IOMSG=msg)
    READ(1,'(2X,F10.3)') u
    CALL MPI_PACK(u,1,MPI_DOUBLE_PRECISION,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(7X,A15)') scheme
    CALL MPI_PACK(scheme,15,MPI_CHARACTER,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(11X,A15)') infunction
    CALL MPI_PACK(infunction,15,MPI_CHARACTER,pack_buf,300,position,MPI_COMM_WORLD,ierr)

    READ(1,'(4X,F7.4)') CFL
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
    !WRITE(*,'(2X,F5.3)') u
    !WRITE(*,'(7X,A15)') scheme
    !WRITE(*,'(11X,A15)') infunction
    !WRITE(*,'(4X,F10.3)') CFL
    !!WRITE(*,'(8X,I6)') npoints
    !WRITE(*,'(3X,F10.3)') xl
    !WRITE(*,'(3X,F10.3)') xr
    !WRITE(*,'(14X,I3)') ncontrolTimes
    !WRITE(*,'(9X,F10.3)') stopTime

    CALL MPI_BCAST(pack_buf,300,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
ELSE
    CALL MPI_BCAST(pack_buf,300,MPI_PACKED,0,MPI_COMM_WORLD,ierr)
    CALL MPI_UNPACK(pack_buf,300,position,u,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    !WRITE(*,*) 'Proc',id,'receives u ',u
    CALL MPI_UNPACK(pack_buf,300,position,scheme,15,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
    !WRITE(*,*) 'Proc',id,'receives scheme ',scheme
    CALL MPI_UNPACK(pack_buf,300,position,infunction,15,MPI_CHARACTER,MPI_COMM_WORLD,ierr)
    !WRITE(*,*) 'Proc',id,'receives infunction ',infunction
    CALL MPI_UNPACK(pack_buf,300,position,CFL,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    !WRITE(*,*) 'Proc',id,'receives CFL ',CFL
    CALL MPI_UNPACK(pack_buf,300,position,npoints,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    !WRITE(*,*) 'Proc',id,'receives npoints ',npoints
    CALL MPI_UNPACK(pack_buf,300,position,xl,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    !WRITE(*,*) 'Proc',id,'receives xl ',xl
    CALL MPI_UNPACK(pack_buf,300,position,xr,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
    !WRITE(*,*) 'Proc',id,'receives xr ',xr
    CALL MPI_UNPACK(pack_buf,300,position,ncontrolTimes,1,MPI_INTEGER,MPI_COMM_WORLD,ierr)
    !WRITE(*,*) 'Proc',id,'receives ncontrolTimes',ncontrolTimes
    CALL MPI_UNPACK(pack_buf,200,position,stopTime,1,MPI_DOUBLE_PRECISION,MPI_COMM_WORLD,ierr)
END IF




ALLOCATE(controlTimes(0:ncontrolTimes));controlTimes(0)=0.0

DO i=1,ncontrolTimes
    controlTimes(i)=stopTime/REAL(ncontrolTimes)*i
END DO

currentTime = controlTimes(0)
iteration = 0
j = 1

!-----------------------------------------------------------
!BLOCK I.I: AUXILIARY CALCULATIONS DUE TO INPUT
!-----------------------------------------------------------

SELECT CASE (scheme) 
    CASE ('upw')
        sst = 1 
        tst = 1
        present = 0
        future = 1
        special_init = .FALSE.
        SCHEME_POINTER => UPWIND 
        !WRITE(*,*) 'Chosen scheme is upwind, stencilSize=', sst, tst
    CASE ('cnt')
        sst = 1
        tst = 1
        present = 0
        future = 1
        special_init = .FALSE.
        SCHEME_POINTER => CENTRAL
        WRITE(*,*) 'Chosen scheme is central, stencilSize=', sst, tst
    CASE ('lax')
        sst = 1
        tst = 1
        present = 0
        future = 1
        special_init = .FALSE.
        SCHEME_POINTER => LAX
        WRITE(*,*) 'Chosen scheme is lax, stencilSize=', sst, tst
    CASE ('lpf')
        sst = 1
        tst = 2
        past = 0
        present = 1
        future = 2
        special_init = .TRUE.
        SCHEME_POINTER => LEAPFROG
        WRITE(*,*) 'Chosen scheme is leapfrog, stencilSize=', sst, tst
    CASE ('lxw')
        sst = 1
        tst = 1
        present = 0
        future = 1
        special_init = .FALSE.
        SCHEME_POINTER => LAXWENDROFF
        WRITE(*,*) 'Chosen scheme is lax, stencilSize=', sst, tst
    CASE ('mcc')
        sst = 1
        tst = 1
        present = 0
        future = 1
        WRITE(*,*) 'Chosen scheme is MacCormack'
        !!pending
    !!! MACCORMACK AND TVD MISSING
END SELECT


SELECT CASE (infunction)
    CASE ('sgn')
        FUNCTION_POINTER => SGN
    CASE ('exp')
        FUNCTION_POINTER => EXPONENTIAL
    CASE ('lin')
        FUNCTION_POINTER => LINEAR
END SELECT



dx = (xr - xl) / REAL(npoints-1) 
dt = CFL * dx / abs(u) 
nperproc = npoints / commsize


!-----------------------------------------------------------!
!BLOCK II: DOMAIN SPLITTING AND INITIALISATION              !
!-----------------------------------------------------------!

IF (id == 0) THEN
!    IF (u > 0.0) THEN
!        send=.TRUE.
!        receive=.FALSE.
!    ELSE IF (u < 0.0) THEN
!        send=.FALSE.
!        receive=.TRUE.
!    END IF
    istart = 0
    iend = nperproc * (id + 1)

    !WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend

ELSE IF (id == commsize -1 ) THEN
!    IF (u < 0.0) THEN
!        send=.TRUE.
!        receive=.FALSE.
!    ELSE IF (u > 0.0) THEN
!        send=.FALSE.
 !       receive=.TRUE.
 !   END IF
    istart = nperproc * id + 1
    iend = npoints - 1
    !WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
ELSE
!    send=.TRUE.
!    receive=.TRUE.
    istart = nperproc * id + 1
    iend = nperproc * (id + 1)
    !WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
    
END IF






ALLOCATE(x(istart-sst:iend+sst));x(:)=0
ALLOCATE(phi(istart-sst:iend+sst, 0:tst));phi(:,:)=0.0

sendsize = SIZE(phi(istart:iend,present))
!PRINT*, 'MY SENDSIZE IS:',sendsize,'PROC',id
CALL MPI_GATHER(sendsize,1,MPI_INTEGER,ssizes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

IF(id == 0) THEN
    PRINT*,ssizes
    !ALLOCATE(x_nodes(0:npoints))
    ALLOCATE(x_tot(0:npoints-1))
    displacements(0)=0
    DO i=1,commsize-1
        displacements(i) = ssizes(i-1) + displacements(i-1)
    END DO
    PRINT*, displacements
    !DO i=0,npoints
    !    x_nodes(i) = xl + (i-0.5D0)*dx
    !END DO
    DO i=0,npoints-1
        x_tot(i) = xl + i*dx
    END DO
END IF

DO i= istart-sst, iend+sst !procs on domain limits will have "ghost elements"
    x(i) = xl + i*dx
    !WRITE(*,'(A2,I2,A2,F8.2,I10)') 'x(',i,')=',x(i), id
END DO



CALL ANALYTICAL(phi(istart-sst:istart+sst,present), istart-sst,iend+sst,FUNCTION_POINTER, &
                present, x, u, currentTime, istart-sst,iend+sst,future,id)
CALL ANALYTICAL(phi(istart-sst:istart+sst,future), istart-sst,iend+sst,FUNCTION_POINTER, &
                present, x, u, currentTime, istart-sst,iend+sst,future,id)



!-----------------------------------------------------------
!BLOCK III: FLUX DETERMINATION AND BOUNDARY EXCHANGE SETUP
!-----------------------------------------------------------
IF ( u >= 0.0 ) THEN
    IF(id == 0) THEN
        pp = MPI_PROC_NULL !THIS WAY there are no IF's inside the time loop to distinguish processors
        np = id +1
        IF(commsize == 1) np= MPI_PROC_NULL
    ELSE IF(id == commsize - 1) THEN
        pp = id - 1
        np = MPI_PROC_NULL
    ELSE
        pp = id - 1 ! rank-1 sends info to this process
        np = id + 1 ! this process sends info to np
    END IF
    
    imsbs = iend - sst + 1    ! i-1 sent buffer start index
    imsbe = iend              ! i-1 sent buffer end index
    imrbs = istart - sst      ! i-1 received buffer start index
    imrbe = istart - 1        ! i-1 received buffer end index

    ipsbs = istart            ! i+1 sent buffer start index
    ipsbe = istart + sst - 1  ! i+1 sent buffer end index
    iprbs = iend + 1          ! i+1 received buffer start index
    iprbe = iend + sst        ! i+1 received buffer end index

    !WRITE(*,100) 'u+',' Processor',id,'sends indices',imsbs,':',imsbe
    !WRITE(*,100) 'u+',' Processor',id,'recvs indices',imrbs,':',imrbe
    
    bstart = istart
    bend = istart + sst -1

!    intstart = istart + sst     !!!THESE ARE IF I USE BUFFERED SEND
!    intend = iend
ELSE IF ( u <= 0.0 ) THEN
    IF(id == 0) THEN
        np = MPI_PROC_NULL
        pp = id +1
    ELSE IF(id == commsize - 1) THEN
        np = id - 1
        pp = MPI_PROC_NULL
    ELSE
        pp = id + 1 ! rank-1 sends info to this process
        np = id - 1 ! this process sends info to np
    END IF
    
    WRITE(*,*) 'Flux goes from right to left'
    imsbs = istart
    imsbe = istart + sst - 1
    imrbs = iend + sst
    imrbe = iend + 1
    !WRITE(*,100) 'u-',' Processor',id,'sends indices',imsbs,':',imsbe
    !WRITE(*,100) 'u-',' Processor',id,'recvs indices',imrbs,':',imrbe
    bstart = iend - sst +1
    bend = iend
!    intstart = istart                    !!!THESE ARE IF I USE BUFFERED SEND
!    intend = iend - sst
ELSE
    WRITE(*,*) 'Flux is neither positive nor negative. It is either zero or an error has ocurred'
    CALL MPI_FINALIZE(ierr)
    STOP
END IF

100 FORMAT(A3,A10,I3,A15,I4,A1,I4)

intstart = istart + sst
intend = iend - sst

!-----------------------------------------------------------
!BLOCK IV: SIMULATION EXECUTION
!-----------------------------------------------------------
010 FORMAT(A9,I2,A25)
005 FORMAT(I3,A20,I5,2F10.2)
020 FORMAT(I4, A25, 2I4, F7.2,A10,I3)


WRITE(*,010) 'Processor',id,'waiting for the barrier'

ALLOCATE(error_array(istart:iend,ncontrolTimes));error_array(:,:) = 3.0D0 !! ALLOCATIONS HAVE "RANDOM" VALUES TO HELP FINDING BUGS
ALLOCATE(analytical_res(istart:iend,ncontrolTimes));analytical_res(:,:) = 8.0D0
ALLOCATE(L1(ncontrolTimes));L1(:) = 5.0D0
ALLOCATE(L2(ncontrolTimes));L2(:) = 6.0D0
ALLOCATE(LINF(ncontrolTimes));LINF(:) = 70.D0
ALLOCATE(saved_results(istart:iend,ncontrolTimes)); saved_results(:,:) = 211.0D0
call MPI_BARRIER(MPI_COMM_WORLD,ierr)



SELECT CASE (scheme)
    CASE('upw')
        GOTO 500
    CASE('mcc')
        GOTO 510
    CASE DEFAULT
        GOTO 505
END SELECT


500 CONTINUE !!UPWIND TIME LOOP STRUCTURE. COMMUNICATION ONLY TO GET THE "I-1" GRID POINT
!SENDS ARE ONLY DOWNSTREAM
!RECEIVES ARE ONLY UPSTREAM
DO 
    dt=min(dt,controlTimes(j)-currentTime)
    !WRITE(*,005) id,'INIT it time dt', iteration, currentTime, dt

    CALL MPI_ISEND(phi(imsbs:imsbe,present),sst, &
        MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,send_req,ierr)

    !WRITE(*,020) id, 'Sent buffer index is',imsbs,imsbe, & 
        !!phi(imsbs:imsbe,present),'iteration',iteration
        
    CALL MPI_IRECV(phi(imrbs:imrbe,present),sst, &
        MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,recv_req,ierr)
        
    CALL SCHEME_POINTER(phi, dx, dt, u, intstart, intend,id, &
        istart-sst, iend+sst,present,future) !compute interior
        
    CALL MPI_WAIT(recv_req,MPI_STATUS_IGNORE,ierr)

    !WRITE(*,'(A15,F8.2)') 'recvbuf is', phi(imrbs:imrbe,present)
        
    !WRITE(*,020) id, 'Received buffer index is',imrbs,imrbe, & 
        !!phi(imrbs:imrbe,present),'iteration',iteration
        
    CALL SCHEME_POINTER(phi, dx, dt, u, bstart, bend,id, &
        istart-sst, iend+sst,present,future)

    CALL MPI_WAIT(send_req,MPI_STATUS_IGNORE,ierr)

    CALL SCHEME_POINTER(phi, dx, dt, u, imsbs, imsbe,id, &
        istart-sst, iend+sst,present,future)

    phi(:,present)=phi(:,future)
    iteration = iteration + 1;
    currentTime = currentTime + dt
    
    IF (currentTime >= controlTimes(j) )THEN
        
        !WRITE(*,*) 'Control point reached, calculating error and norms...'
        
        ! SEGREGATED ERROR AND ANALYTICAL CALCULATION
        CALL ERROR(phi(istart:iend,present), istart,iend,FUNCTION_POINTER, present, x, u, currentTime, & 
                istart,iend,future,id,error_array(istart:iend,j),analytical_res(istart:iend,j))
        
        CALL NORMS(error_array(istart:iend,j),L1(j),L2(j),LINF(j))
        
        PRINT*,'PROC',id,'L1',L1(j),'L2',L2(j),'LINF',LINF(j)
        WRITE(*,'(A5,I1,2(A6,I2), /,(2I3,3(A1,F5.1)))') 'PROC:',id,'istart',istart, 'iend', iend,&
             (id, i, 'a',analytical_res(i,j),'e',error_array(i,j),'f',phi(i,present), i=istart,iend)
        
        saved_results(istart:iend,j) = phi (istart:iend,present)
        j = j+1
        IF(j >= SIZE(controlTimes)) EXIT
        dt = CFL*dx/u
    END IF

        
END DO

GOTO 530

505 CONTINUE !DEFAULT TIME LOOP STRUCTURE INVOLVES COMM TO GET THE I+1 AND I-1 GRID POINTS 
!SENDS ARE: 
!            -DOWNSTREAM SO THAT THE NEXT PROC KNOWS THE I-1 POINT
!            -UPSTREAM SO THAT THE PREVIOUS PROC KNOWS THE I+1 POINT
!RECEIVES ARE:
!            -UPSTREAM SO THAT I KNOW MY I-1 POINT
!            -DOWNSTREAM SO THAT I KNOW THE I+1 POINT

!ALLOCATE(reqs(4))

DO 
    dt=min(dt,controlTimes(j)-currentTime)
    !WRITE(*,005) id,'INIT it time dt', iteration, currentTime, dt

    CALL MPI_ISEND(phi(imsbs:imsbe,present),sst, MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,reqs(1),ierr) !downstream
    CALL MPI_ISEND(phi(ipsbs:ipsbe,present),sst,MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,reqs(2),ierr)  !upstream

    !WRITE(*,020) id, 'Sent buffer index is',imsbs,imsbe, & 
        !!phi(imsbs:imsbe,present),'iteration',iteration
        
    CALL MPI_IRECV(phi(imrbs:imrbe,present),sst, MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,reqs(3),ierr) !upstream
    CALL MPI_IRECV(phi(iprbs:iprbe,present),sst, MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,reqs(4),ierr) !downstream

    CALL SCHEME_POINTER(phi, dx, dt, u, intstart, intend,id, istart-sst, iend+sst,present,future) !compute interior
        
    CALL MPI_WAITALL(4, reqs, MPI_STATUSES_IGNORE,ierr)

    !WRITE(*,'(A15,F8.2)') 'recvbuf is', phi(imrbs:imrbe,present)
        
    !WRITE(*,020) id, 'Received buffer index is',imrbs,imrbe, & 
        !!phi(imrbs:imrbe,present),'iteration',iteration
        
    CALL SCHEME_POINTER(phi, dx, dt, u, ipsbs, ipsbe,id, istart-sst, iend+sst,present,future)

    CALL SCHEME_POINTER(phi, dx, dt, u, imsbs, imsbe,id, istart-sst, iend+sst,present,future)
    IF(id == 0) phi(istart,future)=phi(istart-1,future)  !!Simple way to ensure that f(-40,t) = original BC
    IF(id == commsize - 1) phi(iend,future)=phi(iend+1,future) !!Same for f(40,t) ghost cells out of comp domain are never modified from the initial conditions

    phi(:,present)=phi(:,future)
    iteration = iteration + 1;
    currentTime = currentTime + dt
    
    IF (currentTime >= controlTimes(j) )THEN
        
        !WRITE(*,*) 'Control point reached, calculating error and norms...'
        
        ! SEGREGATED ERROR AND ANALYTICAL CALCULATION
        CALL ERROR(phi(istart:iend,present), istart,iend,FUNCTION_POINTER, present, x, u, currentTime, & 
                istart,iend,future,id,error_array(istart:iend,j),analytical_res(istart:iend,j))
        
        CALL NORMS(error_array(istart:iend,j),L1(j),L2(j),LINF(j))
        
        PRINT*,'PROC',id,'L1',L1(j),'L2',L2(j),'LINF',LINF(j)
        WRITE(*,'(A5,I1,2(A6,I2), /,3(F5.1))') 'PROC:',id,'istart',istart, 'iend', iend,&
             (analytical_res(i,j),error_array(i,j),phi(i,present), i=istart,iend)
        
        saved_results(istart:iend,j) = phi (istart:iend,present)
        j = j+1
        IF(j >= SIZE(controlTimes)) EXIT
        dt = CFL*dx/u
    END IF

END DO

GOTO 530


510 CONTINUE !TIME LOOP STRUCTURE FOR THE MAC-CORMACK METHOD
!THE SENDS AND RECEIVES ARGUMENTS COULD BE SIMPLIFIED HERE BECAUSE OF THE SPECIFICITY OF THE TIME LOOP
!HOWEVER, THEY ARE LEFT THE SAME AS OTHER METHODS TO AVOID PERFORMANCE IMPROVEMENTS THAT COULD LEAD TO CONFUSION EVALUATING THE RESULTS

DO 
    dt=min(dt,controlTimes(j)-currentTime)
    !WRITE(*,005) id,'INIT it time dt', iteration, currentTime, dt
    CALL MPI_ISEND(phi(ipsbs:ipsbe,present),sst,MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,reqs(1),ierr)  !upstream

    !WRITE(*,020) id, 'Sent buffer index is',imsbs,imsbe, & 
        !!phi(imsbs:imsbe,present),'iteration',iteration
    CALL MPI_IRECV(phi(iprbs:iprbe,present),sst, MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,reqs(2),ierr) !downstream
    CALL MPI_WAITALL(2, reqs(1:2), MPI_STATUSES_IGNORE,ierr)

    DO i=istart, iend !PREDICTOR STEP INTO FUTURE COLUMN
        phi(i,future) = phi(i,present) - u*dt/dx*(phi(i+1,present) - phi(i,present))
    END DO

    IF(id == 0) phi(istart,future)=phi(istart-1,future)  !!Simple way to ensure that f(-40,t) = original BC
    IF(id == commsize - 1) phi(iend,future)=phi(iend+1,future) !!Same for f(40,t) ghost cells out of comp domain are never modified from the initial conditions

    CALL MPI_ISEND(phi(imsbs:imsbe,future),sst,MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,reqs(3),ierr)
    CALL MPI_IRECV(phi(imrbs:imrbe,future),sst,MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,reqs(4),ierr)
    CALL MPI_WAITALL(2,reqs(3:4),MPI_STATUSES_IGNORE,ierr)
    
    DO i=iend,istart,-1 !!this doesn't work with reverse flow, just as the rest of the loop
        phi(i,future) = 0.5D0*(phi(i,present) + phi(i,future)) - u*dt/dx*(phi(i,future) - phi(i-1,future))
    END DO !Thanks to computing the subarray in reverse direction we can avoid using an extra array to store the predictor

    IF(id == 0) phi(istart,future)=phi(istart-1,future)  !!Simple way to ensure that f(-40,t) = original BC
    IF(id == commsize - 1) phi(iend,future)=phi(iend+1,future) !!Same for f(40,t) ghost cells out of comp domain are never modified from the initial conditions
  
    phi(:,present)=phi(:,future)
    iteration = iteration + 1;
    currentTime = currentTime + dt
    
    IF (currentTime >= controlTimes(j) )THEN
        
        !WRITE(*,*) 'Control point reached, calculating error and norms...'
        
        ! SEGREGATED ERROR AND ANALYTICAL CALCULATION
        CALL ERROR(phi(istart:iend,present), istart,iend,FUNCTION_POINTER, present, x, u, currentTime, & 
                istart,iend,future,id,error_array(istart:iend,j),analytical_res(istart:iend,j))
        
        CALL NORMS(error_array(istart:iend,j),L1(j),L2(j),LINF(j))
        
        PRINT*,'PROC',id,'L1',L1(j),'L2',L2(j),'LINF',LINF(j)
        WRITE(*,'(A5,I1,2(A6,I2), /,3(F5.1))') 'PROC:',id,'istart',istart, 'iend', iend,&
             (analytical_res(i,j),error_array(i,j),phi(i,present), i=istart,iend)
        
        saved_results(istart:iend,j) = phi (istart:iend,present)
        j = j+1
        IF(j >= SIZE(controlTimes)) EXIT
        dt = CFL*dx/u
    END IF

END DO

GOTO 530






530 CONTINUE
! GATHERING DATA AND WRITING OUTPUT

IF (id == 0) THEN
    ALLOCATE(receive_arr(0:npoints-1,ncontrolTimes))
    ALLOCATE(r_err(0:npoints-1,ncontrolTimes))
    ALLOCATE(r_ana(0:npoints-1,ncontrolTimes))
    ALLOCATE(L1_tot(ncontrolTimes))
    ALLOCATE(L2_tot(ncontrolTimes))
    ALLOCATE(LINF_overall(ncontrolTimes))
END IF

!!! LOS GATHERV HAY QUE HACERLOS EN UN DO O HACER UN NUEVO TIPO DE DATO QUE SEA LA MATRIZ COMPLETA
DO j = 1, ncontrolTimes
    CALL MPI_GATHERV(analytical_res(istart:iend,j),sendsize,MPI_DOUBLE_PRECISION, & 
        r_ana(:,j), ssizes, displacements,MPI_DOUBLE_PRECISION,& 
        0,MPI_COMM_WORLD,ierr)

    CALL MPI_GATHERV(phi(istart:iend,present),sendsize,MPI_DOUBLE_PRECISION, & 
        receive_arr(:,j), ssizes, displacements,MPI_DOUBLE_PRECISION,& 
        0,MPI_COMM_WORLD,ierr)

    CALL MPI_GATHERV(error_array(istart:iend,j),sendsize,MPI_DOUBLE_PRECISION, & 
        r_err(:,j), ssizes, displacements,MPI_DOUBLE_PRECISION,& 
        0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(L1(j),L1_tot(j),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(L2(j),L2_tot(j),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(LINF(j),LINF_overall(j),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
        
END DO

IF (id == 0) THEN

    WRITE(*,*) 'Writing output files'
    !WRITE(*,*) 'L1=', L1(j),'L2=',L2(j),'LINF=',LINF_overall(j)
    WRITE(*,*) (i, receive_arr(i,1), i=0,npoints-1)
    WRITE(CFL_str,'(F8.1)') CFL
    print*, 'CFL_str:', CFL_str
    WRITE(time_str,'(F8.2)') currentTime
    print*, 'time_str:',time_str

    WRITE(npoints_str,'(I15)') npoints
    print*, 'npoints_str:', npoints_str

    filename = TRIM(infunction) // '_' // TRIM(scheme) // '_' // TRIM(ADJUSTL(CFL_str)) // '_'&
     //  TRIM(ADJUSTL(npoints_str)) // '_' // TRIM(ADJUSTL(time_str)) // '.dat'
    
    OPEN(UNIT=100,FILE=filename,STATUS='NEW',ACTION='WRITE', IOSTAT=status,IOMSG=msg)

    WRITE(aux_str,'(A15)') npoints_str
    aux_str=ADJUSTL(aux_str)
    i1 = index(aux_str,' ') - 1
    WRITE(i1_str,*) i1
    print*,'i1', i1, 'aux_str', aux_str
    aux_str= "A" // TRIM(ADJUSTL(i1_str))
    print*, 'aux_str posconcat', aux_str

    format_str = "(A7,A3,A1,A3,A1,F4.2,A1," // aux_str(1:i1) // ",A1,F5.2,A1)"
    format_str = TRIM(format_str)
    print*, 'format_str', format_str
    print*, 'aux_str', aux_str
    WRITE(100,format_str) 'TITLE="', infunction, '_', scheme, '_',  CFL, '_', &
                          TRIM(ADJUSTL(npoints_str)),'_',currentTime, '"'
    WRITE(100,*) 'VARIABLES = "x-nodes", "Analytical", "Numerical", "Error"'
    WRITE(100,*) 'ZONE'
    WRITE(100,format_str) 'T=   "', infunction, '_', scheme, '_',  CFL, '_', &
                         TRIM(ADJUSTL(npoints_str)), '"'
    PRINT*, 'HEADER WRITTEN WITHOUT ERRORS'

    npoints_str= TRIM(ADJUSTL(npoints_str))
    format_str= "(A2," // aux_str(1:i1) // ",A20)"
    WRITE(100,format_str) 'I=', npoints_str, ', DATAPACKING=POINT'
    
    444 FORMAT(4(ES14.7,1X))

    WRITE(100,444) (x_tot(i),r_ana(i,1), receive_arr(i,1), r_err(i,1), i = 0,npoints-1)

    CLOSE(100)
END IF

!GOTO 550

!550 CONTINUE

WRITE(*,121) (id,i, phi(i,present), i=istart,iend)
121 FORMAT(6(I3,I3,F8.2),/)

CALL MPI_FINALIZE(ierr)
END PROGRAM parallel_linear_advection

!!PENDING TO DEFINE NEW MPI DATATYPE AND DO THE ISEND/IRECV USING IT. IF YOU DON'T THEN THE SPATIALST LARGER THAN 1 WILL NOT BE SENT

!------------------!
!     PENDING      !
!------------------!

!-ordenar los allocate y comentar los significados de las cosas mínimamente
!-pasar los case select a un módulo misc_subroutines.f90
!-Cambiar lso nombres de los esquemas a tres letras  OK
!-simplificar los nombres de las variables? Quizás pon un mensaje en el foro   OK A FALTA DE REPASO
!-implementar el Maccormack OK y el TVD
!-sacar el write del time loop añadiendo las columnas necesarias a los arrays de resultado y haciéndolo todo fuera OK SIN TESTEAR
!-archivo de salida con el informe en tiempo de las normas, mira a ver si lo haces tecplot-readable
!-hacer una opción de escupir un paraview-readable (un csv con los títulos de las variables y todo escupido por columnas)
!-limpiar de variables sobrantes tipo x_nodes (que ya la has quitado)
!-hacer la escritura de los tecplot-readable como una serie de datos en tiempo, pendiente que te leas esa parte