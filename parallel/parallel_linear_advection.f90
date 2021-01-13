PROGRAM parallel_linear_advection
USE MPI
USE numerical_schemes !, only: UPWIND, CENTRAL, LAX, LEAPFROG, LAXWENDROFF!, ANALYTICAL
USE input_functions !, only: SGN, EXPONENTIAL, LINEAR, ANALYTICAL, ERROR, NORMS
USE misc_subroutines
!-----------------------------------------------------------------------
! Purpose: 1D code for parallel resolution of the linear advection equation
! using several numerical schemes. The analytical solution is used to verify the numerical results.
!
! Licence: This code is distributed under GNU GPL Licence
! Author: Sócrates Fernández Fernández, s(dot)fernaferna(at)gmail(dot)com
! GitHub: socrates-ferna/parallel-linear-advection
! LAST MOD: 13/01/2021
!-----------------------------------------------------------------------
! COMMENTS
! -I leave several WRITEs throughout the code to facilitate the user the retrieval of info when investigating how the code works
! -
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! BLOCK 0: DECLARATIONS AND MPI INIT
!-----------------------------------------------------------
IMPLICIT NONE
INTEGER :: stat(MPI_STATUS_SIZE), commsize, id, ierr, send_req, recv_req, nbstat, position, reqs(4)
INTEGER :: i, j, npoints, nperproc, sst, tst, istart, iend, iteration, &
             pp, np, imsbs, imsbe, imrbs,ipsbs, ipsbe, iprbs, iprbe, &
             imrbe, past, present, future, sendsize,ncontrolTimes, i1, increment
INTEGER :: bstart, bend, intstart,intend, status, strandid, niters_guess, correctorstart, correctorend
REAL(KIND=8) :: u, CFL, xl, xr, dx, dt, currentTime,stopTime, absu
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: controlTimes, x,  x_tot, L1_tot,L2_tot,LINF_overall, L1, L2, LINF,&
                                           commtime,ovavgcommtime,ovavgtimeloop
DOUBLE PRECISION :: executiontime,extime1,extime2,timeloop1,timeloop2, &
                    writetime1, writetime2, avgtimeloop, commtime2,commtime1, cumcommtime, avgcomm,&
                    commavg,commstdev,timeloopavg,timeloopstdev
DOUBLE PRECISION :: iptime1,iptime2,imtime1,imtime2, gathertime1,gathertime2, gathertime
INTEGER, DIMENSION(:), ALLOCATABLE :: displacements, ssizes
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: phi, error_array, analytical_res, receive_arr, r_err, r_ana, saved_results
LOGICAL :: send, receive, special_init, fileexist
CHARACTER(15) :: scheme, infunction,CFL_str,time_str,npoints_str, aux_str, i1_str
CHARACTER(300) :: msg,pack_buf,filename, format_str, format2_str
PROCEDURE(UPWIND), POINTER :: SCHEME_POINTER
PROCEDURE(SGN), POINTER :: FUNCTION_POINTER
!NAMELIST / input_list / u, scheme, infunction, CFL, npoints, xl, xr, ncontrolTimes, stopTime, strandid


CALL MPI_INIT(ierr)
CALL MPI_COMM_SIZE(MPI_COMM_WORLD, commsize, ierr)
CALL MPI_COMM_RANK(MPI_COMM_WORLD, id, ierr)

IF(id == 0) extime1 = MPI_WTIME()

002 FORMAT(A2,I2,2(A10,F8.3))

!----------------------------------------------------------------!
!BLOCK 0.I: INPUTS. READ FROM A FILE, PACK & BCAST TO ALL PROCS  !
!----------------------------------------------------------------!
position=0  !! Integer that MPI needs to pack & unpack the buffer

CALL INITIAL_READ(id,u,scheme,infunction,CFL,npoints,xl,xr,ncontrolTimes,stopTime,strandid,pack_buf,position)

!-----------------------------------------------------------!
!BLOCK I: SIMULATION CONFIG WITH INPUT DATA                 !
!-----------------------------------------------------------!
ALLOCATE(controlTimes(0:ncontrolTimes));controlTimes(0)=0.0

DO i=1,ncontrolTimes   !USER DEFINED SAVED INSTANTS 
    controlTimes(i)=stopTime/REAL(ncontrolTimes)*i
END DO

currentTime = controlTimes(0)
iteration = 0
j = 1 ! COUNTER FOR SAVING DATA AT REQUIRED TIMES
dx = (xr - xl) / REAL(npoints-1) !UNIFORM SPACING
dt = CFL * dx / abs(u) 
nperproc = npoints / commsize
niters_guess = INT(stopTime/dt + REAL(ncontrolTimes, KIND=8)) + 1
ALLOCATE(commtime(niters_guess)); commtime(:) = 0.0D0 !! to measure communication cost
!DEFAULT VALUES
sst = 1  !SPATIAL STENCIL. !IT WOULD BE MORE PRECISE TO DEFINE AN UPSTREAM STENCIL & A DOWNSTREAM STENCIL
tst = 1  !TIME STENCIL
present = 0 !n time index
future = 1  !n+1 time index (also used as n+1/2 in maccormack)
special_init = .FALSE.   !LEAPFROG REQUIRES INITIALISATION WITH OTHER SCHEME

CALL SCHEMESELECTION(scheme,SCHEME_POINTER,sst,tst,present,future,past,special_init)
CALL FUNCTIONSELECTION(infunction,FUNCTION_POINTER)



!-----------------------------------------------------------!
!BLOCK II: DOMAIN SPLITTING AND INITIALISATION              !
!-----------------------------------------------------------!
!Division as suggested by lecturer, the imbalance of points across the processors will be
!max(2,MOD(npoints/commsize)-1), which is negligible

001 FORMAT(3(A7,I5))
IF (id == 0) THEN 
    istart = 0
    iend = nperproc * (id + 1)
    !WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
    ALLOCATE(ssizes(0:commsize-1))        ! FOR MPI_GATHERV
ELSE IF (id == commsize -1 ) THEN
    istart = nperproc * id + 1
    iend = npoints - 1
    !WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
ELSE
    istart = nperproc * id + 1
    iend = nperproc * (id + 1)
    !WRITE(*,001) 'RANK:', id, 'istart=', istart, 'iend=', iend
END IF

ALLOCATE(x(istart-sst:iend+sst));x(:)=0                !procs on subdomain limits will have "ghost elements"
ALLOCATE(phi(istart-sst:iend+sst, 0:tst));phi(:,:)=0.0


sendsize = SIZE(phi(istart:iend,present)) !Local number of points,  needed for ssizes array, in turn needed for MPI_GATHERV
!PRINT*, 'MY SENDSIZE IS:',sendsize,'PROC',id

CALL MPI_GATHER(sendsize,1,MPI_INTEGER,ssizes,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

IF(id == 0) THEN
    !PRINT*,ssizes
    ALLOCATE(displacements(0:commsize-1)) ! FOR MPI_GATHERV
    ALLOCATE(x_tot(0:npoints-1))

    displacements(0)=0
    DO i=1,commsize-1
        displacements(i) = ssizes(i-1) + displacements(i-1)
    END DO

    !PRINT*, displacements

    DO i=0,npoints-1
        x_tot(i) = xl + i*dx !Proc 0 needs to initialise the whole domain to receive and write results
    END DO
END IF

DO i= istart-sst, iend+sst 
    x(i) = xl + i*dx
    !WRITE(*,'(A2,I2,A2,F8.2,I10)') 'x(',i,')=',x(i), id
END DO



CALL ANALYTICAL(phi(istart-sst:istart+sst,present), istart-sst,iend+sst,FUNCTION_POINTER, &
                present, x, u, currentTime, istart-sst,iend+sst,future,id)
CALL ANALYTICAL(phi(istart-sst:istart+sst,future), istart-sst,iend+sst,FUNCTION_POINTER, &
                present, x, u, currentTime, istart-sst,iend+sst,future,id)
IF( tst > 1) THEN
    CALL ANALYTICAL(phi(:,past), istart,iend,FUNCTION_POINTER, &
                    present, x, u, currentTime, istart,iend,future,id)
END IF  
!--------------------------------------------------------------------!
!BLOCK III: FLUX DIRECTION DETERMINATION AND BOUNDARY EXCHANGE SETUP !
!--------------------------------------------------------------------!
!We only need to do this once because of the equation we are solving (u=const.)

CALL FLUXDETERMINATION(u,id,pp,np,imsbs,imsbe,imrbs,imrbe,ipsbs,ipsbe,iprbs,iprbe,bstart,bend,sst,istart,iend,commsize,increment)
intstart = istart + sst
intend = iend - sst
absu = ABS(u)

!-----------------------------------------------------------
!BLOCK IV: SIMULATION EXECUTION
!-----------------------------------------------------------

010 FORMAT(A9,I2,A25)
005 FORMAT(I3,A20,I5,2F10.2)
020 FORMAT(I4, A25, 2I4, F7.2,A10,I3)

!! ALLOCATIONS HAVE "RANDOM" VALUES TO HELP FINDING BUGS
ALLOCATE(error_array(istart:iend,ncontrolTimes));error_array(:,:) = 3.0D0 
ALLOCATE(analytical_res(istart:iend,ncontrolTimes));analytical_res(:,:) = 8.0D0  
ALLOCATE(L1(ncontrolTimes));L1(:) = 5.0D0   
ALLOCATE(L2(ncontrolTimes));L2(:) = 6.0D0
ALLOCATE(LINF(ncontrolTimes));LINF(:) = 70.D0
ALLOCATE(saved_results(istart:iend,ncontrolTimes)); saved_results(:,:) = 211.0D0


!WRITE(*,010) 'Processor',id,'waiting for the barrier'
call MPI_BARRIER(MPI_COMM_WORLD,ierr)

timeloop1 = MPI_WTIME()
SELECT CASE (scheme)
    CASE('upw','tow')
        GOTO 500
    CASE('mcc')
        GOTO 510
    CASE('lpf')
        SCHEME_POINTER => UPWIND
        DEALLOCATE(controlTimes)
        ALLOCATE(controlTimes(1)); controlTimes(1) = dt
        !!! FALTA DESASIGNAR TODOS LOS ARRAYS DEPENDIENTES DEL CONTADOR J Y
        !!! HACER LA DETERMINACIÓN DE SI SOY LEAPFROG SALTO  A LA DEFAULT STRUCTURE HACIENDO IF(special_init)
        GOTO 500
        !501 CONTINUE

    CASE DEFAULT
        GOTO 505
END SELECT

529 CONTINUE
timeloop2 = MPI_WTIME()
!PRINT*, id,'Out of the time loop in',iteration,'iterations'
avgtimeloop = (timeloop2 - timeloop1) / iteration
avgcomm = SUM(commtime) / iteration !!all procs sum their contribs in REDUCE
GOTO 530


500 CONTINUE !!UPWIND TIME LOOP STRUCTURE. COMMUNICATION ONLY TO GET THE "I-1" (UPSTREAM) GRID POINT
!SENDS ARE ONLY DOWNSTREAM
!RECEIVES ARE ONLY UPSTREAM

DO
    dt=min(dt,controlTimes(j)-currentTime)
    !WRITE(*,005) id,'INIT it time dt', iteration, currentTime, dt
    commtime1 = MPI_WTIME()
    CALL MPI_ISEND(phi(imsbs:imsbe,present),sst, MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,send_req,ierr)

    !WRITE(*,020) id, 'Sent buffer index is',imsbs,imsbe, & 
        !!phi(imsbs:imsbe,present),'iteration',iteration
        
    CALL MPI_IRECV(phi(imrbs:imrbe,present),sst, MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,recv_req,ierr)
    commtime2 = MPI_WTIME()
    cumcommtime = commtime2 - commtime1
    CALL SCHEME_POINTER(phi, dx, dt, absu, intstart, intend,id, istart-sst, iend+sst,present,future,increment) !compute interior
    commtime1 = MPI_WTIME()
    CALL MPI_WAIT(recv_req,MPI_STATUS_IGNORE,ierr)
    commtime2 = MPI_WTIME()
    cumcommtime = cumcommtime + commtime2 - commtime1
    !WRITE(*,'(A15,F8.2)') 'recvbuf is', phi(imrbs:imrbe,present)
        
    !WRITE(*,020) id, 'Received buffer index is',imrbs,imrbe, & 
        !!phi(imrbs:imrbe,present),'iteration',iteration
        
    CALL SCHEME_POINTER(phi, dx, dt, absu, bstart, bend,id, istart-sst, iend+sst,present,future,increment)

    commtime1= MPI_WTIME()
    CALL MPI_WAIT(send_req,MPI_STATUS_IGNORE,ierr)
    commtime2 = MPI_WTIME()
    cumcommtime = cumcommtime + commtime2 - commtime1
    CALL SCHEME_POINTER(phi, dx, dt, absu, imsbs, imsbe,id, istart-sst, iend+sst,present,future,increment)

    phi(:,present)=phi(:,future)
    iteration = iteration + 1;
    currentTime = currentTime + dt
    
    IF (currentTime >= controlTimes(j) )THEN
        saved_results(istart:iend,j) = phi (istart:iend,present)
        j = j+1
        IF(j >= SIZE(controlTimes)) EXIT
        dt = CFL*dx/abs(u)
    END IF

    !IF( iteration == SIZE(commtime) ) THEN !! In case number of iterations exceeds initial guess
    !    CALL DOUBLESIZEARRAY(commtime,aux_array)
    !END IF

    commtime(iteration) = cumcommtime
END DO



GOTO 529

505 CONTINUE !DEFAULT TIME LOOP STRUCTURE INVOLVES COMM TO GET THE UPSTREAM AND DOWNSTREAM STENCIL GRID POINTS 
!SENDS ARE: 
!            -DOWNSTREAM SO THAT THE NEXT PROC KNOWS THE I-1 POINT
!            -UPSTREAM SO THAT THE PREVIOUS PROC KNOWS THE I+1 POINT
!RECEIVES ARE:
!            -UPSTREAM SO THAT I KNOW MY I-1 POINT
!            -DOWNSTREAM SO THAT I KNOW MY I+1 POINT

!ALLOCATE(reqs(4))

DO 
    dt=min(dt,controlTimes(j)-currentTime)
    !WRITE(*,005) id,'INIT it time dt', iteration, currentTime, dt
    commtime1 = MPI_WTIME()
    CALL MPI_ISEND(phi(imsbs:imsbe,present),sst, MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,reqs(1),ierr) !downstream
    CALL MPI_ISEND(phi(ipsbs:ipsbe,present),sst,MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,reqs(2),ierr)  !upstream

    !WRITE(*,020) id, 'Sent buffer index is',imsbs,imsbe, & 
        !!phi(imsbs:imsbe,present),'iteration',iteration
        
    CALL MPI_IRECV(phi(imrbs:imrbe,present),sst, MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,reqs(3),ierr) !upstream
    CALL MPI_IRECV(phi(iprbs:iprbe,present),sst, MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,reqs(4),ierr) !downstream
    commtime2 = MPI_WTIME()
    cumcommtime = commtime2 - commtime1
    CALL SCHEME_POINTER(phi, dx, dt, absu, intstart, intend,id, istart-sst, iend+sst,present,future,increment) !compute interior
    
    commtime1 = MPI_WTIME()
    CALL MPI_WAITALL(4, reqs, MPI_STATUSES_IGNORE,ierr)
    commtime2 = MPI_WTIME()
    cumcommtime = cumcommtime + commtime2 - commtime1
    !WRITE(*,'(A15,F8.2)') 'recvbuf is', phi(imrbs:imrbe,present)
        
    !WRITE(*,020) id, 'Received buffer index is',imrbs,imrbe, & 
        !!phi(imrbs:imrbe,present),'iteration',iteration
        
    CALL SCHEME_POINTER(phi, dx, dt, absu, ipsbs, ipsbe,id, istart-sst, iend+sst,present,future,increment)

    CALL SCHEME_POINTER(phi, dx, dt, absu, imsbs, imsbe,id, istart-sst, iend+sst,present,future,increment)
    IF(id == 0) phi(istart,future)=phi(istart-1,future)  !!Simple way to ensure that f(-40,t) = original BC
    IF(id == commsize - 1) phi(iend,future)=phi(iend+1,future) !!Same for f(40,t) ghost cells out of comp domain are never modified from the initial conditions

    phi(:,present)=phi(:,future)
    iteration = iteration + 1;
    currentTime = currentTime + dt
    
    IF (currentTime >= controlTimes(j) )THEN

        saved_results(istart:iend,j) = phi (istart:iend,present)
        j = j+1
        IF(j >= SIZE(controlTimes)) EXIT
        dt = CFL*dx/absu
    END IF

    commtime(iteration) = cumcommtime

END DO

GOTO 529


510 CONTINUE !TIME LOOP STRUCTURE FOR THE MAC-CORMACK METHOD
!THE SENDS AND RECEIVES ARGUMENTS COULD BE SIMPLIFIED HERE BECAUSE OF THE SPECIFICITY OF THE TIME LOOP
!HOWEVER, THEY ARE LEFT THE SAME AS OTHER METHODS TO AVOID PERFORMANCE IMPROVEMENTS THAT COULD LEAD TO CONFUSION EVALUATING THE RESULTS

IF(u > 0.0) THEN
    correctorstart = iend   !in serial it's =intend because we directly exclude phi(0) & phi(npoints-1)
    correctorend = istart
ELSE IF(u < 0.0) THEN
        correctorstart = istart
        correctorend = iend
END IF

DO 
    dt=min(dt,controlTimes(j)-currentTime)
    !WRITE(*,005) id,'INIT it time dt', iteration, currentTime, dt
    commtime1 = MPI_WTIME()
    CALL MPI_ISEND(phi(ipsbs:ipsbe,present),sst,MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,reqs(1),ierr)  !upstream

    !WRITE(*,020) id, 'Sent buffer index is',imsbs,imsbe, & 
        !!phi(imsbs:imsbe,present),'iteration',iteration
    CALL MPI_IRECV(phi(iprbs:iprbe,present),sst, MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,reqs(2),ierr) !downstream
    CALL MPI_WAITALL(2, reqs(1:2), MPI_STATUSES_IGNORE,ierr)
    commtime2 = MPI_WTIME()
    cumcommtime = commtime2 - commtime1
    !! we do whole interior at once because the theoretical best coordinationb is not worth the overhead of building the extra do-loops and all the stuff
    
    DO i=istart, iend !PREDICTOR STEP INTO FUTURE COLUMN
        phi(i,future) = phi(i,present) - absu*dt/dx*(phi(i+increment,present) - phi(i,present))
        !phi(i,future) = phi(i,present) - absu*dt/dx*(phi(i+increment,present) - phi(i,present)) !!TESTEA QUE ESTO SOPORTA REVERSE FLUX
    END DO

    IF(id == 0) phi(istart,future)=phi(istart-1,future)  !!Simple way to ensure that f(-40,t) = original BC
    IF(id == commsize - 1) phi(iend,future)=phi(iend+1,future) !!Same for f(40,t) ghost cells out of comp domain are never modified from the initial conditions

    commtime1 = MPI_WTIME()
    CALL MPI_ISEND(phi(imsbs:imsbe,future),sst,MPI_DOUBLE_PRECISION,np,0,MPI_COMM_WORLD,reqs(3),ierr)
    CALL MPI_IRECV(phi(imrbs:imrbe,future),sst,MPI_DOUBLE_PRECISION,pp,0,MPI_COMM_WORLD,reqs(4),ierr)
    CALL MPI_WAITALL(2,reqs(3:4),MPI_STATUSES_IGNORE,ierr)
    commtime2 = MPI_WTIME()
    cumcommtime = cumcommtime + commtime2 - commtime1
    DO i=correctorstart,correctorend,-increment !!this doesn't work with reverse flow, just as the rest of the loop
        phi(i,future) = 0.5D0*(phi(i,present) + phi(i,future) - absu*dt/dx*(phi(i,future) - phi(i-increment,future)))
    END DO !Thanks to computing the subarray in reverse direction we can avoid using an extra array to store the predictor

    IF(id == 0) phi(istart,future)=phi(istart-1,future)  !!Simple way to ensure that f(-40,t) = original BC
    IF(id == commsize - 1) phi(iend,future)=phi(iend+1,future) !!Same for f(40,t) ghost cells out of comp domain are never modified from the initial conditions
  
    phi(:,present)=phi(:,future)
    iteration = iteration + 1;
    currentTime = currentTime + dt
    
    IF (currentTime >= controlTimes(j) )THEN
        saved_results(istart:iend,j) = phi (istart:iend,present)
        j = j+1
        IF(j >= SIZE(controlTimes)) EXIT
        dt = CFL*dx/absu
    END IF
    commtime(iteration) = cumcommtime
END DO

GOTO 529

!-----------------------------------------------------------
!BLOCK V: INFORMATION GATHERING AND REPORT FILES WRITING
!-----------------------------------------------------------

530 CONTINUE
! GATHERING DATA AND WRITING OUTPUT


!print*, 'started writing'
IF (id == 0) THEN
    ALLOCATE(receive_arr(0:npoints-1,ncontrolTimes))
    ALLOCATE(r_err(0:npoints-1,ncontrolTimes))
    ALLOCATE(r_ana(0:npoints-1,ncontrolTimes))
    ALLOCATE(L1_tot(ncontrolTimes))
    ALLOCATE(L2_tot(ncontrolTimes))
    ALLOCATE(LINF_overall(ncontrolTimes))
    ALLOCATE(ovavgcommtime(commsize))
    ALLOCATE(ovavgtimeloop(commsize))
END IF

IF(id == 0) gathertime1 = MPI_WTIME()

CALL MPI_GATHER(avgcomm,1,MPI_DOUBLE_PRECISION,ovavgcommtime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
CALL MPI_GATHER(avgtimeloop,1,MPI_DOUBLE_PRECISION,ovavgtimeloop,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)

!!! LOS GATHERV HAY QUE HACERLOS EN UN DO O HACER UN NUEVO TIPO DE DATO QUE SEA LA MATRIZ COMPLETA
DO j = 1, ncontrolTimes

    CALL ERROR(saved_results(istart:iend,j), istart,iend,FUNCTION_POINTER, present, x(istart:iend), u, controlTimes(j), & 
    istart,iend,future,id,error_array(istart:iend,j),analytical_res(istart:iend,j))

    CALL NORMS(error_array(istart:iend,j),L1(j),L2(j),LINF(j))

    !PRINT*,'PROC',id,'L1',L1(j),'L2',L2(j),'LINF',LINF(j)
    !WRITE(*,'(A5,I1,2(A6,I2), /,(2I3,3(A1,F5.1)))') 'PROC:',id,'istart',istart, 'iend', iend,&
    !        (id, i, 'a',analytical_res(i,j),'e',error_array(i,j),'f',phi(i,present), i=istart,iend)

    CALL MPI_GATHERV(analytical_res(istart:iend,j),sendsize,MPI_DOUBLE_PRECISION, & 
        r_ana(:,j), ssizes, displacements,MPI_DOUBLE_PRECISION,& 
        0,MPI_COMM_WORLD,ierr)
    !print*, 'analytical gatherv ok'
    CALL MPI_GATHERV(saved_results(istart:iend,j),sendsize,MPI_DOUBLE_PRECISION, & 
        receive_arr(:,j), ssizes, displacements,MPI_DOUBLE_PRECISION,& 
        0,MPI_COMM_WORLD,ierr)
    !print*, 'phi gatherv ok'
    CALL MPI_GATHERV(error_array(istart:iend,j),sendsize,MPI_DOUBLE_PRECISION, & 
        r_err(:,j), ssizes, displacements,MPI_DOUBLE_PRECISION,& 
        0,MPI_COMM_WORLD,ierr)
    !print*, 'error gatherv ok'
    CALL MPI_REDUCE(L1(j),L1_tot(j),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(L2(j),L2_tot(j),1,MPI_DOUBLE_PRECISION,MPI_SUM,0,MPI_COMM_WORLD,ierr)
    CALL MPI_REDUCE(LINF(j),LINF_overall(j),1,MPI_DOUBLE_PRECISION,MPI_MAX,0,MPI_COMM_WORLD,ierr)
    !print*, 'Reduces ok'
        
END DO

IF (id == 0) THEN
    gathertime2 = MPI_WTIME()
    gathertime = gathertime2 - gathertime1
    commavg = SUM(ovavgcommtime) / commsize
    commstdev = SQRT((commsize*SUM(ovavgcommtime**2) - SUM(ovavgcommtime)**2)/(commsize*(commsize-1)))
    timeloopavg = SUM(ovavgtimeloop) / commsize
    timeloopstdev = SQRT((commsize*SUM(ovavgtimeloop**2) - SUM(ovavgtimeloop)**2)/(commsize*(commsize-1)))
    writetime1 = MPI_WTIME()
    !WRITE(*,*) 'Writing output files'
    !WRITE(*,*) 'L1=', L1(:),'L2=',L2(:),'LINF=',LINF_overall(:)
    !WRITE(*,*) (i, receive_arr(i,1), i=0,npoints-1)
    WRITE(CFL_str,'(F8.1)') CFL
    !print*, 'CFL_str:', CFL_str
    
    !print*, 'time_str:',time_str

    WRITE(npoints_str,'(I15)') npoints
    !print*, 'npoints_str:', npoints_str
    WRITE(aux_str,'(A15)') npoints_str
    aux_str=ADJUSTL(aux_str)
    i1 = INDEX(aux_str,' ') - 1
    WRITE(i1_str,*) i1
    aux_str= "A" // TRIM(ADJUSTL(i1_str))
    format_str = "(A7,A3,A1,A3,A1,F4.2,A1," // aux_str(1:i1) // ",A1,F5.2,A1)"
    format_str = TRIM(format_str)

    DO j=1,ncontrolTimes !WE WRITE ONE TECPLOT-ASCII FILE FOR EVERY SAVED INSTANT
        WRITE(time_str,'(F8.2)') controlTimes(j)
        filename = TRIM(infunction) // '_' // TRIM(scheme) // '_' // TRIM(ADJUSTL(CFL_str)) // '_'&
        //  TRIM(ADJUSTL(npoints_str)) // '_' // TRIM(ADJUSTL(time_str)) // '.dat'
    
        OPEN(UNIT=100,FILE=filename,STATUS='NEW',ACTION='WRITE', IOSTAT=status,IOMSG=msg)

        WRITE(100,format_str) 'TITLE="', infunction, '_', scheme, '_',  CFL, '_', &
                            TRIM(ADJUSTL(npoints_str)),'_',controlTimes(j), '"'
        WRITE(100,*) 'VARIABLES = "x-nodes", "Analytical", "Numerical", "Error"'
        WRITE(100,*) 'ZONE'
        WRITE(100,format_str) 'T=   "', infunction, '_', scheme, '_',  CFL, '_', &
                            TRIM(ADJUSTL(npoints_str)),'_',controlTimes(j), '"'
        

        npoints_str= TRIM(ADJUSTL(npoints_str))
        format2_str= "(A2," // aux_str(1:i1) // ",A30,I0,A15,F6.2)"
        WRITE(100,format2_str) 'I=', npoints_str, ', DATAPACKING=POINT, STRANDID=',strandid,', SOLUTIONTIME=',controlTimes(j)
        !PRINT*, 'HEADER WRITTEN WITHOUT ERRORS'
        444 FORMAT(4(ES14.7,1X))

        WRITE(100,444) (x_tot(i),r_ana(i,j), receive_arr(i,j), r_err(i,j), i = 0,npoints-1)

        CLOSE(100)

    END DO

    i1 = INDEX(filename,' ') - 1 - 4
    filename = 'norms_' // filename(1:i1) // '.csv'
    OPEN(UNIT=101,FILE=filename,STATUS='NEW', ACTION='WRITE',IOSTAT=status,IOMSG=msg)
    WRITE(101,'(A12)') 't,L1,L2,LINF' !File header

    102 FORMAT(4(ES13.5,','))
    DO j = 1, ncontrolTimes
        WRITE(101,102) controlTimes(j),L1_tot(j),L2_tot(j),LINF_overall(j)
    END DO
    CLOSE(101)

    writetime2 = MPI_WTIME()

    INQUIRE(FILE='times.csv', EXIST=fileexist)

    IF(fileexist) THEN
        OPEN(UNIT=200,FILE='times.csv',STATUS='OLD',ACTION='WRITE',POSITION='APPEND',IOSTAT=status,IOMSG=msg)
    ELSE
        OPEN(UNIT=200,FILE='times.csv',STATUS='NEW',ACTION='WRITE',POSITION='APPEND',IOSTAT=status,IOMSG=msg)
        
        WRITE(200,'(A28,A110)') "scheme,function,npoints,CFL,",&
     "iterations,nprocs,do-time,itavgtime,itstdev,commavgtime,commstdev,extime,writetime,L1,L2,LINF,L1_m,L2_m,LINF_m"
    END IF
    extime2 = MPI_WTIME()
    201 FORMAT(2(A3,','),I0,',',F4.2,',',I0,',',I0,13(',',ES15.8))
    WRITE(200,201) scheme,infunction,npoints,CFL,iteration,commsize,timeloopavg*iteration,&
        timeloopavg,timeloopstdev,commavg,commstdev,extime2-extime1, writetime2-writetime1,&
        L1_tot(ncontrolTimes),L2_tot(ncontrolTimes),LINF_overall(ncontrolTimes),L1_tot(ncontrolTimes/2),&
        L2_tot(ncontrolTimes/2),LINF_overall(ncontrolTimes/2)
    CLOSE(200)


END IF

!GOTO 550

!550 CONTINUE

!WRITE(*,121) (id,i, phi(i,present), i=istart,iend)
!121 FORMAT(6(I3,I3,F8.2),/)
!IF(id == 0) THEN
!    extime2 = MPI_WTIME()
!    executiontime = extime2 - extime1
!END IF
CALL MPI_FINALIZE(ierr)
END PROGRAM parallel_linear_advection

!!PENDING TO DEFINE NEW MPI DATATYPE AND DO THE ISEND/IRECV USING IT. IF YOU DON'T THEN THE SPATIALST LARGER THAN 1 WILL NOT BE SENT

!------------------!
!     PENDING      !
!------------------!

!!!-ordenar los allocate y comentar los significados de las cosas mínimamente OK~
!!!-implementar el Maccormack OK, algún TVD NOT YET, THIRD ORDER UPWIND OK AUNQUE NO SE SI ES CORRECTO IGUAL QUE CON EL RESTO ME HACE FALTA CHEQUEAR CON EL SERIAL CODE
!-archivo de salida con el informe en tiempo de las normas, mira a ver si lo haces tecplot-readable OK LO HAGO CSV PARA LEER EN PYTHON
!-hacer una opción de escupir un paraview-readable (un csv con los títulos de las variables y todo escupido por columnas) A FUTURO
!-mírate el capítulo de estabilidad del hoffmann para la numerical viscosity, tienes que hacer el assessment a priori HECHO
!-mírate todos los tutoriales de tecplot incluido pytecplot
!-WARN THE USER ABOUT NUMERICAL DIFFUSION FOR MODIFYING DT WHEN CFL=1.0       NOT NECESSARY, BASTA CON COMENTAR, ES UN CASO PARTICULAR SEGÚN NUMERICAL DIFFUSION
!-revisa de nuevo el tema de los bstart y bend y cómo haces los boundaries en cada time loop structure   OK DECIDO DEJAR EL MACCORMACK COMO ESTÁ
!-testear si la barrera es necesaria   LA DEJO, NO PASA NADA
!-comprobar el lax ahora que la x no se va loca? PENDING
!-maccormack no funciona en reverse por lo de no usar un cacharro externo, arréglalo HECHO SIN TESTEAR
!-cambiar los nombres a los argumentos de las subrutinas, puedes llamarlos igual que en el main sin problema ALGUNOS SÍ
!-hacer el input con una NAMELIST. CAP 14 CHAPMAN   OK
!-intentarás enchufar un WENO de alto orden, por tus cojone
!-errores varios de configuración deben ser alertados al usuario. NO SEAS FLOJO
!-vuélvete a hacer el esquema de la nonblocking comm cuando tengas el serial y ves si puedes eliminar un wait o hacer un buffered send y no llamar dos veces a actualizar el boundary
!-LA NORMA L1 TIENE QUE SER VALOR ABSOLUTO