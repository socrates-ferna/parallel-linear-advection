PROGRAM serial_linear_advection
USE numerical_schemes
USE input_functions
USE misc_subroutines
!-----------------------------------------------------------------------
! Purpose: 1D code for serial resolution of the linear advection equation (wave equation)
! using several numerical schemes. The analytical solution is used to verify the numerical results.
!
! Licence: This code is distributed under GNU GPL Licence
! Author: Sócrates Fernández Fernández, s(dot)fernaferna(at)gmail(dot)com
! GitHub: socrates-ferna
! LAST MOD: 16/12/2020
!-----------------------------------------------------------------------
! A SHORT COMMMENT ON THE CODE:
! I developed this code after the parallel one, therefore, some variables,
! arrays and tweaks are not elegant but needed to reuse subroutines from
! the parallel modules. e.g. x_int array is needed because the arrays in 
! ERROR and ANALYTICAL subroutines are not assumed-shape, so I need this
! extra array to fit the dimensions 0:npoints-1 which in the parallel code
! corresponded to each processor's interval start and end
!-----------------------------------------------------------------------
!! LISTA para actualizar en parallel
!- niters_guess OK
!- arreglo reverse flux mccormack OK NOT TESTED
!- todos los cputime copia la estructura OK NOT TESTED
!- línea de execution time variables entera + iterationtimes array NO
!- doublesizearray subroutine and call inside loop NO
!- avg and whole time report OK
!-fileexist, writetime1 y 2 SI
!- CFL_str dos decimales SI
!- NAMELIST INPUT OK
IMPLICIT NONE
INTEGER :: i, j, npoints, iteration, past, present, future, ncontrolTimes, i1, increment, status, strandid, niters_guess
INTEGER :: correctorstart, correctorend, mstart, mend
REAL(KIND=8) :: u, CFL, xl, xr, dx, dt, currentTime,stopTime, absu
REAL(KIND=8) :: extime1, extime2, timeloop1, timeloop2, iterationtime1, iterationtime2, writetime1, writetime2 !!EXECUTION TIME VARIABLES
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: iterationtimes, aux_array
REAL(KIND=8), DIMENSION(:), ALLOCATABLE :: controlTimes, x, L1, L2, LINF, x_int
REAL(KIND=8), DIMENSION(:,:), ALLOCATABLE :: phi, error_array, analytical_res, saved_results
LOGICAL :: special_init, fileexist
CHARACTER(15) :: scheme, infunction,CFL_str,time_str,npoints_str, aux_str, i1_str
CHARACTER(300) :: msg,filename, format_str, format2_str
PROCEDURE(UPWIND), POINTER :: SCHEME_POINTER
PROCEDURE(SGN), POINTER :: FUNCTION_POINTER
INTEGER :: sst = 1,tst = 1, id = 0, intstart, intend !!!MAKES SUBROUTINES WRITTEN FOR PARALLEL PROGRAM COMPATIBLE WITH SERIAL
NAMELIST / input_list / u, scheme, infunction, CFL, npoints, xl, xr, ncontrolTimes, stopTime, strandid

CALL CPU_TIME(extime1)

OPEN(UNIT=100, FILE='input.nml', DELIM='APOSTROPHE')
READ(UNIT=100,NML=input_list)
CLOSE(100)

WRITE(UNIT=*,NML=input_list)

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
niters_guess = NINT(stopTime/dt + REAL(ncontrolTimes, KIND=8))
ALLOCATE(iterationtimes(niters_guess)); iterationtimes(:) = 0.0D0
present = 0 !n time index
future = 1  !n+1 time index (also used as n+1/2 in maccormack)
special_init = .FALSE.   !LEAPFROG REQUIRES INITIALISATION WITH OTHER SCHEME
mstart = 0 - sst !!memory-allocated index start
mend = npoints -1 + sst !!memory-allocated index end
intstart = 1
intend = npoints - 2

CALL SCHEMESELECTION(scheme,SCHEME_POINTER,sst,tst,present,future,past,special_init)
CALL FUNCTIONSELECTION(infunction,FUNCTION_POINTER)

ALLOCATE(x(mstart:mend)); x(:) = 0.D0
ALLOCATE(x_int(0:npoints-1));
ALLOCATE(phi(mstart:mend,0:tst));phi(:,:) = 0.0D0

DO i = mstart, mend
    x(i) = xl + i*dx
END DO
x_int(:) = x(0:npoints-1) !! ERROR AND ANALYTICAL SUBROUTINES need an array with this dimensions

CALL ANALYTICAL(phi(:,present), mstart,mend,FUNCTION_POINTER, &
                present, x, u, currentTime, mstart,mend,future,id)
CALL ANALYTICAL(phi(:,future), mstart,mend,FUNCTION_POINTER, &
                present, x, u, currentTime, mstart,mend,future,id)
IF( tst > 1) THEN
    CALL ANALYTICAL(phi(:,past), mstart,mend,FUNCTION_POINTER, &
                present, x, u, currentTime, mstart,mend,future,id)
END IF

IF(u > 0.0D0) THEN
    increment = 1
ELSE IF(u  < 0.0D0) THEN
    increment = -1
ELSE IF (u == 0.0) THEN
    PRINT*, 'u=0, no solution. Aborting'
    STOP
END IF

absu = abs(u)

ALLOCATE(error_array(0:npoints-1,ncontrolTimes));error_array(:,:) = 3.0D0 
ALLOCATE(analytical_res(0:npoints-1,ncontrolTimes));analytical_res(:,:) = 8.0D0  
ALLOCATE(L1(ncontrolTimes));L1(:) = 5.0D0   
ALLOCATE(L2(ncontrolTimes));L2(:) = 6.0D0
ALLOCATE(LINF(ncontrolTimes));LINF(:) = 70.D0
ALLOCATE(saved_results(0:npoints-1,ncontrolTimes)); saved_results(:,:) = 211.0D0

SELECT CASE (scheme)
CASE('mcc')
    GOTO 510
CASE('lpf')
    SCHEME_POINTER => UPWIND
    DEALLOCATE(controlTimes)
    ALLOCATE(controlTimes(1)); controlTimes(1) = dt
    !!! FALTA DESASIGNAR TODOS LOS ARRAYS DEPENDIENTES DEL CONTADOR J Y
    !!! HACER LA DETERMINACIÓN DE SI SOY LEAPFROG SALTO  A LA DEFAULT STRUCTURE HACmendO IF(special_init)
    !GOTO 500
    !501 CONTINUE
CASE DEFAULT
    GOTO 505
END SELECT

505 CONTINUE

CALL CPU_TIME(timeloop1)
DO 
    CALL CPU_TIME(iterationtime1)
    dt=min(dt,controlTimes(j)-currentTime)

    CALL SCHEME_POINTER(phi, dx, dt, absu, intstart, intend,id, mstart,mend,present,future,increment) !compute interior

    phi(intstart:intend,present)=phi(intstart:intend,future)
    iteration = iteration + 1
    currentTime = currentTime + dt
    
    IF ( currentTime >= controlTimes(j) )THEN
        saved_results(0:npoints-1,j) = phi (0:npoints-1,present)
        j = j+1
        IF( j >= SIZE(controlTimes) ) EXIT
        dt = CFL*dx/absu
    END IF

    IF( iteration == SIZE(iterationtimes) ) THEN !! In case number of iterations exceeds initial guess
        CALL DOUBLESIZEARRAY(iterationtimes,aux_array)
    END IF

    CALL CPU_TIME(iterationtime2)
    iterationtimes(iteration) = iterationtime2 - iterationtime2
END DO

CALL CPU_TIME(timeloop2)
GOTO 530



510 CONTINUE

IF(u > 0.0) THEN
correctorstart = intend
correctorend = intstart
ELSE IF(u < 0.0) THEN
    correctorstart = intstart
    correctorend = intend
END IF

CALL CPU_TIME(timeloop1)

DO 
    CALL CPU_TIME(iterationtime1)
    dt=min(dt,controlTimes(j)-currentTime)
    DO i=intstart, intend !PREDICTOR STEP INTO FUTURE COLUMN
        phi(i,future) = phi(i,present) - absu*dt/dx*(phi(i+increment,present) - phi(i,present))
    END DO

    DO i=correctorstart,correctorend,-increment !!this doesn't work with reverse flow, just as the rest of the loop
        phi(i,future) = 0.5D0*(phi(i,present) + phi(i,future)) - absu*dt/dx*(phi(i,future) - phi(i-increment,future))
    END DO !Thanks to computing the subarray in reverse direction we can avoid using an extra array to store the predictor
 
    phi(intstart:intend,present)=phi(intstart:intend,future)
    iteration = iteration + 1;
    currentTime = currentTime + dt

    IF (currentTime >= controlTimes(j) )THEN
        saved_results(0:npoints-1,j) = phi (0:npoints-1,present)
        j = j+1
        IF(j >= SIZE(controlTimes)) EXIT
        dt = CFL*dx/absu
    END IF

    IF( iteration == SIZE(iterationtimes) ) THEN !! In case number of iterations exceeds initial guess
        CALL DOUBLESIZEARRAY(iterationtimes,aux_array)
    END IF

    CALL CPU_TIME(iterationtime2)
    iterationtimes(iteration) = iterationtime2 - iterationtime2
END DO
CALL CPU_TIME(timeloop2)
GOTO 530



530 CONTINUE

DO j = 1, ncontrolTimes

    CALL ERROR(saved_results(:,j), 0,npoints-1,FUNCTION_POINTER, present, x_int(:), u, controlTimes(j), & 
    0,npoints-1,future,id,error_array(:,j),analytical_res(:,j))

    CALL NORMS(error_array(:,j),L1(j),L2(j),LINF(j))

    PRINT*,'PROC',id,'L1',L1(j),'L2',L2(j),'LINF',LINF(j)
    WRITE(*,'(A5,I1, /,(2I3,3(A1,F5.1)))') 'PROC:',id,&
            (id, i, 'a',analytical_res(i,j),'e',error_array(i,j),'f',phi(i,present), i=0,npoints-1)
        
END DO

CALL CPU_TIME(writetime1)
WRITE(*,*) 'Writing output files'
WRITE(*,*) 'L1=', L1(:),'L2=',L2(:),'LINF=',LINF(:)
WRITE(*,*) (i, saved_results(i,1), i=0,npoints-1)
WRITE(CFL_str,'(F8.2)') CFL
print*, 'CFL_str:', CFL_str

print*, 'time_str:',time_str

WRITE(npoints_str,'(I15)') npoints
print*, 'npoints_str:', npoints_str
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

    
    !print*,'i1', i1, 'aux_str', aux_str
    
    !print*, 'aux_str posconcat', aux_str


    !print*, 'format_str', format_str
    !print*, 'aux_str', aux_str
    WRITE(100,format_str) 'TITLE="', infunction, '_', scheme, '_',  CFL, '_', &
                        TRIM(ADJUSTL(npoints_str)),'_',controlTimes(j), '"'
    WRITE(100,*) 'VARIABLES = "x-nodes", "Analytical", "Numerical", "Error"'
    WRITE(100,*) 'ZONE'
    WRITE(100,format_str) 'T=   "', infunction, '_', scheme, '_',  CFL, '_', &
                        TRIM(ADJUSTL(npoints_str)),'_',controlTimes(j), '"'
    

    npoints_str= TRIM(ADJUSTL(npoints_str))
    format2_str= "(A2," // aux_str(1:i1) // ",A30,I0,A15,F6.2)"
    WRITE(100,format2_str) 'I=', npoints_str, ', DATAPACKING=POINT, STRANDID=',strandid,', SOLUTIONTIME=',controlTimes(j)
    PRINT*, 'HEADER WRITTEN WITHOUT ERRORS'
    444 FORMAT(4(ES14.7,1X))

    WRITE(100,444) (x(i),analytical_res(i,j), saved_results(i,j), error_array(i,j), i = 0,npoints-1)

    CLOSE(100)

END DO

i1 = INDEX(filename,' ') - 1 - 4
filename = 'norms_' // filename(1:i1) // '.csv'
OPEN(UNIT=101,FILE=filename,STATUS='NEW', ACTION='WRITE',IOSTAT=status,IOMSG=msg)
WRITE(101,'(A12)') 't,L1,L2,LINF' !File header for norms and control times

102 FORMAT(4(F9.5,','))
DO j = 1, ncontrolTimes
    WRITE(101,102) controlTimes(j),L1(j),L2(j),LINF(j)
END DO
CLOSE(101)

CALL CPU_TIME(writetime2)
!!Time reports

INQUIRE(FILE='times.csv', EXIST=fileexist)

IF(fileexist) THEN
    OPEN(UNIT=200,FILE='times.csv',STATUS='OLD',ACTION='WRITE',POSITION='APPEND',IOSTAT=status,IOMSG=msg)
ELSE
    OPEN(UNIT=200,FILE='times.csv',STATUS='NEW',ACTION='WRITE',POSITION='APPEND',IOSTAT=status,IOMSG=msg)
    WRITE(200,*) "scheme,function,npoints,CFL,do-time,itavgtime,extime,writetime,L1,L2,LINF"
END IF
                                
CALL CPU_TIME(extime2)

201 FORMAT(2(A3,','),I0,',',F4.2,7(',',ES15.8))
WRITE(200,201) scheme,infunction,npoints,CFL,timeloop2-timeloop1,&
    (timeloop2-timeloop1)/iteration,extime2-extime1,writetime2-writetime1,L1(ncontrolTimes),L2(ncontrolTimes),&
    LINF(ncontrolTimes)
CLOSE(200)

CALL CPU_TIME(extime2)
PRINT*,'Total execution time* is: ',extime2-extime1
END PROGRAM serial_linear_advection
