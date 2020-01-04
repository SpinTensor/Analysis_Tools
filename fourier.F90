!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NAME: 
!!   fourier
!! DESCRIPTION:
!!   calculates (non) priodic functions of arbitrary
!data
!! COMPILATION:
!!   make
!! USAGE:
!!   cat data.dat | ./fourier.x -p --periodic (periodic?)   (optional)
!!                                 -dt timestep
!!                                 -h --help (show help message)
!!                                 > spectrum.out
!! AUTHOR:
!!   Felix Uhl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM fourier
USE, INTRINSIC :: ISO_FORTRAN_ENV
USE kinds, ONLY : dp, &
                  cdp, &
                  default_string_length
USE re_de_alloc, ONLY : realloc
USE cmd_opt_parser

IMPLICIT NONE

! load FFT3 definitions
#include "fftw3.f"

REAL(KIND=dp), PARAMETER :: pi = 2.0_dp*ACOS(0.0_dp)

TYPE(cmd_opt_type), DIMENSION(:), POINTER :: cmd_opt

! Input constants
REAL(KIND=dp) :: dt = 1.0_dp
LOGICAL :: periodic = .FALSE.
CHARACTER(LEN=default_string_length) :: window = "rectangle"
LOGICAL :: showhelp = .FALSE.

INTEGER :: maxpoints, npoints
REAL(KIND=dp), DIMENSION(:), POINTER :: indata
REAL(KIND=dp) :: tmpread

! FFTW variables
INTEGER(KIND=dp) :: FWD_FFTW3_PLAN
!REAL(KIND=dp), DIMENSION(:), POINTER :: outdata
COMPLEX(KIND=cdp), DIMENSION(:), POINTER :: outdata
REAL(KIND=dp) :: norm

REAL(KIND=dp) :: w, nr
REAL(KIND=dp) :: freq, rep, imp, amp

INTEGER :: i

INTEGER :: readerr = 0

NULLIFY(cmd_opt)

! parse command line arguments
CALL get_cmd_opt("-dt", "--deltat", "Timestep for input data", cmd_opt, 1.0_dp, dt)
CALL get_cmd_opt("-p", "--periodic", "For periodically defined data", cmd_opt, .FALSE., periodic)
CALL get_cmd_opt("-w", "--window", "Window function to be used on the data.           "//&
                                   "Possible values:                                  "//&
                                   "   rectangle, triangle, sin, cos, nutall", cmd_opt, "rectangle", window)
CALL get_cmd_opt("-h", "--help", "Display this help", cmd_opt, .FALSE., showhelp)

IF (showhelp) THEN
   CALL show_help_msg(cmd_opt)
   STOP
END IF

CALL check_invalid_cmd_opt(cmd_opt, readerr)

CALL end_opt_parsing(cmd_opt)

!prepare input data arrays
NULLIFY(indata)

npoints=0
maxpoints=10

ALLOCATE(indata(maxpoints))

! read in the data
readerr = 0
DO WHILE (readerr == 0)
   READ(UNIT=INPUT_UNIT, FMT=*, IOSTAT=readerr) tmpread
   IF (readerr /= 0) EXIT
   npoints = npoints + 1
   !reallocate on demand
   IF (npoints > maxpoints) THEN
      maxpoints = (14*maxpoints)/10
      CALL realloc(indata, maxpoints)
   END IF
   indata(npoints) = tmpread
END DO

WRITE(UNIT=OUTPUT_UNIT, FMT=*) "#Fourier Transform"
!abort if no point was read
IF (npoints == 0) THEN
   WRITE(UNIT=ERROR_UNIT, FMT=*) "No input data provided"
   STOP
END IF

WRITE(UNIT=OUTPUT_UNIT, FMT=*) "#Number of datapoints: ", npoints
WRITE(UNIT=OUTPUT_UNIT, FMT=*) "#Timestep: ", dt

! apply window function
WRITE(UNIT=OUTPUT_UNIT, FMT=*) "#Window function: ", TRIM(window)
SELECT CASE(TRIM(window))
CASE("rectangle")
   CONTINUE
CASE("triangle")
   nr = REAL(npoints-1,dp)
   DO i = 1, npoints
      w = 1.0_dp-ABS(((i-1)-0.5_dp*nr)/(0.5_dp*nr))
      indata(i) = indata(i)*w
   END DO
CASE("sin")
   nr = 1.0_dp/REAL(npoints-1,dp)
   DO i = 1, npoints
      w = SIN(pi*(i-1)*nr)
      indata(i) = indata(i)*w
   END DO
CASE("cos")
   nr = 1.0_dp/REAL(npoints-1,dp)
   DO i = 1, npoints
      w = 0.5_dp - 0.5_dp*COS(2.0_dp*pi*(i-1)*nr)
      indata(i) = indata(i)*w
   END DO
CASE("nuttall")
   nr = 1.0_dp/REAL(npoints-1,dp)
   DO i = 1, npoints
      w =   0.355768_dp &
          - 0.487396_dp*COS(2*pi*(i-1)*nr) &
          + 0.144232_dp*COS(4*pi*(i-1)*nr) &
          - 0.012604_dp*COS(6*pi*(i-1)*nr)
      indata(i) = indata(i)*w
   END DO
CASE DEFAULT
   CONTINUE
END SELECT

! final adjustment for the input data size
WRITE(UNIT=OUTPUT_UNIT, FMT=*) "#Periodic: ", periodic
IF (periodic) THEN
   CALL realloc(indata, npoints)
ELSE
   CALL realloc(indata, 2*npoints)
   ! zero pedding
   indata(npoints+1:2*npoints) = 0.0_dp
   npoints = npoints*2
END IF

! prepare norm for fftw
norm= 1.0_dp/SQRT(REAL(npoints))

NULLIFY(outdata)
ALLOCATE(outdata(npoints/2+1))

CALL dfftw_plan_dft_r2c_1d(FWD_FFTW3_PLAN, & !FFTW plan
                           npoints, &        !number of elements
                           indata, &         !input data
                           outdata, &        !output data
                           FFTW_ESTIMATE)    !FFTW algorithm

CALL dfftw_execute(FWD_FFTW3_PLAN, indata, outdata)

CALL dfftw_destroy_plan(FWD_FFTW3_PLAN)

WRITE(UNIT=OUTPUT_UNIT, FMT=*) "#          Frequency"//&
                               "               Real-coeffs"//&
                               "               Imag-coeffs"//&
                               "                 Amplitude"
outdata = norm*outdata
DO i = 1, npoints/2+1
   freq = REAL(i-1,dp)/(REAL(npoints,dp)*dt)
   rep  = REAL(outdata(i))
   imp  = AIMAG(outdata(i))
   amp  = SQRT(rep*rep+imp*imp)
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) freq, rep, imp, amp
END DO

DEALLOCATE(indata)
NULLIFY(indata)
DEALLOCATE(outdata)
NULLIFY(outdata)

END PROGRAM fourier
