!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NAME: 
!!   autocorr
!! DESCRIPTION:
!!   calculates (non) priodic autocorrelation functions of arbitrary data
!! COMPILATION:
!!   make
!! USAGE:
!!   cat data.cat | ./calccorrfkt.x -p --periodic (periodic?)   (optional)
!!                                  -l --runlength <runlength> 
!!                                  -m --maxlag <max lag>       (optional)
!!                                  --mean <meanvalue>          (optional)
!!                                  -h --help (show help message)
!!                                  > corrfkt.out
!! AUTHOR:
!!    fuhl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM autocorr

USE, INTRINSIC :: ISO_FORTRAN_ENV
USE kinds, ONLY : dp, &
                  idp, &
                  cdp
USE re_de_alloc, ONLY : realloc
USE cmd_opt_parser

IMPLICIT NONE

REAL(KIND=dp), PARAMETER :: sqrt2 = SQRT(2.0_dp)

TYPE(cmd_opt_type), DIMENSION(:), POINTER :: cmd_opt

! Input constants
INTEGER :: runlength = 0
INTEGER :: maxlag = 0
LOGICAL :: laggiven = .FALSE.
LOGICAL :: periodic = .FALSE.
REAL(KIND=dp) :: meanvalue = 0.0_dp
LOGICAL :: readmean = .FALSE.
INTEGER(KIND=idp) :: ntotvals = 0

! FFT method storage
COMPLEX(KIND=cdp), DIMENSION(:), ALLOCATABLE :: ftsq
REAL(KIND=dp) :: tmpread
REAL(KIND=dp), DIMENSION(:), POINTER :: indata
INTEGER :: full_setnum = 0 
INTEGER :: max_setnum = 0
! Rest storage arrays
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: corrfkt
REAL(KIND=dp) :: variance = 0.0_dp
REAL(KIND=dp) :: autocorrtimeold= 0.0_dp, autocorrtimenew = 0.0_dp
REAL(KIND=dp) :: errorbar = 0.0_dp
REAL(KIND=dp) :: runavg = 0.0_dp
INTEGER :: maxintegral = 1
INTEGER, PARAMETER :: maxautcorrestiter = 25 ! Maximum iterations for autocorrelation time estimation

! counter
INTEGER :: i = 0
INTEGER :: icorrf = 0
INTEGER :: iset = 0

! FFT variables
INTEGER(8) :: planfwd,planbwd
COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: datain, dataout

!error variables
INTEGER :: readerr = 0
LOGICAL :: showhelp

#include "fftw3.f"

NULLIFY(cmd_opt)
NULLIFY(indata)

! parse command line arguments
CALL get_cmd_opt("-p", "--periodic", "For periodically defined data", cmd_opt, .FALSE., periodic)
CALL get_cmd_opt("-m", "--maxlag", "Maximum time difference for which the correlation "//&
                 "function should be computet", cmd_opt, 0, maxlag, laggiven)
CALL get_cmd_opt("-l", "--runlength", "Length of one timeseries on stdin", cmd_opt, 1, runlength)
CALL get_cmd_opt("", "--mean", "Average value over complete time data", cmd_opt, 0.0_dp, meanvalue, readmean)
CALL get_cmd_opt("-h", "--help", "Display this help", cmd_opt, .FALSE., showhelp)

IF (showhelp) THEN
   CALL show_help_msg(cmd_opt)
   STOP
END IF

CALL check_invalid_cmd_opt(cmd_opt, readerr)

CALL end_opt_parsing(cmd_opt)


IF (runlength < 1) THEN
   WRITE(UNIT=ERROR_UNIT,FMT=*) "runlength (-l) need to be bigger than 0"
   STOP
END IF

IF (.NOT.laggiven) THEN
   maxlag = runlength / 2
END IF

!Allocate memory for corrfkt
ALLOCATE(corrfkt(0:maxlag))
IF (periodic) THEN
   ALLOCATE(ftsq(0:runlength-1))
   ftsq(:) = CMPLX(0.0_dp,0.0_dp,dp)
   !FFT arrays
   ALLOCATE(datain(0:runlength-1))
   datain(:) = CMPLX(0.0_dp,0.0_dp,dp)
   ALLOCATE(dataout(0:runlength-1))
   dataout(:) = CMPLX(0.0_dp,0.0_dp,dp)
   CALL dfftw_plan_dft_1d(planfwd, &
                          runlength, &
                          datain(:), &
                          dataout(:), &
                          FFTW_FORWARD, &
                          FFTW_ESTIMATE)
   CALL dfftw_plan_dft_1d(planbwd, &
                          runlength, &
                          datain(:), &
                          dataout(:), &
                          FFTW_BACKWARD, &
                          FFTW_ESTIMATE)
ELSE
   ALLOCATE(ftsq(0:2*runlength-1))
   ftsq(:) = CMPLX(0.0_dp,0.0_dp,dp)
   !FFT arrays
   ALLOCATE(datain(0:2*runlength-1))
   datain(:) = CMPLX(0.0_dp,0.0_dp,dp)
   ALLOCATE(dataout(0:2*runlength-1))
   dataout(:) = CMPLX(0.0_dp,0.0_dp,dp)
   CALL dfftw_plan_dft_1d(planfwd, &
                          2*runlength, &
                          datain(:), &
                          dataout(:), &
                          FFTW_FORWARD, &
                          FFTW_ESTIMATE)
   CALL dfftw_plan_dft_1d(planbwd, &
                          2*runlength, &
                          datain(:), &
                          dataout(:), &
                          FFTW_BACKWARD, &
                          FFTW_ESTIMATE)
END IF

full_setnum = 0
max_setnum = 1
CALL realloc(indata, max_setnum*runlength)
readerr = 0
DO WHILE (readerr == 0)
   WRITE(UNIT=ERROR_UNIT,FMT='(x,A,I12)',ADVANCE='NO') &
      "Attempting to Read Data Set Nr. ", full_setnum + 1
   READ(UNIT=INPUT_UNIT,FMT=*,IOSTAT=readerr) tmpread
   IF (readerr == 0) THEN
      full_setnum = full_setnum + 1
      IF (full_setnum > max_setnum) THEN
         max_setnum = CEILING(max_setnum*sqrt2)
         CALL realloc(indata, max_setnum*runlength)
      END IF
      indata((full_setnum - 1)*runlength + 1) = tmpread
      
      DO i = 2, runlength
         READ(UNIT=INPUT_UNIT,FMT=*,IOSTAT=readerr) tmpread
         IF (readerr == 0) THEN
            indata((full_setnum - 1)*runlength + i) = tmpread
         ELSE
            full_setnum = full_setnum - 1
            EXIT
         END IF
      END DO 
   END IF
   IF (readerr /= 0) THEN
      WRITE(UNIT=ERROR_UNIT,FMT='(A)') " [ FAILED ]"
      EXIT
   END IF
   WRITE(UNIT=ERROR_UNIT,FMT='(A)')    " [   OK   ]"
END DO
CALL realloc(indata, full_setnum*runlength)

ntotvals = INT(full_setnum*runlength,idp)
! precalculate mean 
IF (.NOT.readmean) THEN
   readmean = .TRUE.
   meanvalue = SUM(indata) / REAL(ntotvals,dp)
END IF
indata(:) = indata(:) - meanvalue
corrfkt(:) = 0.0_dp

!loop over stdinput until it ends
readerr = 0
DO iset = 0, full_setnum-1
   WRITE(UNIT=ERROR_UNIT, FMT=*) "Processing Data Set Nr.         ", iset + 1
   IF (periodic) THEN
      ! Copy data to fourier array and zeropad
      DO i = 1, runlength
         datain(i-1) = CMPLX(indata(iset*runlength + i),0.0_dp,dp)
      END DO
      ! execute FFTW
      CALL dfftw_execute_dft(planfwd, datain(:), dataout(:))
      ftsq(:) = ftsq(:) + dataout(:)*CONJG(dataout(:))
   ELSE
      ! Copy data to fourier array and zeropad
      DO i = 1, runlength
         datain(i-1) = CMPLX(indata(iset*runlength + i),0.0_dp,dp)
      END DO
      DO i = runlength+1, 2*runlength
         datain(i-1) = CMPLX(0.0_dp,0.0_dp,dp)
      END DO
      ! execute FFTW
      CALL dfftw_execute_dft(planfwd, datain(:), dataout(:))
      ftsq(:) = ftsq(:) + dataout(:)*CONJG(dataout(:))
   END IF
END DO


! summ up data and backtransform
datain = ftsq(:)
CALL dfftw_execute_dft(planbwd, datain(:), dataout(:))
!End normalization
IF (periodic) THEN
   DO i = 0, maxlag
      corrfkt(i) = REAL(dataout(i),dp)/(REAL(ntotvals*runlength*full_setnum,dp))
   END DO
ELSE
   DO i = 0, maxlag
      corrfkt(i) = REAL(dataout(i),dp) / (2.0_dp*REAL(ntotvals*(runlength-i),dp))
   END DO
END IF
variance = DOT_PRODUCT(indata(:), &
                       indata(:)) & 
                       / REAL(ntotvals,dp)

corrfkt(:) = corrfkt(:) / variance

! calculate running average as first guess for integration limit for autocorrelation time guess
DO icorrf = 0,maxlag
   runavg = runavg + corrfkt(icorrf)
   maxintegral = icorrf + 1
   IF (ABS(REAL(runavg,dp) / REAL(maxintegral,dp)) < 0.01_dp) EXIT
END DO
DO icorrf = 0,MIN(maxintegral,maxlag - 1)
   runavg = runavg + 0.5_dp*(corrfkt(icorrf) + corrfkt(icorrf + 1))
END DO

! calculate autocorrelation time selfconsistently as integral from 0 to six times autocorrelation 
autocorrtimenew = 0.0_dp
autocorrtimeold = 1.0_dp
WRITE(UNIT=ERROR_UNIT, FMT=*) "Estimating Autocorr time: ", &
   autocorrtimeold
DO icorrf = 0,MIN(maxintegral,maxlag - 1)
   autocorrtimenew = autocorrtimenew + 0.5_dp*(corrfkt(icorrf) + corrfkt(icorrf + 1))
END DO
autocorrtimenew = autocorrtimenew - 0.5_dp
i = 0
DO WHILE (ABS(ABS(autocorrtimenew/autocorrtimeold) - 1.0_dp) > 0.0001_dp)
   i = i + 1
   WRITE(UNIT=ERROR_UNIT, FMT=*) "Estimating Autocorr time: ", &
      autocorrtimenew, "(",ABS(ABS(autocorrtimenew/autocorrtimeold) - 1.0_dp),")"
   autocorrtimeold = autocorrtimenew
   maxintegral = INT(6.0_dp * autocorrtimenew +0.5_dp)
   autocorrtimenew = 0.0_dp
   DO icorrf = 0,MIN(maxintegral,maxlag - 1)
      autocorrtimenew = autocorrtimenew + 0.5_dp*(corrfkt(icorrf) + corrfkt(icorrf + 1))
   END DO
   autocorrtimenew = autocorrtimenew - 0.5_dp
   autocorrtimenew = 0.5_dp*(autocorrtimenew + autocorrtimeold)
   IF (i >= maxautcorrestiter) THEN
      WRITE(UNIT=ERROR_UNIT, FMT=*) "Reached Maximum Iterations for autocorrelation time estimation"
      IF (ABS(ABS(autocorrtimenew/autocorrtimeold) - 1.0_dp) > 0.0001_dp) THEN
         WRITE(UNIT=ERROR_UNIT, FMT=*) "WARNING: Autocorrelation time estimation did not converge"
         EXIT
      END IF
   END IF
END DO
WRITE(UNIT=ERROR_UNIT, FMT=*) "Estimating Autocorr time: ", &
   autocorrtimenew, "(",ABS(ABS(autocorrtimenew/autocorrtimeold) - 1.0_dp),")"

! Estimation of error bars
errorbar = SQRT((2.0_dp*MAX(0.5_dp,autocorrtimenew))*variance/REAL(ntotvals,dp))

WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# Number of Sets = ", full_setnum
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# Runlength =      ", runlength
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# maxlag =         ", maxlag
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# mean =           ", meanvalue
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# error = +/-      ", errorbar
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# variance =       ", variance
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# stddeviation =   ", SQRT(variance)
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# autocorrtime =   ", autocorrtimenew
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "#"

runavg = 0.0_dp
DO icorrf = 0,maxlag
   WRITE(UNIT=OUTPUT_UNIT,FMT=*) icorrf, corrfkt(icorrf)
END DO
IF ((periodic).AND.(maxlag + 1 == runlength)) THEN
   WRITE(UNIT=OUTPUT_UNIT,FMT=*) maxlag + 1, corrfkt(0)
END IF

DEALLOCATE(ftsq)
DEALLOCATE(corrfkt)
!FFT arrays
DEALLOCATE(datain, dataout)
CALL dfftw_destroy_plan(planfwd)
CALL dfftw_destroy_plan(planbwd)
DEALLOCATE(indata)
NULLIFY(indata)

END PROGRAM autocorr
