!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NAME: 
!!   autocorr
!! DESCRIPTION:
!!   calculates (non) priodic autocorrelation functions of arbitrary data
!! COMPILATION:
!!   make
!! USAGE:
!!   cat data.cat | ./calccorrfkt.x -p --periodic (periodic?)        (optional)
!!                                  -l --runlength <runlength> 
!!                                  -m --maxlag <max lag>            (optional)
!!                                  --mean <meanvalue(1)> <meanvalue(2)> (optional)
!!                                  > corrfkt.out
!! AUTHOR:
!!    fuhl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM crosscorr

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
REAL(KIND=dp), DIMENSION(0) :: defmean
REAL(KIND=dp), DIMENSION(:), POINTER :: meanvalue
LOGICAL :: readmean = .FALSE.
INTEGER(idp) :: ntotvals = 0

! FFT method storage
COMPLEX(cdp), DIMENSION(:), ALLOCATABLE :: ftprod
REAL(KIND=dp) :: tmpreadA, tmpreadB
REAL(KIND=dp), DIMENSION(:), POINTER :: indataA, indataB
INTEGER :: full_setnum = 0 
INTEGER :: max_setnum = 0
! Rest storage arrays
REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: corrfkt
REAL(KIND=dp) :: varianceA = 0.0_dp, varianceB = 0.0_dp

! counter
INTEGER :: i = 0
INTEGER :: icorrf = 0
INTEGER :: iset = 0

! FFT variables
INTEGER(8) :: planfwdA,planbwdA
INTEGER(8) :: planfwdB,planbwdB
COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: datainA, dataoutA
COMPLEX(KIND=8), DIMENSION(:), ALLOCATABLE :: datainB, dataoutB

!error variables
INTEGER :: readerr = 0
LOGICAL :: showhelp

#include "fftw3.f"

NULLIFY(indataA)
NULLIFY(indataB)
NULLIFY(meanvalue)
NULLIFY(cmd_opt)

! parse command line arguments
CALL get_cmd_opt("-p", "--periodic", "For periodically defined data", cmd_opt, .FALSE., periodic)
CALL get_cmd_opt("-m", "--maxlag", "Maximum time difference for which the correlation "//&
                 "function should be computet", cmd_opt, 0, maxlag, laggiven)
CALL get_cmd_opt("-l", "--runlength", "Length of one timeseries on stdin", cmd_opt, 1, runlength)
CALL get_cmd_opt("", "--mean", "Average value over complete time data", cmd_opt, defmean, meanvalue, readmean)
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

IF (readmean) THEN
   IF (SIZE(meanvalue) /= 2) THEN
      WRITE(UNIT=ERROR_UNIT, FMT=*) "Need exactly two meanvalues"
      STOP
   END IF
END IF

IF (.NOT.laggiven) THEN
   maxlag = runlength / 2
END IF

!Allocate memory for corrfkt
ALLOCATE(corrfkt(0:maxlag))
IF (periodic) THEN
   ALLOCATE(ftprod(0:runlength-1))
   ftprod(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   !FFT arrays
   ALLOCATE(datainA(0:runlength-1))
   datainA(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   ALLOCATE(datainB(0:runlength-1))
   datainB(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   ALLOCATE(dataoutA(0:runlength-1))
   dataoutA(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   ALLOCATE(dataoutB(0:runlength-1))
   dataoutB(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   CALL dfftw_plan_dft_1d(planfwdA, &
                          runlength, &
                          datainA(:), &
                          dataoutA(:), &
                          FFTW_FORWARD, &
                          FFTW_ESTIMATE)
   CALL dfftw_plan_dft_1d(planbwdA, &
                          runlength, &
                          datainA(:), &
                          dataoutA(:), &
                          FFTW_BACKWARD, &
                          FFTW_ESTIMATE)

   CALL dfftw_plan_dft_1d(planfwdB, &
                          runlength, &
                          datainB(:), &
                          dataoutB(:), &
                          FFTW_FORWARD, &
                          FFTW_ESTIMATE)
   CALL dfftw_plan_dft_1d(planbwdB, &
                          runlength, &
                          datainB(:), &
                          dataoutB(:), &
                          FFTW_BACKWARD, &
                          FFTW_ESTIMATE)
ELSE
   ALLOCATE(ftprod(0:2*runlength-1))
   ftprod(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   !FFT arrays
   ALLOCATE(datainA(0:2*runlength-1))
   datainA(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   ALLOCATE(datainB(0:2*runlength-1))
   datainB(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   ALLOCATE(dataoutA(0:2*runlength-1))
   dataoutA(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   ALLOCATE(dataoutB(0:2*runlength-1))
   dataoutB(:) = CMPLX(0.0_dp,0.0_dp,cdp)
   CALL dfftw_plan_dft_1d(planfwdA, &
                          2*runlength, &
                          datainA(:), &
                          dataoutA(:), &
                          FFTW_FORWARD, &
                          FFTW_ESTIMATE)
   CALL dfftw_plan_dft_1d(planbwdA, &
                          2*runlength, &
                          datainA(:), &
                          dataoutA(:), &
                          FFTW_BACKWARD, &
                          FFTW_ESTIMATE)
   CALL dfftw_plan_dft_1d(planfwdB, &
                          2*runlength, &
                          datainB(:), &
                          dataoutB(:), &
                          FFTW_FORWARD, &
                          FFTW_ESTIMATE)
   CALL dfftw_plan_dft_1d(planbwdB, &
                          2*runlength, &
                          datainB(:), &
                          dataoutB(:), &
                          FFTW_BACKWARD, &
                          FFTW_ESTIMATE)
END IF

full_setnum = 0
max_setnum = 1
CALL realloc(indataA, max_setnum*runlength)
CALL realloc(indataB, max_setnum*runlength)
readerr = 0
DO WHILE (readerr == 0)
   WRITE(UNIT=ERROR_UNIT,FMT='(x,A,I12)',ADVANCE='NO') &
      "Attempting to Read Data Set Nr. ", full_setnum + 1
   READ(UNIT=INPUT_UNIT,FMT=*,IOSTAT=readerr) tmpreadA, tmpreadB
   IF (readerr == 0) THEN
      full_setnum = full_setnum + 1
      IF (full_setnum > max_setnum) THEN
         max_setnum = CEILING(max_setnum*sqrt2)
         CALL realloc(indataA, max_setnum*runlength)
         CALL realloc(indataB, max_setnum*runlength)
      END IF
      indataA((full_setnum - 1)*runlength + 1) = tmpreadA
      indataB((full_setnum - 1)*runlength + 1) = tmpreadB

      DO i = 2, runlength
         READ(UNIT=INPUT_UNIT,FMT=*,IOSTAT=readerr) tmpreadA, tmpreadB
         IF (readerr == 0) THEN
            indataA((full_setnum - 1)*runlength + i) = tmpreadA
            indataB((full_setnum - 1)*runlength + i) = tmpreadB
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
CALL realloc(indataA, full_setnum*runlength)
CALL realloc(indataB, full_setnum*runlength)

ntotvals = INT(full_setnum*runlength,idp)
! precalculate mean 
IF (.NOT.readmean) THEN
   readmean = .TRUE.
   meanvalue(1) = SUM(indataA) / REAL(ntotvals,dp)
   meanvalue(2) = SUM(indataB) / REAL(ntotvals,dp)
END IF
indataA(:) = indataA(:) - meanvalue(1)
indataB(:) = indataB(:) - meanvalue(2)
corrfkt(:) = 0.0_dp

!loop over stdinput until it ends
readerr = 0
DO iset = 0, full_setnum-1
   WRITE(UNIT=ERROR_UNIT, FMT=*) "Processing Data Set Nr.         ", iset + 1
   IF (periodic) THEN
      ! Copy data to fourier array
      DO i = 1, runlength
         datainA(i-1) = CMPLX(indataA(iset*runlength + i),0.0_dp,cdp)
      END DO
      ! execute FFTW
      CALL dfftw_execute_dft(planfwdA, datainA(:), dataoutA(:))

      ! Copy data to fourier array
      DO i = 1, runlength
         datainB(i-1) = CMPLX(indataB(iset*runlength + i),0.0_dp,cdp)
      END DO
      ! execute FFTW
      CALL dfftw_execute_dft(planfwdB, datainB(:), dataoutB(:))

      ftprod(:) = ftprod(:) + dataoutA(:)*CONJG(dataoutB(:))
   ELSE
      ! Copy data to fourier array and zeropad
      DO i = 1, runlength
         datainA(i-1) = CMPLX(indataA(iset*runlength + i),0.0_dp,cdp)
      END DO
      DO i = runlength+1, 2*runlength
         datainA(i-1) = CMPLX(0.0_dp,0.0_dp,cdp)
      END DO
      ! execute FFTW
      CALL dfftw_execute_dft(planfwdA, datainA(:), dataoutA(:))

      ! Copy data to fourier array and zeropad
      DO i = 1, runlength
         datainB(i-1) = CMPLX(indataB(iset*runlength + i),0.0_dp,dp)
      END DO
      DO i = runlength+1, 2*runlength
         datainB(i-1) = CMPLX(0.0_dp,0.0_dp,cdp)
      END DO
      ! execute FFTW
      CALL dfftw_execute_dft(planfwdB, datainB(:), dataoutB(:))

      ftprod(:) = ftprod(:) + dataoutA(:)*CONJG(dataoutB(:))
   END IF
END DO


! summ up data and backtransform
datainA = ftprod(:)
CALL dfftw_execute_dft(planbwdA, datainA(:), dataoutA(:))
!End normalization
IF (periodic) THEN
   DO i = 0, maxlag
      corrfkt(i) = REAL(dataoutA(i),dp) / (REAL(ntotvals*runlength*full_setnum,dp))
   END DO
ELSE
   DO i = 0, maxlag
      corrfkt(i) = REAL(dataoutA(i),dp) / (2.0_dp*REAL(ntotvals*(runlength-i),dp))
   END DO
END IF
varianceA = DOT_PRODUCT(indataA(:), &
                        indataA(:)) & 
                       / REAL(ntotvals,dp)
varianceB = DOT_PRODUCT(indataB(:), &
                        indataB(:)) & 
                       / REAL(ntotvals,dp)

corrfkt(:) = corrfkt(:) / SQRT(varianceA*varianceB)

WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# Number of Sets = ", full_setnum
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# Runlength =      ", runlength
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# maxlag =         ", maxlag
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# mean =           ", meanvalue(1)
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "#                  ", meanvalue(2)
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# variance =       ", varianceA
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "#                  ",varianceB
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "# stddeviation =   ", SQRT(varianceA)
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "#                  ", SQRT(varianceB)
WRITE(UNIT=OUTPUT_UNIT,FMT=*) "#"

DO icorrf = 0,maxlag
   WRITE(UNIT=OUTPUT_UNIT,FMT=*) icorrf, corrfkt(icorrf)
END DO
IF ((periodic).AND.(maxlag + 1 == runlength)) THEN
   WRITE(UNIT=OUTPUT_UNIT,FMT=*) maxlag + 1, corrfkt(0)
END IF

DEALLOCATE(ftprod)
DEALLOCATE(corrfkt)
!FFT arrays
DEALLOCATE(datainA, dataoutA)
DEALLOCATE(datainB, dataoutB)
CALL dfftw_destroy_plan(planfwdA)
CALL dfftw_destroy_plan(planbwdA)
CALL dfftw_destroy_plan(planfwdB)
CALL dfftw_destroy_plan(planbwdB)
DEALLOCATE(indataA)
NULLIFY(indataA)
DEALLOCATE(indataB)
NULLIFY(indataB)
DEALLOCATE(meanvalue)
NULLIFY(meanvalue)

END PROGRAM crosscorr
