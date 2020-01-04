!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NAME: 
!!   histogram
!! DESCRIPTION:
!!   Program to calculate a histogram from a list of numbers.
!!   It checks for the smallest value and start the intervals from there
!!   negative or zero norm prevents normalizing
!! COMPILATION:
!!   make
!! USAGE:
!!   cat data.dat | ./histogram.x <intervallsize> <normalizationfactor> > histogram.out
!! AUTHOR:
!!   Felix Uhl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM histogram
   USE, INTRINSIC :: ISO_FORTRAN_ENV
   USE kinds,  ONLY : dp, &
                      idp
   USE re_de_alloc, ONLY : realloc, &
                           dealloc
   USE cmd_opt_parser, ONLY: cmd_opt_type, &
                             get_cmd_arg, &
                             show_help_msg, &
                             end_opt_parsing

   IMPLICIT NONE

   TYPE(cmd_opt_type), DIMENSION(:), POINTER :: cmd_opt
   
   INTEGER :: readerr
   REAL(KIND=dp) :: intsize, norm
   REAL(KIND=dp) :: tmpread, smallval, bigval, integral
   REAL(KIND=dp) :: tmpmax, tmpmin
   INTEGER :: histsize = 0
   INTEGER :: histstart = 0
   INTEGER :: histend = 0
   INTEGER(KIND=idp) :: ndata = 0
   INTEGER(KIND=idp), DIMENSION(:), POINTER :: hist
   INTEGER :: i,j, k

   LOGICAL :: opt_present1, opt_present2, show_help

   NULLIFY(hist)
   NULLIFY(cmd_opt)
   show_help = .FALSE.

   CALL get_cmd_arg(1, "Interval size for the histogram", cmd_opt, intsize, opt_present1)
   IF (.NOT.opt_present1) show_help = .TRUE.
   CALL get_cmd_arg(2, "Norm for the histrogram", cmd_opt, norm, opt_present2)
   IF (.NOT.opt_present1) show_help = .TRUE.

   IF (opt_present1) THEN
      IF (ABS(intsize) < EPSILON(1.0_dp)) THEN
         WRITE(UNIT=ERROR_UNIT, FMT=*) "Interval size needs to be positive."
         show_help = .TRUE.
      END IF
   END IF
   IF (opt_present2) THEN
      IF (ABS(norm) < EPSILON(1.0_dp)) THEN
         WRITE(UNIT=ERROR_UNIT, FMT=*) "Norm needs to be positive."
         show_help = .TRUE.
      END IF
   END IF

   IF (show_help) THEN
      CALL show_help_msg(cmd_opt)
      STOP
   END IF

   smallval = 0.0_dp
   bigval = smallval

   readerr = 0
   DO WHILE (readerr == 0)
      READ(UNIT=INPUT_UNIT, FMT=*, IOSTAT=readerr) tmpread
      IF (readerr /= 0) EXIT
      ! for first iteration choose 
      IF (.NOT.ASSOCIATED(hist)) THEN
         bigval = tmpread
         smallval = tmpread
         histstart = FLOOR(tmpread/intsize)
         histend   = FLOOR(tmpread/intsize)+1
      END IF
      k = FLOOR(tmpread/intsize)
      tmpmin = intsize * REAL(k,dp)
      tmpmax = intsize * REAL(k,dp) + intsize
      IF ((tmpmax > bigval).OR.(tmpmin < smallval)) THEN
         bigval = MAX(bigval, tmpmax)
         smallval = MIN(smallval, tmpmin)
         histstart = MIN(k, histstart)
         histend   = MAX(k+1, histend)
         histsize  = histend - histstart + 1
         CALL realloc(hist, histsize, histstart)
      END IF

      ! calculate histogram index
      j = k + 1
      hist(j) = hist(j) + 1
      ndata = ndata + 1
   END DO

   IF (ndata <= 0) STOP

   integral = 0.0_dp
   IF (norm > 0.0_dp) THEN
      DO i = histstart, histend
         integral = integral + REAL(hist(i),dp)
      END DO
      integral = integral * intsize
      norm = integral / norm
   ELSE 
      norm = 1.0_dp
   END IF

   WRITE(UNIT=OUTPUT_UNIT, FMT=*) "# Computet from ", ndata, " samples"
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) smallval - 0.5_dp*intsize, 0.0_dp
   DO i = histstart, histend
      WRITE(UNIT=OUTPUT_UNIT, FMT=*) (REAL(i,dp)-0.5_dp) * intsize, REAL(hist(i),dp)/norm
   END DO
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) (REAL(histend,dp)+0.5_dp) * intsize, 0.0_dp
   
   CALL dealloc(hist)
   
END PROGRAM histogram
