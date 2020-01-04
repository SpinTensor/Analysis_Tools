!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NAME: 
!!   histogram2d
!! DESCRIPTION:
!!   Program to calculate a 2d histogram from a list with two columns of numbers .
!!   It checks for the smallest value and start the intervals from there
!!   negative or zero norm prevents normalizing
!! COMPILATION:
!!   make
!! USAGE:
!!   cat data.dat | ./histogram.x <intervallsize> <normalizationfactor> > histogram.out
!! AUTHOR:
!!   Felix Uhl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM histogram2d
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
   
   INTEGER, PARAMETER :: ndim = 2
   INTEGER :: readerr
   REAL(KIND=dp) :: intsize, norm, integral
   REAL(KIND=dp), DIMENSION(ndim) :: tmpread, smallval, bigval
   REAL(KIND=dp), DIMENSION(ndim) :: tmpmax, tmpmin
   INTEGER, DIMENSION(ndim) :: histsize = 0
   INTEGER, DIMENSION(ndim) :: histstart = 0
   INTEGER, DIMENSION(ndim) :: histend = 0
   INTEGER(KIND=idp), DIMENSION(:,:), POINTER :: hist
   INTEGER(KIND=idp) :: ndata = 0
   INTEGER, DIMENSION(ndim) :: j,k
   INTEGER :: ix, iy

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

   smallval(:) = 0.0_dp
   bigval(:) = smallval(:)

   readerr = 0
   DO WHILE (readerr == 0)
      READ(UNIT=INPUT_UNIT, FMT=*, IOSTAT=readerr) tmpread(:)
      IF (readerr /= 0) EXIT
      ! for first iteration choose 
      IF (.NOT.ASSOCIATED(hist)) THEN
         bigval(:) = tmpread(:)
         smallval(:) = tmpread(:)
         histstart(:) = FLOOR(tmpread(:)/intsize)
         histend(:)   = FLOOR(tmpread(:)/intsize)+1
      END IF
      k(:) = FLOOR(tmpread(:)/intsize)
      tmpmin(:) = intsize * REAL(k(:),dp)
      tmpmax(:) = intsize * REAL(k(:),dp) + intsize
      IF ((tmpmax(1) > bigval(1)).OR.tmpmin(1) < smallval(1).OR. &
          (tmpmax(2) > bigval(2)).OR.tmpmin(2) < smallval(2)) THEN
         bigval(:) = MAX(bigval(:), tmpmax(:))
         smallval(:) = MIN(smallval(:), tmpmin(:))
         histstart(:) = MIN(k(:), histstart(:))
         histend(:)   = MAX(k(:)+1, histend(:))
         histsize(:)  = histend(:) - histstart(:) + 1
         CALL realloc(hist, histsize, histstart)
      END IF
      ! calculate histogram index
      j(:) = k(:) + 1
      hist(j(1),j(2)) = hist(j(1),j(2)) + 1
      ndata = ndata + 1
   END DO

   IF (ndata <= 0) STOP

   integral = 0.0_dp
   IF (norm > 0.0_dp) THEN
      DO iy = histstart(2), histend(2)
         DO ix = histstart(1), histend(1)
            integral = integral + REAL(hist(ix,iy),dp)
         END DO
      END DO
      integral = integral * intsize
      norm = integral / norm
   ELSE 
      norm = 1.0_dp
   END IF

   DO iy = histstart(2), histend(2)
      DO ix = histstart(1), histend(1)
         WRITE(UNIT=OUTPUT_UNIT, FMT=*) (REAL(ix,dp)-0.5_dp) * intsize, &
                                        (REAL(iy,dp)-0.5_dp) * intsize, &
                                        REAL(hist(ix,iy),dp)/norm
      END DO
      WRITE(UNIT=OUTPUT_UNIT, FMT=*) ""
   END DO

   CALL dealloc(hist)
   
END PROGRAM histogram2d
