!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! NAME: 
!!   histogram3d
!! DESCRIPTION:
!!   Program to calculate a 3d histogram from a list with three columns of numbers .
!!   It checks for the smallest value and start the intervals from there
!!   output is a cubefile
!!   negative or zero norm prevents normalizing
!! COMPILATION:
!!   make
!! USAGE:
!!   cat data.dat | ./histogram.x <intervallsize> <normalizationfactor> > histogram.cube
!!   data.dat needs to have the format:
!!   x-coord y-coord z-coord
!!   space seperated coordinates. no atom symbol
!! AUTHOR:
!!   fuhl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
PROGRAM histogram3d
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
   
   INTEGER, PARAMETER :: ndim = 3
   INTEGER :: readerr
   REAL(KIND=dp) :: intsize, norm, integral
   REAL(KIND=dp), DIMENSION(ndim) :: tmpread, smallval, bigval
   REAL(KIND=dp), DIMENSION(ndim) :: tmpmax, tmpmin
   INTEGER, DIMENSION(ndim) :: histsize = 0
   INTEGER, DIMENSION(ndim) :: histstart = 0
   INTEGER, DIMENSION(ndim) :: histend = 0
   INTEGER(KIND=idp) :: ndata = 0
   INTEGER(KIND=idp), DIMENSION(:,:,:), POINTER :: hist
   INTEGER, DIMENSION(ndim) :: idx,k
   INTEGER :: ix, iy, iz
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
      READ(UNIT=INPUT_UNIT, FMT=*, IOSTAT=readerr) tmpread(3), tmpread(2), tmpread(1)
      IF (readerr /= 0) EXIT
      ! for first iteration choose 
      IF (.NOT.ASSOCIATED(hist)) THEN
         bigval(:) = tmpread(:)
         smallval(:) = tmpread(:)
         histstart(:) = FLOOR(tmpread(:)/intsize)
         histend(:) = FLOOR(tmpread(:)/intsize)+1
      END IF
      k(:) = FLOOR(tmpread(:)/intsize)
      tmpmin(:) = intsize * REAL(k(:),dp)
      tmpmax(:) = intsize * REAL(k(:),dp) + intsize
      IF ((tmpmax(1) > bigval(1)).OR.tmpmin(1) < smallval(1).OR. &
          (tmpmax(2) > bigval(2)).OR.tmpmin(2) < smallval(2).OR. &
          (tmpmax(3) > bigval(3)).OR.tmpmin(3) < smallval(3)) THEN
         bigval(:) = MAX(bigval(:), tmpmax(:))
         smallval(:) = MIN(smallval(:), tmpmin(:))
         histstart(:) = MIN(k(:), histstart(:))
         histend(:)   = MAX(k(:)+1, histend(:))
         histsize(:)  = histend(:) - histstart(:) + 1
         CALL realloc(hist, histsize, histstart)
      END IF

      ! calculate histogram index
      idx(:) = k(:) + 1
      hist(idx(1),idx(2),idx(3)) = hist(idx(1),idx(2),idx(3)) + 1
      ndata = ndata + 1
   END DO

   IF (ndata <= 0) STOP

   integral = 0.0_dp
   IF (norm > 0.0_dp) THEN
      DO ix = histstart(3), histend(3)
         DO iy = histstart(2), histend(2)
            DO iz = histstart(1), histend(1)
               integral = integral + REAL(hist(iz,iy,ix),dp)
            END DO
         END DO
      END DO
      integral = integral * intsize
      norm = integral / norm
   ELSE 
      norm = 1.0_dp
   END IF

   !output as cube-file
   !The first two lines of the header are comments, 
   !they are generally ignored by parsing packages or used as two default labels. 
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Gaussian cube file by histogram3d.x (Felix Uhl)"
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) "Created from ", ndata, " datapoints."

   !The third line has the number of atoms included in the file 
   !followed by the position of the origin of the volumetric data.  
   !Atoms will be places at the box edges
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) 8, smallval(3)-intsize, smallval(2)-intsize, smallval(1)-intsize
   
   !The next three lines give the number of voxels along each axis (x, y, z) 
   !followed by the axis vector.
   !If the sign of the number of voxels in a dimension is positive then the units are Bohr, 
   !if negative then Angstroms. 
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) histsize(3), intsize, 0.0_dp , 0.0_dp
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) histsize(2), 0.0_dp , intsize, 0.0_dp
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) histsize(1), 0.0_dp , 0.0_dp , intsize

   !The last section in the header is one line for each atom consisting of 5 numbers, 
   !the first is the atom number, 
   !second (?), 
   !the last three are the x,y,z coordinates of the atom center. 
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) 1, 0.0_dp , smallval(1), smallval(2), smallval(3)
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) 1, 0.0_dp , smallval(1), smallval(2), bigval(3)
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) 1, 0.0_dp , smallval(1), bigval(2),   smallval(3)
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) 1, 0.0_dp , smallval(1), bigval(2),   bigval(3)
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) 1, 0.0_dp , bigval(1),   smallval(2), smallval(3)
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) 1, 0.0_dp , bigval(1),   smallval(2), bigval(3)
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) 1, 0.0_dp , bigval(1),   bigval(2),   smallval(3)
   WRITE(UNIT=OUTPUT_UNIT, FMT=*) 1, 0.0_dp , bigval(1),   bigval(2),   bigval(3)
   
   !The volumetric data is straightforward, one floating point number for each volumetric element. 
   !The original Gaussian format arranged the values in the format shown below in the examplehite space separated format. Traditionally the grid is arranged with the x axis as the outer loop and the z axis as the inner loop
   DO ix = histstart(3), histend(3)
      DO iy = histstart(2), histend(2)
         DO iz = histstart(1), histend(1)
            WRITE(UNIT=OUTPUT_UNIT, FMT='(E15.7,1X)', ADVANCE='no') &
               (REAL(hist(iz,iy,ix),dp)/norm)
         END DO
         WRITE(UNIT=OUTPUT_UNIT, FMT='(A)', ADVANCE='no') NEW_LINE('A')
      END DO
   END DO

   CALL dealloc(hist)
   
END PROGRAM histogram3d
