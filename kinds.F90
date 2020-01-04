!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Description:
!!    contains kinds for standardized precision
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE kinds

   USE, INTRINSIC :: ISO_FORTRAN_ENV

   IMPLICIT NONE

   PRIVATE

#ifdef __HAS_QP
   PUBLIC :: sp, dp, qp
   PUBLIC :: csp, cdp, cqp
#else
   PUBLIC :: csp, cdp
   PUBLIC :: sp, dp
#endif
#ifdef __HAS_IQP
   PUBLIC :: isp, idp, iqp
#else
   PUBLIC :: isp, idp
#endif
   PUBLIC :: default_path_length
   PUBLIC :: default_string_length
   PUBLIC :: default_symbol_length
   PUBLIC :: default_path_ioformat
   PUBLIC :: default_string_ioformat
   PUBLIC :: default_symbol_ioformat

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Real kinds
   !! Variables:
   !!    sp: single precision
   !!    dp: double precision
   !!    qp: quadruple precision
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTEGER, PARAMETER :: sp = KIND(1.0)
   INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(2*PRECISION(1.0_sp))
#ifdef __HAS_QP
   INTEGER, PARAMETER :: qp = SELECTED_REAL_KIND(2*PRECISION(1.0_dp))
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Integer kinds
   !! Variables:
   !!    isp: single precision
   !!    idp: double precision
   !!    iqp: quadruple precision
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTEGER, PARAMETER :: isp = KIND(1)
   INTEGER, PARAMETER :: idp = 2*isp
#ifdef __HAS_IQP
   INTEGER, PARAMETER :: iqp = 2*idp
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    Complex kinds
   !! Variables:
   !!    csp: single precision
   !!    cdp: double precision
   !!    cqp: quadruple precision
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTEGER, PARAMETER :: csp = KIND(1.0)
   INTEGER, PARAMETER :: cdp = SELECTED_REAL_KIND(2*PRECISION(1.0_sp))
#ifdef __HAS_QP
   INTEGER, PARAMETER :: cqp = SELECTED_REAL_KIND(2*PRECISION(1.0_dp))
#endif

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    universal string lengths
   !! Variables:
   !!    default_path_length: length of a file path
   !!    default_string_length: length of a string
   !!    default_symbol_length: length of an atom symbol
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   INTEGER, PARAMETER :: default_path_length = 1024
   INTEGER, PARAMETER :: default_string_length = 1024
   INTEGER, PARAMETER :: default_symbol_length = 3

   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   !! Description:
   !!    universal formats for defined string length
   !! Variables:
   !!    default_path_length: length of a file path
   !!    default_string_length: length of a string
   !!    default_symbol_length: length of an atom symbol
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

   CHARACTER(LEN=5), PARAMETER :: default_path_ioformat = "A1024"
   CHARACTER(LEN=5), PARAMETER :: default_string_ioformat = "A1024"
   CHARACTER(LEN=2), PARAMETER :: default_symbol_ioformat = "A3"

END MODULE kinds
