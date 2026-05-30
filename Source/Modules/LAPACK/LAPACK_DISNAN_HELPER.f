! ##################################################################################################################################

      MODULE LAPACK_DISNAN_HELPER

      CONTAINS

! ##################################################################################################################################

      LOGICAL FUNCTION DLAISNAN_HELPER(DIN1,DIN2)

      DOUBLE PRECISION DIN1,DIN2

      DLAISNAN_HELPER = (DIN1.NE.DIN2)
      RETURN
      END FUNCTION DLAISNAN_HELPER

! ##################################################################################################################################

      LOGICAL FUNCTION DISNAN_HELPER(DIN)

      DOUBLE PRECISION DIN

      DISNAN_HELPER = DLAISNAN_HELPER(DIN,DIN)
      RETURN
      END FUNCTION DISNAN_HELPER

      END MODULE LAPACK_DISNAN_HELPER
