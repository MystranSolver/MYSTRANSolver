! #################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)

! Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
! associated documentation files (the "Software"), to deal in the Software without restriction, including
! without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
! copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to
! the following conditions:

! The above copyright notice and this permission notice shall be included in all copies or substantial
! portions of the Software and documentation.

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
! OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
! FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
! AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.

   MODULE MASS_DOF_EXPANSION

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: EXPAND_MASS_DOFS

   CONTAINS
      SUBROUTINE EXPAND_MASS_DOFS ( M_1DOF )

! Copies per-node mass values to all 3 of each node's translational DOFs in the element mass matrix ME.

      USE PENTIUM_II_KIND, ONLY       :  DOUBLE, LONG
      USE MODEL_STUF, ONLY            :  ELGP, ME
      USE CONSTANTS_1, ONLY           :  ZERO
      USE SCONTR, ONLY                :  SOL_NAME

      IMPLICIT NONE

      INTEGER(LONG)                   :: I,J,L
      INTEGER(LONG)                   :: KI, KJ            ! For converting grid point number to element DOF number

      REAL(DOUBLE), INTENT(INOUT)     :: M_1DOF(ELGP,ELGP)
      REAL(DOUBLE)                    :: ROW_SUM

! **********************************************************************************************************************************

      ! SOL 101, 104, 105, 31(?) must use lumped mass matrix because GRAV and RFORCE assume 6x6 block diagonal structure for MGG.
      ! SOL 103 should use consistent mass matrix because quadratic elements won't solve with this row sum mass lumping. They
      !         will solve with equal-mass-per-node lumping but with slightly worse results.
      IF ( SOL_NAME(1:5) /= 'MODES') THEN

         ! Row sum to convert consistent to lumped mass matrix
         DO I=1,ELGP
            ROW_SUM = ZERO
            DO J=1,ELGP
               ROW_SUM = ROW_SUM + M_1DOF(I,J)
               M_1DOF(I,J) = ZERO
            ENDDO
            M_1DOF(I,I) = ROW_SUM
         ENDDO

      ENDIF

      ME(:,:) = ZERO

                                                           ! Copy to all 3 translational DOFs of each node.
      DO I=1,ELGP
         DO J=1,ELGP
            KI = (I-1) * 6
            KJ = (J-1) * 6
            DO L=1,3
               ME(KI + L, KJ + L) = M_1DOF(I,J)
            ENDDO
         ENDDO
      ENDDO

      RETURN


! **********************************************************************************************************************************

      END SUBROUTINE EXPAND_MASS_DOFS

   END MODULE MASS_DOF_EXPANSION
