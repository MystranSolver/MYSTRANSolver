! ##################################################################################################################################
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

   MODULE SPRING_ELEMENTS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BUSH, ELAS1

   CONTAINS

      SUBROUTINE BUSH ( INT_ELEM_ID, OPT, WRITE_WARN )

! Calculates, for BUSH element

!  1) KE        = element stiffness matrix if OPT(4) = 'Y'
!  2) SE1, STE1 = element stress data recovery matrices if OPT(3) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  BE1, BE2, BUSH_DXA, BUSH_DXB, BUSH_DY, BUSH_DZ, EPROP, KE, OFFDIS_GA_GB, SE1, SE2

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF
      USE RESULT_FORMATTING, ONLY     :  WRT_REAL_TO_CHAR_VAR

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BUSH'
      CHARACTER(LEN=*) , INTENT(IN)   :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices


      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID
      INTEGER(LONG)                   :: I,J               ! DO loop indices

      REAL(DOUBLE)                    :: DUM0(6,12)        ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: DUM1(6,12)        ! Intermediate matrix in a calc
      REAL(DOUBLE)                    :: STIFF(12,12)      ! 12x12 partition from KE (complete stiff matrix for this elem)
      REAL(DOUBLE)                    :: STRE_RCV(2)       ! Two stress recovery values
      REAL(DOUBLE)                    :: STRN_RCV(2)       ! Two strain recovery values



! **********************************************************************************************************************************
! Initialize:

      BUSH_DXA = OFFDIS_GA_GB(1,1)
      BUSH_DY  = OFFDIS_GA_GB(1,2)
      BUSH_DZ  = OFFDIS_GA_GB(1,3)
      BUSH_DXB = OFFDIS_GA_GB(2,1)

! Calc stiffness matrix

      IF ((OPT(3) == 'Y') .OR. (OPT(4) == 'Y')) THEN

         DO I=1,12
            DO J=1,12
               KE(I,J) = ZERO
            ENDDO
         ENDDO

         KE( 1, 1) =  EPROP(1)                             ! NB *** new 09/10/21. Big change. KE (before offsets): just EPROP vals
         KE( 2, 2) =  EPROP(2)
         KE( 3, 3) =  EPROP(3)
         KE( 4, 4) =  EPROP(4)
         KE( 5, 5) =  EPROP(5)
         KE( 6, 6) =  EPROP(6)

         KE( 7, 7) =  KE( 1, 1)
         KE( 8, 8) =  KE( 2, 2)
         KE( 9, 9) =  KE( 3, 3)
         KE(10,10) =  KE( 4, 4)
         KE(11,11) =  KE( 5, 5)
         KE(12,12) =  KE( 6, 6)

         KE( 1, 7) = -KE( 1, 1)
         KE( 2, 8) = -KE( 2, 2)
         KE( 3, 9) = -KE( 3, 3)
         KE( 4,10) = -KE( 4, 4)
         KE( 5,11) = -KE( 5, 5)
         KE( 6,12) = -KE( 6, 6)

         DO I=1,12
            DO J=I,12
               KE(J,I) = KE(I,J)
            ENDDO
         ENDDO

      ENDIF

! Calc stress and strain recovery matrices. The stress recovery coefficients on the Bulk Data entry are intended to be
! multiplied times the element ENGINEERING (NOT NODE) forces to get stresses. In order to calc stresses in the MYSTRAN fasion
! of multiplying SEi times displacements, the SEi matrices are generated with the aid of the matrix that converts element
! nodal forces to element engineering forces. The DUM0 matrix (below), if all of the STRE_RCV numbers were 1.0 is that matrix
! (see conversion code in subr OFP3_ELFE_1D. With the STRE_RCV values in DUM0, then DUM0 times KE would be the SE matrices
! that convert displacements to stresses

      IF (OPT(3) == 'Y') THEN
                                                           ! Calc SEi stress recovery matrices
         STRE_RCV(1) = EPROP(19)                           ! --- Stress/strain recovery coefficients from PBUSH entry
         STRE_RCV(2) = EPROP(20)
         STRN_RCV(1) = EPROP(21)
         STRN_RCV(2) = EPROP(22)

         DO I=1,6
            DO J=1,12
               DUM0(I,J) = ZERO
               DUM1(I,J) = ZERO
            ENDDO
         ENDDO

         DUM0(1, 7) =  STRE_RCV(1)
         DUM0(2, 8) =  STRE_RCV(1)
         DUM0(3, 9) =  STRE_RCV(1)

         DUM0(4,10) =  STRE_RCV(2)
         DUM0(5,11) =  STRE_RCV(2)
         DUM0(6,12) =  STRE_RCV(2)

         DUM0(4, 8) =  STRE_RCV(2)*BUSH_DZ
         DUM0(4, 9) = -STRE_RCV(2)*BUSH_DY
         DUM0(5, 7) = -STRE_RCV(2)*BUSH_DZ
         DUM0(5, 9) = -STRE_RCV(2)*BUSH_DXB
         DUM0(6, 7) =  STRE_RCV(2)*BUSH_DY
         DUM0(6, 8) =  STRE_RCV(2)*BUSH_DXB

         DO I=1,12
            DO J=1,12
               STIFF(I,J) = KE(I,J)
            ENDDO
         ENDDO

         CALL MATMULT_FFF ( DUM0, STIFF, 6, 12, 12, DUM1 ) ! --- Dummy 6x12 matrix that has the SEi matrix data

         DO    I=1,3                                                                                                                                                                        !    ---    SEi
            DO    J=1,12
               SE1(I,J,1) = DUM1(I  ,J)
               SE2(I,J,1) = DUM1(I+3,J)
            ENDDO
         ENDDO

         BE1(1, 1,1) = -STRN_RCV(1)                        ! Calc BEi strain recovery matrices
         BE1(1, 7,1) =  STRN_RCV(1)

         BE1(2, 2,1) = -STRN_RCV(1)
         BE1(2, 8,1) =  STRN_RCV(1)

         BE1(3, 3,1) = -STRN_RCV(1)
         BE1(3, 9,1) =  STRN_RCV(1)

         BE2(1, 4,1) = -STRN_RCV(2)
         BE2(1,10,1) =  STRN_RCV(2)

         BE2(2, 5,1) = -STRN_RCV(2)
         BE2(2,11,1) =  STRN_RCV(2)

         BE2(3, 6,1) = -STRN_RCV(2)
         BE2(3,12,1) =  STRN_RCV(2)

      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE BUSH


      SUBROUTINE ELAS1 ( OPT, WRITE_WARN )

! Calculates 1-D scalar spring matrices in global coordinate system

!  1) SEi  = Element stress recovery matrices if OPT(3) = 'Y'
!  2) KE   = Element stiffness matrix in element coord's if OPT(4) = 'Y'

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE MODEL_STUF, ONLY            :  BGRID, ELAS_COMP, EPROP, FCONV, KE, SE1, TYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_GRID_NUM_COMPS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELAS1'
      CHARACTER(1*BYTE), INTENT(IN)   :: OPT(6)            ! 'Y'/'N' flags for whether to calc certain elem matrices
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not

      INTEGER(LONG)                   :: I1                ! The component no (1-6) at end A that this elem is connected to
      INTEGER(LONG)                   :: I2                ! The component no (1-6) at end B that this elem is connected to
      INTEGER(LONG)                   :: NUM_COMPS_GRID_1  ! No. displ components (1 for SPOINT, 6 for actual grid) for 1st grid


      REAL(DOUBLE)                    :: K                 ! Spring stiffness


! **********************************************************************************************************************************
! Set element property and material constants

      K        = EPROP(1)
      IF ((TYPE(1:5) == 'ELAS1') .OR. (TYPE(1:5) == 'ELAS2') .OR. (TYPE(1:5) == 'ELAS3')) THEN
         FCONV(1) = EPROP(3)
      ELSE
         FCONV(1) = 0.0
      ENDIF
      I1       = ELAS_COMP(1)
      CALL GET_GRID_NUM_COMPS ( BGRID(1), NUM_COMPS_GRID_1, SUBR_NAME )
      I2       = NUM_COMPS_GRID_1 + ELAS_COMP(2)

! **********************************************************************************************************************************
! Calculate the element stiffness matrix in global coordinates.

      IF (OPT(4) == 'Y') THEN
         KE(I1,I1) =  K
         KE(I1,I2) = -KE(I1,I1)
         KE(I2,I1) = -KE(I1,I1)
         KE(I2,I2) =  KE(I1,I1)
      ENDIF

! **********************************************************************************************************************************
! Calculate SE1 matrix for force recovery.

      IF (OPT(3) == 'Y') THEN
         SE1(1,I1,1) =  K
         SE1(1,I2,1) = -K
      ENDIF



      RETURN

! *********************************************************************************************************************************

      END SUBROUTINE ELAS1

   END MODULE SPRING_ELEMENTS
