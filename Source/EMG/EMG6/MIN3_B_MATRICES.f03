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

   MODULE MIN3_B_MATRICES

   USE QUADRILATERAL_B_MATRICES, ONLY :  BCHECK_2D

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BBMIN3, BSMIN3

   CONTAINS

      SUBROUTINE BBMIN3 ( A, B, AREA, MESSAG, WRT_BUG_THIS_TIME, BB )

! Calculate BB bending strain/displacement matrix for MIN3 triangle. Called by subr TPLT2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, WRT_BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_BMAT_BIT, ELDT_BUG_BCHK_BIT, MIN4T_QUAD4_TRIA_NO
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO
      USE MODEL_STUF, ONLY            :  EID, TYPE, XTB, XTL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BBMIN3'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to print out if BCHECK is run
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: JJ                ! A computed index into array BB
      INTEGER(LONG)                   :: ID(9)             ! An input to subr BCHECK, called herein if
      INTEGER(LONG), PARAMETER        :: NR        = 3     ! An input to subr BCHECK, called herein if
      INTEGER(LONG), PARAMETER        :: NC        = 9     ! An input to subr BCHECK, called herein if


      REAL(DOUBLE) , INTENT(IN)       :: A(3)              ! Diffs in x coords of elem
      REAL(DOUBLE) , INTENT(IN)       :: B(3)              ! Diffs in y coords of elem
      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Elem area
      REAL(DOUBLE) , INTENT(OUT)      :: BB(3,9)           ! Output strain-displ matrix for this elem
      REAL(DOUBLE)                    :: BW(3,14)          ! Output from subr BCHECK (matrix of 3 elem strains for 14 various elem
!                                                            rigid body motions/constant strain distortions)

! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,3
         DO J=1,9
            BB(I,J) = ZERO
         ENDDO
      ENDDO

! Calc outputs

      JJ = 0
      DO J=4,6
         JJ = JJ + 1
         BB(1,J) = ZERO
         BB(2,J) = -A(JJ)/(TWO*AREA)
         BB(3,J) = -B(JJ)/(TWO*AREA)
      ENDDO

      JJ = 0
      DO J=7,9
         JJ = JJ + 1
         BB(1,J) = B(JJ)/(TWO*AREA)
         BB(2,J) = ZERO
         BB(3,J) = A(JJ)/(TWO*AREA)
      ENDDO

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(8) > 0)) THEN

         WRITE(BUG,1101) ELDT_BUG_BMAT_BIT, TYPE, EID, MIN4T_QUAD4_TRIA_NO
         WRITE(BUG,8901) SUBR_NAME
         DO I=1,3
            WRITE(BUG,8902) I,(BB(I,J),J=1,9)
            WRITE(BUG,*)
         ENDDO
         WRITE(BUG,*)

      ENDIF

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! As of 03/07/2020 there is some error in the calculation of the check on strain displ matrices for constant strain modes
! of displacement. Therefore, this check is temporarily suspended

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(9) > 0)) THEN
        IF (DEBUG(202) > 0) THEN
         ID(1) =  3
         ID(2) =  9
         ID(3) = 15
         ID(4) =  4
         ID(5) = 10
         ID(6) = 16
         ID(7) =  5
         ID(8) = 11
         ID(9) = 17

         WRITE(BUG,1101) ELDT_BUG_BCHK_BIT, TYPE, EID, MIN4T_QUAD4_TRIA_NO
         WRITE(BUG,9100)
         WRITE(BUG,9101) MESSAG
         WRITE(BUG,9102)
         WRITE(BUG,9103)
         WRITE(BUG,9104)
         CALL BCHECK_2D ( BB, 'B', ID, NR, NC, 3, XTL, XTB, BW )
        ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1101 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------',/,                                                                                                &
             ' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8, ' (triangle number ',I1,' of 4 for the MIN4T QUAD4)'/, &
             ' ==============================================================',/)

 8901 FORMAT(' Strain-displacement matrix BB for bending portion of element in subr ',A,/)

 8902 FORMAT(' Row ',I2,/,9(1ES14.6))

 9100 FORMAT(26X,'Check on strain-displacement matrix BB for bending portion of the element in subr BCHECK'/)

 9101 FORMAT(63X,'S T R A I N S'/,61X,A)

 9102 FORMAT(53X,'Cxx            Cyy            Cxy')

 9103 FORMAT(7X,'Element displacements consistent with:')

 9104 FORMAT(7X,'---------------------------------------')

! **********************************************************************************************************************************

      END SUBROUTINE BBMIN3


      SUBROUTINE BSMIN3 ( XI, A, B, AREA, MESSAG, WRT_BUG_THIS_TIME, BS )

! Calculate BS shear strain/displacement matrix for MIN3 triangle. Called by subr TPLT2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, WRT_BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_BMAT_BIT, ELDT_BUG_BCHK_BIT, MIN4T_QUAD4_TRIA_NO
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO
      USE MODEL_STUF, ONLY            :  EID, TYPE, XTB, XTL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BSMIN3'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to print out if BCHECK is run
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: ID(9)             ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NR        = 2     ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NC        = 9     ! An input to subr BCHECK, called herein


      REAL(DOUBLE) , INTENT(IN)       :: A(3)              ! Vector of x coord differences
      REAL(DOUBLE) , INTENT(IN)       :: AREA              ! Elem area
      REAL(DOUBLE) , INTENT(IN)       :: B(3)              ! Vector of y coord differences
      REAL(DOUBLE) , INTENT(IN)       :: XI(3)             ! Vector of area coordinates
      REAL(DOUBLE) , INTENT(OUT)      :: BS(2,9)           ! Output strain-displ matrix for this elem
      REAL(DOUBLE)                    :: A4                ! Constant for this elem
      REAL(DOUBLE)                    :: BW(2,14)          ! Output from subr BCHECK (matrix of 2 elem strains for 14 various elem
                                                             ! rigid body motions/constant strain distortions)

! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,2
         DO J=1,9
            BS(I,J) = ZERO
         ENDDO
      ENDDO

! Calc outputs

      A4 = 4*AREA

      DO J=1,3
         BS(1,J) = B(J)/(TWO*AREA)
         BS(2,J) = A(J)/(TWO*AREA)
      ENDDO

      BS(1,4) = -(B(1)*( XI(2)*B(3) - XI(3)*B(2))/A4)
      BS(1,5) = -(B(2)*(-XI(1)*B(3) + XI(3)*B(1))/A4)
      BS(1,6) = -(B(3)*( XI(1)*B(2) - XI(2)*B(1))/A4)

      BS(1,7) =  (XI(1)*(A4 - B(2)*A(3) + B(3)*A(2)) + B(1)*(-XI(2)*A(3) + XI(3)*A(2)))/A4
      BS(1,8) =  (XI(2)*(A4 + B(1)*A(3) - B(3)*A(1)) + B(2)*( XI(1)*A(3) - XI(3)*A(1)))/A4
      BS(1,9) =  (XI(3)*(A4 + B(2)*A(1) - B(1)*A(2)) + B(3)*(-XI(1)*A(2) + XI(2)*A(1)))/A4

      BS(2,4) = -(XI(1)*(A4 + A(2)*B(3) - A(3)*B(2)) + A(1)*( XI(2)*B(3) - XI(3)*B(2)))/A4
      BS(2,5) = -(XI(2)*(A4 - A(1)*B(3) + A(3)*B(1)) + A(2)*(-XI(1)*B(3) + XI(3)*B(1)))/A4
      BS(2,6) = -(XI(3)*(A4 - A(2)*B(1) + A(1)*B(2)) + A(3)*( XI(1)*B(2) - XI(2)*B(1)))/A4

      BS(2,7) = A(1)*(-XI(2)*A(3) + XI(3)*A(2))/A4
      BS(2,8) = A(2)*( XI(1)*A(3) - XI(3)*A(1))/A4
      BS(2,9) = A(3)*(-XI(1)*A(2) + XI(2)*A(1))/A4

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(8) > 0)) THEN

         WRITE(BUG,1101) ELDT_BUG_BMAT_BIT, TYPE, EID, MIN4T_QUAD4_TRIA_NO
         WRITE(BUG,8901) SUBR_NAME
         DO I=1,2
            WRITE(BUG,8902) I,(BS(I,J),J=1,9)
            WRITE(BUG,*)
         ENDDO
         WRITE(BUG,*)

      ENDIF

! xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
! As of 03/07/2020 there is some error in the calculation of the check on strain displ matrices for constant strain modes
! of displacement. Therefore, this check is temporarily suspended

     IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(9) > 0)) THEN
        IF (DEBUG(202) > 0) THEN
         ID(1) =  3
         ID(2) =  9
         ID(3) = 15
         ID(4) =  4
         ID(5) = 10
         ID(6) = 16
         ID(7) =  5
         ID(8) = 11
         ID(9) = 17

         WRITE(BUG,1101) ELDT_BUG_BCHK_BIT, TYPE, EID, MIN4T_QUAD4_TRIA_NO
         WRITE(BUG,9100)
         WRITE(BUG,9101) MESSAG
         WRITE(BUG,9102)
         WRITE(BUG,9103)
         WRITE(BUG,9104)
         CALL BCHECK_2D ( BS, 'S', ID, NR, NC, 3, XTL, XTB, BW )
        ENDIF
     ENDIF



      RETURN

! **********************************************************************************************************************************
 1101 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------',/,                                                                                                &
             ' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8, ' (triangle number ',I1,' of 4 for the QUAD4 elem)'/,  &
             ' ==============================================================',/)

 8901 FORMAT(' Strain-displacement matrix BS for transverse shear portion of element in subr ',A,/)

 8902 FORMAT(' Row ',I2,/,9(1ES14.6))

 9100 FORMAT(14X,'Check on strain-displacement matrix BS for transverse shear portion of the element in subr BCHECK'/)

 9101 FORMAT(56X,'S T R A I N S'/,50X,A)

 9102 FORMAT(54X,'Gxz            Gyz')

 9103 FORMAT(7X,'Element displacements consistent with:')

 9104 FORMAT(7X,'---------------------------------------')

90003 FORMAT(1X,3(1ES19.11))

! **********************************************************************************************************************************

      END SUBROUTINE BSMIN3

   END MODULE MIN3_B_MATRICES
