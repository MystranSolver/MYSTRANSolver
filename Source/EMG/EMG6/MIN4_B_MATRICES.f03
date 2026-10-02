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

   MODULE MIN4_B_MATRICES

   USE QUADRILATERAL_B_MATRICES, ONLY :  BCHECK_2D

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BBMIN4, BSMIN4

   CONTAINS

      SUBROUTINE BBMIN4 ( DPSHX, IGAUS, JGAUS, MESSAG, WRT_BUG_THIS_TIME, BB )

! Calculate BB bending strain/displacement matrix for MIN4 quad. Called by subr QPLT2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, WRT_BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_BMAT_BIT, ELDT_BUG_BCHK_BIT
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  EID, TYPE, XEB, XEL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BBMIN4'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Messag to print out if BCHECK is run
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG), INTENT(IN)       :: IGAUS             ! I index of Gaus point (needed for some optional output)
      INTEGER(LONG), INTENT(IN)       :: JGAUS             ! J index of Gaus point (needed for some optional output)
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: JJ                ! A computed index into array BB
      INTEGER(LONG)                   :: ID(8)             ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NR        = 3     ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NC        = 8     ! An input to subr BCHECK, called herein


      REAL(DOUBLE) , INTENT(IN)       :: DPSHX(2,4)        ! Derivatives of the 4 node bilinear isopar interps wrt elem x and y
      REAL(DOUBLE) , INTENT(OUT)      :: BB(3,8)           ! Output strain-displ matrix for this elem
      REAL(DOUBLE)                    :: BW(3,14)          ! Output from subr BCHECK (matrix of 3 elem strains for 14 various elem
                                                           ! rigid body motions/constant strain distortions)

      REAL(DOUBLE)                    :: XB(4,3)           ! First 4 rows of XEB
      REAL(DOUBLE)                    :: XL(4,3)           ! First 4 rows of XEL


! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,3
         DO J=1,8
            BB(I,J) = ZERO
         ENDDO
      ENDDO

! Calc outputs

      JJ = 0

      DO J=1,4

         JJ = JJ + 1
         BB(1,JJ) = ZERO
         BB(2,JJ) = -DPSHX(2,J)
         BB(3,JJ) = -DPSHX(1,J)

         JJ = JJ + 1
         BB(1,JJ) = DPSHX(1,J)
         BB(2,JJ) = ZERO
         BB(3,JJ) = DPSHX(2,J)

      ENDDO

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(8) > 0)) THEN

         WRITE(BUG,1101) ELDT_BUG_BMAT_BIT, TYPE, EID
         WRITE(BUG,8901) IGAUS, JGAUS, SUBR_NAME
         DO I=1,3
            WRITE(BUG,8902) I,(BB(I,J),J=1,8)
            WRITE(BUG,*)
         ENDDO
         WRITE(BUG,*)

      ENDIF

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(9) > 0)) THEN
        IF (DEBUG(202) > 0) THEN
         ID( 1) =  4
         ID( 2) =  5
         ID( 3) = 10
         ID( 4) = 11
         ID( 5) = 16
         ID( 6) = 17
         ID( 7) = 22
         ID( 8) = 23

         DO I=1,4
            DO J=1,3
               XB(I,J) = XEB(I,J)
               XL(I,J) = XEL(I,J)
            ENDDO
         ENDDO

         WRITE(BUG,1101) ELDT_BUG_BCHK_BIT, TYPE, EID
         WRITE(BUG,9100)
         WRITE(BUG,9101) MESSAG, IGAUS, JGAUS
         WRITE(BUG,9102)
         WRITE(BUG,9103)
         WRITE(BUG,9104)
         CALL BCHECK_2D ( BB, 'B', ID, NR, NC, 4, XB, XL, BW )
        ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1101 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------',/,                                                                                                &
             ' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8,/,                                                      &
             ' ==============================================================',/)

 8901 FORMAT(' Strain-displacement matrix BB for Gauss point: I = ',I3,', J = ',I3,' for bending portion of element in subr '      &
             ,A,/)

 8902 FORMAT(' Row ',I2,/,8(1ES14.6))

 9100 FORMAT('                          Check on strain-displacement matrix BB for bending portion of the element in subr BCHECK'/)

 9101 FORMAT('                                                               S T R A I N S'/,                                      &
             '                                           (',A,' for Gauss point: I = ',I3,', J = ',I3,')')

 9102 FORMAT('                                                     Cxx            Cyy            Cxy')

 9103 FORMAT(1X,'      Element displacements consistent with:')

 9104 FORMAT(1X,'      ---------------------------------------')

! **********************************************************************************************************************************

      END SUBROUTINE BBMIN4


      SUBROUTINE BSMIN4 ( PSH, DPSHX, DNXSHX, DNYSHX, IGAUS, JGAUS, MESSAG, WRT_BUG_THIS_TIME, BS )

! Calculate BS shear starin/displacement matrix for MIN4 quad. Called by subr QPLT2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, WRT_BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_BMAT_BIT, ELDT_BUG_BCHK_BIT
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  EID, TYPE, XEB, XEL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BSMIN4'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Messag to print out if BCHECK is run
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG), INTENT(IN)       :: IGAUS             ! I index of Gaus point (needed for some optional output)
      INTEGER(LONG), INTENT(IN)       :: JGAUS             ! J index of Gaus point (needed for some optional output)
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: JJ                ! A computed index into array BS
      INTEGER(LONG)                   :: ID(12)            ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NR        = 2     ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NC        = 12    ! An input to subr BCHECK, called herein


      REAL(DOUBLE) , INTENT(IN)       :: PSH(4)            ! 4 node bilinear isopar interp functions (used for bending)
      REAL(DOUBLE) , INTENT(IN)       :: DPSHX(2,4)        ! Derivatives of PSH shape functions wrt x and y
      REAL(DOUBLE) , INTENT(IN)       :: DNXSHX(2,4)       ! Derivatives of constrained interpolations NX wrt x, y
      REAL(DOUBLE) , INTENT(IN)       :: DNYSHX(2,4)       ! Derivatives of constrained interpolations NY wrt x, y
      REAL(DOUBLE) , INTENT(OUT)      :: BS(2,12)          ! Output strain-displ matrix for this elem
      REAL(DOUBLE)                    :: BW(2,14)          ! Output from subr BCHECK (matrix of 2 elem strains for 14 various elem
      REAL(DOUBLE)                    :: XB(4,3)           ! First 4 rows of XEB
      REAL(DOUBLE)                    :: XL(4,3)           ! First 4 rows of XEL
                                                           ! rigid body motions/constant strain distortions)

! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,2
         DO J=1,12
            BS(I,J) = ZERO
         ENDDO
      ENDDO

! Calc outputs

      JJ = 0

      DO J=1,4

         JJ = JJ + 1
         BS(1,JJ) = DPSHX(1,J)
         BS(2,JJ) = DPSHX(2,J)

         JJ = JJ + 1
         BS(1,JJ) = -DNXSHX(1,J)
         BS(2,JJ) = -DNXSHX(2,J) - PSH(J)

         JJ = JJ + 1
         BS(1,JJ) = DNYSHX(1,J) + PSH(J)
         BS(2,JJ) = DNYSHX(2,J)

      ENDDO

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(8) > 0)) THEN

         WRITE(BUG,1101) ELDT_BUG_BMAT_BIT, TYPE, EID
         WRITE(BUG,8901) IGAUS, JGAUS, SUBR_NAME
         DO I=1,2
            WRITE(BUG,8902) I,(BS(I,J),J=1,12)
            WRITE(BUG,*)
         ENDDO
         WRITE(BUG,*)

      ENDIF

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(9) > 0)) THEN
        IF (DEBUG(202) > 0) THEN
         ID( 1) =  3
         ID( 2) =  4
         ID( 3) =  5
         ID( 4) =  9
         ID( 5) = 10
         ID( 6) = 11
         ID( 7) = 15
         ID( 8) = 16
         ID( 9) = 17
         ID(10) = 21
         ID(11) = 22
         ID(12) = 23

         DO I=1,4
            DO J=1,3
               XB(I,J) = XEB(I,J)
               XL(I,J) = XEL(I,J)
            ENDDO
         ENDDO

         WRITE(BUG,1101) ELDT_BUG_BCHK_BIT, TYPE, EID
         WRITE(BUG,9100)
         WRITE(BUG,9101) MESSAG, IGAUS, JGAUS
         WRITE(BUG,9102)
         WRITE(BUG,9103)
         WRITE(BUG,9104)
         CALL BCHECK_2D ( BS, 'S', ID, NR, NC, 4, XB, XL, BW )
        ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1101 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------',/,                                                                                                &
             ' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8,/,                                                      &
             ' ==============================================================',/)

 8901 FORMAT(' Strain-displacement matrix BS for Gauss point: I = ',I3,', J = ',I3,' for shear portion of element in subr '        &
             ,A,/)

 8902 FORMAT(' Row ',I2,/,6(1ES14.6))


 9100 FORMAT('                  Check on strain-displacement matrix BS for transverse shear portion of the element in subr BCHECK'/)

 9101 FORMAT('                                                        S T R A I N S'/,                                             &
             '                                 (',A,' for Gauss point: I = ',I3,', J = ',I3,')')

 9102 FORMAT('                                                     Gxz            Gyz')

 9103 FORMAT(1X,'      Element displacements consistent with:')

 9104 FORMAT(1X,'      ---------------------------------------')

90004 FORMAT(1X,4(1ES19.11))

! **********************************************************************************************************************************

      END SUBROUTINE BSMIN4

   END MODULE MIN4_B_MATRICES
