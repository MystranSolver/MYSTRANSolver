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

   MODULE QUADRILATERAL_B_MATRICES

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: BMQMEM, BBDKQ, BCHECK_2D

   CONTAINS

      SUBROUTINE BMQMEM ( DPSHX, IGAUS, JGAUS, MESSAG, WRT_BUG_THIS_TIME, BM )

! Calculate BM strain/displ matrix for 4 node membrane isoparametric element (quadratic). Called by subrs QMEM1, QSHEAR

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, WRT_BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_BMAT_BIT, ELDT_BUG_BCHK_BIT
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  BMEANT, EID, HBAR, MXWARP, TYPE, XEB, XEL

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BMQMEM'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Messag to print out if BCHECK is run
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG), INTENT(IN)       :: IGAUS             ! I index of Gauss point (needed for some optional output)
      INTEGER(LONG), INTENT(IN)       :: JGAUS             ! J index of Gauss point (needed for some optional output)
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: JJ                ! A computed index into array BM
      INTEGER(LONG), PARAMETER        :: ID1( 8) = (/ 1, & ! ID1(1) =  1
                                                      2, & ! ID1(2) =  2
                                                      7, & ! ID1(3) =  7
                                                      8, & ! ID1(4) =  8
                                                     13, & ! ID1(5) = 13
                                                     14, & ! ID1(6) = 14
                                                     19, & ! ID1(7) = 19
                                                     20 /) ! ID1(8) = 20

      INTEGER(LONG), PARAMETER        :: ID2(12) = (/ 1, & ! ID2( 1)=  1
                                                      2, & ! ID2( 2)=  2
                                                      3, & ! ID2( 3)=  3
                                                      7, & ! ID2( 4)=  7
                                                      8, & ! ID2( 5)=  8
                                                      9, & ! ID2( 6)=  9
                                                     13, & ! ID2( 7)= 13
                                                     14, & ! ID2( 8)= 14
                                                     15, & ! ID2( 9)= 15
                                                     19, & ! ID2(10)= 19
                                                     20, & ! ID2(11)= 20
                                                     21 /) ! ID2(12)= 21


      REAL(DOUBLE) , INTENT(IN)       :: DPSHX(2,4)        ! Derivatives of the 4 node bilinear isopar interps wrt elem x and y
      REAL(DOUBLE) , INTENT(OUT)      :: BM(3,8)           ! Output strain-displ matrix for this elem
      REAL(DOUBLE)                    :: BM_BMEANT(3,12)   ! Product of BM and BMEANT to be sent to subr BCHECK_2D
      REAL(DOUBLE)                    :: BW(3,14)          ! Output from subr BCHECK (matrix of 3 elem strains for 14 various elem
!                                                            rigid body motions/constant strain distortions)

      REAL(DOUBLE)                    :: XB(4,3)           ! First 4 rows of XEB
      REAL(DOUBLE)                    :: XL(4,3)           ! First 4 rows of XEL


! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,3
         DO J=1,8
            BM(I,J) = ZERO
         ENDDO
      ENDDO

! Calc outputs

      JJ = 0
      DO J=1,4

         JJ = JJ + 1
         BM(1,JJ) = DPSHX(1,J)
         BM(2,JJ) = ZERO
         BM(3,JJ) = DPSHX(2,J)

         JJ = JJ + 1
         BM(1,JJ) = ZERO
         BM(2,JJ) = DPSHX(2,J)
         BM(3,JJ) = DPSHX(1,J)

      ENDDO

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(8) > 0)) THEN

         WRITE(BUG,1101) ELDT_BUG_BMAT_BIT, TYPE, EID
         WRITE(BUG,8901) IGAUS, JGAUS, SUBR_NAME
         DO I=1,3
            WRITE(BUG,8902) I,(BM(I,J),J=1,8)
            WRITE(BUG,*)
         ENDDO
         WRITE(BUG,*)

      ENDIF

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(9) > 0)) THEN
        IF (DEBUG(202) > 0) THEN
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
         IF ((DABS(HBAR) > MXWARP) .AND. (DEBUG(4) ==  0)) THEN
            CALL MATMULT_FFF (BM, BMEANT, 3, 8, 12, BM_BMEANT )
            CALL BCHECK_2D ( BM_BMEANT, 'M', ID2, 3, 12, 4, XL, XB, BW )
         ELSE
            CALL BCHECK_2D ( BM, 'M', ID1, 3, 8, 4, XB, XL, BW )
         ENDIF
        ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1101 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------',/,                                                                                                &
             ' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8,/,                                                      &
             ' ==============================================================',/)

 8901 FORMAT(' Strain-displacement matrix BM for Gauss point: I = ',I3,', J = ',I3,' for membrane portion of element in subr '     &
             ,A,/)

 8902 FORMAT(' Row ',I2,/,8(1ES14.6))

 9100 FORMAT('                          Check on strain-displacement matrix BM for membrane portion of the element in subr BCHECK'/)

 9101 FORMAT('                                                               S T R A I N S'/,                                      &
             '                                             (',A,' for Gauss point: I = ',I3,', J = ',I3,')')

 9102 FORMAT('                                                     Exx            Eyy            Exy')

 9103 FORMAT(1X,'      Element displacements consistent with:')

 9104 FORMAT(1X,'      ---------------------------------------')

! **********************************************************************************************************************************

      END SUBROUTINE BMQMEM


      SUBROUTINE BBDKQ ( DPSHX, XSD, YSD, SLN, IGAUS, JGAUS, MESSAG, WRT_BUG_THIS_TIME, BB )

! Calculate BB strain/displacement matrix for DKQ bending quadrilateral element. Called by subr QPLT1

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, WRT_BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_BMAT_BIT, ELDT_BUG_BCHK_BIT
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, THREE, FOUR
      USE MODEL_STUF, ONLY            :  EID, TYPE, XEB, XEL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BBDKQ'
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Messag to print out if BCHECK is run
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG), INTENT(IN)       :: IGAUS             ! I index of Gaus point (needed for some optional output)
      INTEGER(LONG), INTENT(IN)       :: JGAUS             ! J index of Gaus point (needed for some optional output)
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: ID(12)            ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NR        = 3     ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NC        = 12    ! An input to subr BCHECK, called herein


      REAL(DOUBLE) , INTENT(IN)       :: SLN(4)            ! Quad side lengths
      REAL(DOUBLE) , INTENT(IN)       :: XSD(4)            ! Array of 4 diffs of X dim. of sides
      REAL(DOUBLE) , INTENT(IN)       :: YSD(4)            ! Array of 4 diffs of Y dim. of sides
      REAL(DOUBLE) , INTENT(IN)       :: DPSHX(2,8)        ! Derivatives of the 8 node biquadratic isopar interps wrt elem x and y
      REAL(DOUBLE) , INTENT(OUT)      :: BB(3,12)          ! Output strain-displ matrix for the DKQ elem
      REAL(DOUBLE)                    :: BW(3,14)          ! Output from subr BCHECK (matrix of 3 elem strains for 14 various elem
                                                             ! rigid body motions/constant strain distortions)
      REAL(DOUBLE) , PARAMETER        :: C15 = THREE/TWO   ! Constant = 1.5
      REAL(DOUBLE)                    :: A(4)              ! Intermediate variables used in calculating outputs
      REAL(DOUBLE)                    :: B(4)              ! Intermediate variables used in calculating outputs
      REAL(DOUBLE)                    :: C(4)              ! Intermediate variables used in calculating outputs
      REAL(DOUBLE)                    :: D(4)              ! Intermediate variables used in calculating outputs
      REAL(DOUBLE)                    :: E(4)              ! Intermediate variables used in calculating outputs
      REAL(DOUBLE)                    :: DHXSHX(2,12)      ! Derivatives of Hx with respect to x and y (Hx is a fcn of DPSHX)
      REAL(DOUBLE)                    :: DHYSHX(2,12)      ! Derivatives of Hy with respect to x and y (Hy is a fcn of DPSHX)
      REAL(DOUBLE)                    :: SL2               ! The squares of elem side lengths
      REAL(DOUBLE)                    :: XB(4,3)           ! First 4 rows of XEB
      REAL(DOUBLE)                    :: XL(4,3)           ! First 4 rows of XEL

! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,3
         DO J=1,12
            BB(I,J) = ZERO
         ENDDO
      ENDDO

! Calculate parameters needed

      DO I=1,4
         SL2 = SLN(I)*SLN(I)
         A(I) = -XSD(I)/SL2
         B(I) = THREE*XSD(I)*YSD(I)/(FOUR*SL2)
         C(I) = (XSD(I)*XSD(I)/FOUR - YSD(I)*YSD(I)/TWO)/SL2
         D(I) = -YSD(I)/SL2
         E(I) = (-XSD(I)*XSD(I)/TWO + YSD(I)*YSD(I)/FOUR)/SL2
      ENDDO

! Derivatives of Hx with respect to x

      DHXSHX(1, 1) = C15*(A(1)*DPSHX(1,5) - A(4)*DPSHX(1,8))
      DHXSHX(1, 4) = C15*(A(2)*DPSHX(1,6) - A(1)*DPSHX(1,5))
      DHXSHX(1, 7) = C15*(A(3)*DPSHX(1,7) - A(2)*DPSHX(1,6))
      DHXSHX(1,10) = C15*(A(4)*DPSHX(1,8) - A(3)*DPSHX(1,7))

      DHXSHX(1, 2) = B(1)*DPSHX(1,5) + B(4)*DPSHX(1,8)
      DHXSHX(1, 5) = B(2)*DPSHX(1,6) + B(1)*DPSHX(1,5)
      DHXSHX(1, 8) = B(3)*DPSHX(1,7) + B(2)*DPSHX(1,6)
      DHXSHX(1,11) = B(4)*DPSHX(1,8) + B(3)*DPSHX(1,7)

      DHXSHX(1, 3) = DPSHX(1,1) - C(1)*DPSHX(1,5) - C(4)*DPSHX(1,8)
      DHXSHX(1, 6) = DPSHX(1,2) - C(2)*DPSHX(1,6) - C(1)*DPSHX(1,5)
      DHXSHX(1, 9) = DPSHX(1,3) - C(3)*DPSHX(1,7) - C(2)*DPSHX(1,6)
      DHXSHX(1,12) = DPSHX(1,4) - C(4)*DPSHX(1,8) - C(3)*DPSHX(1,7)

! Derivatives of Hx with respect to y

      DHXSHX(2, 1) = C15*(A(1)*DPSHX(2,5) - A(4)*DPSHX(2,8))
      DHXSHX(2, 4) = C15*(A(2)*DPSHX(2,6) - A(1)*DPSHX(2,5))
      DHXSHX(2, 7) = C15*(A(3)*DPSHX(2,7) - A(2)*DPSHX(2,6))
      DHXSHX(2,10) = C15*(A(4)*DPSHX(2,8) - A(3)*DPSHX(2,7))

      DHXSHX(2, 2) = B(1)*DPSHX(2,5) + B(4)*DPSHX(2,8)
      DHXSHX(2, 5) = B(2)*DPSHX(2,6) + B(1)*DPSHX(2,5)
      DHXSHX(2, 8) = B(3)*DPSHX(2,7) + B(2)*DPSHX(2,6)
      DHXSHX(2,11) = B(4)*DPSHX(2,8) + B(3)*DPSHX(2,7)

      DHXSHX(2, 3) = DPSHX(2,1) - C(1)*DPSHX(2,5) - C(4)*DPSHX(2,8)
      DHXSHX(2, 6) = DPSHX(2,2) - C(2)*DPSHX(2,6) - C(1)*DPSHX(2,5)
      DHXSHX(2, 9) = DPSHX(2,3) - C(3)*DPSHX(2,7) - C(2)*DPSHX(2,6)
      DHXSHX(2,12) = DPSHX(2,4) - C(4)*DPSHX(2,8) - C(3)*DPSHX(2,7)

! Derivatives of Hy with respect to x

      DHYSHX(1, 1) = C15*(D(1)*DPSHX(1,5) - D(4)*DPSHX(1,8))
      DHYSHX(1, 4) = C15*(D(2)*DPSHX(1,6) - D(1)*DPSHX(1,5))
      DHYSHX(1, 7) = C15*(D(3)*DPSHX(1,7) - D(2)*DPSHX(1,6))
      DHYSHX(1,10) = C15*(D(4)*DPSHX(1,8) - D(3)*DPSHX(1,7))

      DHYSHX(1, 2) = -DPSHX(1,1) + E(1)*DPSHX(1,5) + E(4)*DPSHX(1,8)
      DHYSHX(1, 5) = -DPSHX(1,2) + E(2)*DPSHX(1,6) + E(1)*DPSHX(1,5)
      DHYSHX(1, 8) = -DPSHX(1,3) + E(3)*DPSHX(1,7) + E(2)*DPSHX(1,6)
      DHYSHX(1,11) = -DPSHX(1,4) + E(4)*DPSHX(1,8) + E(3)*DPSHX(1,7)

      DHYSHX(1, 3) = -B(1)*DPSHX(1,5) - B(4)*DPSHX(1,8)
      DHYSHX(1, 6) = -B(2)*DPSHX(1,6) - B(1)*DPSHX(1,5)
      DHYSHX(1, 9) = -B(3)*DPSHX(1,7) - B(2)*DPSHX(1,6)
      DHYSHX(1,12) = -B(4)*DPSHX(1,8) - B(3)*DPSHX(1,7)

! Derivatives of Hy with respect to y

      DHYSHX(2, 1) = C15*(D(1)*DPSHX(2,5) - D(4)*DPSHX(2,8))
      DHYSHX(2, 4) = C15*(D(2)*DPSHX(2,6) - D(1)*DPSHX(2,5))
      DHYSHX(2, 7) = C15*(D(3)*DPSHX(2,7) - D(2)*DPSHX(2,6))
      DHYSHX(2,10) = C15*(D(4)*DPSHX(2,8) - D(3)*DPSHX(2,7))

      DHYSHX(2, 2) = -DPSHX(2,1) + E(1)*DPSHX(2,5) + E(4)*DPSHX(2,8)
      DHYSHX(2, 5) = -DPSHX(2,2) + E(2)*DPSHX(2,6) + E(1)*DPSHX(2,5)
      DHYSHX(2, 8) = -DPSHX(2,3) + E(3)*DPSHX(2,7) + E(2)*DPSHX(2,6)
      DHYSHX(2,11) = -DPSHX(2,4) + E(4)*DPSHX(2,8) + E(3)*DPSHX(2,7)

      DHYSHX(2, 3) = -B(1)*DPSHX(2,5) - B(4)*DPSHX(2,8)
      DHYSHX(2, 6) = -B(2)*DPSHX(2,6) - B(1)*DPSHX(2,5)
      DHYSHX(2, 9) = -B(3)*DPSHX(2,7) - B(2)*DPSHX(2,6)
      DHYSHX(2,12) = -B(4)*DPSHX(2,8) - B(3)*DPSHX(2,7)

! Now formulate the BB strain/displacement matrix

      DO J=1,12
         BB(1,J) = DHXSHX(1,J)
         BB(2,J) = DHYSHX(2,J)
         BB(3,J) = DHXSHX(2,J) + DHYSHX(1,J)
      ENDDO

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(8) > 0)) THEN

         WRITE(BUG,1101) ELDT_BUG_BMAT_BIT, TYPE, EID
         WRITE(BUG,8901) IGAUS, JGAUS, SUBR_NAME
         DO I=1,3
            WRITE(BUG,8902) I,(BB(I,J),J=1,12)
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
         CALL BCHECK_2D ( BB, 'B', ID, NR, NC, 4, XB, XL, BW )
        ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1101 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------',/,                                                                                                &
             ' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8, ' (from subr ', A, ')'/,                               &
             ' ==============================================================',/)

 8901 FORMAT(' Strain-displacement matrix BB for Gauss point: I = ',I3,', J = ',I3,' for bending portion of element in subr '      &
             ,A,/)

 8902 FORMAT(' Row ',I2,/,6(1ES14.6))

 9100 FORMAT('                          Check on strain-displacement matrix BB for bending portion of the element in subr BCHECK'/)

 9101 FORMAT('                                                               S T R A I N S'/,                                      &
             '                                           (',A,' for Gauss point: I = ',I3,', J = ',I3,')')

 9102 FORMAT('                                                     Cxx            Cyy            Cxy')

 9103 FORMAT(1X,'      Element displacements consistent with:')

 9104 FORMAT(1X,'      ---------------------------------------')

! **********************************************************************************************************************************

      END SUBROUTINE BBDKQ


      SUBROUTINE BCHECK_2D ( B, BTYPE, ID, NROWB, NCOLB, NUM_GRIDS, XB, XL, BW )

! Checks strain-displacement matrices for rigid body motion and constant strain for 2-D shell elements
! (6 DOF per grid point and up to 4 grid points)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO
      USE MODEL_STUF, ONLY            :  ELDOF, NELGP, TE
      USE MODEL_STUF, ONLY            :  AGRID, ELGP

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE VECTOR_GEOMETRY, ONLY       :  RIGID_BODY_DISP_MAT

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BCHECK_2D'
      CHARACTER(LEN=*), INTENT(IN)    :: BTYPE             ! Type of B matrix ('M' for membrane, 'B' for bending, 'S' for shear)
      CHARACTER(45*BYTE)              :: MESSAG(14)        ! Output messages for the 14 modes of deformation (6 RB + 8 const strain)

      INTEGER(LONG), INTENT(IN)       :: NROWB             ! Number of rows in the input B matrix
      INTEGER(LONG), INTENT(IN)       :: NCOLB             ! Number of cols in the input B matrix
      INTEGER(LONG), INTENT(IN)       :: NUM_GRIDS         ! Number of grids for the input B matrix
      INTEGER(LONG), INTENT(IN)       :: ID(NCOLB)         ! List of elem DOF's for each of the elem grids (e.g 3,4,5 for each of
!                                                            4 grids for a 4 node plate bending elem
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: KK                ! A computed index into array W


      REAL(DOUBLE) , INTENT(IN)       :: B(NROWB,NCOLB)    ! Strain-displ matrix
      REAL(DOUBLE) , INTENT(IN)       :: XB(NUM_GRIDS,3)   ! Basic coords of elem grids (diff than XEB for TPLT2's in a MIN4T QUAD4)
      REAL(DOUBLE) , INTENT(IN)       :: XL(NUM_GRIDS,3)   ! Local coords of elem grids (diff than XEL for TPLT2's in a MIN4T QUAD4)
      REAL(DOUBLE) , INTENT(OUT)      :: BW(NROWB,14)      ! Output from subr BCHECK_2D (matrix of NROWB elem strains for various
!                                                            elem rigid body motions/constant strain distortions)
      REAL(DOUBLE)                    :: GRD_COORDS(3)     ! 3 coords from XEB for one node of the element
      REAL(DOUBLE)                    :: REF_COORDS(3)     ! 3 coords from XEB for node 1
      REAL(DOUBLE)                    :: RB_DISP(6,6)      ! 6 x 6 RB matrix for one grid for this element
      REAL(DOUBLE)                    :: W(24,14)          ! Displs for the 14 modes of elem deformation (6 RB + 8 constant strain)



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NROWB
         DO J=1,14
            BW(I,J) = ZERO
         ENDDO
      ENDDO

! Init W

      DO I=1,24
         DO J=1,14
            W(I,J) = ZERO
         ENDDO
      ENDDO

      MESSAG( 1) = '      Rigid Body Displacement,      T1  = 1.0'
      MESSAG( 2) = '      Rigid Body Displacement,      T2  = 1.0'
      MESSAG( 3) = '      Rigid Body Displacement,      T3  = 1.0'
      MESSAG( 4) = '      Rigid Body Rotation,          R1  = 1.0'
      MESSAG( 5) = '      Rigid Body Rotation,          R2  = 1.0'
      MESSAG( 6) = '      Rigid Body Rotation,          R3  = 1.0'
      MESSAG( 7) = '      Constant in-plane strain,     Exx = 1.0'
      MESSAG( 8) = '      Constant in-plane strain,     Eyy = 1.0'
      MESSAG( 9) = '      Constant in-plane strain,     Exy = 1.0'
      MESSAG(10) = '      Constant curvature,           Cxx = 1.0'
      MESSAG(11) = '      Constant curvature,           Cyy = 1.0'
      MESSAG(12) = '      Constant curvature,           Cxy = 1.0'
      MESSAG(13) = '      Constant transv shear strain, Gxz = 1.0'
      MESSAG(14) = '      Constant transv shear strain, Gyz = 1.0'

! Calc RB modes (cols 1-6 of W)

      DO I=1,NUM_GRIDS
         DO J=1,3
            GRD_COORDS(J) = XB(I,J)
            REF_COORDS(J) = XB(1,J)
         ENDDO
         CALL RIGID_BODY_DISP_MAT ( GRD_COORDS, REF_COORDS, RB_DISP )
         DO J=1,6
            DO K=1,6
               W(6*(I-1)+J,K) = RB_DISP(J,K)
            ENDDO
         ENDDO
      ENDDO

! Calc constant strain modes (cols 7-14 of W)

      DO K=1,NUM_GRIDS

         KK = 6*(K-1)

         W(KK+1, 7) =  XL(K,1)                            ! This gives constant direct x strain

         W(KK+2, 8) =  XL(K,2)                            ! This gives constant direct y strain

         W(KK+1, 9) =  XL(K,2)/TWO                        ! The next 2 give constant shear xy strain
         W(KK+2, 9) =  XL(K,1)/TWO

         W(KK+3,10) = -XL(K,1)*XL(K,1)/TWO                ! This gives constant curvature, Cxx
         W(KK+5,10) =  XL(K,1)

         W(KK+3,11) = -XL(K,2)*XL(K,2)/TWO                ! This gives constant curvature, Cyy
         W(KK+4,11) = -XL(K,2)

         W(KK+3,12) = -XL(K,1)*XL(K,2)/TWO                ! This gives constant curvature, Cxy
         W(KK+4,12) = -XL(K,1)/TWO
         W(KK+5,12) =  XL(K,2)/TWO

         W(KK+3,13) =  XL(K,1)                            ! This gives constant transverse shear strain

         W(KK+3,14) =  XL(K,2)                            ! This gives constant transverse shear strain

      ENDDO

! Calc BW = B*W (but B has fewer cols since element has no stiffness for some DOF's)

      DO I=1,NROWB
         DO J=1,14
            BW(I,J) = ZERO
            DO K=1,NCOLB
               BW(I,J) = BW(I,J) + B(I,K)*W(ID(K),J)
            ENDDO
         ENDDO
      ENDDO

! Write results

      IF ((BTYPE == 'B') .OR. (BTYPE == 'M')) THEN
         DO J=1,6
            WRITE(BUG,9101) MESSAG(J),(BW(I,J),I=1,NROWB)
         ENDDO
      ELSE IF (BTYPE == 'S') THEN
         DO J=1,6
            WRITE(BUG,9201) MESSAG(J),(BW(I,J),I=1,NROWB)
         ENDDO
      ENDIF
      WRITE(BUG,*)

      IF (BTYPE == 'B') THEN
         DO J=7,9
            WRITE(BUG,9101) MESSAG(J),(BW(I,J),I=1,NROWB)
         ENDDO
      ELSE IF (BTYPE == 'M') THEN
         WRITE(BUG,9102) MESSAG(7),(BW(I,7),I=1,NROWB)
         WRITE(BUG,9103) MESSAG(8),(BW(I,8),I=1,NROWB)
         WRITE(BUG,9104) MESSAG(9),(BW(I,9),I=1,NROWB)
      ELSE IF (BTYPE == 'S') THEN
         DO J=7,9
            WRITE(BUG,9201) MESSAG(J),(BW(I,J),I=1,NROWB)
         ENDDO
      ENDIF
      WRITE(BUG,*)

      IF (BTYPE == 'B') THEN
         WRITE(BUG,9102) MESSAG(10),(BW(I,10),I=1,NROWB)
         WRITE(BUG,9103) MESSAG(11),(BW(I,11),I=1,NROWB)
         WRITE(BUG,9104) MESSAG(12),(BW(I,12),I=1,NROWB)
      ELSE IF (BTYPE == 'M') THEN
         DO J=10,12
            WRITE(BUG,9101) MESSAG(J),(BW(I,J),I=1,NROWB)
         ENDDO
      ELSE IF (BTYPE == 'S') THEN
         DO J=10,12
            WRITE(BUG,9201) MESSAG(J),(BW(I,J),I=1,NROWB)
         ENDDO
      ENDIF
      WRITE(BUG,*)

      IF ((BTYPE == 'B') .OR. (BTYPE == 'M')) THEN
         DO J=13,14
           WRITE(BUG,9101) MESSAG(J),(BW(I,J),I=1,NROWB)
         ENDDO
      ELSE IF (BTYPE == 'S') THEN
         WRITE(BUG,9202) MESSAG(13),(BW(I,13),I=1,NROWB)
         WRITE(BUG,9203) MESSAG(14),(BW(I,14),I=1,NROWB)
      ENDIF
      WRITE(BUG,*)



      RETURN

! **********************************************************************************************************************************
 9101 FORMAT(1X,A,3(1ES15.6),'  (should be 0, 0, 0)')

 9102 FORMAT(1X,A,3(1ES15.6),'  (should be 1, 0, 0)')

 9103 FORMAT(1X,A,3(1ES15.6),'  (should be 0, 1, 0)')

 9104 FORMAT(1X,A,3(1ES15.6),'  (should be 0, 0, 1)')

 9201 FORMAT(1X,A,2(1ES15.6),'  (should be 0, 0)')

 9202 FORMAT(1X,A,2(1ES15.6),'  (should be 1, 0)')

 9203 FORMAT(1X,A,2(1ES15.6),'  (should be 0, 1)')


! **********************************************************************************************************************************

      END SUBROUTINE BCHECK_2D

   END MODULE QUADRILATERAL_B_MATRICES
