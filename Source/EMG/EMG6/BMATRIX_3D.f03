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

   MODULE BMATRIX_3D

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: B3D_ISOPARAMETRIC

   CONTAINS

      SUBROUTINE B3D_ISOPARAMETRIC ( DPSHX, GAUSS_PT, IGAUS, JGAUS, KGAUS, MESSAG, WRT_BUG_THIS_TIME, BMAT )

! Generates strain/displ matrix BMAT for solid 3D elements. Called by HEXA, PENTA and TETRA subroutines

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELDT_BUG_BMAT_BIT, ELDT_BUG_BCHK_BIT
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  EID, ELGP, TYPE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'B3D_ISOPARAMETRIC'
      CHARACTER(LEN=*)  , INTENT(IN)  :: MESSAG            ! Messag to print out if BCHECK is run
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME ! If 'Y' then write to BUG file if WRT_BUG array says to

      INTEGER(LONG), INTENT(IN)       :: GAUSS_PT          ! Gauss point number (used for DEBUG output)
      INTEGER(LONG), INTENT(IN)       :: IGAUS             ! I index of Gaus point (needed for some optional output)
      INTEGER(LONG), INTENT(IN)       :: JGAUS             ! J index of Gaus point (needed for some optional output)
      INTEGER(LONG), INTENT(IN)       :: KGAUS             ! K index of Gaus point (needed for some optional output)
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: II,JJ             ! Counters
      INTEGER(LONG)                   :: ID(3*ELGP)        ! An input to subr BCHECK, called herein
      INTEGER(LONG), PARAMETER        :: NR      = 6       ! An input to subr BCHECK, called herein


      REAL(DOUBLE) , INTENT(IN)       :: DPSHX(3,ELGP)     ! Derivatives of the 4 node bilinear isopar interps wrt elem x and y
      REAL(DOUBLE) , INTENT(OUT)      :: BMAT(6,3*ELGP)    ! Output strain-displ matrix for this elem
      REAL(DOUBLE)                    :: BW(6,12)          ! Output from subr BCHECK (matrix of 3 elem strains for 14 various elem
                                                           ! rigid body motions/constant strain distortions)

! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,6
         DO J=1,3*ELGP
            BMAT(I,J) = ZERO
         ENDDO
      ENDDO

! Calc outputs

      JJ = 0
      DO J=1,ELGP

         JJ = JJ + 1
         BMAT(1,JJ) = DPSHX(1,J)
         BMAT(2,JJ) = ZERO
         BMAT(3,JJ) = ZERO
         BMAT(4,JJ) = DPSHX(2,J)
         BMAT(5,JJ) = ZERO
         BMAT(6,JJ) = DPSHX(3,J)

         JJ = JJ + 1
         BMAT(1,JJ) = ZERO
         BMAT(2,JJ) = DPSHX(2,J)
         BMAT(3,JJ) = ZERO
         BMAT(4,JJ) = DPSHX(1,J)
         BMAT(5,JJ) = DPSHX(3,J)
         BMAT(6,JJ) = ZERO

         JJ = JJ + 1
         BMAT(1,JJ) = ZERO
         BMAT(2,JJ) = ZERO
         BMAT(3,JJ) = DPSHX(3,J)
         BMAT(4,JJ) = ZERO
         BMAT(5,JJ) = DPSHX(2,J)
         BMAT(6,JJ) = DPSHX(1,J)

      ENDDO

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(8) > 0)) THEN

         WRITE(BUG,1101) ELDT_BUG_BMAT_BIT, TYPE, EID
         WRITE(BUG,8901) IGAUS, JGAUS, KGAUS, SUBR_NAME
         DO I=1,6
            WRITE(BUG,8902) I,(BMAT(I,J),J=1,3*ELGP)
            WRITE(BUG,*)
         ENDDO
         WRITE(BUG,*)

      ENDIF

      IF ((WRT_BUG_THIS_TIME == 'Y') .AND. (WRT_BUG(9) > 0)) THEN
        IF (DEBUG(202) > 0) THEN
         JJ = 0                                            ! Calculate ID array
         DO I=1,ELGP
            II = 6*(I - 1)
            DO J=1,3
               II = II + 1
               JJ = JJ + 1
               ID(JJ) = II                                 ! Want JJ to range 1 to 3*ELGP and II to range 1 to 45 in increments of 3
            ENDDO
         ENDDO

         WRITE(BUG,1101) ELDT_BUG_BCHK_BIT, TYPE, EID
         WRITE(BUG,9100)
         WRITE(BUG,9101) MESSAG, IGAUS, JGAUS, KGAUS
         WRITE(BUG,9102)
         WRITE(BUG,9103)
         WRITE(BUG,9104)
         CALL BCHECK_3D ( BMAT, ELGP, ID, NR, 3*ELGP, BW )
        ENDIF
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1101 FORMAT(' ------------------------------------------------------------------------------------------------------------------',&
             '-----------------',/,                                                                                                &
             ' ELDATA(',I2,',PRINT) requests for ',A,' element number ',I8, ' (from subr ', A, ')'/,                               &
             ' ==============================================================',/)

 8901 FORMAT(' Strain-displacement matrix B for Gauss point: I = ',I3,', J = ',I3,', K = ',I3,' for element in subr ',A,/)

 8902 FORMAT(' Row ',I2,/,6(1ES14.6))

 9100 FORMAT(47X,'   Check on strain-displacement matrix B for the element in subr BCHECK,',/)

 9101 FORMAT(47X,'                               S T R A I N S'/,                       &
             47X,'        (',A,' for Gauss point: I = ',I3,', J = ',I3,', K = ',I3,')')

 9102 FORMAT(47X,' Exx           Eyy           Ezz           Gxy           Gyz           Gzx      Note:')

 9103 FORMAT(1X,'      Element displacements consistent with:')

 9104 FORMAT(1X,'      ---------------------------------------')

! **********************************************************************************************************************************

      END SUBROUTINE B3D_ISOPARAMETRIC


      SUBROUTINE BCHECK_3D ( B, NUM_GRIDS, ID, NROWB, NCOLB, BW )

! Checks strain-displacement matrices for rigid body motion and constant strain for 3-D solid elements

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, BUG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO
      USE MODEL_STUF, ONLY            :  AGRID, TE, XEB, XEL

      USE VECTOR_GEOMETRY, ONLY       :  RIGID_BODY_DISP_MAT
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'BCHECK_3D'
      CHARACTER(41*BYTE)              :: MESSAG(12)        ! Output messages for the 14 modes of deformation (6 RB + 8 const strain)
      CHARACTER( 3*BYTE)              :: NOTE(12)          ! Output message

      INTEGER(LONG), INTENT(IN)       :: NCOLB             ! Number of cols in the input B matrix
      INTEGER(LONG), INTENT(IN)       :: NROWB             ! Number of rows in the input B matrix
      INTEGER(LONG), INTENT(IN)       :: NUM_GRIDS         ! Number of grids that this solid element has.
      INTEGER(LONG), INTENT(IN)       :: ID(NCOLB)         ! List of elem DOF's for each of the elem grids (e.g 3,4,5 for each of
!                                                            4 grids for a 4 node plate bending elem
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: KK                ! A computed index into array W


      REAL(DOUBLE) , INTENT(IN)       :: B(NROWB,NCOLB)    ! Strain-displ matrix
      REAL(DOUBLE) , INTENT(OUT)      :: BW(NROWB,12)      ! Output from subr BCHECK_3D (matrix of NROWB elem strains for various
!                                                            elem rigid body motions/constant strain distortions)
      REAL(DOUBLE)                    :: GRD_COORDS(3)     ! 3 coords from XEB for one node of the element
      REAL(DOUBLE)                    :: REF_COORDS(3)     ! 3 coords from XEB for node 1
      REAL(DOUBLE)                    :: RB_DISP_0(6,6)    ! 6 x 6 RB matrix for one grid for this element in basic coords
      REAL(DOUBLE)                    :: RB_DISP_L(6,6)    ! 6 x 6 RB matrix for one grid for this element in elenent local coords
      REAL(DOUBLE)                    :: TE6(6,6)          ! 6 x 6 transformation matrix with TE as 2 diagonal 3 x 3 matrices
      REAL(DOUBLE)                    :: W(6*NUM_GRIDS,12) ! Displs for the 14 modes of elem deformation (6 RB + 6 constant strain)



! **********************************************************************************************************************************
! Initialize outputs

      DO I=1,NROWB
         DO J=1,12
            BW(I,J) = ZERO
         ENDDO
      ENDDO

! Form TE6 from 2 diagonal TE matrices

      DO I=1,6
         DO J=1,6
            TE6(I,J) = ZERO
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,3
            TE6(I,J)     = TE(I,J)
            TE6(I+3,J+3) = TE(I,J)
         ENDDO
      ENDDO

! Init W

      DO I=1,6*NUM_GRIDS
         DO J=1,12
            W(I,J) = ZERO
         ENDDO
      ENDDO

      MESSAG( 1) = '      Rigid Body Displacement,  T1  = 1.0:'
      MESSAG( 2) = '      Rigid Body Displacement,  T2  = 1.0:'
      MESSAG( 3) = '      Rigid Body Displacement,  T3  = 1.0:'
      MESSAG( 4) = '      Rigid Body Rotation,      R1  = 1.0:'
      MESSAG( 5) = '      Rigid Body Rotation,      R2  = 1.0:'
      MESSAG( 6) = '      Rigid Body Rotation,      R3  = 1.0:'
      MESSAG( 7) = '      Constant direct x strain, Exx = 1.0:'
      MESSAG( 8) = '      Constant direct y strain, Eyy = 1.0:'
      MESSAG( 9) = '      Constant direct z strain, Ezz = 1.0:'
      MESSAG(10) = '      Constant shear xy strain, Gxy = 1.0:'
      MESSAG(11) = '      Constant shear yz strain, Gxy = 1.0:'
      MESSAG(12) = '      Constant shear zx strain, Gxy = 1.0:'

      NOTE( 1)   = '(a)'
      NOTE( 2)   = '(a)'
      NOTE( 3)   = '(a)'
      NOTE( 4)   = '(a)'
      NOTE( 5)   = '(a)'
      NOTE( 6)   = '(a)'
      NOTE( 7)   = '(b)'
      NOTE( 8)   = '(c)'
      NOTE( 9)   = '(d)'
      NOTE(10)   = '(e)'
      NOTE(11)   = '(f)'
      NOTE(12)   = '(g)'

! Calc RB modes (cols 1-6 of W)

      DO I=1,NUM_GRIDS
         DO J=1,3
            GRD_COORDS(J) = XEB(I,J)
            REF_COORDS(J) = XEB(1,J)
         ENDDO
         CALL RIGID_BODY_DISP_MAT ( GRD_COORDS, REF_COORDS, RB_DISP_0 )
         CALL MATMULT_FFF ( TE6, RB_DISP_0, 6, 6, 6, RB_DISP_L )
         DO J=1,6
            DO K=1,6
               W(6*(I-1)+J,K) = RB_DISP_L(J,K)
            ENDDO
         ENDDO
      ENDDO

! Calc constant strain modes (cols 7-14 of W)

      DO K=1,NUM_GRIDS

         KK = 6*(K-1)

         W(KK+1, 7) =  XEL(K,1)                            ! This gives constant direct x strain

         W(KK+2, 8) =  XEL(K,2)                            ! This gives constant direct y strain

         W(KK+3, 9) =  XEL(K,3)                            ! This gives constant direct z strain

         W(KK+1,10) =  XEL(K,2)/TWO                        ! The next 2 give constant shear xy strain
         W(KK+2,10) =  XEL(K,1)/TWO

         W(KK+2,11) =  XEL(K,3)/TWO                        ! The next 2 give constant shear yz strain
         W(KK+3,11) =  XEL(K,2)/TWO

         W(KK+3,12) =  XEL(K,1)/TWO                        ! The next 2 give constant shear zx strain
         W(KK+1,12) =  XEL(K,3)/TWO

      ENDDO

! Calc BW = B*W (but B has fewer cols since element has no stiffness for 4,5,6 DOF's)

      DO I=1,NROWB
         DO J=1,12
            BW(I,J) = ZERO
            DO K=1,NCOLB
               BW(I,J) = BW(I,J) + B(I,K)*W(ID(K),J)
            ENDDO
         ENDDO
      ENDDO

! Write results

      DO J=1,3
         WRITE(BUG,9101) MESSAG(J),(BW(I,J),I=1,NROWB), NOTE(J)
      ENDDO
      WRITE(BUG,*)

      DO J=4,6
         WRITE(BUG,9101) MESSAG(J),(BW(I,J),I=1,NROWB), NOTE(J)
      ENDDO
      WRITE(BUG,*)

      DO J=7,9
         WRITE(BUG,9101) MESSAG(J),(BW(I,J),I=1,NROWB), NOTE(J)
      ENDDO
      WRITE(BUG,*)

      DO J=10,12
         WRITE(BUG,9101) MESSAG(J),(BW(I,J),I=1,NROWB), NOTE(J)
      ENDDO
      WRITE(BUG,*)

      WRITE(BUG,*)
      WRITE(BUG,9901)



      RETURN

! **********************************************************************************************************************************
 9101 FORMAT(1X,A,6(1ES14.6),2X,A)

 9901 FORMAT(' Notes: (a) All strains should be zero',/,                                                                           &
             '        (b) Exx should be 1.0, others 0',/,                                                                          &
             '        (c) Eyy should be 1.0, others 0',/,                                                                          &
             '        (d) Ezz should be 1.0, others 0',/,                                                                          &
             '        (e) Gxy should be 1.0, others 0',/,                                                                          &
             '        (f) Gyz should be 1.0, others 0',/,                                                                          &
             '        (g) Gzx should be 1.0, others 0',/)

! **********************************************************************************************************************************

      END SUBROUTINE BCHECK_3D

   END MODULE BMATRIX_3D
