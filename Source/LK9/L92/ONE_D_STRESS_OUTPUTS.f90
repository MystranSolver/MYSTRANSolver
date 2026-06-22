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

      SUBROUTINE ONE_D_STRESS_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

! Calculates element stresses for 1D ROD, BAR, ELAS elements
! Calculates element specific stress output from array STRESS (calc'd in subr ELEM_STRE_STRN_ARRAYS) for 1-D elements (ELAS, ROD,
! BAR, BUSH) and puts results into array OGEL for later output to F06 file.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  STRESS, TYPE, ZS
      USE LINK9_STUFF, ONLY           :  MSPRNT, OGEL
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_VECS
      USE PARAMS, ONLY                :  PRTNEU
! --- cbeam_add begin --- !
      USE DEBUG_PARAMETERS
! --- cbeam_add end--- !

      USE ONE_D_STRESS_OUTPUTS_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ONE_D_STRESS_OUTPUTS'
      CHARACTER(1*BYTE)               :: MSP1,MSP2,MSP3     ! Output from subr MARGIN used to control how margins of safety printed
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_OGEL         ! If 'Y' then write data to array OGEL
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_FEMAP        ! If 'Y' then write data to array FEMAP_EL_VECS

      INTEGER(LONG), INTENT(IN)       :: SIZE_ALLOCATED     ! No. of rows allocated to array that will be written to
      INTEGER(LONG), INTENT(IN)       :: NUM_FEMAP_ROWS     ! Number of rows that will be written to FEMAP arrays
      INTEGER(LONG), INTENT(INOUT)    :: NUM1               ! Cum rows written to OGEL prior to running this subr
      INTEGER(LONG)                   :: ICOL               ! An input to subr MARGIN, called by this subr


      REAL(DOUBLE)                    :: C1,C2              ! Coords of point "C" on cross-section of a CBAR where stress is calc'd
      REAL(DOUBLE)                    :: D1,D2              ! Coords of point "D" on cross-section of a CBAR where stress is calc'd
      REAL(DOUBLE)                    :: E1,E2              ! Coords of point "E" on cross-section of a CBAR where stress is calc'd
      REAL(DOUBLE)                    :: F1,F2              ! Coords of point "F" on cross-section of a CBAR where stress is calc'd
      REAL(DOUBLE)                    :: MS1,MS2,MS3        ! Margins of safety calc'd by subr MARGIN
      REAL(DOUBLE)                    :: SA1,SA2,SA3,SA4    ! Stresses at points "C", "D", "E" and "F" at end-a of a CBAR
      REAL(DOUBLE)                    :: SB1,SB2,SB3,SB4    ! Stresses at points "C", "D", "E" and "F" at end-b of a CBAR
      REAL(DOUBLE)                    :: SA,ST              ! Stresses input to subrs that calc margins of safety
      REAL(DOUBLE)                    :: SAMAX,SAMIN        ! Max & min stresses from points "C", "D", "E" and "F" at end-a of CBAR
      REAL(DOUBLE)                    :: SBMAX,SBMIN        ! Max & min stresses from points "C", "D", "E" and "F" at end-b of CBAR
      LOGICAL                         :: WRITE_NEU

! --- cbeam_add begin --- !
      LOGICAL                         :: HAVE_SECTION_POINTS
      INTRINSIC DABS, DMAX1,DMIN1
! --- cbeam_add end--- !

      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Calc engineering stresses from array STRESS and put into array OGEL

! --- cbeam_add begin --- !
! Phase-1 transition: CBEAM is explicitly admitted here, but still uses the legacy end-A/end-B bar-shaped recovery body.
! NX-style intermediate stations will be split into a dedicated BEAM writer path later without changing CBAR behavior.
      IF      ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ')) THEN ! BAR1/CBEAM transition path

! If section y/z stress points are not defined yet, do not fail. Fall back to zero
! fiber offsets so CBEAM can still output axial/torsional terms while station/F06
! work is being developed.
         HAVE_SECTION_POINTS = .FALSE.
         DO ICOL=1,8
            IF (DABS(ZS(ICOL)) > ZERO) THEN
               HAVE_SECTION_POINTS = .TRUE.
               EXIT
            ENDIF
         ENDDO
         IF (HAVE_SECTION_POINTS) THEN
! --- cbeam_add end --- !
            C1 = ZS(1)
            C2 = ZS(2)
            D1 = ZS(3)
            D2 = ZS(4)
            E1 = ZS(5)
            E2 = ZS(6)
            F1 = ZS(7)
            F2 = ZS(8)
! --- cbeam_add begin --- !
         ELSE
            C1 = ZERO
            C2 = ZERO
            D1 = ZERO
            D2 = ZERO
            E1 = ZERO
            E2 = ZERO
            F1 = ZERO
            F2 = ZERO
         ENDIF
! --- cbeam_add end --- !
         SA    =   STRESS(1)
         SA1   = -(C1*STRESS(2) + C2*STRESS(3))
         SA2   = -(D1*STRESS(2) + D2*STRESS(3))
         SA3   = -(E1*STRESS(2) + E2*STRESS(3))
         SA4   = -(F1*STRESS(2) + F2*STRESS(3))
         SB1   = -(C1*STRESS(4) + C2*STRESS(5))
         SB2   = -(D1*STRESS(4) + D2*STRESS(5))
         SB3   = -(E1*STRESS(4) + E2*STRESS(5))
         SB4   = -(F1*STRESS(4) + F2*STRESS(5))
         ST    =   STRESS(6)
         SAMAX = SA + DMAX1(SA1,SA2,SA3,SA4)
         SBMAX = SA + DMAX1(SB1,SB2,SB3,SB4)
         SAMIN = SA + DMIN1(SA1,SA2,SA3,SA4)
         SBMIN = SA + DMIN1(SB1,SB2,SB3,SB4)
! --- cbeam_add begin --- !
         IF (DEBUG(233) > 0) THEN
            WRITE(ERR,9100) TYPE
            WRITE(ERR,9101) STRESS(1), STRESS(2), STRESS(3), STRESS(4), STRESS(5), STRESS(6)
            WRITE(ERR,9102) SA, SA1, SA2, SA3, SA4, SB1, SB2, SB3, SB4, ST, SAMAX, SBMAX, SAMIN, SBMIN
         ENDIF
! --- cbeam_add end --- !
         ICOL = 1
         CALL BAR_MARGIN ( ICOL, SAMAX, SBMAX, SAMIN, SBMIN, ST, MS1, MS2, MS3, MSP1, MSP2, MSP3 )
         IF (WRITE_OGEL == 'Y') THEN
            NUM1 = NUM1 + 1
            IF (NUM1 > SIZE_ALLOCATED) THEN
               WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
               WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            OGEL(NUM1,1) = SA1
            OGEL(NUM1,2) = SA2
            OGEL(NUM1,3) = SA3
            OGEL(NUM1,4) = SA4
            OGEL(NUM1,5) = SA
            OGEL(NUM1,6) = SAMAX
            OGEL(NUM1,7) = SAMIN
            OGEL(NUM1,8) = MS1
            OGEL(NUM1,9) = ST
            MSPRNT(NUM1,1) = MSP1
            NUM1 = NUM1 + 1
            IF (NUM1 > SIZE_ALLOCATED) THEN
               WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
               WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            OGEL(NUM1,1) = SB1
            OGEL(NUM1,2) = SB2
            OGEL(NUM1,3) = SB3
            OGEL(NUM1,4) = SB4
            OGEL(NUM1,5) = SA
            OGEL(NUM1,6) = SBMAX
            OGEL(NUM1,7) = SBMIN
            OGEL(NUM1,8) = MS2
            OGEL(NUM1,9) = MS3
            MSPRNT(NUM1,2) = MSP2
            MSPRNT(NUM1,3) = MSP3
         ENDIF
         IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 1) = SA + SA1
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 2) = SA + SB1
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 3) = SA + SA2
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 4) = SA + SB2
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 5) = SA + SA3
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 6) = SA + SB3
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 7) = SA + SA4
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 8) = SA + SB4
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 9) = SAMAX
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,10) = SBMAX
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,11) = SAMIN
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,12) = SBMIN
         ENDIF

      ELSE IF (TYPE(1:4) == 'BUSH') THEN                        ! BUSH elements
         NUM1 = NUM1 + 1
         IF (NUM1 > SIZE_ALLOCATED) THEN
            WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
            WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF
         OGEL(NUM1,1) = STRESS(1)
         OGEL(NUM1,2) = STRESS(2)
         OGEL(NUM1,3) = STRESS(3)
         OGEL(NUM1,4) = STRESS(4)
         OGEL(NUM1,5) = STRESS(5)
         OGEL(NUM1,6) = STRESS(6)

         IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 1) = STRESS(1)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 2) = STRESS(2)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 3) = STRESS(3)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 4) = STRESS(4)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 5) = STRESS(5)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 6) = STRESS(6)
         ENDIF
      ELSE IF (TYPE(1:4) == 'ELAS') THEN                   ! Stresses for ELAS elements
         IF (WRITE_OGEL == 'Y') THEN
            NUM1 = NUM1 + 1
            IF (NUM1 > SIZE_ALLOCATED) THEN
               WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
               WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            OGEL(NUM1,1) = STRESS(1)
         ENDIF
         IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,1) = STRESS(1)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,2) = STRESS(1)
         ENDIF

      ELSE IF (TYPE == 'ROD     ') THEN                    ! ROD1 elements
         SA = STRESS(1)
         ST = STRESS(2)
         ICOL = 1
         CALL ROD_MARGIN ( ICOL, SA, ST, MS1, MS2, MSP1, MSP2 )
         IF (WRITE_OGEL == 'Y') THEN
            NUM1 = NUM1 + 1
            IF (NUM1 > SIZE_ALLOCATED) THEN
               WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
               WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            OGEL(NUM1,1)   = SA
            OGEL(NUM1,2)   = MS1
            OGEL(NUM1,3)   = ST
            OGEL(NUM1,4)   = MS2
            MSPRNT(NUM1,1) = MSP1
            MSPRNT(NUM1,2) = MSP2
         ENDIF
         IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,1) = SA
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,2) = SA
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,3) = ST
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,4) = ST
         ENDIF

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9203) SUBR_NAME, TYPE
         WRITE(F06,9203) SUBR_NAME, TYPE
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN
! --- cbeam_add begin --- !
! **********************************************************************************************************************************
 9100 FORMAT('*** CBEAM ONE_D STRESS DEBUG ***********************************************',/, &
             '  Element type     : "',A8,'"')
 9101 FORMAT('  Raw STRESS(1:6)  :',6(1X,1ES14.6))
 9102 FORMAT('  SA/pts/max/min   :',14(1X,1ES14.6),/, &
             '***************************************************************************')
! --- cbeam_add end --- !
! **********************************************************************************************************************************
 9200 FORMAT(' *ERROR  9200: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ARRAY OGEL WAS ALLOCATED TO HAVE ',I12,' ROWS. ATTEMPT TO WRITE TO OGEL BEYOND THIS')

! **********************************************************************************************************************************
 9203 FORMAT(' *ERROR  9203: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INCORRECT ELEMENT TYPE = "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE ONE_D_STRESS_OUTPUTS
