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

   MODULE ELEMENT_STRAIN_RECOVERY

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: CALC_ELEM_STRAINS

   CONTAINS

      SUBROUTINE CALC_ELEM_STRAINS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

! Calls routines that process the STRAIN array calculated in subr ELEM_STRE_STRN_ARRAYS to obtain element specific strain values
! that will be written to the F06 file.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  TYPE
      USE CONSTANTS_1, ONLY           :  ZERO
      USE ROD_BAR_OUTPUT, ONLY        :  ROD_MARGIN

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_ELEM_STRAINS'
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_OGEL        ! If 'Y' then write data to array OGEL
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_FEMAP       ! If 'Y' then write data to array FEMAP arrays

      INTEGER(LONG), INTENT(IN)       :: SIZE_ALLOCATED    ! No. of rows allocated to array that will be written to in a subr
!                                                            called here so we can check that we don't try to write more rows
      INTEGER(LONG), INTENT(IN)       :: NUM_FEMAP_ROWS    ! Number of rows that will be written to FEMAP arrays
      INTEGER(LONG), INTENT(INOUT)    :: NUM1              ! Cum rows written to OGEL prior to running this subr




! **********************************************************************************************************************************
! Calculate STRAIN for shell and solid elements

      IF      ((TYPE(1:3) == 'BAR'  ) .OR. (TYPE(1:4) == 'BUSH') .OR. (TYPE(1:4) == 'ELAS') .OR. (TYPE(1:3) == 'ROD'  )) THEN
         CALL ONE_D_STRAIN_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8') .OR.                                 &
               (TYPE(1:5) == 'SHEAR') .OR. (TYPE(1:6) == 'USERIN')) THEN
         CALL SHELL_STRAIN_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

      ELSE IF ((TYPE(1:4) == 'HEXA' ) .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         CALL SOLID_STRAIN_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9203) SUBR_NAME, TYPE
         WRITE(F06,9203) SUBR_NAME, TYPE
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
 9203 FORMAT(' *ERROR  9203: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INCORRECT ELEMENT TYPE = "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE CALC_ELEM_STRAINS


      SUBROUTINE ONE_D_STRAIN_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

! Calculates element strains for 1D ROD, BAR, ELAS, BUSH elements
! Calculates element specific strain output from array STRAIN (calc'd in subr ELEM_STRE_STRN_ARRAYS) for 1-D elements (ELAS, ROD,
! BAR, BUSH) and puts results into array OGEL for later output to F06 file.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  STRAIN, TYPE, ZS
      USE LINK9_STUFF, ONLY           :  MSPRNT, OGEL
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_VECS
      USE PARAMS, ONLY                :  PRTNEU

      USE ROD_BAR_OUTPUT, ONLY        :  ROD_MARGIN
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ONE_D_STRAIN_OUTPUTS'
      CHARACTER(1*BYTE)               :: MSP1,MSP2,MSP3     ! Output from subr MARGIN used to control how margins of safety printed
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_OGEL         ! If 'Y' then write data to array OGEL
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_FEMAP        ! If 'Y' then write data to array FEMAP_EL_VECS

      INTEGER(LONG), INTENT(IN)       :: SIZE_ALLOCATED     ! No. of rows allocated to array that will be written to
      INTEGER(LONG), INTENT(IN)       :: NUM_FEMAP_ROWS     ! Number of rows that will be written to FEMAP arrays
      INTEGER(LONG), INTENT(INOUT)    :: NUM1               ! Cum rows written to OGEL prior to running this subr
      INTEGER(LONG)                   :: ICOL               ! An input to subr MARGIN, called by this subr


      REAL(DOUBLE)                    :: C1,C2              ! Coords of point "C" on cross-section of a CBAR where strain is calc'd
      REAL(DOUBLE)                    :: D1,D2              ! Coords of point "D" on cross-section of a CBAR where strain is calc'd
      REAL(DOUBLE)                    :: E1,E2              ! Coords of point "E" on cross-section of a CBAR where strain is calc'd
      REAL(DOUBLE)                    :: F1,F2              ! Coords of point "F" on cross-section of a CBAR where strain is calc'd
      REAL(DOUBLE)                    :: MS1,MS2,MS3        ! Margins of safety calc'd by subr MARGIN
      REAL(DOUBLE)                    :: SA1,SA2,SA3,SA4    ! Strains at points "C", "D", "E" and "F" at end-a of a CBAR
      REAL(DOUBLE)                    :: SB1,SB2,SB3,SB4    ! Strains at points "C", "D", "E" and "F" at end-b of a CBAR
      REAL(DOUBLE)                    :: SA,ST              ! Strains input to subrs that calc margins of safety
      REAL(DOUBLE)                    :: SAMAX,SAMIN        ! Max & min strains from points "C", "D", "E" and "F" at end-a of CBAR
      REAL(DOUBLE)                    :: SBMAX,SBMIN        ! Max & min strains from points "C", "D", "E" and "F" at end-b of CBAR
      LOGICAL                         :: WRITE_NEU

      INTRINSIC DMAX1,DMIN1


      WRITE_NEU = (PRTNEU == 'Y')
! **********************************************************************************************************************************
! Calc engineering strains from array STRAIN and put into array OGEL

      IF      (TYPE == 'BAR     ') THEN                    ! BAR1 elements
         ! TODO: not validated
         C1 = ZS(1)
         C2 = ZS(2)
         D1 = ZS(3)
         D2 = ZS(4)
         E1 = ZS(5)
         E2 = ZS(6)
         F1 = ZS(7)
         F2 = ZS(8)
         SA    =   STRAIN(1)
         SA1   = -(C1*STRAIN(2) + C2*STRAIN(3))
         SA2   = -(D1*STRAIN(2) + D2*STRAIN(3))
         SA3   = -(E1*STRAIN(2) + E2*STRAIN(3))
         SA4   = -(F1*STRAIN(2) + F2*STRAIN(3))
         SB1   = -(C1*STRAIN(4) + C2*STRAIN(5))
         SB2   = -(D1*STRAIN(4) + D2*STRAIN(5))
         SB3   = -(E1*STRAIN(4) + E2*STRAIN(5))
         SB4   = -(F1*STRAIN(4) + F2*STRAIN(5))
         ST    =   STRAIN(6)
         SAMAX = SA + DMAX1(SA1,SA2,SA3,SA4)
         SBMAX = SA + DMAX1(SB1,SB2,SB3,SB4)
         SAMIN = SA + DMIN1(SA1,SA2,SA3,SA4)
         SBMIN = SA + DMIN1(SB1,SB2,SB3,SB4)
         ICOL = 1
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
         IF (WRITE_OGEL == 'Y') THEN
            NUM1 = NUM1 + 1
            IF (NUM1 > SIZE_ALLOCATED) THEN
               WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
               WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            OGEL(NUM1,1) = STRAIN(1)
            OGEL(NUM1,2) = STRAIN(2)
            OGEL(NUM1,3) = STRAIN(3)
            OGEL(NUM1,4) = STRAIN(4)
            OGEL(NUM1,5) = STRAIN(5)
            OGEL(NUM1,6) = STRAIN(6)
         ENDIF
         IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,1) = STRAIN(1)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,2) = STRAIN(2)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,3) = STRAIN(3)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,4) = STRAIN(4)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,5) = STRAIN(5)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,6) = STRAIN(6)
         ENDIF
      ELSE IF (TYPE(1:4) == 'ELAS') THEN                   ! Strains for ELAS elements
         IF (WRITE_OGEL == 'Y') THEN
            NUM1 = NUM1 + 1
            IF (NUM1 > SIZE_ALLOCATED) THEN
               WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
               WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            OGEL(NUM1,1) = STRAIN(1)
         ENDIF
         IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,1) = STRAIN(1)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,2) = STRAIN(1)
         ENDIF

      ELSE IF (TYPE == 'ROD     ') THEN                    ! ROD1 elements
         ! TODO: not validated
         SA = STRAIN(1)
         ST = STRAIN(2)
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

! **********************************************************************************************************************************
 9200 FORMAT(' *ERROR  9200: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ARRAY OGEL WAS ALLOCATED TO HAVE ',I12,' ROWS. ATTEMPT TO WRITE TO OGEL BEYOND THIS')

! **********************************************************************************************************************************
 9203 FORMAT(' *ERROR  9203: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INCORRECT ELEMENT TYPE = "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE ONE_D_STRAIN_OUTPUTS


      SUBROUTINE SHELL_STRAIN_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

! Calculates element specific strain output from array STRAIN (calc'd in subr ELEM_STRE_STRN_ARRAYS) for shell elements (TRIA3,
! QUAD4, SHEAR) and puts results into array OGEL for later output to F06 file.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  ANY_FAILURE_THEORY, FAILURE_THEORY, PCOMP_PROPS, STRAIN, STRESS, TYPE, ZS
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_OPT, STRN_CUR
      USE LINK9_STUFF, ONLY           :  FTNAME, OGEL
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_VECS
      USE PARAMS, ONLY                :  PRTNEU

      USE PRINCIPAL_VALUES, ONLY      :  PRINCIPAL_STRAIN_2D
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: GET_COMP_SHELL_ALLOWS
      USE COMPOSITE_FAILURE, ONLY     :  INDEP_FAILURE_INDEX, POLY_FAILURE_INDEX

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SHELL_STRAIN_OUTPUTS'
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_OGEL         ! If 'Y' then write data to array OGEL
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_FEMAP        ! If 'Y' then write data to array FEMAP_EL_VECS

      INTEGER(LONG), INTENT(IN)       :: SIZE_ALLOCATED     ! No. of rows allocated to array that will be written to
      INTEGER(LONG), INTENT(IN)       :: NUM_FEMAP_ROWS     ! Number of rows that will be written to FEMAP arrays
      INTEGER(LONG), INTENT(INOUT)    :: NUM1               ! Cum rows written to OGEL prior to running this subr
      INTEGER(LONG)                   :: I,J                ! DO loop indices
      INTEGER(LONG)                   :: NUM_ROWS           ! Number of rows of stress for an element (plates have 2 ZS vals)


      REAL(DOUBLE)                    :: ANGLE              ! Angle of prin strains in plate elems
      REAL(DOUBLE)                    :: FAILURE_INDEX      ! Failure index (scalar value)
      REAL(DOUBLE)                    :: MEAN               ! Mean strains
      REAL(DOUBLE)                    :: STREi(6)           ! 6 components of stress
      REAL(DOUBLE)                    :: STRNi(6)           ! 6 components of strain
      REAL(DOUBLE)                    :: SMAJ,SMIN          ! Major/minor prin strains in plate elems
      REAL(DOUBLE)                    :: STRE_ALLOWABLES(9) ! Allowable strains (incl tension and compr for normal strains)
      REAL(DOUBLE)                    :: STRN_ALLOWABLES(9) ! Allowable strains  (incl tension and compr for normal strains)
      REAL(DOUBLE)                    :: SX,SY,SXY          ! In-plane strains in plate elements
      REAL(DOUBLE)                    :: SXZ,SYZ            ! Transverse shear strains in plate elements
      REAL(DOUBLE)                    :: SXYMAX             ! Max shear strain in plate elems
      REAL(DOUBLE)                    :: VONMISES           ! von Mises strain
      REAL(DOUBLE)                    :: FIBER_Z            ! Value for the fiber distance or strain curvature column.
      LOGICAL                         :: WRITE_NEU

      INTRINSIC DMAX1,DMIN1


      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Calculates strain output for shell elements (TRIA3, QUAD4, SHEAR) and puts results into array OGEL for later output to F06 file

      IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8') .OR. (TYPE(1:5) == 'SHEAR')) THEN

         IF (PCOMP_PROPS == 'Y') THEN

            SX  = STRAIN(1)
            SY  = STRAIN(2)

            SXY = STRAIN(3)
            SXZ = STRAIN(7)
            SYZ = STRAIN(8)
            CALL PRINCIPAL_STRAIN_2D ( SX, SY, SXY, ANGLE, SMAJ, SMIN, SXYMAX, MEAN, VONMISES )
            IF (WRITE_OGEL == 'Y') THEN
               NUM1 = NUM1 + 1
               IF (NUM1 > SIZE_ALLOCATED) THEN
                  WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
                  WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
               OGEL(NUM1,1) = SX
               OGEL(NUM1,2) = SY
               OGEL(NUM1,3) = SXY
               OGEL(NUM1,4) = SXZ
               OGEL(NUM1,5) = SYZ
               OGEL(NUM1,6) = ANGLE
               OGEL(NUM1,7) = SMAJ
               OGEL(NUM1,8) = SMIN
               IF (STRN_OPT == 'VONMISES') THEN
                  OGEL(NUM1,9) = VONMISES
               ELSE
                  OGEL(NUM1,9) = SXYMAX
               ENDIF
               CALL GET_COMP_SHELL_ALLOWS ( STRE_ALLOWABLES, STRN_ALLOWABLES )

               FAILURE_INDEX = ZERO
               IF (FAILURE_THEORY == 'NONE') THEN

                  FTNAME(NUM1)  = 'None'

               ELSE

                  STREi(1) = STRESS(1)                      ! STREi stress order for failure theory calc has #1 = x  direct stress
                  STREi(2) = STRESS(2)                      ! STREi stress order for failure theory calc has #2 = y  direct stress
                  STREi(3) = ZERO                           ! STREi stress order for failure theory calc has #3 = z  direct stress
                  STREi(4) = STRESS(8)                      ! STREi stress order for failure theory calc has #4 = yz shear  stress
                  STREi(5) = STRESS(7)                      ! STREi stress order for failure theory calc has #5 = zx shear  stress
                  STREi(6) = STRESS(3)                      ! STREi stress order for failure theory calc has #6 = xy shear  stress

                  STRNi(1) = STRAIN(1)                      ! STRNi strain order for failure theory calc has #1 = x  direct strain
                  STRNi(2) = STRAIN(2)                      ! STRNi strain order for failure theory calc has #2 = y  direct strain
                  STRNi(3) = ZERO                           ! STRNi strain order for failure theory calc has #3 = z  direct strain
                  STRNi(4) = STRAIN(8)                      ! STRNi strain order for failure theory calc has #4 = yz shear  strain
                  STRNi(5) = STRAIN(7)                      ! STRNi strain order for failure theory calc has #5 = zx shear  strain
                  STRNi(6) = STRAIN(3)                      ! STRNi strain order for failure theory calc has #6 = xy shear  strain

                  ANY_FAILURE_THEORY = 'Y'

                  IF ((FAILURE_THEORY == 'HILL') .OR. (FAILURE_THEORY == 'HOFF') .OR. (FAILURE_THEORY == 'TSAI')) THEN
                     CALL POLY_FAILURE_INDEX  ( STREi, STRE_ALLOWABLES, FAILURE_INDEX )
                  ELSE IF ((FAILURE_THEORY == 'STRE') .OR. (FAILURE_THEORY == 'STRN')) THEN
                     CALL INDEP_FAILURE_INDEX ( STREi, STRNi, STRE_ALLOWABLES, STRN_ALLOWABLES, FAILURE_INDEX )
                  ELSE
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,9205) '', FAILURE_THEORY, 'HILL, HOFF or TSAI'
                     WRITE(F06,9205) '', FAILURE_THEORY, 'HILL, HOFF or TSAI'
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF

                  OGEL(NUM1,10) = FAILURE_INDEX
                  FTNAME(NUM1)  = FAILURE_THEORY(1:4)

               ENDIF

               IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,1) = SX
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,2) = SY
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,3) = SXY
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,4) = SXZ
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,5) = SYZ
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,6) = ANGLE
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,7) = SMAJ
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,8) = SMIN
                  IF (STRN_OPT == 'VONMISES') THEN
                     FEMAP_EL_VECS(NUM_FEMAP_ROWS,9) = VONMISES
                  ELSE
                     FEMAP_EL_VECS(NUM_FEMAP_ROWS,9) = SXYMAX
                  ENDIF
               ENDIF

            ENDIF

         ELSE                                              ! Not PCOMP

            NUM_ROWS = 2
            IF (TYPE(1:5) == 'SHEAR') THEN
               NUM_ROWS = 1
            ENDIF

            DO I=1,NUM_ROWS
               IF      (STRN_CUR == 'FIBER') THEN
                  FIBER_Z = ZS(I)
                  SX  = STRAIN(1) + ZS(I)*STRAIN(4)
                  SY  = STRAIN(2) + ZS(I)*STRAIN(5)
                  SXY = STRAIN(3) + ZS(I)*STRAIN(6)
               ELSE
                  FIBER_Z = 1 - I                          ! 0 for membrane strain row and -1 for curvature row.
                  SX  = (FIBER_Z + 1) * STRAIN(1) + FIBER_Z * STRAIN(4)
                  SY  = (FIBER_Z + 1) * STRAIN(2) + FIBER_Z * STRAIN(5)
                  SXY = (FIBER_Z + 1) * STRAIN(3) + FIBER_Z * STRAIN(6)
               ENDIF
               SXZ = STRAIN(7)
               SYZ = STRAIN(8)
               CALL PRINCIPAL_STRAIN_2D ( SX, SY, SXY, ANGLE, SMAJ, SMIN, SXYMAX, MEAN, VONMISES )
               IF (WRITE_OGEL == 'Y') THEN
                  NUM1 = NUM1 + 1
                  IF (NUM1 > SIZE_ALLOCATED) THEN
                     WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
                     WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
                     FATAL_ERR = FATAL_ERR + 1
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
                  IF (TYPE(1:5) == 'SHEAR') THEN
                     OGEL(NUM1, 1) = SX
                     OGEL(NUM1, 2) = SY
                     OGEL(NUM1, 3) = SXY
                     DO J=4,7
                        OGEL(NUM1,J) = ZERO
                     ENDDO
                  ELSE
                     OGEL(NUM1, 1) = FIBER_Z
                     OGEL(NUM1, 2) = SX
                     OGEL(NUM1, 3) = SY
                     OGEL(NUM1, 4) = SXY
                     OGEL(NUM1, 5) = ANGLE
                     OGEL(NUM1, 6) = SMAJ
                     OGEL(NUM1, 7) = SMIN
                  ENDIF
                  IF (STRN_OPT == 'VONMISES') THEN
                     OGEL(NUM1, 8) = VONMISES
                  ELSE
                     OGEL(NUM1, 8) = SXYMAX
                  ENDIF
                  OGEL(NUM1, 9) = SXZ
                  OGEL(NUM1,10) = SYZ
               ENDIF
               IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 1) = SX
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 2) = SY
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 3) = SXY
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 4) = SMAJ
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 5) = SMIN
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 6) = ANGLE
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 7) = MEAN
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 8) = SXYMAX
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 9) = VONMISES
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+10) = SXZ
!                 FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+11) = SYZ

                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 1) = SX
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 2) = SY
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 3) = SXY
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 4) = SMAJ
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 5) = SMIN
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 6) = ANGLE
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 7) = MEAN
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 8) = SXYMAX
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+ 9) = VONMISES
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+10) = SXZ
                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,11*(I-1)+11) = SYZ

!xx               IF (STRN_OPT == 'VONMISES') THEN
!xx                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,9*(I-1)+8) = VONMISES
!xx               ELSE
!xx                  FEMAP_EL_VECS(NUM_FEMAP_ROWS,9*(I-1)+8) = SXYMAX
!xx               ENDIF
               ENDIF
            ENDDO

         ENDIF

      ELSE IF (TYPE(1:6) == 'USERIN') THEN                 ! Stresses for USERIN elements

         IF (WRITE_OGEL == 'Y') THEN
            NUM1 = NUM1 + 1
            IF (NUM1 > SIZE_ALLOCATED) THEN
               WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
               WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )
            ENDIF
            DO J=1,9
               OGEL(NUM1,J) = STRAIN(J)
            ENDDO
         ENDIF
         IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
            DO J=1,9
               FEMAP_EL_VECS(NUM_FEMAP_ROWS,J) = STRAIN(J)
            ENDDO
         ENDIF

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9203) SUBR_NAME, TYPE
         WRITE(F06,9203) SUBR_NAME, TYPE
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
 9200 FORMAT(' *ERROR  9200: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ARRAY OGEL WAS ALLOCATED TO HAVE ',I12,' ROWS. ATTEMPT TO WRITE TO OGEL BEYOND THIS')

 9203 FORMAT(' *ERROR  9203: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INCORRECT ELEMENT TYPE = "',A,'"')

 9205 FORMAT(' *ERROR  9205: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INVALID ',A,' FAILURE THEORY = ',A,'. VALID ONES ARE: ',A)


! **********************************************************************************************************************************

      END SUBROUTINE SHELL_STRAIN_OUTPUTS


      SUBROUTINE SOLID_STRAIN_OUTPUTS ( SIZE_ALLOCATED, NUM1, NUM_FEMAP_ROWS, WRITE_OGEL, WRITE_FEMAP )

! Calculates element specific strain outputs from array STRESS (generated in subr ELEM_STRE_STRN_ARRAYS) for solid elements

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  STRAIN, TYPE
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_OPT
      USE LINK9_STUFF, ONLY           :  OGEL
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_VECS
      USE PARAMS, ONLY                :  PRTNEU

      USE PRINCIPAL_VALUES, ONLY      :  PRINCIPAL_3D
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SOLID_STRAIN_OUTPUTS'
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_OGEL         ! If 'Y' then write data to array OGEL
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRITE_FEMAP        ! If 'Y' then write data to array FEMAP_EL_VECS

      INTEGER(LONG), INTENT(IN)       :: SIZE_ALLOCATED     ! No. of rows allocated to array that will be written to
      INTEGER(LONG), INTENT(IN)       :: NUM_FEMAP_ROWS     ! Number of rows that will be written to FEMAP arrays
      INTEGER(LONG), INTENT(INOUT)    :: NUM1               ! Cum rows written to OGEL prior to running this subr
      INTEGER(LONG)                   :: J                  ! DO loop indices


      REAL(DOUBLE)                    :: MEAN               ! Mean strains
      REAL(DOUBLE)                    :: PRINCIPAL_STRAIN(3)! Principal strains
      REAL(DOUBLE)                    :: SIG_OCT            ! Octrahedral normal strain
      REAL(DOUBLE)                    :: TAU_OCT            ! Octrahedral shear  strain
      REAL(DOUBLE)                    :: VONMISES           ! von Mises strain
      LOGICAL                         :: WRITE_NEU

      INTRINSIC DMAX1,DMIN1


      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Calc engineering strains from array STRAIN and put into array OGEL

      IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         CALL PRINCIPAL_3D ( STRAIN, PRINCIPAL_STRAIN, MEAN, VONMISES, SIG_OCT, TAU_OCT )
         IF (WRITE_OGEL == 'Y') THEN
            NUM1 = NUM1 + 1
            IF (NUM1 > SIZE_ALLOCATED) THEN
               WRITE(ERR,9200) SUBR_NAME,SIZE_ALLOCATED
               WRITE(F06,9200) SUBR_NAME,SIZE_ALLOCATED
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Coding error (dim of array OGEL too small), so quit
            ENDIF
            DO J=1,6
               OGEL(NUM1,J) = STRAIN(J)
            ENDDO
            IF (STRN_OPT == 'VONMISES') THEN
               OGEL(NUM1,7) = VONMISES
               OGEL(NUM1,8) = ZERO
            ELSE
               OGEL(NUM1,7) = SIG_OCT
               OGEL(NUM1,8) = TAU_OCT
            ENDIF
            OGEL(NUM1,9) = PRINCIPAL_STRAIN(1)
            OGEL(NUM1,10) = PRINCIPAL_STRAIN(2)
            OGEL(NUM1,11) = PRINCIPAL_STRAIN(3)
            OGEL(NUM1,12) = MEAN
         ENDIF
         IF (WRITE_NEU .AND. (WRITE_FEMAP == 'Y')) THEN
            DO J=1,6
               FEMAP_EL_VECS(NUM_FEMAP_ROWS,J) = STRAIN(J)
            ENDDO
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 7) = PRINCIPAL_STRAIN(1)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 8) = PRINCIPAL_STRAIN(2)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS, 9) = PRINCIPAL_STRAIN(3)
            FEMAP_EL_VECS(NUM_FEMAP_ROWS,10) = MEAN
            IF (STRN_OPT == 'VONMISES') THEN
               FEMAP_EL_VECS(NUM_FEMAP_ROWS,11) = VONMISES
               FEMAP_EL_VECS(NUM_FEMAP_ROWS,12) = ZERO
            ELSE
               FEMAP_EL_VECS(NUM_FEMAP_ROWS,11) = SIG_OCT
               FEMAP_EL_VECS(NUM_FEMAP_ROWS,12) = TAU_OCT
            ENDIF
         ENDIF

      ELSE

         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,9203) SUBR_NAME, TYPE
         WRITE(F06,9203) SUBR_NAME, TYPE
         CALL OUTA_HERE ( 'Y' )

      ENDIF



      RETURN

! **********************************************************************************************************************************
 9200 FORMAT(' *ERROR  9200: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ARRAY OGEL WAS ALLOCATED TO HAVE ',I12,' ROWS. ATTEMPT TO WRITE TO OGEL BEYOND THIS')

 9203 FORMAT(' *ERROR  9203: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' INCORRECT ELEMENT TYPE = "',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE SOLID_STRAIN_OUTPUTS

   END MODULE ELEMENT_STRAIN_RECOVERY
