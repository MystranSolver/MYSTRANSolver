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

   MODULE MODAL_OUTPUT_WRITERS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: WRITE_MEFFMASS, WRITE_MPFACTOR, WRITE_SUBCASE_EIGENVEC_HEADER

   CONTAINS

      SUBROUTINE WRITE_MEFFMASS

      ! Writes output for modal effective mass
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NVEC
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE, TWO, ONE_HUNDRED, PI
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, MEFFMASS
      USE MODEL_STUF, ONLY            :  MEFM_RB_MASS, LABEL, STITLE, TITLE
      USE PARAMS, ONLY                :  EPSIL, GRDPNT, MEFMCORD, MEFMGRID, MEFMLOC, SUPINFO, WTMASS

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_MEFFMASS'
      CHARACTER(14*BYTE)              :: CHAR_PCT(6)       ! Character representation of MEFFMASS sum percents of total model mass
      CHARACTER(1*BYTE)               :: IHDR   = 'Y'      ! Indicator of whether to write an output header

      INTEGER(LONG)                   :: I,J               ! DO loop indices


      REAL(DOUBLE)                    :: CYCLES            ! Circular frequency of a mode
      REAL(DOUBLE)                    :: EPS1              ! Small number to compare against zero
      REAL(DOUBLE)                    :: MEFM_TOTALS(6)    ! Totals for the 6 modal effective masses over all modes
      REAL(DOUBLE)                    :: MODES_PCT(6)      ! Modal mass as % of total mass
      !LOGICAL                        :: WRITE_F06  ! flag
      !LOGICAL                        :: WRITE_OP2  ! flag
      LOGICAL                         :: IS_LOW_PRECISION  ! Print MPFACTOR, MEFFMASS values with 2 decimal places of accuracy rather than 6




! **********************************************************************************************************************************
      IS_LOW_PRECISION = (DEBUG(174) == 0)
      !--------------------------------------------------

      EPS1 = EPSIL(1)

      ! Write output headers.
      IF (IHDR == 'Y') THEN
         WRITE(F06,900)
         ! There is always a TITLE(1), etc (even if they are blank)
         WRITE(F06,909) TITLE(1)
         WRITE(F06,909) STITLE(1)
         WRITE(F06,909) LABEL(1)
         WRITE(F06,*)
      ENDIF

      ! Write modal effective masses
      WRITE(F06,9102) MEFMCORD

      IF      (MEFMLOC == 'GRDPNT') THEN
         IF (MEFMGRID == 0) THEN
            WRITE(F06,9103)
         ELSE
            WRITE(F06,9104) GRDPNT
         ENDIF
      ELSE IF (MEFMLOC == 'CG    ') THEN
         WRITE(F06,9105)
      ELSE IF (MEFMLOC == 'GRID  ') THEN
         WRITE(F06,9106) MEFMGRID
      ENDIF

      IF (IS_LOW_PRECISION) THEN
         WRITE(F06,9107)
      ELSE
         WRITE(F06,9108)
      ENDIF

      DO J=1,6
         MEFM_TOTALS(J) = ZERO
      ENDDO

      DO I=1,NVEC
         CYCLES = DSQRT(DABS(EIGEN_VAL(I)))/(TWO*PI)

         IF (IS_LOW_PRECISION) THEN ! 6 digits
            WRITE(F06,9110) I, CYCLES, (MEFFMASS(I,J)/WTMASS,J=1,6)
         ELSE ! low precision (2 digits)
            WRITE(F06,9111) I, CYCLES, (MEFFMASS(I,J)/WTMASS,J=1,6)
         ENDIF

         DO J=1,6
            MEFM_TOTALS(J) = MEFM_TOTALS(J) + MEFFMASS(I,J)/WTMASS
         ENDDO

      ENDDO

      IF (IS_LOW_PRECISION) THEN
         WRITE(F06,9112) (MEFM_TOTALS(J),J=1,6)
      ELSE
         WRITE(F06,9113) (MEFM_TOTALS(J),J=1,6)
      ENDIF
                                                           ! MEFM_RB_MASS is in the same units as in the DAT file
      IF (IS_LOW_PRECISION) THEN
         WRITE(F06,9116) (MEFM_RB_MASS(I,I),I=1,6)
      ELSE
         WRITE(F06,9117) (MEFM_RB_MASS(I,I),I=1,6)
      ENDIF

      ! For each of the 6 modal masses, calc % of total mass.
      ! A character variable is used to store the % so that blank percentages can
      ! be printed if zero modal mass exists for a component (T1 - R3) or
      ! if a denominator in the % expression is zero
      IF (DABS(MEFM_RB_MASS(1,1)) > EPS1) THEN
         MODES_PCT(1) = ONE_HUNDRED*MEFM_TOTALS(1)/MEFM_RB_MASS(1,1)
         WRITE(CHAR_PCT(1),1001) MODES_PCT(1), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(1)(1:) = ' '
         IF (MEFM_TOTALS(1) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2001)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2001)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(2,2)) > EPS1) THEN
         MODES_PCT(2) = ONE_HUNDRED*MEFM_TOTALS(2)/MEFM_RB_MASS(2,2)
         WRITE(CHAR_PCT(2),1001) MODES_PCT(2), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(2)(1:) = ' '
         IF (MEFM_TOTALS(2) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2002)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2002)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(3,3)) > EPS1) THEN
         MODES_PCT(3) = ONE_HUNDRED*MEFM_TOTALS(3)/MEFM_RB_MASS(3,3)
         WRITE(CHAR_PCT(3),1001) MODES_PCT(3), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(3)(1:) = ' '
         IF (MEFM_TOTALS(3) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2003)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2003)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(4,4)) > EPS1) THEN
         MODES_PCT(4) = ONE_HUNDRED*MEFM_TOTALS(4)/MEFM_RB_MASS(4,4)
         WRITE(CHAR_PCT(4),1001) MODES_PCT(4), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(4)(1:) = ' '
         IF (MEFM_TOTALS(4) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2004)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2004)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(5,5)) > EPS1) THEN
         MODES_PCT(5) = ONE_HUNDRED*MEFM_TOTALS(5)/MEFM_RB_MASS(5,5)
         WRITE(CHAR_PCT(5),1001) MODES_PCT(5), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(5)(1:) = ' '
         IF (MEFM_TOTALS(5) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2005)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2005)
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(MEFM_RB_MASS(6,6)) > EPS1) THEN
         MODES_PCT(6) = ONE_HUNDRED*MEFM_TOTALS(6)/MEFM_RB_MASS(6,6)
         WRITE(CHAR_PCT(6),1001) MODES_PCT(6), '%'
      ELSE
         ! Denominator is zero so leave % blank
         CHAR_PCT(6)(1:) = ' '
         IF (MEFM_TOTALS(6) > EPS1) THEN
            ! Error: denominator is zero but numerator is not, so write message
            WRITE(ERR,2006)
            IF (SUPINFO == 'N') THEN
               WRITE(F06,2006)
            ENDIF
         ENDIF
      ENDIF

      WRITE(F06,9118) (CHAR_PCT(I),I=1,6)



      RETURN

! **********************************************************************************************************************************
  900 FORMAT('--------------------------------------------------------------------------------------------------------------------'&
            ,'----------------')

  909 FORMAT(1X,A)

 1001 FORMAT(F13.2,A)

 2001 FORMAT(' *INFORMATION: CANNOT CALCULATE T1 MODAL EFFECTIVE MASS PERCENT OF MODEL MASS SO IT IS LEFT BLANK')

 2002 FORMAT(' *INFORMATION: CANNOT CALCULATE T2 MODAL EFFECTIVE MASS PERCENT OF MODEL MASS SO IT IS LEFT BLANK')

 2003 FORMAT(' *INFORMATION: CANNOT CALCULATE T3 MODAL EFFECTIVE MASS PERCENT OF MODEL MASS SO IT IS LEFT BLANK')

 2004 FORMAT(' *INFORMATION: CANNOT CALCULATE R1 MODAL EFFECTIVE MASS PERCENT OF MODEL IXX  SO IT IS LEFT BLANK')

 2005 FORMAT(' *INFORMATION: CANNOT CALCULATE R2 MODAL EFFECTIVE MASS PERCENT OF MODEL IYY  SO IT IS LEFT BLANK')

 2006 FORMAT(' *INFORMATION: CANNOT CALCULATE R3 MODAL EFFECTIVE MASS PERCENT OF MODEL IZZ  SO IT IS LEFT BLANK')

 9102 FORMAT(14X,'                    E F F E C T I V E   M O D A L   M A S S E S   O R   W E I G H T S',/,                        &
             14X,'                                     (in coordinate system ',I8,')',/,                                           &
             14X,'                       Units are same as units for mass input in the Bulk Data Deck')

 9103 FORMAT(14X,'                          Reference point is the basic coordinate system origin',/)

 9104 FORMAT(14X,'                            Reference point is the PARAM GRDPNT grid: ',I8,/)

 9105 FORMAT(14X,'                              Reference point is the model center of gravity',/)

 9106 FORMAT(14X,'                                    Reference point is grid ',I8,/)

 9107 FORMAT(13X,'MODE     CYCLES          T1            T2            T3            R1            R2            R3',/,            &
             13X,' NUM')

 9108 FORMAT(13X,'MODE       CYCLES          T1            T2            T3            R1            R2            R3',/,          &
             13X,' NUM')

 9110 FORMAT(9X,I8,7(1ES14.6))

 9111 FORMAT(9X,I8,7(1ES14.2))

 9112 FORMAT(32X,' ------------  ------------  ------------  ------------  ------------  ------------',/,                          &
             17X,'Sum all modes:',6(1ES14.6))

 9113 FORMAT(32X,'     --------      --------      --------      --------      --------      --------',/,                          &
             17X,'Sum all modes:',6(1ES14.2))

 9116 FORMAT(14X,'Total model mass:',6(1ES14.6))

 9117 FORMAT(14X,'Total model mass:',6(1ES14.2))

 9118 FORMAT(8X,'Modes % of total mass*:',6(A14),//,' *If all modes are calculated the % of total mass should be 100% of the '     &
               ,'free mass (i.e. not counting mass at constrained DOF''s).',/,                                                     &
               '  Percentages are only printed for components that have finite model mass.',/,                                     &
               '                                                               -----')

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_MEFFMASS


      SUBROUTINE WRITE_MPFACTOR                ! ( IHDR )

      ! Writes output for modal participation factors
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NDOFG, NDOFR, NVEC, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, PI
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE EIGEN_MATRICES_1, ONLY      :  EIGEN_VAL, MPFACTOR_NR, MPFACTOR_N6
      USE MODEL_STUF, ONLY            :  LABEL, STITLE, TITLE
      USE PARAMS, ONLY                :  GRDPNT, MEFMCORD, MEFMGRID, MEFMLOC, MPFOUT
      USE DOF_TABLES, ONLY            :  TDOFI

      USE DOF_NUMBERING, ONLY         :  TDOF_COL_NUM

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_MPFACTOR'
!xx   CHARACTER(LEN=*) , INTENT(IN)   :: IHDR              ! Indicator of whether to write an output header
      CHARACTER(1*BYTE)               :: IHDR   = 'Y'      ! Indicator of whether to write an output header

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: K                 ! Counter
      INTEGER(LONG)                   :: R_SET_GRIDS(NDOFR)! Array of grids for the R-set
      INTEGER(LONG)                   :: R_SET_COMPS(NDOFR)! Array of displ components for the R-set
      INTEGER(LONG)                   :: R_SET_COL         ! Col in TDOFI array where R-set exists


      REAL(DOUBLE)                    :: CYCLES            ! Circular frequency of a mode
      !LOGICAL                        :: WRITE_F06  ! flag
      !LOGICAL                        :: WRITE_OP2  ! flag
      LOGICAL                         :: IS_LOW_PRECISION  ! Print MPFACTOR, MEFFMASS values with 2 decimal places of accuracy rather than 6



! **********************************************************************************************************************************
      IS_LOW_PRECISION = (DEBUG(174) == 0)
      !--------------------------------------------------

      CALL TDOF_COL_NUM ( 'R ', R_SET_COL )
      K = 0
      DO I=1,NDOFG
         IF (TDOFI(I,R_SET_COL) /= 0) THEN
            K = K + 1
            R_SET_GRIDS(K) = TDOFI(I,1)
            R_SET_COMPS(K) = TDOFI(I,2)
         ENDIF
      ENDDO

      WRITE(F06,*)
      ! Write output headers.
      IF (IHDR == 'Y') THEN
         WRITE(F06,9000)
         ! There is always a TITLE(1), etc (even if they are blank)
         WRITE(F06,9003) TITLE(1)
         WRITE(F06,9003) STITLE(1)
         WRITE(F06,9003) LABEL(1)
         WRITE(F06,*)
      ENDIF
                                                           ! Write modal participation factors for CB analyses
      IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (MPFOUT == 'R')) THEN

         WRITE(F06,9004) MEFMCORD

         IF (IS_LOW_PRECISION) THEN
            WRITE(F06,9101) (I,I=1,NDOFR)
            WRITE(F06,9102) (R_SET_GRIDS(I), R_SET_COMPS(I),I=1,NDOFR)
            WRITE(F06,9103)
         ELSE
            WRITE(F06,9201) (I,I=1,NDOFR)
            WRITE(F06,9202) (R_SET_GRIDS(I), R_SET_COMPS(I),I=1,NDOFR)
            WRITE(F06,9203)
         ENDIF

         DO I=1,NVEC

            CYCLES = DSQRT(DABS(EIGEN_VAL(I)))/(TWO*PI)

            IF (IS_LOW_PRECISION) THEN
               WRITE(F06,9301) I, CYCLES, (MPFACTOR_NR(I,J),J=1,NDOFR)
            ELSE
               WRITE(F06,9302) I, CYCLES, (MPFACTOR_NR(I,J),J=1,NDOFR)
            ENDIF

         ENDDO

      ELSE

         WRITE(F06,9005) MEFMCORD
         IF      (MEFMLOC == 'GRDPNT') THEN
            IF (MEFMGRID == 0) THEN
               WRITE(F06,9006)
            ELSE
               WRITE(F06,9007) GRDPNT
            ENDIF
         ELSE IF (MEFMLOC == 'CG    ') THEN
            WRITE(F06,9008)
         ELSE IF (MEFMLOC == 'GRID  ') THEN
            WRITE(F06,9009) MEFMGRID
         ENDIF

         IF (IS_LOW_PRECISION) THEN
            WRITE(F06,9501)
         ELSE
            WRITE(F06,9502)
         ENDIF

         DO I=1,NVEC

            CYCLES = DSQRT(DABS(EIGEN_VAL(I)))/(TWO*PI)

            IF (IS_LOW_PRECISION) THEN
               WRITE(F06,9503) I, CYCLES, (MPFACTOR_N6(I,J),J=1,6)
            ELSE
               WRITE(F06,9504) I, CYCLES, (MPFACTOR_N6(I,J),J=1,6)
            ENDIF

         ENDDO

      ENDIF

      WRITE(F06,*)



      RETURN

! **********************************************************************************************************************************
 9000 FORMAT('--------------------------------------------------------------------------------------------------------------------'&
            ,'----------------')

 9003 FORMAT(1X,A)

 9004 FORMAT(13X,'                           M O D A L   P A R T I C I P A T I O N   F A C T O R S',/,                             &
             13X,'              (dimensionless, in coordinate sys ',I8,' with cols marked by R-set grid/comp)',/)

 9005 FORMAT(13X,'                           M O D A L   P A R T I C I P A T I O N   F A C T O R S',/,                             &
             13X,'                                (dimensionless, in coordinate sys ',I8,')')

 9006 FORMAT(14X,'                          Reference point is the basic coordinate system origin',/)

 9007 FORMAT(14X,'                            Reference point is the PARAM GRDPNT grid: ',I8,/)

 9008 FORMAT(14X,'                              Reference point is the model center of gravity',/)

 9009 FORMAT(14X,'                                    Reference point is grid ',I8,/)

 9101 FORMAT(32X,32767(I8,6X))

 9102 FORMAT(13X,'MODE     CYCLES  ',32767(2X,I8,'-',I1,2X))

 9103 FORMAT(13X,' NUM')

 9201 FORMAT(34X,32767(I8,6X))

 9202 FORMAT(13X,'MODE       CYCLES  ',32767(2X,I8,'-',I1,2X))

 9203 FORMAT(13X,' NUM')

 9301 FORMAT(9X,I8,32767(1ES14.6))

 9302 FORMAT(9X,I8,32767(1ES14.2))

 9501 FORMAT(13X,'MODE     CYCLES          T1            T2            T3            R1            R2            R3',/,            &
             13X,' NUM')

 9502 FORMAT(13X,'MODE       CYCLES          T1            T2            T3            R1            R2            R3',/,          &
             13X,' NUM')

 9503 FORMAT(9X,I8,7(1ES14.6))

 9504 FORMAT(9X,I8,7(1ES14.2))

! **********************************************************************************************************************************

      END SUBROUTINE WRITE_MPFACTOR


      SUBROUTINE WRITE_SUBCASE_EIGENVEC_HEADER ( JSUB, WRITE_F06 )

! Writes the complete per-vector block header to F06 for all LINK9 WRITE_* subroutines:
!   - Two blank separator lines
!   - "OUTPUT FOR SUBCASE x"         (all except GEN CB MODEL)
!   - "OUTPUT FOR EIGENVECTOR y"     (MODES and BUCKLING step 2 only)
!   - TITLE / SUBTITLE / LABEL lines  (each written only if non-blank)
!   - One trailing blank line
!
! This is the central single point for the per-vector F06 header emitted by all LINK9 WRITE_* subroutines.
! For GEN CB MODEL the SUBCASE/EIGENVECTOR lines are skipped; the caller writes the CB DOF line after returning.
!
! JSUB  = global solution-vector index (subcase number for STATICS, or global eigenvector index for MODES/BUCKLING)
!
! Module variables consumed (set by LINK9 before each call into the WRITE_* routines):
!   INT_SC_NUM  - owning internal subcase index (for TITLE/STITLE/LABEL lookup in the caller)
!   INT_EIG_NUM - per-subcase local eigenvector counter (1-based, reset per subcase); 0 for non-eigen solutions
!
! Craig-Bampton: the two blank lines are written here; the caller is responsible for the CB DOF line itself,
! since that requires grid/component lookup data not available here.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, INT_EIG_NUM, INT_SC_NUM, NDOFR, NUM_CB_DOFS, NVEC, SOL_NAME
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE MODEL_STUF, ONLY            :  LABEL, SCNUM, STITLE, TITLE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'WRITE_SUBCASE_EIGENVEC_HEADER'

      INTEGER(LONG), INTENT(IN)       :: JSUB              ! Global solution-vector index passed in by the caller
      LOGICAL,       INTENT(IN)       :: WRITE_F06         ! If .FALSE., suppress all F06 output (mirrors caller guard)

! **********************************************************************************************************************************

      IF (.NOT. WRITE_F06) RETURN

      WRITE(F06,*)
      WRITE(F06,*)

      IF    ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN

         WRITE(F06,9011) SCNUM(JSUB)

      ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 1)) THEN

         WRITE(F06,9011) SCNUM(JSUB)

      ELSE IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN

         WRITE(F06,9011) SCNUM(INT_SC_NUM)
         WRITE(F06,9012) INT_EIG_NUM

      ELSE IF (SOL_NAME(1:5) == 'MODES') THEN

         WRITE(F06,9011) SCNUM(INT_SC_NUM)
         WRITE(F06,9012) INT_EIG_NUM

      ! GEN CB MODEL: caller must write the CB DOF line -- just emit the blank lines (done above) and return
      ENDIF

      IF (TITLE(INT_SC_NUM)(1:)   /= ' ') WRITE(F06,9013) TITLE(INT_SC_NUM)
      IF (STITLE(INT_SC_NUM)(1:)  /= ' ') WRITE(F06,9013) STITLE(INT_SC_NUM)
      IF (LABEL(INT_SC_NUM)(1:)   /= ' ') WRITE(F06,9013) LABEL(INT_SC_NUM)
      WRITE(F06,*)

      RETURN

! **********************************************************************************************************************************
 9011 FORMAT(' OUTPUT FOR SUBCASE ',I8)
 9012 FORMAT(' OUTPUT FOR EIGENVECTOR ',I8)
 9013 FORMAT(1X,A)

      END SUBROUTINE WRITE_SUBCASE_EIGENVEC_HEADER

   END MODULE MODAL_OUTPUT_WRITERS
