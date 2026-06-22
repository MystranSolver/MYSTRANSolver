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

      SUBROUTINE MAXREQ_OGEL

! Count number of output requests to determine required leading dimension of array OGEL so memory can be allocated to it

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, IBIT, LSUB, NDOFG, NELE, NGRID, METYPE, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE MODEL_STUF, ONLY            :  ELMTYP, ELOUT, ESORT2, ETYPE, GROUT, MEFFMASS_CALC, MPFACTOR_CALC, NELGP, NUM_PLIES,      &
                                         PCOMP_PROPS, SCNUM, TYPE
      ! --- cbeam_add begin --- !
      USE MODEL_STUF, ONLY            :  EDAT, EPNT, PBEAM_NSTATIONS
      ! --- cbeam_add end --- !
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_LOC, STRE_LOC, FORC_LOC
      USE LINK9_STUFF, ONLY           :  MAXREQ
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE MAXREQ_OGEL_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'MAXREQ_OGEL'
      CHARACTER(25*BYTE)              :: GROUT_NAME(0:15)
      CHARACTER(25*BYTE)              :: ELOUT_NAME(0:15)
      CHARACTER(46*BYTE)              :: MPF_MEFM_MSG         ! Message written if MAXREQ impacted by req for output of MPF/MEFM
      CHARACTER( 1*BYTE)              :: SKIP_ELEM_TYPE

      INTEGER(LONG)                   :: NUMBER_ROWS(0:15)    ! For elem output, this accounts for more than 1 row written to OGEL
      INTEGER(LONG)                   :: I,J,K,L              ! DO loop indices
      INTEGER(LONG)                   :: IB                   ! Result of IAND to determine if there is an output request
      INTEGER(LONG)                   :: INT_ELEM_ID          ! Internal element number
      INTEGER(LONG)                   :: MAXGROUT_SC          ! Max no. grid pt output requests in array GROUT for any 1 subcase
      INTEGER(LONG)                   :: MAXELOUT_SC          ! Max no. elem output requests in array ELOUT for any 1 subcase
      INTEGER(LONG)                   :: MAXGROUT             ! Max of MAXGROUT_SC for all subcases
      INTEGER(LONG)                   :: MAXELOUT             ! Max of MAXELOUT_SC for all subcases
      INTEGER(LONG)                   :: NREQ_EL(METYPE,0:15) ! No. of requests in ELOUT for each bit of ELOUT for 1 subcase
      INTEGER(LONG)                   :: NREQ_GR(0:15)        ! No. of requests in GROUT for each bit of GROUT for 1 subcase
      ! --- cbeam_add begin --- !
      INTEGER(LONG)                   :: BEAM_FORCE_ROWS_SC   ! Conservative per-subcase count of beam engr-force rows
      INTEGER(LONG)                   :: BEAM_STRESS_ROWS_SC  ! Conservative per-subcase count of beam stress rows
      INTEGER(LONG)                   :: BEAM_STRAIN_ROWS_SC  ! Conservative per-subcase count of beam strain rows
      INTEGER(LONG)                   :: BEAM_NSTA            ! Conservative beam station count used for LINK9 allocations
      ! --- cbeam_add end --- !


      INTRINSIC                       :: IAND, MAX



! **********************************************************************************************************************************
! Initialize outputs

      MAXREQ = ZERO

      DO I=0,15
         GROUT_NAME(I) = '******** No Name ********'
         ELOUT_NAME(I) = '******** No Name ********'
      ENDDO

      GROUT_NAME(0) = 'grid displacement'
      GROUT_NAME(1) = 'grid applied load'
      GROUT_NAME(2) = 'grid SPC Force'
      GROUT_NAME(3) = 'grid MPC Force'
      GROUT_NAME(4) = 'grid point force balance'
      ELOUT_NAME(0) = 'element nodal force'
      ELOUT_NAME(1) = 'element engineering force'
      ELOUT_NAME(2) = 'element stress'
      ELOUT_NAME(3) = 'element strain'

! Count MAXGROUT, the max number of requests in GROUT.

      IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '11' )

      MAXGROUT = 0
      DO I=1,LSUB

         DO K=0,15
            NREQ_GR(K) = 0
         ENDDO

         DO J=1,NGRID
            DO K=0,15
               IB = IAND(GROUT(J,I),IBIT(K))
               IF (IB > 0) THEN
                  NREQ_GR(K) = NREQ_GR(K) + 1
               ENDIF
            ENDDO
         ENDDO

         MPF_MEFM_MSG(1:) = ' '
         IF ((MEFFMASS_CALC == 'Y') .OR. (MPFACTOR_CALC == 'Y')) THEN! Need to make sure MAXGROUT can cover MPF, MEFM
            IF (SOL_NAME /= 'GEN CB MODEL') THEN
               IF (NREQ_GR(2) < NDOFG) THEN
                  NREQ_GR(2) = NDOFG
                  MPF_MEFM_MSG = '(required for calculation of MPFACTOR/MEFMASS)'
               ENDIF
            ENDIF
         ENDIF

         IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '12' )

         MAXGROUT_SC = 0
         DO K=0,15
            IF (NREQ_GR(K) > 0) THEN
               IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '13' )
               IF (NREQ_GR(K) > MAXGROUT_SC) THEN
                  MAXGROUT_SC = NREQ_GR(K)
               ENDIF
               IF (NREQ_GR(K) > MAXGROUT) THEN
                  MAXGROUT = NREQ_GR(K)
               ENDIF
            ENDIF
         ENDDO
         IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '14' )
      ENDDO

      IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '15' )

! **********************************************************************************************************************************
! Count MAXELOUT, the max number of requests in ELOUT.

      DO I=0,15
         NUMBER_ROWS(I) = 0
      ENDDO

      MAXELOUT = 0
      DO I=1,LSUB

         DO L=1,METYPE
            DO K=0,15
               NREQ_EL(L,K) = 0
            ENDDO
         ENDDO

         DO L=1,METYPE
            DO J=1,NELE
               INT_ELEM_ID = ESORT2(J)
               DO K=0,15
                  TYPE = ETYPE(INT_ELEM_ID)
                  IF (TYPE == ELMTYP(L)) THEN
                     IB = IAND(ELOUT(INT_ELEM_ID,I),IBIT(K))
                     IF (IB > 0) THEN
                        CALL GET_NUM_ROWS_OUTPUT( L )      ! Determines the num rows of output for this element and request type
                        NREQ_EL(L,K) = NREQ_EL(L,K) + NUMBER_ROWS(K)
                     ENDIF
                  ENDIF
               ENDDO
            ENDDO
         ENDDO

         IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '21' )

         DO L=1,METYPE
            SKIP_ELEM_TYPE = 'Y'
            DO K=0,15
               IF (NREQ_EL(L,K) > 0) THEN
                  SKIP_ELEM_TYPE = 'N'
                  EXIT
               ENDIF
            ENDDO
            MAXELOUT_SC = 0
            IF (SKIP_ELEM_TYPE == 'N') THEN
               IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '22' )
               DO K=0,15
                  IF (NREQ_EL(L,K) > 0) THEN
                     IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '23' )
                     IF (NREQ_EL(L,K) > MAXELOUT_SC) THEN
                        MAXELOUT_SC = NREQ_EL(L,K)
                     ENDIF
                     IF (NREQ_EL(L,K) > MAXELOUT) THEN
                        MAXELOUT = NREQ_EL(L,K)
                     ENDIF
                  ENDIF
               ENDDO
            ENDIF
            IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '24' )
         ENDDO

      ENDDO

      IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '25' )
! --- cbeam_add begin --- !
! Conservative beam-row pass: some simple 2-station PBEAM definitions are stored with
! PBEAM_NSTATIONS=1 even though LINK9 beam force/stress/strain output writes both ends.
      DO I=1,LSUB
         BEAM_FORCE_ROWS_SC  = 0
         BEAM_STRESS_ROWS_SC = 0
         BEAM_STRAIN_ROWS_SC = 0
         DO J=1,NELE
            INT_ELEM_ID = ESORT2(J)
            TYPE = ETYPE(INT_ELEM_ID)
            IF (TYPE(1:4) /= 'BEAM') CYCLE
            BEAM_NSTA = PBEAM_NSTATIONS(EDAT(EPNT(INT_ELEM_ID)+1))
            IF (BEAM_NSTA <= 1) BEAM_NSTA = 2
            IF (IAND(ELOUT(INT_ELEM_ID,I),IBIT(1)) > 0) BEAM_FORCE_ROWS_SC  = BEAM_FORCE_ROWS_SC  + BEAM_NSTA
            IF (IAND(ELOUT(INT_ELEM_ID,I),IBIT(2)) > 0) BEAM_STRESS_ROWS_SC = BEAM_STRESS_ROWS_SC + 2*BEAM_NSTA
            IF (IAND(ELOUT(INT_ELEM_ID,I),IBIT(3)) > 0) BEAM_STRAIN_ROWS_SC = BEAM_STRAIN_ROWS_SC + 2*BEAM_NSTA
         ENDDO
         MAXELOUT = MAX(MAXELOUT, BEAM_FORCE_ROWS_SC)
         MAXELOUT = MAX(MAXELOUT, BEAM_STRESS_ROWS_SC)
         MAXELOUT = MAX(MAXELOUT, BEAM_STRAIN_ROWS_SC)
      ENDDO
! --- cbeam_add end --- !
      MAXREQ = MAX(MAXGROUT,MAXELOUT)

      IF (DEBUG(91) == 1) CALL MAXREQ_OGEL_DEB ( '31' )



      RETURN

! **********************************************************************************************************************************




! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_NUM_ROWS_OUTPUT ( LETYPE )

      USE MODEL_STUF, ONLY            :  NUM_SEi

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: LETYPE            ! Index identifying which row to use from NUM_SEi array

! **********************************************************************************************************************************
      IF      (K == 0) THEN                                ! K = 0 is elem node force requests. (NELGP num of rows of output/elem)
!                                                            -----
         NUMBER_ROWS(K) = NELGP(L)

      ELSE IF (K == 1) THEN                                ! K = 1 is elem engr force output requests. (only 1 row of output/elem)
!                                                            -----
         NUMBER_ROWS(K) = 1
! --- cbeam_add begin --- !
         IF (TYPE(1:4) == 'BEAM') THEN
            NUMBER_ROWS(K) = PBEAM_NSTATIONS(EDAT(EPNT(INT_ELEM_ID)+1))
            IF (NUMBER_ROWS(K) <= 0) NUMBER_ROWS(K) = 5
         ENDIF
! --- cbeam_add end --- !
         IF (TYPE(1:5) == 'QUAD4') THEN
            IF (FORC_LOC == 'CENTER  ') THEN            !    PSHELL requires 2 rows of output/elem for FORC_LOC = 'CENTER'
               NUMBER_ROWS(K) = 1
            ELSE                                        !    PSHELL requires more lines of output for other FORC_LOC
               NUMBER_ROWS(K) = NUM_SEi(LETYPE)
            ENDIF
         ENDIF

      ELSE IF (K == 2) THEN                                ! K = 2 is elem stress output requests
!                                                            -----
         NUMBER_ROWS(K) = 1                                !    1 row of stress output/elem unless elem is BAR or shell

         CALL IS_ELEM_PCOMP_PROPS ( INT_ELEM_ID )          !    See if this is a PCOMP elem and get NUM_PLIES
         IF (PCOMP_PROPS == 'Y') THEN
            CALL GET_ELEM_NUM_PLIES ( INT_ELEM_ID )
         ENDIF

         IF  (TYPE(1:3) == 'BAR  ') THEN
               NUMBER_ROWS(K) = 2                          !    BAR stresses require 2 rows of output/elem
! --- cbeam_add begin --- !
         ELSE IF  (TYPE == 'BEAM    ') THEN
               NUMBER_ROWS(K) = 2*PBEAM_NSTATIONS(EDAT(EPNT(INT_ELEM_ID)+1))
               IF (NUMBER_ROWS(K) <= 0) NUMBER_ROWS(K) = 10
! --- cbeam_add end --- !
         ELSE IF ((TYPE(1:5) == 'TRIA3' ) .OR. (TYPE(1:5) == 'QUAD4')) THEN
            IF (PCOMP_PROPS == 'Y') THEN
               NUMBER_ROWS(K) = NUM_PLIES                  !    PCOMP requires NUM_PLIES rows of output/elem
            ELSE
               IF (STRE_LOC == 'CENTER  ') THEN            !    PSHELL requires 2 rows of output/elem for STRE_LOC = 'CENTER'
                  NUMBER_ROWS(K) = 2
               ELSE                                        !    PSHELL requires more lines of output for other STRE_LOC
                  NUMBER_ROWS(K) = 2*NUM_SEi(LETYPE)
               ENDIF
            ENDIF
         ELSE IF (TYPE(1:5) == 'QUAD8' ) THEN
            IF (PCOMP_PROPS == 'Y') THEN
               NUMBER_ROWS(K) = NUM_PLIES                  !    PCOMP requires NUM_PLIES rows of output/elem
            ELSE
               NUMBER_ROWS(K) = 2*NUM_SEi(LETYPE)          !    CQUAD8 stress output is CORNER even if CENTER is specified.
            ENDIF
         ELSE IF ((TYPE(1:4) == 'HEXA' ) .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
            NUMBER_ROWS(K) = NUM_SEi(LETYPE)
         ENDIF

      ELSE IF (K == 3) THEN                                ! K = 3 is elem strain output requests
!                                                            -----
         NUMBER_ROWS(K) = 1                                !    1 row of strain output/elem unless elem is BAR or shell

         CALL IS_ELEM_PCOMP_PROPS ( INT_ELEM_ID )          !    See if this is a PCOMP elem and get NUM_PLIES
         IF (PCOMP_PROPS == 'Y') THEN
            CALL GET_ELEM_NUM_PLIES ( INT_ELEM_ID )
         ENDIF


         IF ((TYPE(1:5) == 'TRIA3' ) .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'SHEAR')) THEN
            IF (PCOMP_PROPS == 'Y') THEN
               NUMBER_ROWS(K) = NUM_PLIES                  !    PCOMP requires NUM_PLIES rows of output/elem
            ELSE
               IF (STRN_LOC == 'CENTER  ') THEN            !    PSHELL requires 2 rows of output/elem for STRN_LOC = 'CENTER'
                  NUMBER_ROWS(K) = 2
               ELSE                                        !    PSHELL requires more lines of output for other STRN_LOC
                  NUMBER_ROWS(K) = 2*NUM_SEi(LETYPE)
               ENDIF
            ENDIF
         ELSE IF ((TYPE(1:4) == 'HEXA' ) .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
            NUMBER_ROWS(K) = NUM_SEi(LETYPE)
! --- cbeam_add begin --- !
         ELSE IF (TYPE(1:4) == 'BEAM') THEN
            NUMBER_ROWS(K) = 2*PBEAM_NSTATIONS(EDAT(EPNT(INT_ELEM_ID)+1))
            IF (NUMBER_ROWS(K) <= 0) NUMBER_ROWS(K) = 10
! --- cbeam_add end--- !
         ENDIF

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE GET_NUM_ROWS_OUTPUT

! ##################################################################################################################################

      SUBROUTINE MAXREQ_OGEL_DEB ( WHICH )

      CHARACTER( 2*BYTE)              :: WHICH             ! Decides what to print out for this call to this subr

! **********************************************************************************************************************************
      IF      (WHICH == '11') THEN

         WRITE(F06,1101)
         WRITE(F06,1102)
         WRITE(F06,1103)

      ELSE IF (WHICH == '12') THEN

         IF  ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            WRITE(F06,1201) SCNUM(I)
         ELSE
            WRITE(F06,1202)
         ENDIF

      ELSE IF (WHICH == '13') THEN

         WRITE(F06,1301) GROUT_NAME(K), NREQ_GR(K), MPF_MEFM_MSG

      ELSE IF (WHICH == '14') THEN

         IF (MAXGROUT_SC > 0) THEN
            IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
               WRITE(F06,1401) MAXGROUT_SC
            ENDIF
         ENDIF

      ELSE IF (WHICH == '15') THEN

         IF  ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            WRITE(F06,1501) MAXGROUT
         ELSE
            WRITE(F06,1502) MAXGROUT
         ENDIF
         WRITE(F06,1103)

! ---------------------------------------------------------------------------------------------------------------------------------
      ELSE IF (WHICH == '21') THEN

         IF  ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            WRITE(F06,2101) SCNUM(I)
         ELSE
            WRITE(F06,2102)
         ENDIF

      ELSE IF (WHICH == '22') THEN

         WRITE(F06,2201) ELMTYP(L)

      ELSE IF (WHICH == '23') THEN

         WRITE(F06,1301) ELOUT_NAME(K), NREQ_EL(L,K)

      ELSE IF (WHICH == '24') THEN

         IF (MAXELOUT_SC > 0) THEN
            WRITE(F06,2401) ELMTYP(L), MAXELOUT_SC
         ENDIF

      ELSE IF (WHICH == '25') THEN

         WRITE(F06,2501) MAXELOUT
         WRITE(F06,1103)

      ELSE IF (WHICH == '31') THEN

         WRITE(F06,3101) MAXREQ
         WRITE(F06,3102)

      ENDIF

! **********************************************************************************************************************************
 1101 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' :::::::::::::::::::::::::::::::::::::::START DEBUG(91) OUTPUT FROM SUBROUTINE MAXREQ_OGEL:::::::::::::::::::::::::',&
              ':::::::::::::::::',/)

 1102 FORMAT(10X,'C A L C U L A T E   N U M B E R   O F   R O W S   O F   O U T P U T   N E E D E D   F O R   A R R A Y   O G E L'/)

 1103 FORMAT(1X,'****************************************************************************************************************',&
                '*******************')

 1201 FORMAT(1X,/,' Grid point related output requests for subcase ',I8                                                          ,/&
                  ' -------------------------------------------------------',/)

 1202 FORMAT(1X,/,' Grid point related output requests'                                                                          ,/&
                  ' ----------------------------------',/)

 1301 FORMAT('   Number of rows of output needed for ',A,':',I8,2X,A)

 1401 FORMAT(1X,/,'   The maximum number of rows needed for grid related outputs for this subcase is    : ',I8,/)

 1501 FORMAT(1X,/,'   The maximum number of rows needed for grid related outputs for all subcases is    : ',I8,/)

 1502 FORMAT(1X,/,'   The maximum number of rows needed for grid related outputs is                     : ',I8,/)

 2101 FORMAT(1X,/,' Element related output requests for subcase ',I8                                                             ,/&
             ' ----------------------------------------------------',/)

 2102 FORMAT(1X,/,' Element related output requests'                                                                             ,/&
                  ' -------------------------------',/)

 2201 FORMAT('   For element type: ',A                                                                                          ,/,&
             '   -------------------------')

 2401 FORMAT(1X,/,'   The maximum number of rows needed for ',A,' element outputs for this subcase is: ',I8,/)

 2402 FORMAT(1X,/,'   The maximum number of rows needed for ',A,' element outputs is                    : ',I8,/)

 2501 FORMAT(1X,/,'   The maximum number of rows needed for all element outputs for all subcases is     : ',I8,/)

 3101 FORMAT(1X,/,'   The maximum number of rows of output needed for any subcase is                    : ',I8,/)

 3102 FORMAT(' ::::::::::::::::::::::::::::::::::::::END DEBUG(91) OUTPUT FROM SUBROUTINE MAXREQ_OGEL::::::::::::::::::::::::::::',&
              ':::::::::::::::::'                                                                                               ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE MAXREQ_OGEL_DEB

      END SUBROUTINE MAXREQ_OGEL



