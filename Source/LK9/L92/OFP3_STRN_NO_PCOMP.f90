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

      SUBROUTINE OFP3_STRN_NO_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element strain output requests for non PCOMP elements for one subcase. Also write Output Transformation Matrices (OTM's)
! for strains for Craig-Bampton models)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_STRN_BIT, FATAL_ERR, IBIT, INT_SC_NUM,                                &
                                         MAX_STRESS_POINTS, MBUG, MOGEL,                                                           &
                                         NELE, NCBAR, NCBUSH, NCELAS1, NCELAS2, NCELAS3, NCELAS4, NCHEXA8, NCHEXA20, NCPENTA6,     &
                                         NCPENTA15,NCTETRA4, NCTETRA10, NCQUAD4, NCQUAD4K, NCROD, NCSHEAR, NCTRIA3, NCTRIA3K,      &
                                         SOL_NAME


      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, TWO, FOUR
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS
      USE PARAMS, ONLY                :  OTMSKIP, PRTNEU
      USE MODEL_STUF, ONLY            :  AGRID, ANY_STRN_OUTPUT, EDAT, EPNT, ETYPE, EID, ELGP, ELMTYP, ELOUT,                      &
                                         METYPE, NUM_SEi, NUM_EMG_FATAL_ERRS, PCOMP_PROPS, PLY_NUM, STRAIN, TYPE, SHELL_STR_ANGLE 
      USE CC_OUTPUT_DESCRIBERS, ONLY  :  STRN_LOC, STRN_OPT
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, GID_OUT_ARRAY, MAXREQ, OGEL, POLY_FIT_ERR, POLY_FIT_ERR_INDEX
      USE OUTPUT4_MATRICES, ONLY      :  OTM_STRN, TXT_STRN
! --- cbeam_add begin --- !
      USE SCONTR, ONLY                :  NCBEAM
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_VECS
      USE MODEL_STUF, ONLY            :  CBEAM_ACTIVE_NSTATIONS, CBEAM_ACTIVE_XL, PBEAM_NSTATIONS, ZS
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE LINK9_STUFF, ONLY           :  CBEAM_XL_OUT
! --- cbeam_add end ---

      USE PLANE_COORD_TRANS_21_Interface
      USE TRANSFORM_SHELL_STR_Interface
      USE OFP3_STRN_NO_PCOMP_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_STRN_NO_PCOMP'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message
      CHARACTER(20*BYTE)              :: STRAIN_ITEM(22)   ! Char description of element strains

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: ELOUT_STRN        ! If > 0, there are STRAIN   requests for some elems
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: IERROR    = 0     ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: NDUM              ! Dummy valye needed in call to CALC_ELEM_ENFR_FORCES
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
! --- cbeam_add begin --- !
      INTEGER(LONG)                   :: NUM_PTS_ELEM      ! Num strain stations/points for current element
      INTEGER(LONG)                   :: NUM_PTS_CUR       ! Actual number of strain points for the current element
! --- cbeam_add end --- !
      INTEGER(LONG)                   :: NUM_OGEL_ROWS     ! No. elems processed prior to writing results to F06 file
      INTEGER(LONG)                   :: NUM_FROWS         ! No. elems processed for FEMAP
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_OGEL_ROWS since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)
                                                           ! Indicator for output of elem data to BUG file
      INTEGER(LONG)                   :: NUM_OTM_ENTRIES   ! Number of entries in OGEL for a particular element type
      INTEGER(LONG)                   :: NUM_PTS(METYPE)   ! Num diff strain points for one element (3rd dim in arrays SEi, STEi)

                                                           ! Strain index (1 through 9) where poly fit err is max
      INTEGER(LONG)                   :: STRAIN_OUT_ERR_INDEX(MAX_STRESS_POINTS+1) ! add+1 for cbeam
      INTEGER(LONG)                   :: NUM_STR_POINTS      


                                                           ! Array of %errs from subr POLYNOM_FIT_STRE_STRN (only NUM_PTS vals used)
      REAL(DOUBLE)                    :: STRAIN_OUT_PCT_ERR(MAX_STRESS_POINTS+1) ! add+1 for cbeam

      REAL(DOUBLE)                    :: PCT_ERR_MAX       ! Max value from array STRAIN_OUT_PCT_ERR
! --- cbeam_add begin --- !
      REAL(DOUBLE)                    :: C1,C2,D1,D2,E1,E2,F1,F2
      REAL(DOUBLE)                    :: EA0,EA1,EA2,EA3,EA4,EAMAX,EAMIN
      REAL(DOUBLE)                    :: EB0,EB1,EB2,EB3,EB4,EBMAX,EBMIN
! --- cbeam_add end --- !
                                                           ! Array of values from array STRAIN for all stress points
      REAL(DOUBLE)                    :: STRAIN_RAW(9,MAX_STRESS_POINTS+1)

                                                           ! Array of output stress values after surface fit
      REAL(DOUBLE)                    :: STRAIN_OUT(9,MAX_STRESS_POINTS+1)
      REAL(DOUBLE)                    :: TEL(3,3)          ! Transformation matrix from cartesian local (L) to element (E) coordinates.

      ! OP2 stuff
      CHARACTER(8*BYTE)               :: TABLE_NAME   ! name of the op2 table name
      INTEGER(LONG)                   :: ITABLE       ! the subtable
      LOGICAL                         :: WRITE_NEU
! --- cbeam_add begin --- !
      LOGICAL                         :: HAVE_SECTION_POINTS

      INTRINSIC DABS, DMAX1, DMIN1, IAND
! --- cbeam_add end --- !
      ITABLE = 0
      TABLE_NAME = "OES ERR "

      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Process element strain output (STRAIN) requests for all elems except composite shells

      OPT(1) = 'N'                                         ! OPT(1) is for calc of ME
      OPT(2) = 'N'                                         ! OPT(2) is for calc of PTE
      OPT(3) = 'Y'                                         ! OPT(3) is for calc of SEi, STEi
      OPT(4) = 'N'                                         ! OPT(4) is for calc of KE-linear
      OPT(5) = 'N'                                         ! OPT(5) is for calc of PPE
      OPT(6) = 'N'                                         ! OPT(6) is for calc of KE-diff stiff

! Find out how many output requests were made for each element type.

      DO I=1,METYPE                                        ! Initialize the array containing the no. requests/elem.
         NELREQ(I) = 0
      ENDDO

      DO I=1,METYPE  !metype
                                       ! Only count requests for elem types that can have strain output
! --- cbeam_add begin --- !
         IF((ELMTYP(I)(1:3) == 'BAR'  ) .OR. (ELMTYP(I)(1:4) == 'BEAM' ) .OR. (ELMTYP(I)(1:3) == 'ROD'  ) .OR.                    &
! --- cbeam_add end --- !
            (ELMTYP(I)(1:4) == 'ELAS' ) .OR. (ELMTYP(I)(1:4) == 'BUSH' ) .OR. (ELMTYP(I)(1:5) == 'TRIA3') .OR.                    &
            (ELMTYP(I)(1:5) == 'QUAD4') .OR. (ELMTYP(I)(1:5) == 'SHEAR') .OR. (ELMTYP(I)(1:4) == 'HEXA' ) .OR.                    &
            (ELMTYP(I)(1:5) == 'PENTA') .OR. (ELMTYP(I)(1:5) == 'TETRA') .OR. (ELMTYP(I)(1:5) == 'QUAD8')) THEN
            DO J=1,NELE !nele
               CALL IS_ELEM_PCOMP_PROPS ( J )
               IF (PCOMP_PROPS == 'N') THEN !pcomp_props
                  IF (ETYPE(J) == ELMTYP(I)) THEN !etype j

! --- cbeam_add begin --- !
                     IF (ETYPE(J) == 'BEAM    ') THEN
                        NUM_PTS_ELEM = PBEAM_NSTATIONS(EDAT(EPNT(J)+1))
                        IF (NUM_PTS_ELEM <= 0) NUM_PTS_ELEM = 5
                        IF (NUM_PTS_ELEM > NUM_PTS(I)) NUM_PTS(I) = NUM_PTS_ELEM
                     ELSE !not beam
! --- cbeam_add end --- !
                        IF ((STRN_LOC == 'CORNER  ') .OR.                                                                         &
                            (STRN_LOC == 'GAUSS   ') .OR.                                                                         &
                            (ETYPE(J)(1:4) == 'HEXA') .OR.                                                                        &
                            (ETYPE(J)(1:5) == 'PENTA') .OR.                                                                       &
                            (ETYPE(J)(1:5) == 'TETRA') .OR.                                                                       &
                            (ETYPE(J)(1:5) == 'QUAD8')) THEN
                            NUM_PTS(I) = NUM_SEi(I)! 
                            NUM_PTS_ELEM = NUM_SEi(I) 
                        ELSE ! bush maybe 
                           NUM_PTS_ELEM = 1
                           NUM_PTS(I) = 1
                        ENDIF ! corner
                     ENDIF ! beam or not
                     
                     ELOUT_STRN = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRN_BIT))
                     IF (ELOUT_STRN > 0) THEN
                        NELREQ(I) = NELREQ(I) + NUM_PTS(I) ! NELREQ(I) = NELREQ(I) + NUM_PTS_ELEM  !Checkagain
                     ENDIF
                  ENDIF ! tyoe j
               ENDIF ! pcomp
            ENDDO ! nele
         ENDIF !bar or something
      ENDDO ! metype

      OGEL = ZERO

!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element strain'
reqs7:DO I=1,METYPE
         IF (NELREQ(I) == 0) CYCLE reqs7
         NUM_OGEL_ROWS = 0
         NUM_OGEL = 0
elems_7: DO J = 1,NELE

            EID   = EDAT(EPNT(J))
            TYPE  = ETYPE(J)
            IF (ETYPE(J) == ELMTYP(I)) THEN
               ELOUT_STRN = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRN_BIT))
               IF (ELOUT_STRN > 0) THEN
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )! Calc SEi matrices
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE elems_7
                  ENDIF
                  CALL ELMDIS

! --- cbeam_add begin --- !
                  NUM_STR_POINTS =  NUM_PTS(I) 
                  NUM_PTS_CUR = NUM_PTS(I)
                  IF (TYPE == 'BEAM    ') THEN
                     NUM_PTS_CUR = CBEAM_ACTIVE_NSTATIONS
                     IF (NUM_PTS_CUR <= 0) NUM_PTS_CUR = 1
                     NUM_STR_POINTS  = NUM_PTS_CUR
                  ENDIF
! --- cbeam_add end --- !
                  DO M=1,NUM_STR_POINTS  
                     CALL ELEM_STRE_STRN_ARRAYS ( M )
                     DO K=1,9
                        STRAIN_RAW(K,M) = STRAIN(K)
                     ENDDO
                  ENDDO

! --- cbeam_stations begin --- !
                  IF (TYPE == 'BEAM    ') THEN
                     STRAIN_OUT(:,:) = STRAIN_RAW(:,:)
                  ELSE
                     STRAIN_OUT(:,1) = STRAIN(:)              ! Set STRAIN_OUT for NUM_PTS(I) = 1
                  ENDIF
! --- cbeam_stations end --- !

                  IF ((STRN_LOC == 'CORNER  ') .OR.                                                                                &
                      (STRN_LOC == 'GAUSS   ') .OR.                                                                                &
                      (TYPE(1:4) == 'HEXA') .OR.                                                                                   &
                      (TYPE(1:5) == 'PENTA') .OR.                                                                                  &
                      (TYPE(1:5) == 'TETRA') .OR.                                                                                  &
                      (TYPE(1:5) == 'QUAD8')) THEN

                     IF (TYPE(1:5) == 'QUAD4') THEN
                        CALL POLYNOM_FIT_STRE_STRN ( STRAIN_RAW, 9, NUM_PTS_CUR, STRAIN_OUT, STRAIN_OUT_PCT_ERR,                    &
                                                     STRAIN_OUT_ERR_INDEX, PCT_ERR_MAX )

                     ELSE IF (TYPE(1:5) == 'QUAD8') THEN
                        CALL POLYNOM_FIT_STRE_STRN ( STRAIN_RAW, 9, NUM_PTS_CUR, STRAIN_OUT, STRAIN_OUT_PCT_ERR,                    &
                                                     STRAIN_OUT_ERR_INDEX, PCT_ERR_MAX )

                                                           ! Transform strain from the cartesian local coordinate system to
                                                           ! the element coordinate system
                        DO M=1,NUM_PTS(I)
                           CALL PLANE_COORD_TRANS_21( SHELL_STR_ANGLE( M ), TEL, '')
                           CALL TRANSFORM_SHELL_STR( TEL, STRAIN_OUT(:,M), TWO)
                        ENDDO

                                                           ! Center strain is the average of corner strain in element coordinates.
                                                           ! This is how MSC does it.
                        STRAIN_OUT(:,1) = (STRAIN_OUT(:,2) + STRAIN_OUT(:,3) + STRAIN_OUT(:,4) + STRAIN_OUT(:,5)) / FOUR

                     ELSE IF ((TYPE(1:4) == 'HEXA') .OR.                                                                           &
                              (TYPE(1:5) == 'PENTA') .OR.                                                                          &
                              (TYPE(1:5) == 'TETRA')) THEN
! Strains are directly evaluated at the corner grid points. If they are going to be evaluated at Gauss points
! then extrapolated to grid points, that should be done here, in POLYNOM_FIT_STRE_STRN, or in an equivalent subroutine.
                        STRAIN_OUT(:,:) = STRAIN_RAW(:,:)

                     ENDIF

                  ENDIF
                  IF (TYPE == 'BEAM    ') THEN
			NUM_STR_POINTS = NUM_PTS_CUR
		  ELSE
			NUM_STR_POINTS = NUM_PTS(I)
		  ENDIF

do_strain_pts:    DO M=1,NUM_STR_POINTS 

                     DO K=1,9
                        STRAIN(K) = STRAIN_OUT(K,M)
                     ENDDO
                     CALL CALC_ELEM_STRAINS ( MAXREQ, NUM_OGEL, J, 'Y', 'N' )
                                                              ! If CB soln, write rows of OGEL, from CALC_ELEM_STRAINS, to OTM_STRN
                     IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

                        CALL GET_STRAIN_ITEM_DATA

                        IF ((TYPE == 'TRIA3   ') .OR. (TYPE == 'QUAD4   ') .OR. (TYPE == 'SHEAR   ')) THEN
                           DO L=1,2
                              DO K=1,NUM_OTM_ENTRIES
                                 OT4_EROW = OT4_EROW + 1
                                 OTM_STRN(OT4_EROW,JVEC) = OGEL(NUM_OGEL-2+L,K)
                                 IF (JVEC == 1) THEN
                                    IF ((STRN_LOC == 'CORNER  ') .OR. (STRN_LOC == 'GAUSS   ')) THEN
                                       IF (M == 1) THEN
                                          IF (TYPE(1:5) == 'QUAD4') THEN
                                             WRITE(TXT_STRN(OT4_EROW),9190) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                   &
                                                                           STRAIN_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                          ELSE
                                             WRITE(TXT_STRN(OT4_EROW), 9193) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                  &
                                                                             STRAIN_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                          ENDIF
                                       ELSE
                                          WRITE(TXT_STRN(OT4_EROW),9191) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, AGRID(M-1),          &
                                                                         STRAIN_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                       ENDIF
                                    ELSE
                                       WRITE(TXT_STRN(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                        &
                                                                       STRAIN_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                    ENDIF
                                 ENDIF
                              ENDDO
                           ENDDO
                        ELSE
                           DO K=1,NUM_OTM_ENTRIES
                              OT4_EROW = OT4_EROW + 1
                              OTM_STRN(OT4_EROW,JVEC) = OGEL(NUM_OGEL,K)
                              IF (JVEC == 1) THEN
                                 WRITE(TXT_STRN(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID, STRAIN_ITEM(K)
                              ENDIF
                           ENDDO
                        ENDIF
                     ENDIF

                     IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                        DO K=1,OTMSKIP                        ! Write OTMSKIP blank separator lines
                           OT4_EROW = OT4_EROW + 1
                           WRITE(TXT_STRN(OT4_EROW), 9199)
                        ENDDO
                     ENDIF

                     NUM_OGEL_ROWS = NUM_OGEL_ROWS + 1
                     EID_OUT_ARRAY(NUM_OGEL_ROWS,1) = EID
! --- cbeam_add begin --- !
                     IF (TYPE == 'BEAM    ') THEN
                        CBEAM_XL_OUT(NUM_OGEL_ROWS) = CBEAM_ACTIVE_XL(M)
                     ELSE
                        CBEAM_XL_OUT(NUM_OGEL_ROWS) = ZERO
                     ENDIF
! --- cbeam_add end --- !
                     GID_OUT_ARRAY(NUM_OGEL_ROWS,1) = 0
                     IF ((STRN_LOC == 'CORNER  ') .OR. (STRN_LOC == 'GAUSS   ')) THEN
                        IF (TYPE(1:5) == 'QUAD4') THEN
                           POLY_FIT_ERR(NUM_OGEL_ROWS)       = STRAIN_OUT_PCT_ERR(M)
                           POLY_FIT_ERR_INDEX(NUM_OGEL_ROWS) = STRAIN_OUT_ERR_INDEX(M)
                        ENDIF
                     ENDIF
                     DO K=1,ELGP
                        GID_OUT_ARRAY(NUM_OGEL_ROWS,K+1) = AGRID(K)
                     ENDDO

                  ENDDO do_strain_pts

                  IF (ETYPE(J)(1:5) /='USER1') THEN
                     IF (NUM_OGEL_ROWS == NELREQ(I)) THEN
                        CALL CHK_OGEL_ZEROS ( NUM_OGEL )
 100                    FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE=",A,"; TABLE_NAME=",A,"; ITABLE=",I8)
                        WRITE(ERR,100) "A",TYPE,TABLE_NAME,ITABLE
                        CALL SET_OST_TABLE_NAME(TYPE, TABLE_NAME, ITABLE)
                        WRITE(ERR,100) "B",TYPE,TABLE_NAME,ITABLE
                        CALL WRITE_ELEM_STRAINS ( JVEC, NUM_OGEL_ROWS, IHDR, NUM_STR_POINTS, ITABLE )
                        EXIT
                     ENDIF
                  ENDIF

               ENDIF

            ENDIF

         ENDDO elems_7

      ENDDO reqs7

      IF ((TABLE_NAME .NE. "OES ERR ") .AND. (ITABLE < 0)) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF
!===========================
      IF (WRITE_NEU .AND. (ANY_STRN_OUTPUT > 0)) THEN

         NDUM = 0
! --- cbeam_add begin --- !
         NUM_FROWS= 0                                      ! Write out BEAM strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBEAM, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:4) == 'BEAM') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=1,12
                     FEMAP_EL_VECS(NUM_FROWS,K) = ZERO
                  ENDDO
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  HAVE_SECTION_POINTS = .FALSE.
                  DO K=1,8
                     IF (DABS(ZS(K)) > ZERO) THEN
                        HAVE_SECTION_POINTS = .TRUE.
                        EXIT
                     ENDIF
                  ENDDO
                  IF (HAVE_SECTION_POINTS) THEN
                     C1 = ZS(1); C2 = ZS(2)
                     D1 = ZS(3); D2 = ZS(4)
                     E1 = ZS(5); E2 = ZS(6)
                     F1 = ZS(7); F2 = ZS(8)
                  ELSE
                     C1 = ZERO; C2 = ZERO
                     D1 = ZERO; D2 = ZERO
                     E1 = ZERO; E2 = ZERO
                     F1 = ZERO; F2 = ZERO
                  ENDIF
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  EA0   = STRAIN(1)
                  EA1   = EA0 - (C1*STRAIN(2) + C2*STRAIN(3))
                  EA2   = EA0 - (D1*STRAIN(2) + D2*STRAIN(3))
                  EA3   = EA0 - (E1*STRAIN(2) + E2*STRAIN(3))
                  EA4   = EA0 - (F1*STRAIN(2) + F2*STRAIN(3))
                  EAMAX = DMAX1(EA1,EA2,EA3,EA4)
                  EAMIN = DMIN1(EA1,EA2,EA3,EA4)
                  FEMAP_EL_VECS(NUM_FROWS, 1) = EA1
                  FEMAP_EL_VECS(NUM_FROWS, 3) = EA2
                  FEMAP_EL_VECS(NUM_FROWS, 5) = EA3
                  FEMAP_EL_VECS(NUM_FROWS, 7) = EA4
                  FEMAP_EL_VECS(NUM_FROWS, 9) = EAMAX
                  FEMAP_EL_VECS(NUM_FROWS,11) = EAMIN
                  IF (CBEAM_ACTIVE_NSTATIONS > 1) THEN
                     CALL ELEM_STRE_STRN_ARRAYS ( CBEAM_ACTIVE_NSTATIONS )
                  ELSE
                     CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  ENDIF
                  EB0   = STRAIN(1)
                  EB1   = EB0 - (C1*STRAIN(2) + C2*STRAIN(3))
                  EB2   = EB0 - (D1*STRAIN(2) + D2*STRAIN(3))
                  EB3   = EB0 - (E1*STRAIN(2) + E2*STRAIN(3))
                  EB4   = EB0 - (F1*STRAIN(2) + F2*STRAIN(3))
                  EBMAX = DMAX1(EB1,EB2,EB3,EB4)
                  EBMIN = DMIN1(EB1,EB2,EB3,EB4)
                  FEMAP_EL_VECS(NUM_FROWS, 2) = EB1
                  FEMAP_EL_VECS(NUM_FROWS, 4) = EB2
                  FEMAP_EL_VECS(NUM_FROWS, 6) = EB3
                  FEMAP_EL_VECS(NUM_FROWS, 8) = EB4
                  FEMAP_EL_VECS(NUM_FROWS,10) = EBMAX
                  FEMAP_EL_VECS(NUM_FROWS,12) = EBMIN
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'BEAM    ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
! --- cbeam_add end --- !
         NUM_FROWS= 0                                      ! Write out BUSH strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCBUSH, 6, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:4) == 'BUSH') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCBUSH, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'BUSH  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3K strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3K, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'TRIA3K') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCTRIA3K, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'TRIA3K  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTRIA3, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'TRIA3 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCTRIA3, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'TRIA3   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4K strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4K, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'QUAD4K') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCQUAD4K, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'QUAD4K  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCQUAD4, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'QUAD4 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCQUAD4, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'QUAD4   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out HEXA8 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCHEXA8, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'HEXA8 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCHEXA8, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'HEXA8   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out HEXA20 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCHEXA20, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:6) == 'HEXA20') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCHEXA20, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'HEXA20  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PENTA6 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCPENTA6, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PENTA6 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCPENTA6, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'PENTA6  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out PENTA15 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCPENTA15, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'PENTA15') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCPENTA15, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'PENTA15 ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TETRA4 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTETRA4, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'TETRA4 ') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCTETRA4, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'TETRA4  ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TETRA10 strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCTETRA10, 12, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:7) == 'TETRA10') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCTETRA10, NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'TETRA10 ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

         NUM_FROWS= 0                                      ! Write out SHEAR strains
         CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NCSHEAR, 22, SUBR_NAME )
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'N') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'SHEAR') THEN
                  NUM_FROWS= NUM_FROWS+ 1
                  DO K=0,MBUG-1
                     WRT_BUG(K) = 0
                  ENDDO
                  PLY_NUM = 1                              ! 'N' in call to EMG means do not write to BUG file
                  CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                  FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                  IF (NUM_EMG_FATAL_ERRS > 0) THEN
                     IERROR = IERROR + 1
                     CYCLE
                  ENDIF
                  CALL ELMDIS
                  CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                  CALL CALC_ELEM_STRAINS ( NCSHEAR , NDUM, NUM_FROWS, 'N', 'Y' )
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL WRITE_FEMAP_STRN_VECS ( 'SHEAR   ', 'N', NUM_FROWS, FEMAP_SET_ID )
         ENDIF
         CALL DEALLOCATE_FEMAP_DATA

      ENDIF

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT STRAIN'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9190 FORMAT(I8,1X,A,A8,I8,4X,'CENTER',9X,A20)

 9191 FORMAT(I8,1X,A,A8,I8,4X,'GRID',I8,3X,A20)

 9192 FORMAT(I8,1X,A,A8,I8,4X,A20)

 9193 FORMAT(I8,1X,A,A8,I8,19X,A20)

 9199 FORMAT(' ')

 9201 FORMAT(' *ERROR  9201: DUE TO ABOVE LISTED ERRORS, CANNOT CALCULATE ',A,' REQUESTS FOR ',A,' ELEMENT ID = ',I8)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_STRAIN_ITEM_DATA

      IMPLICIT NONE

      INTEGER(LONG)                   :: II               ! DO loop index

! **********************************************************************************************************************************
      DO II=1,18
         STRAIN_ITEM(II)(1:) = ' '
      ENDDO

      IF  (TYPE(1:3) == 'ROD'  ) THEN
         NUM_OTM_ENTRIES = 4
         STRAIN_ITEM( 1) = 'Axial Strain        '
         STRAIN_ITEM( 2) = 'MS - Axial          '
         STRAIN_ITEM( 3) = 'Torsional Strain    '
         STRAIN_ITEM( 4) = 'MS - Torsion        '

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4')) THEN
         NUM_OTM_ENTRIES = 10
         STRAIN_ITEM( 1) = 'Fibre Dist      -Z1 '  ;  STRAIN_ITEM(11) = 'Fibre Dist      +Z1 '
         STRAIN_ITEM( 2) = 'Normal X Strain -Z1 '  ;  STRAIN_ITEM(12) = 'Normal X Strain +Z1 '
         STRAIN_ITEM( 3) = 'Normal Y Strain -Z1 '  ;  STRAIN_ITEM(13) = 'Normal Y Strain +Z1 '
         STRAIN_ITEM( 4) = 'Shear XY Strain -Z1 '  ;  STRAIN_ITEM(14) = 'Shear XY Strain +Z1 '
         STRAIN_ITEM( 5) = 'Princ Angle  at -Z1 '  ;  STRAIN_ITEM(15) = 'Princ Angle  at +Z1 '
         STRAIN_ITEM( 6) = 'Major Strain at -Z1 '  ;  STRAIN_ITEM(16) = 'Major Strain at +Z1 '
         STRAIN_ITEM( 7) = 'Minor Strain at -Z1 '  ;  STRAIN_ITEM(17) = 'Minor Strain at +Z1 '
         IF      (STRN_OPT(1:8) == 'VONMISES') THEN
            STRAIN_ITEM( 8) = 'von Mises XY at -Z1 '  ;  STRAIN_ITEM(18) = 'von Mises XY at +Z1 '
         ELSE IF (STRN_OPT(1:4) == 'MAXS'    ) THEN
            STRAIN_ITEM( 8) = 'Max Shear XY at -Z1 '  ;  STRAIN_ITEM(18) = 'Max Shear XY at +Z1 '
         ELSE
            STRAIN_ITEM( 8) = '**** undefined **** '  ;  STRAIN_ITEM(18) = '**** undefined **** '
         ENDIF
         STRAIN_ITEM( 9) = 'Transv Shear XZ avg '  ;  STRAIN_ITEM(19) = 'Transv Shear XZ avg '
         STRAIN_ITEM(10) = 'Transv Shear YZ avg '  ;  STRAIN_ITEM(20) = 'Transv Shear YZ avg '

      ELSE IF  (TYPE(1:6) == 'USERIN') THEN
         NUM_OTM_ENTRIES = 9
         STRAIN_ITEM( 1) = '*** Not Defined ****'
         STRAIN_ITEM( 2) = '*** Not Defined ****'
         STRAIN_ITEM( 3) = '*** Not Defined ****'
         STRAIN_ITEM( 4) = '*** Not Defined ****'
         STRAIN_ITEM( 5) = '*** Not Defined ****'
         STRAIN_ITEM( 6) = '*** Not Defined ****'
         STRAIN_ITEM( 7) = '*** Not Defined ****'
         STRAIN_ITEM( 8) = '*** Not Defined ****'
         STRAIN_ITEM( 9) = '*** Not Defined ****'

      ELSE IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         NUM_OTM_ENTRIES = 8
         STRAIN_ITEM( 1) = 'Normal x Strain     '
         STRAIN_ITEM( 2) = 'Normal y Strain     '
         STRAIN_ITEM( 3) = 'Normal z Strain     '
         STRAIN_ITEM( 4) = 'Shear xy Strain     '
         STRAIN_ITEM( 5) = 'Shear yz Strain     '
         STRAIN_ITEM( 6) = 'Shear zx Strain     '
         STRAIN_ITEM( 7) = 'Oct Direct Strain   '
         STRAIN_ITEM( 8) = 'Oct Shear Strain    '

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE GET_STRAIN_ITEM_DATA

      END SUBROUTINE OFP3_STRN_NO_PCOMP

