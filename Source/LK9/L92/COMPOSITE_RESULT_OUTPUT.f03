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

   MODULE COMPOSITE_RESULT_OUTPUT

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: OFP3_STRE_PCOMP, OFP3_STRN_PCOMP

   CONTAINS

      SUBROUTINE OFP3_STRE_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element ply stress output requests for PCOMP elements (TRIA3, QUAD4) for one subcase

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_STRE_BIT, FATAL_ERR, IBIT, INT_SC_NUM, MBUG, MOGEL,                   &
                                         NELE, NCQUAD4, NCSHEAR,  NCTRIA3, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS
      USE PARAMS, ONLY                :  OTMSKIP, PRTNEU
      USE MODEL_STUF, ONLY            :  ANY_STRE_OUTPUT, EDAT, EPNT, ETYPE, EID, ELMTYP, ELOUT, METYPE, NUM_EMG_FATAL_ERRS,       &
                                         NUM_PLIES, PCOMP_PROPS, PLY_NUM, THETA_PLY, TYPE
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, MAXREQ, OGEL
      USE OUTPUT4_MATRICES, ONLY      :  OTM_STRE, TXT_STRE

      USE COMPOSITE_SHELL_PREPARATION, ONLY:  GET_ELEM_NUM_PLIES, IS_ELEM_PCOMP_PROPS, ROT_COMP_ELEM_AXES
      USE EMG_MOD, ONLY               :  EMG
      USE RESULT_COORDINATES, ONLY    :  ELMDIS, ELMDIS_PLY
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: ELEM_STRE_STRN_ARRAYS
      USE ELEMENT_STRESS_RECOVERY, ONLY: CALC_ELEM_STRESSES
      USE SUSER1_RECOVERY, ONLY       :  SUSER1
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_OGEL_ZEROS
      USE COMPOSITE_PLY_OUTPUT, ONLY  :  WRITE_PLY_STRESSES
      USE LINK9_WORKSPACE, ONLY       :  ALLOCATE_FEMAP_DATA, DEALLOCATE_FEMAP_DATA
      USE FEMAP_OUTPUT_WRITERS, ONLY  :  WRITE_FEMAP_STRE_VECS

      USE OP2_GEOMETRY_OUTPUT, ONLY   :  END_OP2_TABLE
      USE OP2_STRESS_OUTPUT, ONLY     :  SET_OESC_TABLE_NAME
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_STRE_PCOMP'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message
      CHARACTER(20*BYTE)              :: STRESS_ITEM(18)   ! Char description of element strains

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: ELOUT_STRE        ! If > 0, there are STRESS   requests for some elems
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: IERROR       = 0  ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: NDUM              ! Dummy valye needed in call to CALC_ELEM_ENFR_FORCES
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
      INTEGER(LONG)                   :: NUM_FROWS         ! No. elems processed for FEMAP
      INTEGER(LONG)                   :: NUM_LINES         ! No. lines of output for ply stresses
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_ELEM since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)
                                                           ! Indicator for output of elem data to BUG file
      INTEGER(LONG)                   :: NUM_OTM_ENTRIES   ! Number of entries in OGEL for a particular element type


      INTEGER(LONG)                   :: ITABLE           ! the op2 subtable number
      CHARACTER(8*BYTE)               :: TABLE_NAME       ! the op2 table name
      LOGICAL                         :: IS_RESULT        ! is there a result
      LOGICAL                         :: WRITE_NEU

      INTRINSIC IAND
! **********************************************************************************************************************************
      ! initial values
      IS_RESULT = .FALSE.
      TABLE_NAME = "OES1C"
      ITABLE = 0

      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Process element stress output (STRESS) requests for composite shell elements

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

      DO I=1,METYPE
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'Y') THEN
               IF ((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4')) THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                     ELOUT_STRE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRE_BIT))
                     IF (ELOUT_STRE > 0) THEN
                        CALL GET_ELEM_NUM_PLIES ( J )
                        NELREQ(I) = NELREQ(I) + NUM_PLIES
                        IS_RESULT = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      OGEL = ZERO

!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element stress'
reqs4:DO I=1,METYPE
         IF (NELREQ(I) == 0) CYCLE reqs4
         NUM_LINES = 0
         NUM_OGEL  = 0

elems_4: DO J = 1,NELE

            CALL IS_ELEM_PCOMP_PROPS ( J )
pcomp:      IF (PCOMP_PROPS == 'Y') THEN

               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  ELOUT_STRE = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRE_BIT))
                  IF (ELOUT_STRE > 0) THEN
                     DO K=0,MBUG-1
                        WRT_BUG(K) = 0
                     ENDDO

                     CALL GET_ELEM_NUM_PLIES ( J )
do_plies_4:          DO M=1,NUM_PLIES                         ! Cycle over number of plies, processing stress for each ply
                        PLY_NUM = M                           ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE elems_4
                        ENDIF
                        CALL ELMDIS                           ! Get elem displs (at mid plane of the elem) in elem coords
                        IF (NUM_PLIES > 1) THEN
                           CALL ELMDIS_PLY                    ! If there is more than 1 ply, get displs from elem mid plane displs
                        ENDIF

                        IF (ETYPE(J)(1:5) == 'USER1') THEN
                           CALL SUSER1
                        ELSE                                  ! Get ply stresses in ply coords (also need strains for failure index)
                           CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                           CALL ROT_COMP_ELEM_AXES ( J, M, THETA_PLY, '2-1' )
                           CALL CALC_ELEM_STRESSES ( MAXREQ, NUM_OGEL, J, 'Y', 'N' )
                        ENDIF
                                                              ! If CB soln, write rows of OGEL, from CALC_ELEM_STRESSES, to OTM_STRE
                        IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

                           CALL GET_STRESS_ITEM_DATA

                           DO L=1,2
                              DO K=1,NUM_OTM_ENTRIES
                                 OT4_EROW = OT4_EROW + 1
                                 OTM_STRE(OT4_EROW,JVEC) = OGEL(NUM_OGEL-2+L,K)
                                 IF (JVEC == 1) THEN
                                    WRITE(TXT_STRE(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                           &
                                                                     STRESS_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDIF

                        IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                           DO K=1,OTMSKIP                        ! Write OTMSKIP blank separator lines
                              OT4_EROW = OT4_EROW + 1
                              WRITE(TXT_STRE(OT4_EROW), 9199)
                           ENDDO
                        ENDIF

                        NUM_LINES = NUM_LINES + 1
                        EID_OUT_ARRAY(NUM_LINES,1) = EID
                        EID_OUT_ARRAY(NUM_LINES,2) = M        ! Ply number for EID
                        IF ((ETYPE(J)(1:5) =='TRIA3') .OR. (ETYPE(J)(1:5) == "QUAD4")) THEN
                           IF (NUM_LINES == NELREQ(I)) THEN
                              CALL CHK_OGEL_ZEROS ( NUM_OGEL )

 100                          FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE=",A,"; TABLE_NAME=",A,"; ITABLE=",I8)
                              WRITE(ERR,100) "OES_PCOMP",TYPE,TABLE_NAME,ITABLE
                              CALL SET_OESC_TABLE_NAME(TABLE_NAME, ITABLE)
                              WRITE(ERR,100) "OES_PCOMP",ETYPE(J)(1:8),TABLE_NAME,ITABLE

                              CALL WRITE_PLY_STRESSES ( JVEC, NUM_LINES, IHDR, ETYPE(J)(1:8), ITABLE )
                              EXIT
                           ENDIF
                        ENDIF

                     ENDDO do_plies_4

                  ENDIF

               ENDIF

            ENDIF pcomp

         ENDDO elems_4

      ENDDO reqs4

  10   FORMAT("*DEBUG:      OES_PCOMP_END:    TABLE_NAME",A)
      WRITE(ERR,10) TABLE_NAME
      IF (ITABLE < 0) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF

      IF (WRITE_NEU .AND. (ANY_STRE_OUTPUT > 0)) THEN
         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3 stresses
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'Y') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'TRIA3') THEN
                  CALL GET_ELEM_NUM_PLIES ( J )
                  NUM_FROWS = NUM_FROWS + NUM_PLIES
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NUM_FROWS, 9, SUBR_NAME )
            NUM_FROWS = 0
            DO J=1,NELE
               CALL IS_ELEM_PCOMP_PROPS ( J )
               IF (PCOMP_PROPS == 'Y') THEN
                  EID   = EDAT(EPNT(J))
                  TYPE  = ETYPE(J)
                  IF (ETYPE(J)(1:5) == 'TRIA3') THEN
                     CALL GET_ELEM_NUM_PLIES ( J )
                     DO K=0,MBUG-1
                        WRT_BUG(K) = 0
                     ENDDO
                     DO M=1,NUM_PLIES
                        PLY_NUM = M                        ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        NUM_FROWS = NUM_FROWS+ 1
                        FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                        FEMAP_EL_NUMS(NUM_FROWS,2) = M
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE
                        ENDIF
                        CALL ELMDIS
                        IF (NUM_PLIES > 1) THEN
                           CALL ELMDIS_PLY
                        ENDIF
                        CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                        CALL ROT_COMP_ELEM_AXES ( J, M, THETA_PLY, '2-1' )
                        CALL CALC_ELEM_STRESSES ( NCTRIA3, NDUM, NUM_FROWS, 'N', 'Y' )
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            CALL WRITE_FEMAP_STRE_VECS ( 'TRIA3   ', 'Y', NUM_FROWS, FEMAP_SET_ID )
            CALL DEALLOCATE_FEMAP_DATA
         ENDIF

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4 stresses
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'Y') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'QUAD4') THEN
                  CALL GET_ELEM_NUM_PLIES ( J )
                  NUM_FROWS = NUM_FROWS + NUM_PLIES
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NUM_FROWS, 9, SUBR_NAME )
            NUM_FROWS = 0
            DO J=1,NELE
               CALL IS_ELEM_PCOMP_PROPS ( J )
               IF (PCOMP_PROPS == 'Y') THEN
                  EID   = EDAT(EPNT(J))
                  TYPE  = ETYPE(J)
                  IF (ETYPE(J)(1:5) == 'QUAD4') THEN
                     CALL GET_ELEM_NUM_PLIES ( J )
                     NUM_FROWS = NUM_FROWS+ 1
                     DO K=0,MBUG-1
                        WRT_BUG(K) = 0
                     ENDDO
                     DO M=1,NUM_PLIES
                        PLY_NUM = M                        ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                        FEMAP_EL_NUMS(NUM_FROWS,2) = M
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE
                        ENDIF
                        CALL ELMDIS
                        IF (NUM_PLIES > 1) THEN
                           CALL ELMDIS_PLY
                        ENDIF
                        CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                        CALL ROT_COMP_ELEM_AXES ( J, M, THETA_PLY, '2-1' )
                        CALL CALC_ELEM_STRESSES ( NCQUAD4, NDUM, NUM_FROWS, 'N', 'Y' )
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            CALL WRITE_FEMAP_STRE_VECS ( 'QUAD4   ', 'Y', NUM_FROWS, FEMAP_SET_ID )
            CALL DEALLOCATE_FEMAP_DATA
         ENDIF

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out SHEAR stresses
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'Y') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'SHEAR') THEN
                  CALL GET_ELEM_NUM_PLIES ( J )
                  NUM_FROWS = NUM_FROWS + NUM_PLIES
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NUM_FROWS, 9, SUBR_NAME )
            NUM_FROWS = 0
            DO J=1,NELE
               CALL IS_ELEM_PCOMP_PROPS ( J )
               IF (PCOMP_PROPS == 'Y') THEN
                  EID   = EDAT(EPNT(J))
                  TYPE  = ETYPE(J)
                  IF (ETYPE(J)(1:5) == 'SHEAR') THEN
                     CALL GET_ELEM_NUM_PLIES ( J )
                     DO K=0,MBUG-1
                        WRT_BUG(K) = 0
                     ENDDO
                     DO M=1,NUM_PLIES
                        PLY_NUM = M                        ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        NUM_FROWS = NUM_FROWS+ 1
                        FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                        FEMAP_EL_NUMS(NUM_FROWS,2) = M
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE
                        ENDIF
                        CALL ELMDIS
                        IF (NUM_PLIES > 1) THEN
                           CALL ELMDIS_PLY
                        ENDIF
                        CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                        CALL ROT_COMP_ELEM_AXES ( J, M, THETA_PLY, '2-1' )
                        CALL CALC_ELEM_STRESSES ( NCSHEAR, NDUM, NUM_FROWS, 'N', 'Y' )
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            CALL WRITE_FEMAP_STRE_VECS ( 'SHEAR   ', 'Y', NUM_FROWS, FEMAP_SET_ID )
            CALL DEALLOCATE_FEMAP_DATA
         ENDIF

      ENDIF

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT PLY STRESSES'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9192 FORMAT(I8,1X,A,A8,I8,4X,A20)

 9199 FORMAT(' ')

 9201 FORMAT(' *ERROR  9201: DUE TO ABOVE LISTED ERRORS, CANNOT CALCULATE ',A,' REQUESTS FOR ',A,' ELEMENT ID = ',I8)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE GET_STRESS_ITEM_DATA

      IMPLICIT NONE

      INTEGER(LONG)                   :: II               ! DO loop index

! **********************************************************************************************************************************
      DO II=1,18
         STRESS_ITEM(II)(1:) = ' '
      ENDDO

      IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE == 'QUAD4(1:5)')) THEN
         NUM_OTM_ENTRIES = 8
         STRESS_ITEM( 1) = 'Fibre Dist      -Z1 '  ;  STRESS_ITEM( 9) = 'Fibre Dist      +Z1 '
         STRESS_ITEM( 2) = 'Normal x Stress -Z1 '  ;  STRESS_ITEM(10) = 'Normal x Stress +Z1 '
         STRESS_ITEM( 3) = 'Normal y Stress -Z1 '  ;  STRESS_ITEM(11) = 'Normal y Stress +Z1 '
         STRESS_ITEM( 4) = 'Shear xy Stress -Z1 '  ;  STRESS_ITEM(12) = 'Shear xy Stress +Z1 '
         STRESS_ITEM( 5) = 'Princ Angle  at -Z1 '  ;  STRESS_ITEM(13) = 'Princ Angle  at +Z1 '
         STRESS_ITEM( 6) = 'Major Stress at -Z1 '  ;  STRESS_ITEM(14) = 'Major Stress at +Z1 '
         STRESS_ITEM( 7) = 'Minor Stress at -Z1 '  ;  STRESS_ITEM(15) = 'Minor Stress at +Z1 '
         STRESS_ITEM( 8) = 'Max Shear    at -Z1 '  ;  STRESS_ITEM(16) = 'Max Shear    at +Z1 '

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE GET_STRESS_ITEM_DATA

      END SUBROUTINE OFP3_STRE_PCOMP


      SUBROUTINE OFP3_STRN_PCOMP ( JVEC, FEMAP_SET_ID, ITE, OT4_EROW )

! Processes element ply strain output requests for PCOMP elements (TRIA3, QUAD4) for one subcase

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ELOUT_STRN_BIT, FATAL_ERR, IBIT, INT_SC_NUM, MBUG, MOGEL,                   &
                                         NELE, NCQUAD4, NCSHEAR, NCTRIA3, SOL_NAME, WARN_ERR, SOL_NAME
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO
      USE FEMAP_ARRAYS, ONLY          :  FEMAP_EL_NUMS
      USE PARAMS, ONLY                :  OTMSKIP, PRTNEU
      USE MODEL_STUF, ONLY            :  ANY_STRN_OUTPUT, EDAT, EPNT, ETYPE, EID, ELMTYP, ELOUT, METYPE, NUM_EMG_FATAL_ERRS,       &
                                         NUM_PLIES, PCOMP_PROPS, PLY_NUM, THETA_PLY, TYPE
      USE LINK9_STUFF, ONLY           :  EID_OUT_ARRAY, MAXREQ, OGEL
      USE OUTPUT4_MATRICES, ONLY      :  OTM_STRN, TXT_STRN

      USE COMPOSITE_SHELL_PREPARATION, ONLY:  GET_ELEM_NUM_PLIES, IS_ELEM_PCOMP_PROPS, ROT_COMP_ELEM_AXES
      USE EMG_MOD, ONLY               :  EMG
      USE RESULT_COORDINATES, ONLY    :  ELMDIS, ELMDIS_PLY
      USE ELEMENT_RECOVERY_SUPPORT, ONLY: ELEM_STRE_STRN_ARRAYS
      USE ELEMENT_STRAIN_RECOVERY, ONLY: CALC_ELEM_STRAINS
      USE SUSER1_RECOVERY, ONLY       :  SUSER1
      USE DIAGNOSTICS_MEMORY_REPORTING, ONLY:  CHK_OGEL_ZEROS
      USE COMPOSITE_PLY_OUTPUT, ONLY  :  WRITE_PLY_STRAINS
      USE LINK9_WORKSPACE, ONLY       :  ALLOCATE_FEMAP_DATA, DEALLOCATE_FEMAP_DATA
      USE FEMAP_OUTPUT_WRITERS, ONLY  :  WRITE_FEMAP_STRN_VECS

      USE OP2_GEOMETRY_OUTPUT, ONLY   :  END_OP2_TABLE
      USE OP2_STRESS_OUTPUT, ONLY     :  SET_OESC_TABLE_NAME
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'OFP3_STRN_PCOMP'
      CHARACTER( 1*BYTE), PARAMETER   :: IHDR      = 'Y'   ! An input to subr WRITE_GRID_OUTPUTS, called herein
      CHARACTER( 1*BYTE)              :: OPT(6)            ! Option indicators for subr EMG, called herein
      CHARACTER(31*BYTE)              :: OT4_DESCRIPTOR    ! Descriptor for rows of OT4 file
      CHARACTER(30*BYTE)              :: REQUEST           ! Text for error message
      CHARACTER(20*BYTE)              :: STRAIN_ITEM(18)   ! Char description of element strains

      INTEGER(LONG), INTENT(IN)       :: FEMAP_SET_ID      ! Set ID for FEMAP output
      INTEGER(LONG), INTENT(IN)       :: ITE               ! Unit number for text files for OTM row descriptors
      INTEGER(LONG), INTENT(IN)       :: JVEC              ! Solution vector number
      INTEGER(LONG), INTENT(INOUT)    :: OT4_EROW          ! Row number in OT4 file for elem related OTM descriptors
      INTEGER(LONG)                   :: ELOUT_STRN        ! If > 0, there are STRAIN   requests for some elems
      INTEGER(LONG)                   :: I,J,K,L,M         ! DO loop indices
      INTEGER(LONG)                   :: IERROR       = 0  ! Local error count
!xx   INTEGER(LONG)                   :: IROW_MAT          ! Row number in OTM's
!xx   INTEGER(LONG)                   :: IROW_TXT          ! Row number in OTM text file
      INTEGER(LONG)                   :: NDUM              ! Dummy valye needed in call to CALC_ELEM_ENFR_FORCES
      INTEGER(LONG)                   :: NELREQ(METYPE)    ! Count of the no. of requests for ELFORCE(NODE or ENGR) or STRESS
      INTEGER(LONG)                   :: NUM_FROWS         ! No. elems processed for FEMAP
      INTEGER(LONG)                   :: NUM_LINES         ! No. lines of output for ply stresses
      INTEGER(LONG)                   :: NUM_OGEL          ! No. rows written to array OGEL prior to writing results to F06 file
!                                                            (this can be > NUM_ELEM since more than 1 row is written to OGEL
!                                                            for ELFORCE(NODE) - elem nodal forces)
                                                           ! Indicator for output of elem data to BUG file
      INTEGER(LONG)                   :: NUM_OTM_ENTRIES   ! Number of entries in OGEL for a particular element type


      INTEGER(LONG)                   :: ITABLE           ! the op2 subtable number
      CHARACTER(8*BYTE)               :: TABLE_NAME       ! the op2 table name
      LOGICAL                         :: IS_RESULT        ! is there a result
      LOGICAL                         :: WRITE_NEU

      INTRINSIC IAND
! **********************************************************************************************************************************
      ! initial values
      IS_RESULT = .FALSE.
      TABLE_NAME = "OSTR1C"
      ITABLE = 0

      WRITE_NEU = (PRTNEU == 'Y')

! **********************************************************************************************************************************
! Process element strain output (STRAIN) requests for composite shell elements

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

      DO I=1,METYPE
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'Y') THEN
               IF ((ETYPE(J)(1:5) == 'TRIA3') .OR. (ETYPE(J)(1:5) == 'QUAD4') .OR. (ETYPE(J)(1:5) == 'SHEAR')) THEN
                  IF (ETYPE(J) == ELMTYP(I)) THEN
                     ELOUT_STRN = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRN_BIT))
                     IF (ELOUT_STRN > 0) THEN
                        CALL GET_ELEM_NUM_PLIES ( J )
                        NELREQ(I) = NELREQ(I) + NUM_PLIES
                        IS_RESULT = .TRUE.
                     ENDIF
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

      OGEL = ZERO

!xx   IROW_MAT = 0
!xx   IROW_TXT = 0
      OT4_DESCRIPTOR = 'Element strain'
reqs6:DO I=1,METYPE
         IF (NELREQ(I) == 0) CYCLE reqs6
         NUM_LINES = 0
         NUM_OGEL  = 0

elems_6: DO J = 1,NELE

            CALL IS_ELEM_PCOMP_PROPS ( J )
pcomp:      IF (PCOMP_PROPS == 'Y') THEN

               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J) == ELMTYP(I)) THEN
                  ELOUT_STRN = IAND(ELOUT(J,INT_SC_NUM),IBIT(ELOUT_STRN_BIT))
                  IF (ELOUT_STRN > 0) THEN
                     DO K=0,MBUG-1
                        WRT_BUG(K) = 0
                     ENDDO

                     CALL GET_ELEM_NUM_PLIES ( J )
do_plies_6:          DO M=1,NUM_PLIES                         ! Cycle over number of plies, processing strain for each ply
                        PLY_NUM = M                           ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE elems_6
                        ENDIF
                        CALL ELMDIS                           ! Get elem displs (at mid plane of the elem) in elem coords
                        IF (NUM_PLIES > 1) THEN
                           CALL ELMDIS_PLY                    ! If there is more than 1 ply, get displs from elem mid plane displs
                        ENDIF

                        IF (ETYPE(J)(1:5) == 'USER1') THEN
                           CALL SUSER1
                        ELSE                                  ! Get ply strains inply coords
                           CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                           CALL ROT_COMP_ELEM_AXES ( J, M, THETA_PLY, '2-1' )
                           CALL CALC_ELEM_STRAINS ( MAXREQ, NUM_OGEL, J, 'Y', 'N' )
                        ENDIF
                                                              ! If CB soln, write rows of OGEL, from CALC_ELEM_STRAINS, to OTM_STRN
                        IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN

                           CALL GET_STRAIN_ITEM_DATA

                           DO L=1,2
                              DO K=1,NUM_OTM_ENTRIES
                                 OT4_EROW = OT4_EROW + 1
                                 OTM_STRN(OT4_EROW,JVEC) = OGEL(NUM_OGEL-2+L,K)
                                 IF (JVEC == 1) THEN
                                    WRITE(TXT_STRN(OT4_EROW), 9192) OT4_EROW, OT4_DESCRIPTOR, TYPE, EID,                           &
                                                                     STRAIN_ITEM(K+(L-1)*NUM_OTM_ENTRIES)
                                 ENDIF
                              ENDDO
                           ENDDO
                        ENDIF

                        IF ((SOL_NAME(1:12) == 'GEN CB MODEL') .AND. (JVEC == 1) .AND. (OT4_EROW >= 1)) THEN
                           DO K=1,OTMSKIP                        ! Write OTMSKIP blank separator lines
                              OT4_EROW = OT4_EROW + 1
                              WRITE(TXT_STRN(OT4_EROW), 9199)
                           ENDDO
                        ENDIF

                        NUM_LINES = NUM_LINES + 1
                        EID_OUT_ARRAY(NUM_LINES,1) = EID
                        EID_OUT_ARRAY(NUM_LINES,2) = M        ! Ply number for EID
                        IF ((ETYPE(J)(1:5) =='TRIA3') .OR. (ETYPE(J)(1:5) == "QUAD4")) THEN
                           IF (NUM_LINES == NELREQ(I)) THEN
                              CALL CHK_OGEL_ZEROS ( NUM_OGEL )

 100                          FORMAT("*DEBUG:      ",A,"; ELEMENT_TYPE=",A,"; TABLE_NAME=",A,"; ITABLE=",I8)
                              WRITE(ERR,100) "OSTR1C_PCOMP",TYPE,TABLE_NAME,ITABLE
                              CALL SET_OESC_TABLE_NAME(TABLE_NAME, ITABLE)
                              WRITE(ERR,100) "OSTR1C_PCOMP",ETYPE(J)(1:8),TABLE_NAME,ITABLE

                              CALL WRITE_PLY_STRAINS ( JVEC, NUM_LINES, IHDR, ETYPE(J)(1:8), ITABLE  )
                              EXIT
                           ENDIF
                        ENDIF

                     ENDDO do_plies_6

                  ENDIF

               ENDIF

            ENDIF pcomp

         ENDDO elems_6

      ENDDO reqs6

  10   FORMAT("*DEBUG:      OSTR_PCOMP_END:    TABLE_NAME",A)
      WRITE(ERR,10) TABLE_NAME
      IF (ITABLE < 0) THEN
        CALL END_OP2_TABLE(ITABLE)
      ENDIF

      IF (WRITE_NEU .AND. (ANY_STRN_OUTPUT > 0)) THEN

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out TRIA3 strains
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'Y') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'TRIA3') THEN
                  CALL GET_ELEM_NUM_PLIES ( J )
                  NUM_FROWS = NUM_FROWS + NUM_PLIES
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NUM_FROWS, 18, SUBR_NAME )
            NUM_FROWS = 0
            DO J=1,NELE
               CALL IS_ELEM_PCOMP_PROPS ( J )
               IF (PCOMP_PROPS == 'Y') THEN
                  EID   = EDAT(EPNT(J))
                  TYPE  = ETYPE(J)
                  IF (ETYPE(J)(1:5) == 'TRIA3') THEN
                     CALL GET_ELEM_NUM_PLIES ( J )
                     DO K=0,MBUG-1
                        WRT_BUG(K) = 0
                     ENDDO
                     DO M=1,NUM_PLIES
                        PLY_NUM = M                        ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        NUM_FROWS = NUM_FROWS+ 1
                        FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                        FEMAP_EL_NUMS(NUM_FROWS,2) = M
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE
                        ENDIF
                        CALL ELMDIS
                        IF (NUM_PLIES > 1) THEN
                           CALL ELMDIS_PLY
                        ENDIF
                        CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                        CALL ROT_COMP_ELEM_AXES ( J, M, THETA_PLY, '2-1' )
                        CALL CALC_ELEM_STRAINS ( NCTRIA3, NDUM, NUM_FROWS, 'N', 'Y' )
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            CALL WRITE_FEMAP_STRN_VECS ( 'TRIA3   ', 'Y', NUM_FROWS, FEMAP_SET_ID )
            CALL DEALLOCATE_FEMAP_DATA
         ENDIF

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out QUAD4 strains
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'Y') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'QUAD4') THEN
                  CALL GET_ELEM_NUM_PLIES ( J )
                  NUM_FROWS = NUM_FROWS + NUM_PLIES
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NUM_FROWS, 18, SUBR_NAME )
            NUM_FROWS = 0
            DO J=1,NELE
               CALL IS_ELEM_PCOMP_PROPS ( J )
               IF (PCOMP_PROPS == 'Y') THEN
                  EID   = EDAT(EPNT(J))
                  TYPE  = ETYPE(J)
                  IF (ETYPE(J)(1:5) == 'QUAD4') THEN
                     CALL GET_ELEM_NUM_PLIES ( J )
                     DO K=0,MBUG-1
                        WRT_BUG(K) = 0
                     ENDDO
                     DO M=1,NUM_PLIES
                        PLY_NUM = M                        ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        NUM_FROWS= NUM_FROWS+ 1
                        FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                        FEMAP_EL_NUMS(NUM_FROWS,2) = M
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE
                        ENDIF
                        CALL ELMDIS
                        IF (NUM_PLIES > 1) THEN
                           CALL ELMDIS_PLY
                        ENDIF
                        CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                        CALL ROT_COMP_ELEM_AXES ( J, M, THETA_PLY, '2-1' )
                        CALL CALC_ELEM_STRAINS ( NCQUAD4, NDUM, NUM_FROWS, 'N', 'Y' )
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            CALL WRITE_FEMAP_STRN_VECS ( 'QUAD4   ', 'Y', NUM_FROWS, FEMAP_SET_ID )
            CALL DEALLOCATE_FEMAP_DATA
         ENDIF

         NDUM = 0
         NUM_FROWS= 0                                      ! Write out SHEAR strains
         DO J=1,NELE
            CALL IS_ELEM_PCOMP_PROPS ( J )
            IF (PCOMP_PROPS == 'Y') THEN
               EID   = EDAT(EPNT(J))
               TYPE  = ETYPE(J)
               IF (ETYPE(J)(1:5) == 'SHEAR') THEN
                  CALL GET_ELEM_NUM_PLIES ( J )
                  NUM_FROWS = NUM_FROWS + NUM_PLIES
               ENDIF
            ENDIF
         ENDDO
         IF (NUM_FROWS > 0) THEN
            CALL ALLOCATE_FEMAP_DATA ( 'FEMAP ELEM ARRAYS', NUM_FROWS, 18, SUBR_NAME )
            NUM_FROWS = 0
            DO J=1,NELE
               CALL IS_ELEM_PCOMP_PROPS ( J )
               IF (PCOMP_PROPS == 'Y') THEN
                  EID   = EDAT(EPNT(J))
                  TYPE  = ETYPE(J)
                  IF (ETYPE(J)(1:5) == 'SHEAR') THEN
                     CALL GET_ELEM_NUM_PLIES ( J )
                     DO K=0,MBUG-1
                        WRT_BUG(K) = 0
                     ENDDO
                     DO M=1,NUM_PLIES
                        PLY_NUM = M                        ! 'N' in call to EMG means do not write to BUG file
                        CALL EMG ( J   , OPT, 'N', SUBR_NAME, 'N' )
                        NUM_FROWS = NUM_FROWS+ 1
                        FEMAP_EL_NUMS(NUM_FROWS,1) = EID
                        FEMAP_EL_NUMS(NUM_FROWS,2) = M
                        IF (NUM_EMG_FATAL_ERRS > 0) THEN
                           IERROR = IERROR + 1
                           CYCLE
                        ENDIF
                        CALL ELMDIS
                        IF (NUM_PLIES > 1) THEN
                           CALL ELMDIS_PLY
                        ENDIF
                        CALL ELEM_STRE_STRN_ARRAYS ( 1 )
                        CALL ROT_COMP_ELEM_AXES ( J, M, THETA_PLY, '2-1' )
                        CALL CALC_ELEM_STRAINS ( NCSHEAR, NDUM, NUM_FROWS, 'N', 'Y' )
                     ENDDO
                  ENDIF
               ENDIF
            ENDDO
            CALL WRITE_FEMAP_STRN_VECS ( 'SHEAR   ', 'Y', NUM_FROWS, FEMAP_SET_ID )
            CALL DEALLOCATE_FEMAP_DATA
         ENDIF

      ENDIF

      IF (IERROR > 0) THEN
         REQUEST = 'ELEMENT PLY STRAINS'
         WRITE(ERR,9201) TYPE, REQUEST, EID
         WRITE(F06,9201) TYPE, REQUEST, EID
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9192 FORMAT(I8,1X,A,A8,I8,4X,A20)

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

      IF      ((TYPE(1:5) == 'TRIA3') .OR. (TYPE == 'QUAD4(1:5)') .OR. (ETYPE(J)(1:5) == 'SHEAR')) THEN
         NUM_OTM_ENTRIES = 8
         STRAIN_ITEM( 1) = 'Fibre Dist      -Z1 '  ;  STRAIN_ITEM( 9) = 'Fibre Dist      +Z1 '
         STRAIN_ITEM( 2) = 'Normal x Strain -Z1 '  ;  STRAIN_ITEM(10) = 'Normal x Strain +Z1 '
         STRAIN_ITEM( 3) = 'Normal y Strain -Z1 '  ;  STRAIN_ITEM(11) = 'Normal y Strain +Z1 '
         STRAIN_ITEM( 4) = 'Shear xy Strain -Z1 '  ;  STRAIN_ITEM(12) = 'Shear xy Strain +Z1 '
         STRAIN_ITEM( 5) = 'Princ Angle  at -Z1 '  ;  STRAIN_ITEM(13) = 'Princ Angle  at +Z1 '
         STRAIN_ITEM( 6) = 'Major Strain at -Z1 '  ;  STRAIN_ITEM(14) = 'Major Strain at +Z1 '
         STRAIN_ITEM( 7) = 'Minor Strain at -Z1 '  ;  STRAIN_ITEM(15) = 'Minor Strain at +Z1 '
         STRAIN_ITEM( 8) = 'Max Shear    at -Z1 '  ;  STRAIN_ITEM(16) = 'Max Shear    at +Z1 '

      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE GET_STRAIN_ITEM_DATA

      ENDSUBROUTINE OFP3_STRN_PCOMP

   END MODULE COMPOSITE_RESULT_OUTPUT
