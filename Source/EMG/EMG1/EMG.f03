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

   MODULE EMG_MOD

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: EMG

   CONTAINS

      SUBROUTINE EMG ( INT_ELEM_ID, OPT, WRITE_WARN, CALLING_SUBR, WRT_BUG_THIS_TIME )

! Main driver routine for calculation of matrices for all elements. This routine initializes appropriate arrays and calls
! other routines to calculate element matrices:

!  1) ME             = element mass matrix                             , if OPT(1) = 'Y'
!  2) PTE            = element thermal load vectors                    , if OPT(2) = 'Y'
!  3) SEi, STEi, BEi = element stress and strain data recovery matrices, if OPT(3) = 'Y'
!  4) KE             = element linea stiffness matrix                  , if OPT(4) = 'Y'
!  5) PPE            = element pressure load matrix                    , if OPT(5) = 'Y'
!  6) KED            = element differen stiff matrix calc              , if OPT(6) = 'Y'

! Also, calculate:

!  TE        = Coord transformation matrix (basic to elem)
!  ZS        = Stress recovery coeff's
!  FCONV     = Constants to convert stress to engineering force
!              NOTE: may need to calc KE to get this

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MBUG, MEDAT0_CUSERIN, MELDOF, MEMATC, MOFFSET, NSUB, NTSUB
      USE SCONTR, ONLY                :  DEDAT_Q4_MATANG_KEY, DEDAT_Q4_POFFS_KEY, DEDAT_Q4_SHELL_KEY, DEDAT_Q4_THICK_KEY,          &
                                         DEDAT_T3_MATANG_KEY, DEDAT_T3_POFFS_KEY, DEDAT_T3_SHELL_KEY, DEDAT_T3_THICK_KEY,          &
                                         WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  SUPINFO, SUPWARN, QUAD4TYP
      USE CONSTANTS_1, ONLY           :  CONV_DEG_RAD, CONV_RAD_DEG, ZERO, ONE
      USE MODEL_STUF, ONLY            :  CAN_ELEM_TYPE_OFFSET, EDAT, EID, EPNT, ETYPE, ISOLID, MATANGLE, NUM_EMG_FATAL_ERRS,       &
                                         PCOMP_PROPS, PLY_NUM, TE_IDENT, THETAM, TYPE, XEL, TE

      USE COMPOSITE_SHELL_PREPARATION, ONLY:  IS_ELEM_PCOMP_PROPS, SHELL_ABD_MATRICES
      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE ELEMENT_DATA_GATHERING, ONLY:  ELMDAT1, ELMDAT2
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE ELEMENT_GEOMETRY_PREPARATION, ONLY:  ELMGM1, ELMGM2, ELMGM3
      USE ELEMENT_LOOKUPS, ONLY       :  GET_MATANGLE_FROM_CID
      USE MATERIAL_PROPERTIES, ONLY   :  MATERIAL_PROPS_2D, MATERIAL_PROPS_3D
      USE MATERIAL_TRANSFORMATIONS, ONLY:  ROT_AXES_MATL_TO_LOC
      USE ELEMENT_DIAGNOSTICS, ONLY   :  ELMOFF, ELMOUT
      USE SPRING_ELEMENTS, ONLY       :  BUSH, ELAS1
      USE LINE_ELEMENTS, ONLY         :  BREL1
      USE TRIANGULAR_SHELL_ELEMENTS, ONLY:  TREL1
      USE QUADRILATERAL_ELEMENT_DISPATCH, ONLY:  QDEL1
      USE SOLID_ELEMENTS, ONLY        :  HEXA, PENTA, TETRA
      USE USER_DEFINED_ELEMENTS, ONLY :  KUSER1, USERIN
      USE MITC8_MOD, ONLY             :  MITC8
      USE MITC4_MOD, ONLY:  MITC4
      USE ELEMENT_GEOMETRY_PREPARATION, ONLY:  ELMGM1_BUSH
      USE MITC8_MOD, ONLY             :  MITC8
      USE SHELL_ROTATION_SUPPORT, ONLY:  CALC_K6ROT
      USE MITC4_MOD, ONLY             :  MITC4
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'EMG'
      CHARACTER(LEN=*)  , INTENT(IN)  :: CALLING_SUBR       ! Name of subr that called this one
      CHARACTER( 1*BYTE), INTENT(IN)  :: OPT(6)             ! Array of EMG option indicators explained above
      CHARACTER(LEN=*)  , INTENT(IN)  :: WRITE_WARN         ! If 'Y" write warning messages, otherwise do not
      CHARACTER( 1*BYTE), INTENT(IN)  :: WRT_BUG_THIS_TIME  ! If 'Y' then write to BUG file if WRT_BUG array says to
      CHARACTER( 2*BYTE)              :: LOC                ! Location where THETAM is calculated (for DEBUG output purposes)
      CHARACTER( 1*BYTE)              :: FIX_EDAT      = 'N'! If 'Y', run code to change order of grids in EDAT for 3D elems
      CHARACTER( 1*BYTE)              :: RED_INT_SHEAR = 'N'! If 'Y', use Gaussian weighted average of B matrices for shear terms

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID        ! Internal element ID for which
      INTEGER(LONG)                   :: CASE_NUM    = 0    ! Can be subcase number (e.g. for UEL, PEL outout)
      INTEGER(LONG)                   :: DUM_BUG(0:MBUG-1)  ! Values from WRT_BUG sent to subr ELMOUT in a particular call
      INTEGER(LONG)                   :: EPNTK              ! Value from array EPNT at the row for this internal elem ID. It is the
!                                                             row number in array EDAT where data begins for this element.
      INTEGER(LONG)                   :: I                  ! DO loop index
      INTEGER(LONG)                   :: INT_ORDER   = 0    ! Gaussian integration order for element
      INTEGER(LONG)                   :: IORD_IJ            ! Integration order in the triangular plane for PENTA elements
      INTEGER(LONG)                   :: IORD_K             ! Integration order in Z direction for PENTA elements
      INTEGER(LONG)                   :: INT41,INT42        ! An integer used in getting MATANGLE


! **********************************************************************************************************************************
      EPNTK = EPNT(INT_ELEM_ID)
      EID   = EDAT(EPNTK)
      TYPE  = ETYPE(INT_ELEM_ID)
      CALL IS_ELEM_PCOMP_PROPS ( INT_ELEM_ID )

      NUM_EMG_FATAL_ERRS = 0


      CALL EMG_INIT

! **********************************************************************************************************************************
! Call ELMDAT1 subr to get some of the data needed for this elem.

      IF ((TYPE == 'ELAS1   ') .OR. (TYPE == 'ELAS2   ') .OR. (TYPE == 'ELAS3   ') .OR. (TYPE == 'ELAS4   ') .OR.                  &
          (TYPE == 'ROD     ') .OR. (TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'BUSH    ') .OR.                  &
          (TYPE == 'HEXA8   ') .OR. (TYPE == 'HEXA20  ') .OR.                                                                      &
          (TYPE == 'PENTA6  ') .OR. (TYPE == 'PENTA15 ') .OR.                                                                      &
          (TYPE == 'TETRA4  ') .OR. (TYPE == 'TETRA10 ') .OR.                                                                      &
          (TYPE == 'USER1   ') .OR. (TYPE == 'USERIN  ') .OR. (TYPE == 'PLOTEL  ') .OR.                                            &
          (TYPE == 'SHEAR   ') .OR. (TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4'   ) .OR. (TYPE(1:5) == 'QUAD8'   )) THEN
         CALL ELMDAT1 ( INT_ELEM_ID, WRITE_WARN )
      ELSE
         WRITE(ERR,1916) SUBR_NAME,EID,TYPE
         WRITE(F06,1916) SUBR_NAME,EID,TYPE
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Coding error, so quit
      ENDIF

      IF (NUM_EMG_FATAL_ERRS > 0)   CALL EMG_QUIT

! **********************************************************************************************************************************
! For all elements but ELASi, USERIN, get routine to check element geometry, get element grid point coords in local elem
! coord system, and generate the coord transformation matrix to be used in transforming the stiffness matrices in
! elem system to basic sys. If elem is 'ELAS' the calcd stiff matrix IS in global coord's, so TE matrix is not needed.

!xx   IF      ((TYPE      == 'ROD     ') .OR. (TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'BUSH    ') .OR.        &
      IF      ((TYPE      == 'ROD     ') .OR. (TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR.                                  &
               (TYPE(1:5) == 'TRIA3'   ) .OR.                                                                                      &
               (TYPE      == 'PENTA6  ') .OR. (TYPE == 'PENTA15 ') .OR.                                                            &
               (TYPE      == 'TETRA4  ') .OR. (TYPE == 'TETRA10 ')) THEN
         CALL ELMGM1 ( INT_ELEM_ID, WRITE_WARN )
         FIX_EDAT = 'N'

         IF (TYPE(1:5) == 'PENTA') THEN
            IF (XEL(4,3) < ZERO) THEN
               FIX_EDAT = 'Y'
            ENDIF
         ELSE IF (TYPE(1:5) == 'TETRA') THEN
            IF (XEL(4,3) < ZERO) THEN
               FIX_EDAT = 'Y'
            ENDIF
         ENDIF

         IF (FIX_EDAT == 'Y') THEN
            CALL EDAT_FIXUP
            CALL ELMDAT1 ( INT_ELEM_ID, WRITE_WARN )
            CALL ELMGM1 ( INT_ELEM_ID, WRITE_WARN )
         ENDIF

      ELSE IF (TYPE == 'BUSH    ') THEN
         CALL ELMGM1_BUSH ( INT_ELEM_ID, WRITE_WARN )

      ELSE IF ((TYPE == 'QUAD4   ') .OR. (TYPE == 'QUAD4K  ') .OR. (TYPE == 'QUAD8   ') .OR. (TYPE == 'SHEAR   ')) THEN
         CALL ELMGM2 ( WRITE_WARN )

      ELSE IF ((TYPE == 'HEXA8   ') .OR. (TYPE == 'HEXA20  ')) THEN
         CALL ELMGM3 ( WRITE_WARN )
         FIX_EDAT = 'N'

         IF      (TYPE(1:4) == 'HEXA') THEN                ! See if we need to reorder grids to get local z in proper direction
            IF (XEL(5,3) < ZERO) THEN
               FIX_EDAT = 'Y'
            ENDIF
         ENDIF

         IF (FIX_EDAT == 'Y') THEN
            CALL EDAT_FIXUP
            CALL ELMDAT1 ( INT_ELEM_ID, WRITE_WARN )
            CALL ELMGM3 ( WRITE_WARN )
         ENDIF

      ELSE IF  (TYPE(1:4) == 'ELAS') THEN
         TE_IDENT = 'Y'

      ELSE IF ((TYPE == 'USER1   ') .OR. (TYPE == 'USERIN  ')) THEN
         CONTINUE

      ENDIF

      IF (NUM_EMG_FATAL_ERRS > 0)   CALL EMG_QUIT

! **********************************************************************************************************************************
! Generate element material matrices. Material matrices for shell elements with PCOMP props are generated elsewhere.
! Matrices of material props are not generated for 1-D elements
! --------

      IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE == 'SHEAR   ')) THEN
         IF (PCOMP_PROPS == 'N') THEN                      ! SHEAR elem does not use PCOMP props

            THETAM = ZERO

            IF      (TYPE(1:5) == 'QUAD4') THEN
               INT41 = EDAT(EPNTK+DEDAT_Q4_MATANG_KEY)     ! Key to say whether matl angle is an actual angle or a coord ID
               INT42 = EDAT(EPNTK+DEDAT_Q4_MATANG_KEY+1)   ! Key to say (if INT41 is neg) whether coord ID is basic or otherwise
               IF      (INT41 >  0) THEN                   ! Angle is defined in array MATANGLE at row INT41
                  LOC = '#1'
                  THETAM = CONV_DEG_RAD*MATANGLE( INT41 )
               ELSE IF (INT41 <  0) THEN                   ! Angle is defined by a coord sys ID whose value is -INT41
                  LOC = '#2'
                  CALL GET_MATANGLE_FROM_CID ( -INT41 )
               ELSE IF (INT41 ==  0) THEN                  ! Angle is either specified as defined by basic coord sys or angle is 0.
                  IF      (INT42 == 1) THEN
                     LOC = '#3'
                     CALL GET_MATANGLE_FROM_CID ( 0 )
                  ELSE IF (INT42 == 0) THEN
                     LOC = '#4'
                     THETAM = ZERO
                  ENDIF
               ENDIF

            ELSE IF (TYPE(1:5) == 'TRIA3') THEN
               INT41 = EDAT(EPNTK+DEDAT_T3_MATANG_KEY)     ! Key to say whether matl angle is an actual angle or a coord ID
               INT42 = EDAT(EPNTK+DEDAT_T3_MATANG_KEY+1)   ! Key to say (if INT41 is neg) whether coord ID is basic or otherwise
               IF      (INT41 >  0) THEN                   ! Angle is defined in array MATANGLE at row INT41
                  LOC = '#1'
                  THETAM = CONV_DEG_RAD*MATANGLE( INT41 )
               ELSE IF (INT41 <  0) THEN                   ! Angle is defined by a coord sys ID whose value is -INT41
                  LOC = '#2'
                  CALL GET_MATANGLE_FROM_CID ( -INT41 )
               ELSE IF (INT41 ==  0) THEN                  ! Angle is either specified as defined by basic coord sys or angle is 0.
                  IF      (INT42 == 1) THEN
                     LOC = '#3'
                     CALL GET_MATANGLE_FROM_CID ( 0 )
                  ELSE IF (INT42 == 0) THEN
                     LOC = '#4'
                     THETAM = ZERO
                  ENDIF
               ENDIF

            ENDIF
                                                           ! Use WRITE_WARN even though the following is not a warning message
!                                                            this will allow THETAM to be printed out in only 1 call to EMG
            IF ((DEBUG(112) > 0) .AND. (WRITE_WARN == 'Y')) THEN
               IF (INT41 > 0) THEN
                  WRITE(F06,1001) TYPE, EID, CONV_RAD_DEG*THETAM, LOC, INT41, INT42
               ELSE
                  WRITE(F06,1002) TYPE, EID, CONV_RAD_DEG*THETAM, LOC, INT41, INT42
               ENDIF
            ENDIF

            CALL MATERIAL_PROPS_2D ( WRITE_WARN )
            CALL ROT_AXES_MATL_TO_LOC ( WRITE_WARN )

         ELSE

            IF (TYPE(1:5) == 'QUAD4') THEN
               INT41 = EDAT(EPNTK+DEDAT_Q4_MATANG_KEY)
               INT42 = EDAT(EPNTK+DEDAT_Q4_MATANG_KEY+1)
            ENDIF

            IF (TYPE(1:5) == 'TRIA3') THEN
               INT41 = EDAT(EPNTK+DEDAT_T3_MATANG_KEY)
               INT42 = EDAT(EPNTK+DEDAT_T3_MATANG_KEY+1)
            ENDIF

            IF (WRITE_WARN == 'Y') THEN
               IF ((INT41 /= 0) .OR. (INT42 /= 0)) THEN

                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,1003) TYPE, EID
                  IF (SUPWARN == 'N') THEN
                     WRITE(F06,1003) TYPE, EID
                  ENDIF
               ENDIF
            ENDIF

         ENDIF

      ENDIF

      IF (TYPE == 'QUAD8   ') THEN

        CALL MATERIAL_PROPS_2D ( WRITE_WARN )
        ! Don't transform the material properties here because the transformation is different at each point in the element.

      ENDIF

      IF ((TYPE == 'HEXA8   ') .OR. (TYPE == 'HEXA20  ') .OR.                                                                      &
          (TYPE == 'PENTA6  ') .OR. (TYPE == 'PENTA15 ') .OR.                                                                      &
          (TYPE == 'TETRA4  ') .OR. (TYPE == 'TETRA10 ')) THEN
         CALL MATERIAL_PROPS_3D ( WRITE_WARN )
         CALL ROT_AXES_MATL_TO_LOC ( WRITE_WARN )
      ENDIF

! Call ELMOUT to output element data for item 0, 1

      IF (WRT_BUG(0) > 0) THEN
         CASE_NUM = 0
         DO I=0,MBUG-1
            DUM_BUG(I) = 0
         ENDDO
         DUM_BUG(0) = 1
         CALL ELMOUT ( INT_ELEM_ID, DUM_BUG, CASE_NUM, OPT )
      ENDIF

      IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4') .OR. (TYPE(1:5) == 'QUAD8') .OR. (TYPE(1:6) == 'SHEAR') .OR.          &
          (TYPE == 'USER1   ')) THEN
         CALL SHELL_ABD_MATRICES ( INT_ELEM_ID, WRITE_WARN )
      ENDIF

      IF (WRT_BUG(1) > 0) THEN
         CASE_NUM = 0
         DO I=0,MBUG-1
            DUM_BUG(I) = 0
         ENDDO
         DUM_BUG(1) = 1
         CALL ELMOUT ( INT_ELEM_ID, DUM_BUG, CASE_NUM, OPT )
      ENDIF

! Quick return if all OPT are 'N' (if we only need ELMDAT, ELMGMi called, and not all the other subr's which generate elem matrices

      IF ((OPT(1) == 'N') .AND. (OPT(2) == 'N') .AND. (OPT(3) == 'N') .AND. (OPT(4) == 'N') .AND. (OPT(5) == 'N') .AND.            &
          (OPT(6) == 'N')) THEN
         RETURN
      ENDIF

! **********************************************************************************************************************************
! For all but USERIN elem, call ELMDAT2 subr to get the rest of the data needed to calculate the matrices for this element.

      IF ((TYPE(1:4) == 'ELAS'    ) .OR. (TYPE      == 'ROD     ') .OR. (TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR.        &
          (TYPE(1:5) == 'TRIA3'   ) .OR. (TYPE(1:5) == 'QUAD4'   ) .OR. (TYPE == 'SHEAR   ') .OR. (TYPE == 'USER1   ') .OR.        &
          (TYPE      == 'HEXA8   ') .OR. (TYPE      == 'HEXA20  ') .OR.                                                            &
          (TYPE      == 'PENTA6  ') .OR. (TYPE      == 'PENTA15 ') .OR.                                                            &
          (TYPE      == 'TETRA4  ') .OR. (TYPE      == 'TETRA10 ')) THEN
         CALL ELMDAT2 ( INT_ELEM_ID, OPT, WRITE_WARN )
      ENDIF

      IF (NUM_EMG_FATAL_ERRS > 0)   CALL EMG_QUIT

! **********************************************************************************************************************************
! Now get the individual elem routines to calc the required elem matrices: ME and/or PTE and/or (SE1, SE2, STE1,STE2)
! and/or KE).

      IF (TYPE(1:4) == 'ELAS') THEN
         CALL ELAS1 ( OPT, WRITE_WARN )

      ELSE IF ((TYPE == 'ROD     ') .OR. (TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ')) THEN
         CALL BREL1 ( OPT, WRITE_WARN )
         IF (NUM_EMG_FATAL_ERRS > 0)   CALL EMG_QUIT

      ELSE IF (TYPE(1:4) == 'BUSH') THEN
         CALL BUSH ( INT_ELEM_ID, OPT, WRITE_WARN )
         IF (NUM_EMG_FATAL_ERRS > 0)   CALL EMG_QUIT

      ELSE IF (TYPE(1:5) == 'TRIA3') THEN
         CALL TREL1 ( OPT, WRITE_WARN )
         IF (NUM_EMG_FATAL_ERRS > 0)   CALL EMG_QUIT

      ELSE IF (((TYPE == 'QUAD4   ') .AND. ((QUAD4TYP == 'MIN4  ') .OR. (QUAD4TYP == 'MIN4T '))) .OR.                              &
                (TYPE == 'QUAD4K  ') .OR.                                                                                          &
                (TYPE == 'SHEAR   ')) THEN
         CALL QDEL1 ( OPT, INT_ELEM_ID, WRITE_WARN )
         IF (NUM_EMG_FATAL_ERRS > 0)   CALL EMG_QUIT

      ELSE IF ((TYPE == 'QUAD4   ') .AND. ((QUAD4TYP == 'MITC4 ') .OR. (QUAD4TYP == 'MITC4+'))) THEN
         CALL MITC4 ( OPT, INT_ELEM_ID )
         IF (NUM_EMG_FATAL_ERRS > 0)   CALL EMG_QUIT

      ELSE IF (TYPE(1:5) == 'QUAD8') THEN
         CALL MITC8 ( OPT, INT_ELEM_ID )
         IF (NUM_EMG_FATAL_ERRS > 0)   CALL EMG_QUIT

      ELSE IF ((TYPE == 'HEXA8   ') .OR. (TYPE == 'HEXA20  ') .OR.                                                                 &
               (TYPE == 'PENTA6  ') .OR. (TYPE == 'PENTA15 ') .OR.                                                                 &
               (TYPE == 'TETRA4  ') .OR. (TYPE == 'TETRA10 ')) THEN

         IF (ISOLID(6) == 0) THEN                          ! Integration scheme
            RED_INT_SHEAR = 'Y'
         ELSE
            RED_INT_SHEAR = 'N'
         ENDIF

         IF ((TYPE == 'HEXA8   ') .OR. (TYPE == 'HEXA20  ')) THEN
            INT_ORDER = ISOLID(4)                          ! Integration order
            CALL HEXA ( OPT, INT_ELEM_ID, INT_ORDER, RED_INT_SHEAR, WRITE_WARN )
         ELSE IF ((TYPE == 'PENTA6  ') .OR. (TYPE == 'PENTA15 ')) THEN
            IF (ISOLID(4) == 2) THEN                       ! Integration order
               IORD_K  = 2
               IORD_IJ = 3
            ELSE
               IORD_K  = 3
               IORD_IJ = 7
            ENDIF
            CALL PENTA ( OPT, INT_ELEM_ID, IORD_IJ, IORD_K, RED_INT_SHEAR, WRITE_WARN )
         ELSE IF ((TYPE == 'TETRA4  ') .OR. (TYPE == 'TETRA10 ')) THEN
            IF (ISOLID(4) == 2) THEN                       ! Integration order
               INT_ORDER  = 1
            ELSE
               INT_ORDER  = 4
            ENDIF
            CALL TETRA ( OPT, INT_ELEM_ID, INT_ORDER, RED_INT_SHEAR, WRITE_WARN )
         ENDIF


      ELSE IF (TYPE == 'USER1   ') THEN
         CALL KUSER1 ( OPT, WRITE_WARN )

      ELSE IF (TYPE == 'USERIN  ') THEN
         CALL USERIN ( INT_ELEM_ID, OPT, CALLING_SUBR, WRITE_WARN )

      ENDIF


! **********************************************************************************************************************************
! For plate elements, process offsets (since they are specified in local element coordinates)
! Then add drilling stiffness at the the grid points.

      IF ((TYPE(1:5) == 'TRIA3') .OR. (TYPE(1:5) == 'QUAD4')) THEN

         CALL ELMOFF ( OPT, WRITE_WARN )

         IF (OPT(4) == 'Y') THEN
            CALL CALC_K6ROT()
         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Call ELMOUT to output for data items 2-5

      IF (WRT_BUG_THIS_TIME == 'Y') THEN

         IF ((WRT_BUG(2) > 0) .OR. (WRT_BUG(3) > 0) .OR. (WRT_BUG(4) > 0) .OR. (WRT_BUG(5) > 0)) THEN
            CASE_NUM = 0
            DO I=0,MBUG-1
               DUM_BUG(I) = 0
            ENDDO
            DUM_BUG(2) = WRT_BUG(2)
            DUM_BUG(3) = WRT_BUG(3)
            DUM_BUG(4) = WRT_BUG(4)
            DUM_BUG(5) = WRT_BUG(5)
            CALL ELMOUT ( INT_ELEM_ID, DUM_BUG, CASE_NUM, OPT )
         ENDIF

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1916 FORMAT(' *ERROR  1916: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT',I8,' IS TYPE ',A,'. NO CODE WRITTEN FOR THIS TYPE.')

 1001 FORMAT(' THETAM material angle for ',a,i8,' = ',f8.3,' deg (at code location ',a,' in subr EMG: Row number in array MATANGLE'&
            ,' INT42 = ',I8,I3)

 1002 FORMAT(' THETAM material angle for ',a,i8,' = ',f8.3,' deg (at code location ',a,' in subr EMG: Material angle coord system' &
            ,' INT42 = ',I8,I3)

 1003 FORMAT(' WARNING     : MATERIAL ANGLE FOR ',A,I8,' IS NOT USED FOR PLATE ELEMENTS WITH PCOMP PROPERTIES')

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE EMG_INIT

      USE SCONTR, ONLY                :  LSUB, MELDOF, MELGP, MDT, MPRESS, MAX_STRESS_POINTS
      USE MODEL_STUF, ONLY            :  AGRID, BE1, BE2, BE3, BGRID, DOFPIN, DT, ELAS_COMP, FCONV, KE, KED, ME,                   &
                                         OFFDIS, OFFSET, PEB, PEG, PEL, PPE, PRESS, PTE, SE1, SE2, SE3, STE1, STE2, STE3,          &
                                         UEB, UEG, UEL, UGG, XEB, XEL
      IMPLICIT NONE

      INTEGER(LONG)                   :: II,JJ,KK

! **********************************************************************************************************************************
                                                           !  1
      IF (ALLOCATED(AGRID)) THEN
         AGRID=0
      ENDIF
                                                           !  2
      IF (ALLOCATED(BE1)) THEN
         BE1=ZERO
      ENDIF
                                                           !  3
      IF (ALLOCATED(BE2)) THEN
         BE2=ZERO
      ENDIF
                                                           !  4
      IF (ALLOCATED(BE3)) THEN
         BE3=ZERO
      ENDIF
                                                           !  5
      IF (ALLOCATED(BGRID)) THEN
         BGRID=0
      ENDIF
                                                           !  6
      IF (ALLOCATED(DOFPIN)) THEN
         DOFPIN = 0
      ENDIF
                                                           !  7
      IF (ALLOCATED(DT)) THEN
         DT = 0
      ENDIF
                                                           !  8
      IF (ALLOCATED(KE)) THEN
         KE = ZERO
      ENDIF
                                                           !  9
      IF (ALLOCATED(KED)) THEN
               KED = ZERO
      ENDIF
                                                           ! 10
      IF (ALLOCATED(ME)) THEN
               ME = ZERO
      ENDIF
                                                           ! 11
      IF (ALLOCATED(OFFDIS)) THEN
               OFFDIS = ZERO
      ENDIF
                                                           ! 12
      IF (ALLOCATED(OFFSET)) THEN
            OFFSET = ' '
      ENDIF
                                                           ! 13
      IF (ALLOCATED(PEB)) THEN
            PEB = ZERO
      ENDIF
                                                           ! 14
      IF (ALLOCATED(PEG)) THEN
            PEG = ZERO
      ENDIF
                                                           ! 15
      IF (ALLOCATED(PEL)) THEN
            PEL = ZERO
      ENDIF
                                                           ! 16
      IF (ALLOCATED(PPE)) THEN
               PPE = ZERO
      ENDIF
                                                           ! 17
      IF (ALLOCATED(PRESS)) THEN
               PRESS = ZERO
      ENDIF
                                                           ! 18
      IF (ALLOCATED(PTE)) THEN
         PTE = ZERO
      ENDIF
                                                           ! 19
      IF (ALLOCATED(SE1)) THEN
         SE1 = ZERO
      ENDIF
                                                           ! 20
      IF (ALLOCATED(SE2)) THEN
         SE2 = ZERO
      ENDIF
                                                           ! 21
      IF (ALLOCATED(SE3)) THEN
         SE3 = ZERO
      ENDIF
                                                           ! 22
      IF (ALLOCATED(STE1)) THEN
         STE1 = ZERO
      ENDIF
                                                           ! 23
      IF (ALLOCATED(STE2)) THEN
         STE2 = ZERO
      ENDIF
                                                           ! 24
      IF (ALLOCATED(STE3)) THEN
         STE3 = ZERO
      ENDIF
                                                           ! 25
      IF (ALLOCATED(UEB)) THEN
            UEB = ZERO
      ENDIF
                                                           ! 26
      IF (ALLOCATED(UEG)) THEN
            UEG = ZERO
      ENDIF
                                                           ! 27
      IF (ALLOCATED(UEL)) THEN
            UEL = ZERO
      ENDIF
                                                           ! 28
      IF (ALLOCATED(UGG)) THEN
            UGG = ZERO
      ENDIF
                                                           ! 29
      IF (ALLOCATED(XEB)) THEN
               XEB = ZERO
      ENDIF
                                                           ! 30
      IF (ALLOCATED(XEL)) THEN
               XEL = ZERO
      ENDIF

      FCONV = ZERO

! **********************************************************************************************************************************

      END SUBROUTINE EMG_INIT

! ##################################################################################################################################

      SUBROUTINE EDAT_FIXUP

      USE MODEL_STUF, ONLY            :  ELGP

      IMPLICIT NONE

      INTEGER(LONG)                   :: DUM_AGRID(ELGP)   ! List og grids to which an elem is attached

! **********************************************************************************************************************************

      EPNTK = EPNT(INT_ELEM_ID)
      TYPE  = ETYPE(INT_ELEM_ID)
      EID   = EDAT(EPNTK)

      DO I=1,ELGP
         DUM_AGRID(I) = EDAT(EPNTK+I)
      ENDDO

      IF      (TYPE == 'HEXA8   ') THEN

         DUM_AGRID(1)  = EDAT(EPNTK+2)
         DUM_AGRID(2)  = EDAT(EPNTK+3)
         DUM_AGRID(3)  = EDAT(EPNTK+4)
         DUM_AGRID(4)  = EDAT(EPNTK+5)
         DUM_AGRID(5)  = EDAT(EPNTK+6)
         DUM_AGRID(6)  = EDAT(EPNTK+7)
         DUM_AGRID(7)  = EDAT(EPNTK+8)
         DUM_AGRID(8)  = EDAT(EPNTK+9)
         EDAT(EPNTK+2) = DUM_AGRID(5)
         EDAT(EPNTK+3) = DUM_AGRID(6)
         EDAT(EPNTK+4) = DUM_AGRID(7)
         EDAT(EPNTK+5) = DUM_AGRID(8)
         EDAT(EPNTK+6) = DUM_AGRID(1)
         EDAT(EPNTK+7) = DUM_AGRID(2)
         EDAT(EPNTK+8) = DUM_AGRID(3)
         EDAT(EPNTK+9) = DUM_AGRID(4)
         WRITE(ERR,9902) TYPE, EID, (EDAT(EPNTK+I+1),I=1,8)
         IF (SUPINFO == 'N') THEN
            WRITE(F06,9902) TYPE, EID, (EDAT(EPNTK+I+1),I=1,8)
         ENDIF

      ELSE IF (TYPE == 'PENTA6  ') THEN

         DUM_AGRID(1)  = EDAT(EPNTK+2)
         DUM_AGRID(2)  = EDAT(EPNTK+3)
         DUM_AGRID(3)  = EDAT(EPNTK+4)
         DUM_AGRID(4)  = EDAT(EPNTK+5)
         DUM_AGRID(5)  = EDAT(EPNTK+6)
         DUM_AGRID(6)  = EDAT(EPNTK+7)
         EDAT(EPNTK+2) = DUM_AGRID(4)
         EDAT(EPNTK+3) = DUM_AGRID(5)
         EDAT(EPNTK+4) = DUM_AGRID(6)
         EDAT(EPNTK+5) = DUM_AGRID(1)
         EDAT(EPNTK+6) = DUM_AGRID(2)
         EDAT(EPNTK+7) = DUM_AGRID(3)

      ELSE IF (TYPE == 'TETRA4  ') THEN

         DUM_AGRID(1)  = EDAT(EPNTK+2)
         DUM_AGRID(2)  = EDAT(EPNTK+3)
         DUM_AGRID(3)  = EDAT(EPNTK+4)
         DUM_AGRID(4)  = EDAT(EPNTK+5)
         EDAT(EPNTK+2) = DUM_AGRID(1)
         EDAT(EPNTK+3) = DUM_AGRID(3)
         EDAT(EPNTK+4) = DUM_AGRID(2)
         EDAT(EPNTK+5) = DUM_AGRID(4)

      ENDIF

! **********************************************************************************************************************************
 9902 FORMAT(' *INFORMATION: ',A,I8,' GRID ORDER HAS BEEN CHANGED TO: ',20I8)

! **********************************************************************************************************************************

      END SUBROUTINE EDAT_FIXUP

! ##################################################################################################################################

      SUBROUTINE EMG_QUIT

      IMPLICIT NONE

! **********************************************************************************************************************************
      WRITE(ERR,9999) NUM_EMG_FATAL_ERRS
      WRITE(F06,9999) NUM_EMG_FATAL_ERRS
      CALL OUTA_HERE ( 'Y' )

! **********************************************************************************************************************************
 9999 FORMAT(' Stopping processing due to ', I8, ' errors found in subr EMG')

! **********************************************************************************************************************************

      END SUBROUTINE EMG_QUIT

      END SUBROUTINE EMG

   END MODULE EMG_MOD
