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

      SUBROUTINE ELMDAT1 ( INT_ELEM_ID, WRITE_WARN )

! Generates small arrays of elem data, for use by subroutine EMG, one elem at a time for all elems. Arrays generated are:

!   XEB      : basic coords of grids for 1 elem
!   V VECTOR : for some 1-D elems
!   EPROP    : array of elem geometric properties
!   ISOLID   : data for 3-D elems defining options from the PSOLID Bulk Data entry
!   EMAT     : material property data
!   PIN FLAG : Pin flag data for some elems
!   OFFSETS  : offsets for some elems

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  FATAL_ERR, MEDAT0_CUSERIN, MELGP, MEMATC, MEMATR, MEPROP, METYPE, MOFFSET, MRMATLC,       &
                                         MRPBAR, MRPBEAM, MRPBUSH, MRPELAS, MRPROD, MRPSHEAR, MRPUSER1, MPSOLID, BLNK_SUB_NAM,     &
                                         NCORD, NGRID, SOL_NAME
      USE SCONTR, ONLY                :  DEDAT_Q4_MATANG_KEY, DEDAT_Q4_THICK_KEY, DEDAT_Q4_POFFS_KEY,                              &
                                         DEDAT_T3_MATANG_KEY, DEDAT_T3_THICK_KEY, DEDAT_T3_POFFS_KEY,                              &
                                                              DEDAT_Q8_THICK_KEY, DEDAT_Q8_POFFS_KEY
      USE PARAMS, ONLY                :  EPSIL, TSTM_DEF
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONEPM4, ONE, TWO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  AGRID, BAROFF, BUSH_CID, BUSH_OCID, BUSH_VVEC, BUSH_VVEC_OR_CID, BUSHOFF, BGRID,          &
                                         CAN_ELEM_TYPE_OFFSET, CORD, DOFPIN, EDAT, EID, ELAS_COMP, ELDOF, ELEM_LEN_12, ELGP,       &
                                         ELMTYP, EMAT, EOFF, NUM_EMG_FATAL_ERRS, EPROP, EPNT, ETYPE, GRID, RGRID, GRID_ID,         &
                                         INTL_MID, INTL_PID, ISOLID, MATANGLE, MATL, MTRL_TYPE, NUM_SEi, OFFDIS, OFFDIS_O, OFFSET, &
                                         PBAR, PBEAM, PCOMP, PCOMP_PROPS, PLATEOFF, PLATETHICK, PROD, PSHEAR, PSHEL, PSOLID,       &
                                         PUSER1, PUSERIN, RMATL, RPBAR, RPBEAM, RPBUSH, RPELAS, RPROD, RPSHEAR, RPSHEL, RPUSER1,   &
                                         TYPE, VVEC, XEB, ZOFFS
      USE MODEL_STUF, ONLY            :  USERIN_ACT_GRIDS, USERIN_ACT_COMPS, USERIN_CID0, USERIN_IN4_INDEX,                        &
                                         USERIN_MAT_NAMES, USERIN_NUM_BDY_DOF, USERIN_NUM_ACT_GRDS, USERIN_NUM_SPOINTS,            &
                                         USERIN_MASS_MAT_NAME, USERIN_LOAD_MAT_NAME, USERIN_RBM0_MAT_NAME, USERIN_STIF_MAT_NAME
! --- cbeam_add begin --- !
      USE MODEL_STUF, ONLY            :  CBEAM_ACTIVE_NSTATIONS, CBEAM_ACTIVE_XL, CBEAM_ACTIVE_RPROPS, PBEAM_NSTATIONS, PBEAM_XL,  &
                                         PBEAM_RPROPS
      USE SCONTR, ONLY                :  MPBEAM_STATIONS
! --- cbeam_add end --- !
      USE ELMDAT1_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMDAT1'
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN         ! If 'Y' write warning messages, otherwise do not
      CHARACTER( 1*BYTE)              :: GET_VVEC           ! If 'Y' run code to get VVEC
      CHARACTER( 1*BYTE)              :: VVEC_DEFINED       ! If 'Y' then a VVEC was found for elements that require one

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID        ! Internal element ID
      INTEGER(LONG)                   :: ACID               ! Actual coord sys ID

      INTEGER(LONG)                   :: EPNTK              ! Value from array EPNT at the row for this internal elem ID. It is the
!                                                             row number in array EDAT where data begins for this element.

      INTEGER(LONG)                   :: IPNTR              ! Pointer into an array
! --- cbeam_add begin --- !
      INTEGER(LONG)                   :: ISTA               ! Loop index for active CBEAM station metadata
! --- cbeam_add end --- !
      INTEGER(LONG)                   :: VVEC_FLAG          ! Either actual grid ID for V vector or -IVVEC

      INTEGER(LONG)                   :: I,J                ! DO loop indices
      INTEGER(LONG)                   :: IERROR             ! Local error count
      INTEGER(LONG)                   :: ICID               ! Internal coord sys ID
      INTEGER(LONG)                   :: IGID_V             ! Internal grid ID for V vector
      INTEGER(LONG)                   :: USERIN_INDEX_COMP  ! Index in EDAT where displ comp list begins
      INTEGER(LONG)                   :: USERIN_INDEX_GRID  ! Index in EDAT where grid list begins
      INTEGER(LONG)                   :: IPIN(2)            ! Pinflag fields from EDAT for a BAR or BEAM elem
      INTEGER(LONG)                   :: IREM               ! Indicator that a pinflag has digits different than 1,2,3,4,5 or 6
      INTEGER(LONG)                   :: IROW               ! Row no. in a real array where data is found for this element
      INTEGER(LONG)                   :: IVVEC              ! Row number in VVEC where V vector for this elem is found
      INTEGER(LONG)                   :: DELTA              ! Delta EDAT row count to get to plate thickness key for QUAD, TRIA
      INTEGER(LONG)                   :: NFLAG              ! Row number in array DOFPIN
      INTEGER(LONG)                   :: NUM_COMPS          ! No. displ components (1 for SPOINT, 6 for actual grid)
      INTEGER(LONG)                   :: NUMMAT             ! No. matl properties for an element type


      REAL(DOUBLE)                    :: DXI                ! An offset distance in direction 1
      REAL(DOUBLE)                    :: DYI                ! An offset distance in direction 2
      REAL(DOUBLE)                    :: DZI                ! An offset distance in direction 3
      REAL(DOUBLE)                    :: EPS1               ! A small number to compare for real zero
      REAL(DOUBLE)                    :: K_BUSH(6),B_BUSH(6)! Values read from PBUSH for stiffness, damping
      REAL(DOUBLE)                    :: PHID, THETAD       ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: THICK_AVG          ! Average of all THICK(i)
      REAL(DOUBLE)                    :: T0G(3,3)           ! Matrix to transform V vector from global to basic coords
      REAL(DOUBLE)                    :: VV(3)              ! V vector in basic coords for this elem

      INTRINSIC                       :: MOD, FLOOR



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      EPNTK     = EPNT(INT_ELEM_ID)
      TYPE      = ETYPE(INT_ELEM_ID)
      EID       = EDAT(EPNTK)
      INTL_PID  = EDAT(EPNTK+1)

! --- cbeam_add begin --- !
      CBEAM_ACTIVE_NSTATIONS = 0
      DO ISTA=1,MPBEAM_STATIONS
         CBEAM_ACTIVE_XL(ISTA) = ZERO
         CBEAM_ACTIVE_RPROPS(ISTA,1) = ZERO
         CBEAM_ACTIVE_RPROPS(ISTA,2) = ZERO
         CBEAM_ACTIVE_RPROPS(ISTA,3) = ZERO
         CBEAM_ACTIVE_RPROPS(ISTA,4) = ZERO
         CBEAM_ACTIVE_RPROPS(ISTA,5) = ZERO
         CBEAM_ACTIVE_RPROPS(ISTA,6) = ZERO
      ENDDO
! --- cbeam_add end --- !

! ELGP is the number of G.P.'s for this elem. Call GET_ELGP to find out how many grids there are for elem type TYPE

      CALL GET_ELGP ( INT_ELEM_ID )


! Make sure that NUM_SEi < ELGP+1 (for all but USERIN). This required in LINK9 where stresses are output with element corner values

      IF (TYPE(1:6) /= 'USERIN') THEN
         DO J=1,METYPE
            IF (ELMTYP(J) == TYPE) THEN
! --- cbeam_add begin --- !
               IF ((TYPE == 'BEAM    ') .AND. (NUM_SEi(J) > 0)) THEN
                  CONTINUE
               ELSE IF (NUM_SEi(J) > (ELGP + 1)) THEN
! --- cbeam_add end --- !
                  WRITE(ERR,1957) SUBR_NAME, TYPE, NUM_SEi(J), ELGP
                  WRITE(F06,1957) SUBR_NAME, TYPE, NUM_SEi(J), ELGP
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
            ENDIF
         ENDDO
      ENDIF

! Set CAN_ELEM_TYPE_OFFSET and test to make sure we have dimensioned OFFDIS, OFFSET properly

      ! *** NOTE: CHECK CODE FOR 3D ELEMS IF THEY ARE TO HAVE OFFSET. GRID ORDER MAY GET CHANGED IN SUBR EDAT_FIXUP (SEE EMG)
      IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'BUSH    ') .OR. (TYPE(1:5) == 'TRIA3'   ) .OR.             &
          ((TYPE(1:5) == 'QUAD4'   ) ) .OR. (TYPE(1:5) == 'QUAD8'   )) THEN
         CAN_ELEM_TYPE_OFFSET = 'Y'
      ELSE
         CAN_ELEM_TYPE_OFFSET = 'N'
      ENDIF

      IF ((CAN_ELEM_TYPE_OFFSET == 'Y') .AND. (ELGP > MOFFSET)) THEN
         WRITE(ERR,1954) SUBR_NAME, MOFFSET, ELGP, TYPE
         WRITE(F06,1954) SUBR_NAME, MOFFSET, ELGP, TYPE
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

! **********************************************************************************************************************************
! AGRID/BGRID contain the G.P. no's (actual/internal) for points that the elem connects to (not for the v vector)

      CALL GET_ELEM_AGRID_BGRID ( INT_ELEM_ID, 'Y' )

! ELDOF are the number of DOF's for this elem

      ELDOF = 0
      DO I=1,ELGP
         CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
         ELDOF = ELDOF + NUM_COMPS
      ENDDO

! Obtain the coordinates of the G.P.'s for this element


      XEB = ZERO

      DO I=1,ELGP
         XEB(I,1) = RGRID(BGRID(I),1)
         XEB(I,2) = RGRID(BGRID(I),2)
         XEB(I,3) = RGRID(BGRID(I),3)
      ENDDO

! For BUSH, make sure CID (field 9 of CBUSH) was specified if the 2 grids are coincident

      IF (TYPE == 'BUSH    ') THEN
         BUSH_CID = EDAT(EPNTK+6)
         IF (ELGP > 1) THEN
            DXI = XEB(2,1) - XEB(1,1)
            DYI = XEB(2,2) - XEB(1,2)
            DZI = XEB(2,3) - XEB(1,3)
!xx         IF (DXI*DXI + DYI*DYI + DZI*DZI < ONEPM6) THEN
            IF ((DABS(DXI) < ONEPM4) .AND. (DABS(DYI) < ONEPM4) .AND. (DABS(DZI) < ONEPM4)) THEN
               IF (BUSH_CID < 0) THEN
                  FATAL_ERR = FATAL_ERR + 1
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS +1
                  WRITE(ERR,1960) TYPE, EID, BUSH_CID
                  WRITE(F06,1960) TYPE, EID, BUSH_CID
               ENDIF
            ENDIF
         ENDIF
      ENDIF

! Obtain the V vector if this is BAR, BEAM, BUSH or USER1 elem. The V vector specification is in EDAT just after the elem G.P.'s.
! If this is a positive number, then it is the G.P. no. from the connection card giving a point on the V vector.
! If it is a neg number it is in VVEC at the row equal to the positive value of the neg number in EDAT. However we
! will have to transform V to basic coords if AGRID(I) global system is not basic. The V vector in basic coord dirs
! is put into XEB after the grid point locations of the elem grid points. Thus, the last row in XEB is a vector,
! measured from grid a, in the basic coord directions along the v vector.
! NOTE: For BUSH elements with a CID (field 9 on CBUSH) > 0, the element coord system is defined by CID and no V vector is required

      GET_VVEC         = 'Y'
      VVEC_DEFINED     = 'N'
      BUSH_VVEC_OR_CID = 'Y'

      ! Initialize V vector
      VV = ZERO

      IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'BUSH    ') .OR. (TYPE == 'USER1   ')) THEN

         VVEC_FLAG = 0
         IF     ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ')) THEN
            VVEC_FLAG = EDAT(EPNTK+4)
         ELSE IF (TYPE == 'USER1   ') THEN
            VVEC_FLAG = EDAT(EPNTK+6)
         ELSE IF (TYPE == 'BUSH    ') THEN
            VVEC_FLAG = EDAT(EPNTK+5)
            BUSH_VVEC = VVEC_FLAG                          ! BUSH_VVEC is in module MODEL_STUF and is used in subr ELMGM1_BUSH
!xx         ACID = BUSH_CID
            IF (BUSH_CID >= 0) THEN                        ! CID >= 0 for BUSH means VVEC is defined by coord sys in BUSH_CID
               GET_VVEC = 'N'
            ENDIF
         ENDIF
         IF (GET_VVEC == 'Y') THEN

            IF (VVEC_FLAG > 0) THEN                        ! V vector is defined by a grid point

               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, VVEC_FLAG, IGID_V )
               IF (IGID_V == -1) THEN
                  WRITE(ERR,1900) AGRID(I),EID,TYPE
                  WRITE(F06,1900) AGRID(I),EID,TYPE
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR = FATAL_ERR + 1
               ELSE
                  VVEC_DEFINED = 'Y'
                  XEB(ELGP+1,1) = RGRID(IGID_V,1)
                  XEB(ELGP+1,2) = RGRID(IGID_V,2)
                  XEB(ELGP+1,3) = RGRID(IGID_V,3)
               ENDIF

            ELSE IF (VVEC_FLAG < 0) THEN                   ! V vector is an actual vector

               IF (BGRID(1) /= -1) THEN                    ! If grid not defined an error was set in subr GET_ELEM_AGRID_BGRID

                  ACID = GRID(BGRID(1),3)                  ! Global coord sys ID
                  IVVEC = -VVEC_FLAG
                  IF (ACID /= 0) THEN                      ! Need to transform V vector to basic coords
                     DO I=1,NCORD
                        IF (ACID == CORD(I,2)) THEN        ! ACID global coord system exists. It was checked in CORDP_PROC
                           ICID = I
                           EXIT
                        ENDIF
                     ENDDO

                     CALL GEN_T0L ( BGRID(1), ICID, THETAD, PHID, T0G )
                     DO J=1,3
                        VV(J) = T0G(J,1)*VVEC(IVVEC,1) + T0G(J,2)*VVEC(IVVEC,2) + T0G(J,3)*VVEC(IVVEC,3)
                     ENDDO
                     DO J=1,3
                        XEB(ELGP+1,J) = XEB(1,J) + VV(J)   ! Basic coords at end of V vector
                     ENDDO

                  ELSE                                     ! V vector was in basic coords

                     DO J=1,3
                        XEB(ELGP+1,J) = XEB(1,J) + VVEC(IVVEC,J)
                     ENDDO

                  ENDIF

               ENDIF

               VVEC_DEFINED = 'Y'

            ELSE IF (VVEC_FLAG == 0)  THEN

               IF (TYPE == 'BUSH    ') THEN                ! vvec and CID not spec for BUSH so need to test EPROP later
                  BUSH_VVEC_OR_CID = 'N'                   ! (NOTE: already made sure CID was defined if GA, GB coincident)
               ENDIF

            ENDIF

         ENDIF

      ENDIF

      IF ((TYPE == 'BUSH    ') .AND. (DEBUG(110) > 0)) THEN
         write(f06,*) ' ==================================================================='
         write(f06,*) ' DEBUG(110) for subr DEBUG_ELMDAT1_FOR_BUSH needs rewrite 2021-08-28'
         write(f06,*) ' ==================================================================='
         CALL DEBUG_ELMDAT1_FOR_BUSH
      ENDIF

! If BAR, BEAM elem, make sure a V vector was specified.

      IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .AND. (GET_VVEC == 'Y')) THEN
         IF (VVEC_DEFINED /= 'Y') THEN
            WRITE(ERR,1902) TYPE, EID
            WRITE(F06,1902) TYPE, EID
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDIF

! **********************************************************************************************************************************
! For BUSH, get coord systems

      IF (TYPE == 'BUSH    ') THEN
         BUSH_CID  = EDAT(EPNTK+6)
         BUSH_OCID = EDAT(EPNTK+7)
      ENDIF

! **********************************************************************************************************************************
! Generate EPROP elem properties for this element. Shell elems with PCOMP props are handled in subr SHELL_ABD_MATRICES.
! Note solid elements have no real property data

      EPROP = ZERO

      IF      (TYPE == 'BAR     ') THEN
         DO I=1,MRPBAR
            EPROP(I) = RPBAR(INTL_PID,I)
         ENDDO

      ELSE IF (TYPE == 'BEAM    ') THEN
         DO I=1,MRPBEAM
            EPROP(I) = RPBEAM(INTL_PID,I)
         ENDDO
! --- cbeam_add begin --- !
!	 Debug Print what pbeam had
         IF (DEBUG(233) > 0) THEN
            WRITE(F06,'(A)') '*** CBEAM PROPERTY DEBUG *******************************************************'
            WRITE(F06,'(A,I0)') '  Element ID       : ', EID
            WRITE(F06,'(A,I0)') '  Property ID      : ', INTL_PID
            WRITE(F06,'(A,1P,ES14.6)') '  Area            : ', EPROP(1)
            WRITE(F06,'(A,1P,ES14.6)') '  I1              : ', EPROP(2)
            WRITE(F06,'(A,1P,ES14.6)') '  I2              : ', EPROP(3)
            WRITE(F06,'(A,1P,ES14.6)') '  I12             : ', EPROP(4)
            WRITE(F06,'(A,1P,ES14.6)') '  J               : ', EPROP(5)
            WRITE(F06,'(A,1P,ES14.6)') '  NSM             : ', EPROP(6)
            WRITE(F06,'(A,1P,ES14.6)') '  K1              : ', EPROP(30)
            WRITE(F06,'(A,1P,ES14.6)') '  K2              : ', EPROP(31)
         ENDIF
! --- cbeam_stations begin --- !
         CBEAM_ACTIVE_NSTATIONS = PBEAM_NSTATIONS(INTL_PID)
         IF (CBEAM_ACTIVE_NSTATIONS > MPBEAM_STATIONS) CBEAM_ACTIVE_NSTATIONS = MPBEAM_STATIONS
         IF (CBEAM_ACTIVE_NSTATIONS <= 1) THEN
            CBEAM_ACTIVE_NSTATIONS = 1
            CBEAM_ACTIVE_XL(1) = ZERO
            CBEAM_ACTIVE_RPROPS(1,1) = EPROP(1)
            CBEAM_ACTIVE_RPROPS(1,2) = EPROP(2)
            CBEAM_ACTIVE_RPROPS(1,3) = EPROP(3)
            CBEAM_ACTIVE_RPROPS(1,4) = EPROP(4)
            CBEAM_ACTIVE_RPROPS(1,5) = EPROP(5)
            CBEAM_ACTIVE_RPROPS(1,6) = EPROP(6)
         ELSE
            DO ISTA=1,CBEAM_ACTIVE_NSTATIONS
               CBEAM_ACTIVE_XL(ISTA) = PBEAM_XL(INTL_PID,ISTA)
               CBEAM_ACTIVE_RPROPS(ISTA,1) = PBEAM_RPROPS(INTL_PID,ISTA,1)
               CBEAM_ACTIVE_RPROPS(ISTA,2) = PBEAM_RPROPS(INTL_PID,ISTA,2)
               CBEAM_ACTIVE_RPROPS(ISTA,3) = PBEAM_RPROPS(INTL_PID,ISTA,3)
               CBEAM_ACTIVE_RPROPS(ISTA,4) = PBEAM_RPROPS(INTL_PID,ISTA,4)
               CBEAM_ACTIVE_RPROPS(ISTA,5) = PBEAM_RPROPS(INTL_PID,ISTA,5)
               CBEAM_ACTIVE_RPROPS(ISTA,6) = PBEAM_RPROPS(INTL_PID,ISTA,6)
            ENDDO
         ENDIF
         IF (DEBUG(233) > 0) THEN
            WRITE(F06,'(A,I0)') '  Stored stations  : ', PBEAM_NSTATIONS(INTL_PID)
            WRITE(F06,'(A,I0)') '  Active stations  : ', CBEAM_ACTIVE_NSTATIONS
            IF (CBEAM_ACTIVE_NSTATIONS > 0) THEN
               WRITE(F06,'(A,11(1X,ES12.5))') '  Active x/L      :', (CBEAM_ACTIVE_XL(ISTA), ISTA=1,CBEAM_ACTIVE_NSTATIONS)
            ENDIF
            WRITE(F06,'(A)') '***************************************************************************'
         ENDIF
! --- cbeam_stations end --- !
! --- cbeam_add end --- !

      ELSE IF (TYPE == 'BUSH    ') THEN
         DO I=1,MRPBUSH
            EPROP(I) = RPBUSH(INTL_PID,I)
            IF (BUSH_VVEC_OR_CID == 'N') THEN              ! Check K2,3,5,6 and B2,3,5,6 = zero when no vvec or CID on BUSH
               DO J=1,6                                    ! (NOTE: already made sure CID was defined if GA, GB coincident)
                  K_BUSH(J) = EPROP(J)
                  B_BUSH(J) = EPROP(J+6)
               ENDDO
               IF ((ABS(K_BUSH(2)) > 0) .OR. (ABS(K_BUSH(3)) > 0) .OR. (ABS(K_BUSH(5)) > 0) .OR. (ABS(K_BUSH(6)) > 0) .OR.         &
                   (ABS(B_BUSH(2)) > 0) .OR. (ABS(B_BUSH(3)) > 0) .OR. (ABS(B_BUSH(5)) > 0) .OR. (ABS(B_BUSH(6)) > 0)) THEN
                  WRITE(ERR,1902) TYPE, EID
                  WRITE(F06,1902) TYPE, EID
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR = FATAL_ERR + 1
               ENDIF
            ENDIF
         ENDDO

      ELSE IF (TYPE(1:4) == 'ELAS') THEN
         DO I=1,MRPELAS
            EPROP(I) = RPELAS(INTL_PID,I)
         ENDDO

      ELSE IF (TYPE == 'ROD     ') THEN
         DO I=1,MRPROD
            EPROP(I) = RPROD(INTL_PID,I)
         ENDDO

      ELSE IF (TYPE == 'SHEAR   ') THEN
         DO I=1,MRPSHEAR
            EPROP(I) = RPSHEAR(INTL_PID,I)
         ENDDO

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. ((TYPE(1:5) == 'QUAD4') ) .OR. (TYPE(1:5) == 'QUAD8')) THEN

                                                           ! For elems that not composites do EPROP in subr SHELL_ABD_MATRICES)
         IF (PCOMP_PROPS == 'N') THEN                      ! Shell properties are in array PSHELL (except maybe membrane thickness)

            IF(TYPE(1:5) == 'QUAD8') THEN                  ! Features that aren't currently supported by CQUAD8.

               IF(SOL_NAME(1:7) /= 'STATICS') THEN
                  WRITE(ERR,*) ' *ERROR: CQUAD8 IS NOT ALLOWED WITH SOL', SOL_NAME
                  WRITE(F06,*) ' *ERROR: CQUAD8 IS NOT ALLOWED WITH SOL', SOL_NAME
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR          = FATAL_ERR + 1
               ENDIF

                                                              ! MID2 /= MID3
               IF(PSHEL(INTL_PID, 2) /= PSHEL(INTL_PID, 3)) THEN
                  WRITE(ERR,*) ' *ERROR: PSHELL FIELDS MID1 AND MID2 MUST HAVE THE SAME VALUE FOR QUAD8'
                  WRITE(F06,*) ' *ERROR: PSHELL FIELDS MID1 AND MID2 MUST HAVE THE SAME VALUE FOR QUAD8'
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR          = FATAL_ERR + 1
               ENDIF

               IF(PSHEL(INTL_PID, 5) /= 0) THEN               ! MID4
                  WRITE(ERR,*) ' *ERROR: PSHELL FIELD MID4 MUST BE BLANK FOR QUAD8'
                  WRITE(F06,*) ' *ERROR: PSHELL FIELD MID4 MUST BE BLANK FOR QUAD8'
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR          = FATAL_ERR + 1
               ENDIF

               IF(RPSHEL(INTL_PID, 2) /= ONE) THEN
                  WRITE(ERR,*) ' *ERROR: PSHELL FIELD 12I/TM**3 MUST BE BLANK OR 1.0 FOR QUAD8'
                  WRITE(F06,*) ' *ERROR: PSHELL FIELD 12I/TM**3 MUST BE BLANK OR 1.0 FOR QUAD8'
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR          = FATAL_ERR + 1
               ENDIF

               IF(RPSHEL(INTL_PID, 4) /= ZERO) THEN           ! NSM
                  WRITE(ERR,*) ' *ERROR: PSHELL FIELD NSM MUST BE BLANK OR 0.0 FOR QUAD8'
                  WRITE(F06,*) ' *ERROR: PSHELL FIELD NSM MUST BE BLANK OR 0.0 FOR QUAD8'
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR          = FATAL_ERR + 1
               ENDIF
                                                              ! Default Z1 set in BD_PSHEL
               IF (RPSHEL(INTL_PID,5) /= -RPSHEL(INTL_PID,1)/TWO) THEN
                  WRITE(ERR,*) ' *ERROR: PSHELL FIELD Z1 MUST BE BLANK FOR QUAD8'
                  WRITE(F06,*) ' *ERROR: PSHELL FIELD Z1 MUST BE BLANK FOR QUAD8'
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR          = FATAL_ERR + 1
               ENDIF
                                                              ! Default Z2 set in BD_PSHEL
               IF (RPSHEL(INTL_PID,6) /= RPSHEL(INTL_PID,1)/TWO) THEN
                  WRITE(ERR,*) ' *ERROR: PSHELL FIELD Z2 MUST BE BLANK FOR QUAD8'
                  WRITE(F06,*) ' *ERROR: PSHELL FIELD Z2 MUST BE BLANK FOR QUAD8'
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR          = FATAL_ERR + 1
               ENDIF
            ENDIF


            EPROP( 2) = RPSHEL(INTL_PID, 2)                ! 12*IB/(TM^3)

            IF (PSHEL(INTL_PID,6) == 0) THEN               ! Use TS/TM from PSHEL card for EPROP( 3)
               EPROP( 3) = RPSHEL(INTL_PID, 3)
            ELSE                                           ! Use default TS/TM (=5/6 unless PARAM TSTM card reset this)
               EPROP( 3) = TSTM_DEF
            ENDIF

            EPROP( 4) = RPSHEL(INTL_PID, 4)                ! NSM
            EPROP( 5) = RPSHEL(INTL_PID, 5)                ! ZS(1)
            EPROP( 6) = RPSHEL(INTL_PID, 6)                ! ZS(2)

            THICK_AVG = ZERO                               ! DELTA locates where thickness key is in EDAT (rel to EID) for plates
            IF ((TYPE(1:5) == 'QUAD4')) THEN
               DELTA = DEDAT_Q4_THICK_KEY
            ELSE IF (TYPE(1:5) == 'TRIA3') THEN
               DELTA = DEDAT_T3_THICK_KEY
            ELSE IF (TYPE(1:5) == 'QUAD8') THEN
               DELTA = DEDAT_Q8_THICK_KEY
            ENDIF

            IF (EDAT(EPNTK+DELTA) > 0) THEN                ! Membrane thickness was defined as grid thicknesses on connection entry
               IPNTR = EDAT(EPNTK+DELTA)
               DO I=1,ELGP
                  EPROP(6+I) = PLATETHICK(IPNTR+I-1)
                  THICK_AVG  = THICK_AVG + PLATETHICK(IPNTR+I-1)/ELGP
               ENDDO
               IF (DABS(EPROP(5)) == ZERO) THEN            ! Since thick was defined on conn entry, reset Z1,2 if they were zero
                  EPROP(5) = -THICK_AVG/TWO
               ENDIF
               IF (DABS(EPROP(6)) == ZERO) THEN
                  EPROP(6) =  THICK_AVG/TWO
               ENDIF
            ELSE                                           ! Membrane thickness was defined on PSHELL
               THICK_AVG = RPSHEL(INTL_PID,1)
            ENDIF

            EPROP(1) = THICK_AVG
            IF (DABS(THICK_AVG) < DABS(EPS1)) THEN         ! Membrane thickness cannot be zero
               WRITE(ERR,1949) TYPE, EID, THICK_AVG
               WRITE(F06,1949) TYPE, EID, THICK_AVG
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               FATAL_ERR          = FATAL_ERR + 1
            ENDIF

         ELSE                                              ! EPROP will be set in subr SHELL_ABD_MATRICES for PCOMP's

         ENDIF


      ELSE IF (TYPE == 'USER1   ') THEN
         DO I=1,MRPUSER1
            EPROP(I) = RPUSER1(INTL_PID,I)
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! Generate ISOLID array of solid element integer data (matl coord system, integration order, stress location, scheme)

      IF ((TYPE(1:4) == 'HEXA') .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN
         DO I=1,MPSOLID
            ISOLID(I) = PSOLID(INTL_PID,I)
         ENDDO
      ENDIF

! Set default integration order if IN = blank

      IF (ISOLID(4) == -1) THEN

        IF (TYPE == 'HEXA8'  ) THEN
            WRITE(ERR,1964) TYPE
            WRITE(F06,1964) TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        ELSE IF (TYPE == 'HEXA20' ) THEN
            ISOLID(4) = 3
        ELSE IF (TYPE == 'PENTA6' ) THEN
            WRITE(ERR,1964) TYPE
            WRITE(F06,1964) TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        ELSE IF (TYPE == 'PENTA15') THEN
            ISOLID(4) = 3
        ELSE IF (TYPE == 'TETRA4' ) THEN
            ISOLID(4) = 2
        ELSE IF (TYPE == 'TETRA10') THEN
            ISOLID(4) = 3
        ELSE
            WRITE(ERR,1967) SUBR_NAME
            WRITE(F06,1967) SUBR_NAME
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        ENDIF

      ELSE IF (ISOLID(4) == 2) THEN

        IF (TYPE == 'TETRA10') THEN
            WRITE(ERR,1965) "IN","TWO",TYPE
            WRITE(F06,1965) "IN","TWO",TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        END IF

      END IF

! Set reduced/standard integration option if ISOP = blank, 0, or 1

      IF (ISOLID(6) == -3) then ! -3 means ISOP = 1

        IF ((TYPE == 'HEXA8'  ) .OR. (TYPE == 'HEXA20' ) .OR. (TYPE == 'PENTA6' ) .OR. (TYPE == 'PENTA15')) THEN
          ISOLID(6) = 1
        ELSE IF ((TYPE == 'TETRA4' ) .OR. (TYPE == 'TETRA10')) THEN
            WRITE(ERR,1965) "ISOP","1",TYPE
            WRITE(F06,1965) "ISOP","1",TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        ELSE
            WRITE(ERR,1967) SUBR_NAME
            WRITE(F06,1967) SUBR_NAME
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        ENDIF

      ELSE IF (ISOLID(6) == -2) then ! -2 means ISOP = 0

        IF ((TYPE == 'HEXA8'  ) .OR. (TYPE == 'HEXA20' ) .OR. (TYPE == 'PENTA6' ) .OR. (TYPE == 'PENTA15')) THEN
          ISOLID(6) = 0
        ELSE IF ((TYPE == 'TETRA4' ) .OR. (TYPE == 'TETRA10')) THEN
            WRITE(ERR,1965) "ISOP","0",TYPE
            WRITE(F06,1965) "ISOP","0",TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        ELSE
            WRITE(ERR,1967) SUBR_NAME
            WRITE(F06,1967) SUBR_NAME
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        ENDIF

      ELSE IF (ISOLID(6) == -1) then ! -1 means ISOP = blank

        IF ((TYPE == 'HEXA8'  ) .OR. (TYPE == 'HEXA20' ) .OR. (TYPE == 'PENTA6' ) .OR. (TYPE == 'PENTA15')) THEN
          ISOLID(6) = 0
        ELSE IF ((TYPE == 'TETRA4' ) .OR. (TYPE == 'TETRA10')) THEN
          ISOLID(6) = 1
        ELSE
            WRITE(ERR,1967) SUBR_NAME
            WRITE(F06,1967) SUBR_NAME
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        ENDIF

      ELSE IF (ISOLID(6) == 0) then ! 0 means ISOP = REDUCED

        IF ((TYPE == 'TETRA4' ) .OR. (TYPE == 'TETRA10')) THEN
            WRITE(ERR,1965) "ISOP","REDUCED",TYPE
            WRITE(F06,1965) "ISOP","REDUCED",TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
        ENDIF

      END IF


! **********************************************************************************************************************************
! Generate EMAT matl props for this elem.  Shell elems with PCOMP_PROPS are handled in subr SHELL_ABD_MATRICES.
! The internal matl ID is in the PROD, PBAR, PBEAM, PSHEL, etc, data sets. There are no matl props for BUSH

      DO I=1,MEMATR
         DO J=1,MEMATC
            EMAT(I,J) = ZERO
         ENDDO
      ENDDO

      DO I=1,MEMATC
         MTRL_TYPE(I) = 0
      ENDDO

      IF      (TYPE == 'BAR     ') THEN

         INTL_MID(1)  = PBAR(INTL_PID,2)
         MTRL_TYPE(1) = MATL(INTL_MID(1),2)
         IF (MTRL_TYPE(1) /= 1) THEN                       ! Must be MAT1 for BAR
            WRITE(ERR,1915) TYPE, EID
            WRITE(F06,1915) TYPE, EID
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         NUMMAT = 1

      ELSE IF (TYPE == 'BEAM    ') THEN

         INTL_MID(1)  = PBEAM(INTL_PID,2)
         MTRL_TYPE(1) = MATL(INTL_MID(1),2)
         IF (MTRL_TYPE(1) /= 1) THEN                       ! Must be MAT1 for BEAM
            WRITE(ERR,1915) TYPE, EID
            WRITE(F06,1915) TYPE, EID
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         NUMMAT = 1

      ELSE IF ((TYPE(1:4) == 'HEXA')  .OR. (TYPE(1:5) == 'PENTA') .OR. (TYPE(1:5) == 'TETRA')) THEN

         INTL_MID(1) = PSOLID(INTL_PID,2)
         MTRL_TYPE(1) = MATL(INTL_MID(1),2)                ! Must be MAT1 or MAT9 for solids
         IF ((MTRL_TYPE(1) /= 1) .AND. (MTRL_TYPE(1) /= 9)) THEN
            WRITE(ERR,1958) TYPE, EID
            WRITE(F06,1958) TYPE, EID
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         NUMMAT = 1

      ELSE IF (TYPE == 'ROD     ') THEN

         INTL_MID(1) = PROD(INTL_PID,2)
         MTRL_TYPE(1) = MATL(INTL_MID(1),2)
         IF (MTRL_TYPE(1) /= 1) THEN                       ! Must be MAT1 for ROD
            WRITE(ERR,1915) TYPE, EID
            WRITE(F06,1915) TYPE, EID
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         NUMMAT = 1

      ELSE IF (TYPE == 'SHEAR   ') THEN

         INTL_MID(1)  = PSHEAR(INTL_PID,2)
         MTRL_TYPE(1) = MATL(INTL_MID(1),2)
         IF (MTRL_TYPE(1) /= 1) THEN                       ! Must be MAT1 for SHEAR
            WRITE(ERR,1915) TYPE, EID
            WRITE(F06,1915) TYPE, EID
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
         NUMMAT = 1

      ELSE IF ((TYPE(1:5) == 'TRIA3') .OR. ((TYPE(1:5) == 'QUAD4') ) .OR. (TYPE(1:5) == 'QUAD8')) THEN
                                                           ! For elems that are not composites do EMAT in subr SHELL_ABD_MATRICES)
         IF (PCOMP_PROPS == 'N') THEN
            INTL_MID(1) = PSHEL(INTL_PID,2)
            IF (INTL_MID(1) /= 0) THEN
               MTRL_TYPE(1) = MATL(INTL_MID(1),2)
            ENDIF

            INTL_MID(2) = PSHEL(INTL_PID,3)
            IF (INTL_MID(2) /= 0) THEN
               MTRL_TYPE(2) = MATL(INTL_MID(2),2)
            ENDIF

            INTL_MID(3) = PSHEL(INTL_PID,4)
            IF (INTL_MID(3) /= 0) THEN
               MTRL_TYPE(3) = MATL(INTL_MID(3),2)
            ENDIF

            INTL_MID(4) = PSHEL(INTL_PID,5)
            IF (INTL_MID(4) /= 0) THEN
               MTRL_TYPE(4) = MATL(INTL_MID(4),2)
            ENDIF

            NUMMAT = 4

            IF (TYPE(1:5) == 'QUAD8') THEN
               DO I=1,NUMMAT                               ! Must be MAT1 for CQUAD8 elements
                  IF (MTRL_TYPE(I) /= 0) THEN              ! as long as a material was defined
                     IF (MTRL_TYPE(I) /= 1) THEN
                        WRITE(ERR,1921) TYPE, EID
                        WRITE(F06,1921) TYPE, EID
                        NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                        FATAL_ERR = FATAL_ERR + 1
                     ENDIF
                  ENDIF
               ENDDO
                                                           ! Density is not allowed in case it's used for
                                                           ! gravity which doesn't currently work.
               IF(RMATL(INTL_MID(1),4) /= ZERO) THEN
                  WRITE(ERR,*) ' *ERROR: MAT1 DENSITY MUST BE 0.0 OR BLANK FOR QUAD8'
                  WRITE(F06,*) ' *ERROR: MAT1 DENSITY MUST BE 0.0 OR BLANK FOR QUAD8'
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR          = FATAL_ERR + 1
               ENDIF
            ENDIF

            DO I=1,NUMMAT                                  ! Must be MAT1, MAT2 or MAT8 for plate elems
               IF (MTRL_TYPE(I) /= 0) THEN                 ! as long as a material was defined
                  IF ((MTRL_TYPE(I) /= 1) .AND. (MTRL_TYPE(I) /= 2) .AND. (MTRL_TYPE(I) /= 8)) THEN
                     WRITE(ERR,1920) TYPE, EID
                     WRITE(F06,1920) TYPE, EID
                     NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
               ENDIF
            ENDDO

         ENDIF

      ELSE IF (TYPE == 'USER1   ') THEN

         INTL_MID(1)  = PUSER1(INTL_PID,2)
         MTRL_TYPE(1) = MATL(INTL_MID(1),2)

         INTL_MID(2)  = PUSER1(INTL_PID,3)
         MTRL_TYPE(2) = MATL(INTL_MID(2),2)

         INTL_MID(3)  = PUSER1(INTL_PID,4)
         MTRL_TYPE(3) = MATL(INTL_MID(3),2)

         NUMMAT = 3

         DO I=1,NUMMAT                                     ! Must be MAT1 or MAT8 for USER1 elements
            IF ((MTRL_TYPE(I) /= 1) .AND. (MTRL_TYPE(I) /= 8)) THEN
               WRITE(ERR,1920) TYPE, EID
               WRITE(F06,1920) TYPE, EID
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
         ENDDO

      ENDIF

      IF ((TYPE /= 'BUSH    ') .AND. (TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'PLOTEL  ') .AND.                                        &
          (TYPE /= 'USERIN  ') .AND. (PCOMP_PROPS == 'N')) THEN
         DO I=1,MRMATLC
            DO J=1,NUMMAT                                  ! NUMMAT are num of matls for this elem type (e.g. 4 for TRIA's & QUAD's)
               IF (INTL_MID(J) /= 0) THEN
                  EMAT(I,J) = RMATL(INTL_MID(J),I)
               ENDIF
            ENDDO
         ENDDO
      ENDIF

! Set transverse shear alloawbles to same as in-plane shear allowables for non PCOMP shells. The transverse shear allowables go
! in rows 19 and 20 of EMAT

      IF ((TYPE(1:5) == 'TRIA3') .OR. ((TYPE(1:5) == 'QUAD4') )) THEN
         if (PCOMP_PROPS == 'N') THEN
            DO I=1,NUMMAT
               IF      (INTL_MID(I) == 1) THEN
                  EMAT(19,I) = EMAT(10,I)
                  EMAT(20,I) = EMAT(10,I)
               ELSE IF (INTL_MID(I) == 2) THEN
                  EMAT(19,I) = EMAT(15,I)
                  EMAT(20,I) = EMAT(15,I)
               ELSE IF (INTL_MID(I) == 8) THEN
                  EMAT(19,I) = EMAT(15,I)
                  EMAT(20,I) = EMAT(15,I)
               ENDIF
            ENDDO
         ENDIF
      ENDIF

! **********************************************************************************************************************************
! Process pin flag data and create array of elem DOF numbers for the DOF's that are pinned. The input pin flag data
! is integer for each G.P. and is stored in EDAT. Call this IPIN(I) for each G.P. and convert this to a list of (up
! to 6 times the no.of elem G.P.'s) DOF numbers for the pinned DOF's. Call this list, DOFPIN.

      DO I=1,ELDOF
         DOFPIN(I) = 0
      ENDDO

! Process pin flag data in EDAT

      IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ')) THEN

         NFLAG = 0

         IF (ELGP > 2) THEN                                ! IPIN is only dimensioned to 2 so set error

            WRITE(ERR,1953) SUBR_NAME, ELGP,TYPE
            WRITE(F06,1953) SUBR_NAME, ELGP,TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1

         ELSE

            DO I=1,ELGP
               IPIN(I) = EDAT(EPNTK+4+I)
               DO J=1,6
                  IREM    = MOD(IPIN(I),10)                ! e.g., if IPIN(I) = 136 then IREM = 6
                  IPIN(I) = FLOOR(DBLE(IPIN(I)/10))        ! e.g., if IPIN(I) = 136 then IPIN(I) is changed to 13
                  IF ((IREM < 0) .OR. (IREM > 6)) THEN
                     WRITE(ERR,1901) EID,TYPE
                     WRITE(F06,1901) EID,TYPE
                     NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
                  IF (IREM > 0) THEN
                     NFLAG = NFLAG + 1
                     DOFPIN(NFLAG) = IREM + 6*(I-1)
                  ENDIF
               ENDDO
            ENDDO
            IF ((DEBUG(9) > 0) .AND. (NFLAG > 0)) THEN
               WRITE(F06,*) 'In ELMDAT1: DOFPIN for ',type,eid,' is = ',(dofpin(i),i=1,nflag)
               WRITE(F06,*)
            ENDIF

         ENDIF

      ENDIF

! **********************************************************************************************************************************
! Get offset data for this elem. This data is in the BAROFF array for BAR, BEAM elems and in BUSHOFF for BUSH elements
! and in RPSHEL array for shell elems. For BUSH, we need ELEM_LEN_12 since the offsets are fractions of this length

! Distance between grids 1 and 2 is:

      ELEM_LEN_12 = DSQRT( (XEB(2,1) - XEB(1,1))*(XEB(2,1) - XEB(1,1))                                                             &
                         + (XEB(2,2) - XEB(1,2))*(XEB(2,2) - XEB(1,2))                                                             &
                         + (XEB(2,3) - XEB(1,3))*(XEB(2,3) - XEB(1,3)) )
      IF (CAN_ELEM_TYPE_OFFSET == 'Y') THEN

         OFFDIS = ZERO
	
         IF      ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ')) THEN

            IROW = EDAT(EPNTK + 7)
            IF (IROW > 0) THEN
               EOFF(INT_ELEM_ID) = 'Y'
               OFFDIS(1,1) = BAROFF(IROW,1)
               OFFDIS(1,2) = BAROFF(IROW,2)
               OFFDIS(1,3) = BAROFF(IROW,3)
               OFFDIS(2,1) = BAROFF(IROW,4)
               OFFDIS(2,2) = BAROFF(IROW,5)
               OFFDIS(2,3) = BAROFF(IROW,6)
            ELSE
               EOFF(INT_ELEM_ID) = 'N'
            ENDIF

         ELSE IF (TYPE(1:4) == 'BUSH') THEN
            IROW = EDAT(EPNTK + 8)
            IF (IROW > 0) THEN
                  IF(BUSH_OCID == -1) THEN                 ! Offsets are along line GA-GB and use only S
                     IF (ELEM_LEN_12 >= .0001) THEN
                                                           ! Offsets require an axis, so nonzero length
                        EOFF(INT_ELEM_ID) = 'Y'            ! is needed since we can't use GA-GB when OCID=-1
                        OFFDIS(1,1) = ELEM_LEN_12*BUSHOFF(IROW,1)
                        OFFDIS(1,2) = ELEM_LEN_12*BUSHOFF(IROW,2)
                        OFFDIS(1,3) = ELEM_LEN_12*BUSHOFF(IROW,3)

                        OFFDIS(2,1) =-ELEM_LEN_12*(ONE - BUSHOFF(IROW,1))
                        OFFDIS(2,2) = ELEM_LEN_12*BUSHOFF(IROW,2)
                        OFFDIS(2,3) = ELEM_LEN_12*BUSHOFF(IROW,3)
                     ELSE
                        EOFF(INT_ELEM_ID) = 'N'           ! No axis, no OCID, no offset.
                     END IF
                  ELSE
                                               ! Offsets are 3 dimensional and we can't get OFFDIS for grid GB yet
                     EOFF(INT_ELEM_ID) = 'Y'            ! is required since we can't use GA-GB (OCID=-1)

                     OFFDIS(1,1) = BUSHOFF(IROW,1)
                     OFFDIS(1,2) = BUSHOFF(IROW,2)
                     OFFDIS(1,3) = BUSHOFF(IROW,3)
                                                           ! We can't calc OFFDIS for grid GB if OCID >= 0 until subr ELMGM1_BUSH
                     OFFDIS(2,1) = OFFDIS(1,1)             ! So just load grid GA values into OFFDIS for grid GB temporarily
                     OFFDIS(2,2) = OFFDIS(1,2)
                     OFFDIS(2,3) = OFFDIS(1,3)

                  ENDIF


               DO I=1,MOFFSET                              ! Now set OFFDIS_O to be these original OFFDIS values so we can preserve
                  DO J=1,3                                 ! the original values based on input
                     OFFDIS_O(I,J) = OFFDIS(I,J)
                  ENDDO
               ENDDO

            ELSE
               EOFF(INT_ELEM_ID) = 'N'
            ENDIF

         ELSE IF (TYPE(1:5) == 'TRIA3') THEN

            IROW = EDAT(EPNTK + DEDAT_T3_POFFS_KEY)
            IF (IROW > 0) THEN                             ! Elem has offset. IROW > 0 is the row in PLATEOFF where ZOFFS is
               EOFF(INT_ELEM_ID) = 'Y'
               ZOFFS = PLATEOFF(IROW)
               IF (PCOMP_PROPS == 'N') THEN                ! This elem has PSHELL props so OK to have offset
                  OFFDIS(1,3) = ZOFFS
                  OFFDIS(2,3) = ZOFFS
                  OFFDIS(3,3) = ZOFFS
               ELSE                                        ! This elem has PCOMP props so not OK to have offset
                  IF (DABS(ZOFFS) > 0.D0) THEN
                     WRITE(ERR,1945) TYPE, EID, ZOFFS, PCOMP(INTL_PID,1)
                     WRITE(F06,1945) TYPE, EID, ZOFFS, PCOMP(INTL_PID,1)
                     NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
               ENDIF
            ELSE
               EOFF(INT_ELEM_ID) = 'N'
               ZOFFS = ZERO
            ENDIF

         ELSE IF ((TYPE(1:5) == 'QUAD4') ) THEN

            IROW = EDAT(EPNTK + DEDAT_Q4_POFFS_KEY)
            IF (IROW > 0) THEN                             ! Elem has offset. IROW > 0 is the row in PLATEOFF where ZOFFS is
               EOFF(INT_ELEM_ID) = 'Y'
               ZOFFS = PLATEOFF(IROW)
               IF (PCOMP_PROPS == 'N') THEN                ! This elem has PSHELL props so OK to have offset
                  OFFDIS(1,3) = ZOFFS
                  OFFDIS(2,3) = ZOFFS
                  OFFDIS(3,3) = ZOFFS
                  OFFDIS(4,3) = ZOFFS
               ELSE                                        ! This elem has PCOMP props so not OK to have offset
                  IF (DABS(ZOFFS) > 0.D0) THEN
                     WRITE(ERR,1945) TYPE, EID, ZOFFS, PCOMP(INTL_PID,1)
                     WRITE(F06,1945) TYPE, EID, ZOFFS, PCOMP(INTL_PID,1)
                     NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                     FATAL_ERR = FATAL_ERR + 1
                  ENDIF
               ENDIF
            ELSE
               EOFF(INT_ELEM_ID) = 'N'
               ZOFFS = ZERO
            ENDIF

         ELSE IF (TYPE(1:5) == 'QUAD8') THEN

            IROW = EDAT(EPNTK + DEDAT_Q8_POFFS_KEY)
            IF (IROW > 0) THEN                             ! Elem has offset. IROW > 0 is the row in PLATEOFF where ZOFFS is
              EOFF(INT_ELEM_ID) = 'Y'
              ZOFFS = PLATEOFF(IROW)
              IF (DABS(ZOFFS) > 0.D0) THEN                 ! Offset not yet allowed on quad8
                WRITE(ERR,*) ' *ERROR  : OFFSET CANNOT BE USED FOR CQUAD8'
                WRITE(F06,*) ' *ERROR  : OFFSET CANNOT BE USED FOR CQUAD8'
                NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                FATAL_ERR = FATAL_ERR + 1
              ENDIF
            ELSE
               EOFF(INT_ELEM_ID) = 'N'
               ZOFFS = ZERO
            ENDIF

         ENDIF

! Determine which grids have finite offsets for this element

         DO I=1,ELGP
            OFFSET(I) = 'N'
         ENDDO

         DO I=1,ELGP
            DXI = DABS(OFFDIS(I,1))
            DYI = DABS(OFFDIS(I,2))
            DZI = DABS(OFFDIS(I,3))
            IF ((DXI >= EPS1) .OR. (DYI >= EPS1) .OR. (DZI >= EPS1)) THEN
               OFFSET(I) = 'Y'
            ENDIF
         ENDDO

      ENDIF

! **********************************************************************************************************************************
! For ELAS elems form the list of G.P. DOF no's that the elem connects to. This is data from the CELAS card and is stored in EDAT

      IF ((TYPE == 'ELAS1   ') .OR. (TYPE == 'ELAS2   ')) THEN

         DO I=1,2
            ELAS_COMP(I) = 0
         ENDDO

         IERROR = 0
         DO I=1,2                                          ! If displ comps on CELAS1,2 entry were blank or 0, change to 1,2
            ELAS_COMP(I) = EDAT(EPNTK+3+I)                 ! (i.e. ELAS has 2 components of displ)
            IF (ELAS_COMP(I) == 0) THEN
               CALL GET_GRID_NUM_COMPS ( BGRID(I), NUM_COMPS, SUBR_NAME )
               IF (NUM_COMPS > 1) THEN
                  NUM_EMG_FATAL_ERRS  = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR = FATAL_ERR + 1
                  IERROR    = IERROR + 1
                  WRITE(ERR,1950) TYPE, EID
                  WRITE(F06,1950) TYPE, EID
               ELSE
                  ELAS_COMP(I) = 1
               ENDIF
            ENDIF
         ENDDO
         IF (IERROR > 0) THEN
            CALL OUTA_HERE ( 'Y' )
         ENDIF

      ELSE IF ((TYPE == 'ELAS3   ') .OR. (TYPE == 'ELAS4   ')) THEN

         ELAS_COMP(1) = 1                                  ! These are the components at the SPOINT, not the cols of the KE matrix
         ELAS_COMP(2) = 1                                  ! (for cols of the KE matrix, see subr ELAS1, they depend on ELAS_COMP)

      ENDIF

! **********************************************************************************************************************************
! Generate USERIN specific data

      IF (TYPE == 'USERIN') THEN

         USERIN_NUM_ACT_GRDS  = EDAT(EPNTK+2)
         USERIN_NUM_SPOINTS   = EDAT(EPNTK+3)
         USERIN_CID0          = EDAT(EPNTK+4)
         USERIN_NUM_BDY_DOF   = EDAT(EPNTK + MEDAT0_CUSERIN + 2*USERIN_NUM_ACT_GRDS + USERIN_NUM_SPOINTS)

         USERIN_INDEX_GRID    = EPNTK + MEDAT0_CUSERIN
         USERIN_INDEX_COMP    = EPNTK + MEDAT0_CUSERIN + USERIN_NUM_ACT_GRDS + USERIN_NUM_SPOINTS

         USERIN_IN4_INDEX     = PUSERIN(INTL_PID,2)

         USERIN_STIF_MAT_NAME = USERIN_MAT_NAMES(INTL_PID,1)
         USERIN_MASS_MAT_NAME = USERIN_MAT_NAMES(INTL_PID,2)
         USERIN_LOAD_MAT_NAME = USERIN_MAT_NAMES(INTL_PID,3)
         USERIN_RBM0_MAT_NAME = USERIN_MAT_NAMES(INTL_PID,4)
                                                           ! Calc array of grid num, comp num for the DOF's that elem connects to
         CALL DEALLOCATE_MODEL_STUF ( 'USERIN_ACT_GRDS, USERIN_ACT_COMPS' )
         CALL ALLOCATE_MODEL_STUF   ( 'USERIN_ACT_GRDS, USERIN_ACT_COMPS', SUBR_NAME )
         DO I=1,USERIN_NUM_ACT_GRDS
            USERIN_ACT_GRIDS(I) = EDAT(USERIN_INDEX_GRID+I-1)
            USERIN_ACT_COMPS(I) = EDAT(USERIN_INDEX_COMP+I-1)
         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1900 FORMAT(' *ERROR  1900: GRID ',I8,' ON ELEMENT ',I8,' TYPE ',A,' NOT  DEFINED')

 1901 FORMAT(' *ERROR  1901: A PIN FLAG ON ELEMENT ',I8,' TYPE ',A,' IS < 0  OR > 6')

 1902 FORMAT(' *ERROR  1902: NO V VECTOR SPECIFIED FOR ',A,I8)

 1915 FORMAT(' *ERROR  1915: MATERIAL PROPERTIES FOR ',A8,' ID = ',I8,' MUST BE DEFINED ON A MAT1 BULK DATA ENTRY')

 1920 FORMAT(' *ERROR  1920: MATERIAL PROPERTIES FOR ',A8,' ID = ',I8,' MUST BE DEFINED ON A MAT1, MAT2 OR MAT8 BULK DATA ENTRY')

 1921 FORMAT(' *ERROR  1921: MATERIAL PROPERTIES FOR ',A8,' ID = ',I8,' MUST BE DEFINED ON A MAT1 BULK DATA ENTRY')

 1945 FORMAT(' *ERROR  1945: ',A,' ELEMENT ',I8,' HAS OFFSET ZOFFS = ',1ES9.2,' SPECIFIED. HOWEVER, IT HAS PROPERTIES DEFINED ON', &
                             ' PCOMP ',I8,'.'                                                                                      &
                      ,/,14X,' ZOFFS CANNOT BE USED FOR SHELL ELEMENTS WITH PCOMP PROPERTIES. USE THE Z0 OPTION ON THE PCOMP')

 1949 FORMAT(' *ERROR  1949: MEMBRANE THICKNESS FOR ',A,I8,' = ',1ES9.2,'. CHECK TM ON PSHELL OR Ti ON CONNECTION ENTRY TO MAKE'   &
                    ,/,14X,' SURE THAT ELEM AVG THICKNESS IS > 0')

 1950 FORMAT(' *ERROR  1950: ',A,I8,' HAS A DISPLACEMENT COMPONENT = 0 (OR BLANK) ON THE CELAS B.D. ENTRY. FOR A PHYSICAL GRID IT',&
                           ' MUST BE AN INTEGER >= 1 AND <= 6')

 1953 FORMAT(' *ERROR  1953: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' DIMENSION OF ARRAY IPIN IS ONLY 2 BUT MUST BE ',I8,' FOR ELEM TYPE ',A)

 1954 FORMAT(' *ERROR  1954: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' DIMENSION OF ARRAYS OFFDIS, OFFSET ARE ONLY ',I8,' BUT MUST BE ',I8,' FOR ELEM TYPE ',A)

 1957 FORMAT(' *ERROR  1957: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT TYPE "',A,'" HAS NUM_SEi = ',I3,' WHICH MUST BE AT LEAST 1 GREATER THAN THE NUMBER OF GRIDS', &
                           ' = ',I3)

 1958 FORMAT(' *ERROR  1958: MATERIAL PROPERTIES FOR ',A8,' ID = ',I8,' MUST BE DEFINED ON A MAT1 OR MAT9 BULK DATA ENTRY')

 1960 FORMAT(' *ERROR  1960: ',A,I8,' MUST HAVE A CID COORD SYSTEM (FIELD 9 OF CBUSH BULK DATA ENTRY) DEFINED WHEN THE 2 GRIDS ARE'&
                                   ,' COINCIDENT. HOWEVER, CID = ',I8)

 1964 FORMAT(' *ERROR  1964: PSOLID ENTRY HAD FIELD FOR IN BLANK WITH ELEMENT TYPE ',A)

 1965 FORMAT(' *ERROR  1965: PSOLID ENTRY HAD FIELD FOR ',A,' = "',A,'" WITH ELEMENT TYPE ',A)

 1967 FORMAT(' *ERROR  1967: PROGRAMMING ERROR IN SUBROUTINE ',A)








99978 format(' Element grid thicknesses for ',a,i8,' (pointer, avg thickness, grid thicknesses)',/,i8,20(1es14.6))

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEBUG_ELMDAT1_FOR_BUSH

      IMPLICIT NONE

      CHARACTER(17*BYTE)              :: CBUSH_MSG(9)

! **********************************************************************************************************************************
      WRITE(F06,*)
      WRITE(F06,98720)
      WRITE(F06,'(A,A,I8,A,I4,A)') ' In ELMDAT1 for ', type, eid, ' (EPNTK = ',epntk,')'
      WRITE(F06,*) '----------------------------------------------'
      CBUSH_MSG(1) = '  elem ID'
      CBUSH_MSG(2) = '  prop ID'
      CBUSH_MSG(3) = '  grid A coords:'
      CBUSH_MSG(4) = '  grid B coords:'
      CBUSH_MSG(5) = '  VVEC_FLAG'
      CBUSH_MSG(6) = '  CID'
      CBUSH_MSG(7) = '  OCID'
      CBUSH_MSG(8) = '  offset flag'
      CBUSH_MSG(9) = '  Location flag'
      DO I=1,2
         WRITE(F06,'(A,I2,A,I8,A)') '    EDAT(EPNTK+',I-1,') = ',EDAT(EPNTK+I-1), CBUSH_MSG(I)
      ENDDO
      WRITE(F06,'(A,I2,A,I8,A,3(1ES14.6))') '    EDAT(EPNTK+',2,') = ',EDAT(EPNTK+2), CBUSH_MSG(3), (RGRID(BGRID(1),J),J=1,3)
      WRITE(F06,'(A,I2,A,I8,A,3(1ES14.6))') '    EDAT(EPNTK+',3,') = ',EDAT(EPNTK+3), CBUSH_MSG(4), (RGRID(BGRID(2),J),J=1,3)
      DO I=5,8
         WRITE(F06,'(A,I2,A,I8,A)') '    EDAT(EPNTK+',I-1,') = ',EDAT(EPNTK+I-1), CBUSH_MSG(I)
      ENDDO
      WRITE(F06,*)

      WRITE(F06,'(A,2I14,A14)') '    VVEC_FLAG, CID, GET_VVEC                = ', VVEC_FLAG, ACID, GET_VVEC

      IF (VVEC_FLAG < 0) THEN
         WRITE(F06,'(A,I8,A,3(1ES14.6))') '    VVEC in input coord system ',acid,'     = ', (VVEC(IVVEC,J),J=1,3)
         WRITE(F06,'(A,3(1ES14.6))') '    VVEC in basic coords                    = ', (VV(J),J=1,3)
         WRITE(F06,*)
         WRITE(F06,'(A,3(1ES14.6))') '    XEB basic sys coords of origin of VVEC  = ', (XEB(1     ,J),j=1,3)
         WRITE(F06,'(A,3(1ES14.6))') '    XEB basic sys coords of end    of VVEC  = ', (XEB(ELGP+1,J),j=1,3)
         WRITE(F06,*)
      ENDIF

      WRITE(F06,*)
      WRITE(F06,98799)
      WRITE(F06,*)

! **********************************************************************************************************************************
98720 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::::::::START DEBUG(110) OUTPUT FROM SUBROUTINE ELMDAT1:::::::::::::::::::::::::',&
             ':::::::::::::::::',/)

98799 FORMAT(' :::::::::::::::::::::::::::::::::::::::::::END DEBUG(110) OUTPUT FROM SUBROUTINE ELMDAT1::::::::::::::::::::::::::',&
             ':::::::::::::::::'                                                                                                ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_ELMDAT1_FOR_BUSH

      END SUBROUTINE ELMDAT1
