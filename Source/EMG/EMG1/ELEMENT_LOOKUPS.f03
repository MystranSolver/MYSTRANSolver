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

   MODULE ELEMENT_LOOKUPS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: GET_ELEM_AGRID_BGRID, GET_ELEM_ONAME, GET_ELGP, GRID_ELEM_CONN_TABLE, GET_MATANGLE_FROM_CID

   CONTAINS

      SUBROUTINE GET_ELEM_AGRID_BGRID ( INT_ELEM_ID, CHECK_AGRID )

! Gets element actual and internal grid numbers given the element's internal ID

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, medat0_cuserin, MELGP, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  AGRID, BGRID, EDAT, EID, ELGP, EPNT, ETYPE, GRID, GRID_ID, TYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DOF_ARRAY_INDEXING, ONLY    :  GET_ARRAY_ROW_NUM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_ELEM_AGRID_BGRID'
      CHARACTER(LEN=*), INTENT(IN)    :: CHECK_AGRID       ! If 'Y' perform check on AGRID's to see if appropriate type

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID for which
      INTEGER(LONG)                   :: IERR        = 0   ! Local error count for BGRID not defined
      INTEGER(LONG)                   :: EPNTK             ! Value from array EPNT at the row for this internal elem ID. It is the
!                                                            row number in array EDAT where data begins for this element.
      INTEGER(LONG)                   :: GRID_ID_ROW_NUM   ! Row num in GRID_ID where AGRID(I) exists
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: DELTA             ! Offset in EDAT (from 1st record for an elem) where grid no's begin




! **********************************************************************************************************************************
      EPNTK = EPNT(INT_ELEM_ID)
      TYPE  = ETYPE(INT_ELEM_ID)
      EID   = EDAT(EPNTK)

! AGRID/BGRID contain the G.P. no's (actual/internal) for points that the elem connects to (not for the v vector)

      DO I=1,MELGP+1
         AGRID(I) = 0
         BGRID(I) = 0
      ENDDO

      CALL GET_ELGP ( INT_ELEM_ID )

      DO I=1,ELGP
         DELTA = 1
         IF (TYPE == 'BUSH    ') THEN                      ! 1st grid in EDAT for BUSH is at EPNTK+3 since "Num grids" is EPNTK+2
            DELTA = 2
         ENDIF
         IF (TYPE == 'USERIN  ') THEN
            DELTA = MEDAT0_CUSERIN - 1
         ENDIF
         AGRID(I) = EDAT(EPNTK+I+DELTA)
         CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(I), BGRID(I) )
         IF (BGRID(I) == -1) THEN
            WRITE(ERR,1900) AGRID(I), EID, TYPE
            WRITE(F06,1900) AGRID(I), EID, TYPE
            IERR = IERR + 1
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF
      ENDDO

! Test to determine if AGRID's are appropriate for the elem TYPE (do test only if grid exists)

      IF (CHECK_AGRID == 'Y') THEN
         IF ((TYPE(1:4) /= 'ELAS') .AND. (TYPE /= 'USERIN  ')) THEN
            DO I=1,ELGP
               CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(I), GRID_ID_ROW_NUM )
               IF (GRID_ID_ROW_NUM > 0) THEN
                  IF (GRID(GRID_ID_ROW_NUM,6) /= 6) THEN
                     IERR = IERR + 1
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1951) TYPE, EID, AGRID(I)
                     WRITE(F06,1951) TYPE, EID, AGRID(I)
                  ENDIF
               ENDIF
            ENDDO
         ENDIF
      ENDIF

      IF (IERR > 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1900 FORMAT(' *ERROR  1900: GRID ',I8,' ON ELEMENT ',I8,' TYPE ',A,' NOT  DEFINED')

 1951 FORMAT(' *ERROR  1951: ',A,I8,' USES GRID ',I8,' WHICH IS A SCALAR POINT. SCALAR POINTS NOT ALLOWED FOR THIS ELEM TYPE')

! **********************************************************************************************************************************

      END SUBROUTINE GET_ELEM_AGRID_BGRID


      SUBROUTINE GET_ELEM_ONAME ( NAME )

! Gets element output name (used in LINK9 subr's which write elem and/or ply outputs) for a given element type

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, METYPE
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  ELEM_ONAME, ELMTYP, TYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM))            :: SUBR_NAME = 'GET_ELEM_ONAME'
      CHARACTER( 1*BYTE)                          :: FOUND       = 'N' ! 'Y' if we find the requested element tyoe
      CHARACTER(LEN=LEN(ELEM_ONAME)), INTENT(OUT) :: NAME              ! Name of an elem for output purposes in LINK9 WRTELi subr's

      INTEGER(LONG)                               :: I                 ! DO loop index




! **********************************************************************************************************************************
      NAME = ' '

      FOUND = 'N'
      DO I=1,METYPE
         IF (TYPE == ELMTYP(I)) THEN
            NAME  = ELEM_ONAME(I)
            FOUND = 'Y'
         ENDIF
      ENDDO

! Make sure we found a valid element type

      IF (FOUND == 'N') THEN
         WRITE(ERR,1940) SUBR_NAME, TYPE
         WRITE(F06,1940) SUBR_NAME, TYPE
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
 1940 FORMAT(' *ERROR  1940: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT TYPE "',A,'" NOT FOUND IN ARRAY ELMTYP')



      RETURN

      END SUBROUTINE GET_ELEM_ONAME


      SUBROUTINE GET_ELGP ( INT_ELEM_ID )

! Gets number of grid points for a given element based on the element's internal ID

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE IOUNT1, ONLY                :  ERR, F06, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MEFE, MELGP, METYPE
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, EID, ELGP, ELMTYP, etype, EMG_IFE, EPNT, ERR_SUB_NAM, NELGP, NUM_EMG_FATAL_ERRS, TYPE

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_ELGP'
      CHARACTER( 1*BYTE)              :: FOUND       = 'N' ! 'Y' if we find the requested element tyoe

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID       ! Internal element ID
      INTEGER(LONG)                   :: EPNTK             ! Value from array EPNT at the row for this internal elem ID. It is the
!                                                            row number in array EDAT where data begins for this element.
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: NG                ! Number of GRID's for USERIN elem
      INTEGER(LONG)                   :: NS                ! Number of SPOINT's for USERIN elem




! **********************************************************************************************************************************
      EPNTK = EPNT(INT_ELEM_ID)
      FOUND = 'N'
      TYPE  = ETYPE(INT_ELEM_ID)
      EID   = EDAT(EPNTK)

      IF (TYPE == 'USERIN  ') THEN

         NG = EDAT(EPNTK+2)
         NS = EDAT(EPNTK+3)
         ELGP = NG + NS
         FOUND = 'Y'

      ELSE IF (TYPE == 'BUSH    ') THEN

         ELGP  = EDAT(EPNTK+2)
         FOUND = 'Y'

      ELSE

         DO I=1,METYPE
            IF (TYPE == ELMTYP(I)) THEN
               ELGP  = NELGP(I)
               FOUND = 'Y'
            ENDIF
         ENDDO
         IF ((TYPE(1:4) == 'ELAS') .AND. (EDAT(EPNTK+3) <= 0)) THEN
            ELGP = 1                                       ! A grounded spring: point B is 0 or -1 (subr ELAS_GROUND_END)
         ENDIF

      ENDIF

! If we didn't find a valid element type write error and quit

      IF (FOUND == 'N') THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1940) SUBR_NAME, TYPE
            WRITE(F06,1940) SUBR_NAME, TYPE
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1940
            ENDIF
         ENDIF
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Check ELGP against max allowable (for coding error)

      IF ((ELGP < 1) .OR. (ELGP > MELGP)) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1933) SUBR_NAME, TYPE, EID, ELGP, MELGP
            WRITE(F06,1933) SUBR_NAME, TYPE, EID, ELGP, MELGP
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1933
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = ELGP
               EMG_IFE(NUM_EMG_FATAL_ERRS,3) = MELGP
            ENDIF
         ENDIF
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! **********************************************************************************************************************************
 1933 FORMAT(' *ERROR  1933: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE NUMBER OF GRID POINTS FOR ',A,' ELEMENT ',I8,', IS: ',I8                                          &
                    ,/,14X,' THE MINIMUM NUMBER IS 1 AND THE MAXIMUM NUMBER ALLOWED FOR ANY ELEMENT IS: ',I8)

 1940 FORMAT(' *ERROR  1940: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT TYPE "',A,'" NOT FOUND IN ARRAY ELMTYP')



      RETURN

      END SUBROUTINE GET_ELGP


      SUBROUTINE GRID_ELEM_CONN_TABLE

! Calculates an array that has as many rows as there are grids and as many cols as the max number of elements connected to any grid.
! This array is used for 2 purposes:

!   1) If user has requested any grid point force balance output requests, or
!   2) the user has a PARAM PRTCONN Bulk data entry

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MAX_ELEM_DEGREE, NELE, NGRID
      USE IOUNT1, ONLY                :  F06
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  AGRID, ELGP, ETYPE, ESORT1, ESORT2, GRID_ID, GRID_ELEM_CONN_ARRAY
      USE PARAMS, ONLY                :  PRTCONN

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE DOF_ARRAY_INDEXING, ONLY    :  ARRAY_SIZE_ERROR_1, GET_ARRAY_ROW_NUM
      USE MODEL_STORAGE_ALLOCATION, ONLY:  ALLOCATE_MODEL_STUF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GRID_ELEM_CONN_TABLE'

      INTEGER(LONG)                   :: GRD_NUM_ELEM(NGRID)! Array that specifies the number of elements connected to each grid
      INTEGER(LONG)                   :: I,J,K              ! DO loop indices
      INTEGER(LONG)                   :: IGRID              ! Internal grid ID (row in array GRID_ID where an act grid num exists)




! **********************************************************************************************************************************
      DO I=1,NGRID
         GRD_NUM_ELEM(I) = 0
      ENDDO

      DO J=1,NELE
         CALL GET_ELEM_AGRID_BGRID ( J, 'N' )
         DO K=1,ELGP
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(K), IGRID )
            GRD_NUM_ELEM(IGRID) = GRD_NUM_ELEM(IGRID) + 1
         ENDDO
      ENDDO

      MAX_ELEM_DEGREE = 0
      DO I=1,NGRID
         IF (GRD_NUM_ELEM(I) > MAX_ELEM_DEGREE) THEN
            MAX_ELEM_DEGREE = GRD_NUM_ELEM(I)
         ENDIF
      ENDDO

      CALL ALLOCATE_MODEL_STUF ( 'GRID_ELEM_CONN_ARRAY', SUBR_NAME )

      DO I=1,NGRID
         GRID_ELEM_CONN_ARRAY(I,1) = GRID_ID(I)
         GRID_ELEM_CONN_ARRAY(I,2) = GRD_NUM_ELEM(I)
      ENDDO

      DO I=1,NGRID
         GRD_NUM_ELEM(I) = 0
      ENDDO

      DO J=1,NELE
         CALL GET_ELEM_AGRID_BGRID ( J, 'N' )
         DO K=1,ELGP
            CALL GET_ARRAY_ROW_NUM ( 'GRID_ID', SUBR_NAME, NGRID, GRID_ID, AGRID(K), IGRID )
            GRD_NUM_ELEM(IGRID) = GRD_NUM_ELEM(IGRID) + 1
            IF (GRD_NUM_ELEM(IGRID) > MAX_ELEM_DEGREE) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, MAX_ELEM_DEGREE, 'GRID_ELEM_CONN_ARRAY')
            GRID_ELEM_CONN_ARRAY(IGRID,GRD_NUM_ELEM(IGRID)+2) = ESORT1(J)
         ENDDO
      ENDDO

      IF (PRTCONN > 0) THEN
         WRITE(F06,*)
         WRITE(F06,9100)
         DO I=1,NGRID
            WRITE(F06,9101) (GRID_ELEM_CONN_ARRAY(I,J),J=1,GRD_NUM_ELEM(I)+2)
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
 9100 FORMAT(1X,'                     Table of elements connected to each grid ',//,                                               &
             1X,'         Grid    Num elems     ID''s of elements connected to this grid -->',/)

 9101 FORMAT(1X,2I13,32767I9)

! **********************************************************************************************************************************

      END SUBROUTINE GRID_ELEM_CONN_TABLE




      SUBROUTINE GET_MATANGLE_FROM_CID ( ACID )

! Calcs THETAM for plate elements that have the material angle specified via a coord sys ID (ACID here)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NCORD
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  CONV_DEG_RAD, ZERO, ONE
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  CORD, EID, NUM_EMG_FATAL_ERRS, NUM_EMG_FATAL_ERRS, RCORD, TE, THETAM, TYPE, QUAD_DELTA

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE VECTOR_GEOMETRY, ONLY       :  PROJ_VEC_ONTO_PLANE

      USE VECTOR_GEOMETRY, ONLY       :  CROSS
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_MATANGLE_FROM_CID'
      CHARACTER( 1*BYTE)              :: CORD_FND          ! If 'Y', ACID internal coord sys ID was found in array CORD

      INTEGER(LONG), INTENT(IN)       :: ACID              ! Actual coord system ID for the sys that defines the material axes
      INTEGER(LONG)                   :: I                 ! DO loop indices
      INTEGER(LONG)                   :: ICID              ! Internal coord sys ID for ACID


      REAL(DOUBLE)                    :: DOT_XM            ! Dot product of VEC_XE and VEC_ME
      REAL(DOUBLE)                    :: CROSS_XM(3)       ! Cross product of VEC_XE and VEC_ME
      REAL(DOUBLE)                    :: Y                 ! First parameter for ATAN2 function
      REAL(DOUBLE)                    :: EPS1              ! A small number to comapre to zero
      REAL(DOUBLE)                    :: MAG2_ME           ! Magnitude squared of VEC_ME
      REAL(DOUBLE)                    :: MAG_ME            ! Magnitude of VEC_ME
      REAL(DOUBLE)                    :: VEC_XE(3)         ! Vector in x direction in element coord sys
      REAL(DOUBLE)                    :: VEC_XM(3)         ! Vector in x direction in material angle coord sys
      REAL(DOUBLE)                    :: VEC_ZE(3)         ! Vector in z direction in element coord sys
      REAL(DOUBLE)                    :: VEC_ME(3)         ! Vector proj of VEC_XM onto elem plane



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      DO I=1,3
         VEC_XM(I) = ZERO
      ENDDO

      CORD_FND = 'N'
      ICID     = -1
      IF (ACID /= 0) THEN
i_do1:   DO I=1,NCORD
            IF (ACID == CORD(I,2)) THEN
               CORD_FND = 'Y'
               ICID = I
               EXIT i_do1
            ENDIF
         ENDDO i_do1
         IF (CORD_FND == 'N') THEN
            FATAL_ERR = FATAL_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            WRITE(ERR,1822) 'COORD SYSTEM ', ACID, TYPE, EID
            WRITE(F06,1822) 'COORD SYSTEM ', ACID, TYPE, EID
            RETURN
         ENDIF
         DO I=1,3
            VEC_XM(I) = RCORD(ICID,3*I+1)
         ENDDO
      ELSE
         VEC_XM(1) = ONE
      ENDIF

      DO I=1,3
         VEC_XE(I) = TE(1,I)
         VEC_ZE(I) = TE(3,I)
      ENDDO

      CALL PROJ_VEC_ONTO_PLANE ( VEC_XM, VEC_ZE, VEC_ME )

                                                           ! Check for the MCID x direction being normal to the element.
      MAG2_ME = ZERO
      DO I=1,3
         MAG2_ME = MAG2_ME + VEC_ME(I)*VEC_ME(I)
      ENDDO
      MAG_ME = DSQRT(MAG2_ME)

      IF (MAG_ME <= EPS1) THEN
         FATAL_ERR = FATAL_ERR + 1
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         WRITE(ERR,1829) TYPE, EID, ACID, MAG_ME
         WRITE(F06,1829) TYPE, EID, ACID, MAG_ME
         THETAM = ZERO
      ENDIF

                                                           ! Find the angle of ME from XE about ZE.
                                                           ! THETAM = atan2(ZE dot (XE cross ME) , XE dot ME )

      CALL CROSS(VEC_XE, VEC_ME, CROSS_XM)                 ! XE cross ME
      DOT_XM = VEC_XE(1)*VEC_ME(1) + VEC_XE(2)*VEC_ME(2) + VEC_XE(3)*VEC_ME(3)  ! XE dot ME
      Y = VEC_ZE(1)*CROSS_XM(1) + VEC_ZE(2)*CROSS_XM(2) + VEC_ZE(3)*CROSS_XM(3) ! ZE dot (CROSS_XM)
      THETAM = ATAN2(Y, DOT_XM)

      IF (TYPE(1:5) == 'QUAD4') THEN
        THETAM = THETAM + QUAD_DELTA                       ! Convert angle from element x axis to angle from 1-2 edge.
      ENDIF




      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1829 FORMAT(' *ERROR  1829: CANNOT FIND MATERIAL ANGLE FOR ',A,I8,' USING COORD SYSTEM ',I8,'. '                                  &
                    ,/,14X,' THE PROJECTION (VEC_ME) OF THE VECTOR IN THE X DIR OF THE ABOVE COORD SYSTEM IS NULL:'                &
                    ,/,14X,' MAG VEC_ME = ',1ES9.2)

! **********************************************************************************************************************************

      END SUBROUTINE GET_MATANGLE_FROM_CID

   END MODULE ELEMENT_LOOKUPS
