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

   MODULE ELEMENT_GEOMETRY_PREPARATION

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ELMGM1, ELMGM1_BUSH, ELMGM2, ELMGM3

   CONTAINS

      SUBROUTINE ELMGM1 ( INT_ELEM_ID, WRITE_WARN )

! Calculates and checks some elem geometry for ROD, BAR, BEAM, triangles and provides a transformation matrix (TE) to transform the
! element stiffness matrix in the element system to the basic coordinate system. Calculates grid point coords in local coord system.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MELGP, MOFFSET, NCORD, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  BGRID, CAN_ELEM_TYPE_OFFSET, CORD, EID, ELEM_LEN_12, ELEM_LEN_AB, ELGP,NUM_EMG_FATAL_ERRS,&
                                         EOFF, GRID, OFFDIS, OFFDIS_O, OFFDIS_B, OFFDIS_G, RCORD, TE, TE_IDENT, TYPE, XEB, XEL

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE VECTOR_GEOMETRY, ONLY       :  CROSS, GEN_T0L
      USE SORTING, ONLY               :  CALC_VEC_SORT_ORDER

      IMPLICIT NONE

      CHARACTER( 1*BYTE)              :: ID(3)              ! Used in deciding whether TE_IDENT = 'Y' or 'N'
      CHARACTER( 5*BYTE)              :: SORT_ORDER         ! Order in which the VX(i) have been sorted in subr CALC_VEC_SORT_ORDER
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMGM1'
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN         ! If 'Y" write warning messages, otherwise do not
      CHARACTER(1*BYTE)               :: CORD_FND           ! = 'Y' if coord sys ID on CONM2 defined, 'N' otherwise
      CHARACTER(1*BYTE)               :: DO_IT              ! = 'Y' execute code that follows (see where DO_IT is initialized)

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID        ! Internal element ID for which
      INTEGER(LONG)                   :: ACID_G             ! Actual coordinate system ID
      INTEGER(LONG)                   :: I,J,K              ! DO loop indices
      INTEGER(LONG)                   :: I3_IN(3)           ! Integer array used in sorting VX.

      INTEGER(LONG)                   :: I3_OUT(3)          ! Integer array giving order of VX comps. If VX is in the order with
!                                                             comp 2 smallest then comp 3 then comp 1 then I3_OUT is 2, 3, 1

      INTEGER(LONG)                   :: ICID               ! Internal coord sys no. corresponding to an actual coord sys no.
      INTEGER(LONG)                   :: ROWNUM             ! A row number in an array


      REAL(DOUBLE)                    :: DX1(3)             ! Array used in intermediate calc's
      REAL(DOUBLE)                    :: DX2(3)             ! Array used in intermediate calc's
      REAL(DOUBLE)                    :: EPS1               ! A small number to compare to real zero
      REAL(DOUBLE)                    :: LX(3)              ! Distances
      REAL(DOUBLE)                    :: MAGY               ! Magnitude of vector VY
      REAL(DOUBLE)                    :: MAGZ               ! Magnitude of vector VZ
      REAL(DOUBLE)                    :: PHID, THETAD       ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: T0G(3,3)           ! Matrix to transform offsets from global to basic  coords
      REAL(DOUBLE)                    :: TG0(3,3)           ! Matrix to transform offsets from basic  to global coords
      REAL(DOUBLE)                    :: TET(3,3)           ! Transpose of TE: UEL = TE*UEB
      REAL(DOUBLE)                    :: VX(3)              ! A vector in the elem x dir
      REAL(DOUBLE)                    :: VY(3)              ! A vector in the elem y dir
      REAL(DOUBLE)                    :: VZ(3)              ! A vector in the elem z dir
      REAL(DOUBLE)                    :: V13(3)             ! A vector from grid 1 to grid 3 (for BAR, BEAM or USER1 it is V vector)



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Initialize

      DO I=1,MELGP
        DO J=1,3
           XEL(I,J) = ZERO
        ENDDO
      ENDDO

      DO I=1,3
        DO J=1,3
           TE(I,J) = ZERO
        ENDDO
      ENDDO

      DO I=1,3
         VX(I) = ZERO
      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------
! Make sure the number of elem grids is not larger than OFFDIS is dimensioned for

      IF ((CAN_ELEM_TYPE_OFFSET == 'Y') .AND. (ELGP > MOFFSET)) THEN
         WRITE(ERR,1954) SUBR_NAME, MOFFSET, ELGP, TYPE
         WRITE(F06,1954) SUBR_NAME, MOFFSET, ELGP, TYPE
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

! Init OFFDIS_B. Some elements will not require the offsets transformed to basic

      IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN
         DO I=1,ELGP
            DO J=1,3
               OFFDIS_B(I,J) = ZERO
            ENDDO
         ENDDO
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! If there are offsets that are specified in global coords, calculate the value of the offsets in basic
! coords so that they can be used below to find the element axes in basic coords

      IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN

         IF (EOFF(INT_ELEM_ID) == 'Y') THEN
            DO I=1,ELGP
               ACID_G = GRID(BGRID(I),3)                   ! Get global coord sys for this grid
               IF (ACID_G /= 0) THEN                       ! Need to transform offset vector from global to basic coords
                  ICID = 0
                  DO J=1,NCORD
                     IF (ACID_G == CORD(J,2)) THEN         ! ACID_G global coord system exists. It was checked in CORDP_PROC
                        ICID = J
                        EXIT
                     ENDIF
                  ENDDO

                  CALL GEN_T0L ( BGRID(I), ICID, THETAD, PHID, T0G )
                  DO J=1,3
                     OFFDIS_B(I,J) = T0G(J,1)*OFFDIS(I,1) + T0G(J,2)*OFFDIS(I,2) + T0G(J,3)*OFFDIS(I,3)
                  ENDDO
               ELSE                                        ! Offset was in basic coords
                  DO J=1,3
                     OFFDIS_B(I,J) = OFFDIS(I,J)
                  ENDDO
               ENDIF
            ENDDO
         ELSE                                              ! There are no offsets so set OFFDIS_B to zero
            DO I=1,ELGP
               DO J=1,3
                  OFFDIS_B(I,J) = ZERO
               ENDDO
            ENDDO
         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Calculate a vector between ends of the element in basic coords (not between grids if there are offsets).

      IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'ROD     ')) THEN
         VX(1) = ( XEB(2,1) + OFFDIS_B(2,1) ) - ( XEB(1,1) + OFFDIS_B(1,1) )
         VX(2) = ( XEB(2,2) + OFFDIS_B(2,2) ) - ( XEB(1,2) + OFFDIS_B(1,2) )
         VX(3) = ( XEB(2,3) + OFFDIS_B(2,3) ) - ( XEB(1,3) + OFFDIS_B(1,3) )
      ELSE
         VX(1) = XEB(2,1) - XEB(1,1)
         VX(2) = XEB(2,2) - XEB(1,2)
         VX(3) = XEB(2,3) - XEB(1,3)
      ENDIF
      LX(1) = VX(1)
      LX(2) = VX(2)
      LX(3) = VX(3)

! Length of element between ends is:

      ELEM_LEN_AB = DSQRT( LX(1)*LX(1) + LX(2)*LX(2) + LX(3)*LX(3) )

! If ELEM_LEN_AB is equal to zero write error and return.

      IF (ELEM_LEN_AB <= EPS1) THEN
         WRITE(ERR,1904) TYPE, EID, ELEM_LEN_AB
         WRITE(F06,1904) TYPE, EID, ELEM_LEN_AB
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Unit vector in element X direction

      DO I=1,3
         TE(1,I) = VX(I)/ELEM_LEN_AB
      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------
! Calculate y and z axes for ROD (since it has no v-vector to help). The y-elem and z-elem axes will be calculated based on the
! procedure referenced below from the internet ("Some Basic Vector Operations In IDL")

      IF (TYPE == 'ROD     ') THEN                         ! NB *** new 09/13/21

         DO I=1,3
            I3_IN(I)   = I
            I3_OUT(I)  = I3_IN(I)
         ENDDO
         CALL CALC_VEC_SORT_ORDER ( VX, SORT_ORDER, I3_OUT)! Use this rather than SORT_INT1_REAL1 - didn't work for vec 10., 0., 0.
         IF (SORT_ORDER == '     ') THEN                   ! Subr CALC_VEC_SORT_ORDER did not find a sort order
            FATAL_ERR = FATAL_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            WRITE(ERR,1944) SUBR_NAME, TYPE, EID
            WRITE(F06,1944) SUBR_NAME, TYPE, EID
            RETURN
         ENDIF
                                                           ! See notes on "Some Basic Vector Operations In IDL" on web site:
                                                           ! http://fermi.jhuapl.edu/s1r/idl/s1rlib/vectors/v_basic.html
         VY(I3_OUT(1)) =  ZERO                             !  (a) Component of VY in direction of min VX is set to zero
         VY(I3_OUT(2)) =  VX(I3_OUT(3))                    !  (b) Other 2 VY(i) are corresponding VX(i) switched with one x(-1)
         VY(I3_OUT(3)) = -VX(I3_OUT(2))
         MAGY  = DSQRT(VY(1)*VY(1) + VY(2)*VY(2) + VY(3)*VY(3))


         IF (DABS(MAGY) < EPS1) THEN
            FATAL_ERR = FATAL_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            WRITE(ERR,1938) SUBR_NAME, 'Y', TYPE, EID, (VY(I),I=1,3)
            WRITE(F06,1938) SUBR_NAME, 'Y', TYPE, EID, (VY(I),I=1,3)
            RETURN
         ENDIF

         DO I=1,3
            TE(2,I) = VY(I)/MAGY
         ENDDO

         CALL CROSS ( VX, VY, VZ )

         MAGZ  = DSQRT(VZ(1)*VZ(1) + VZ(2)*VZ(2) + VZ(3)*VZ(3))

         IF (DABS(MAGZ) < EPS1) THEN
            FATAL_ERR = FATAL_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            WRITE(ERR,1938) SUBR_NAME, 'Z', TYPE, EID, (VZ(I),I=1,3)
            WRITE(F06,1938) SUBR_NAME, 'Z', TYPE, EID, (VZ(I),I=1,3)
            RETURN
         ENDIF

         DO I=1,3
            TE(3,I) = VZ(I)/MAGZ
         ENDDO

! Check if TE is the identity matrix and set a flag

         TE_IDENT = 'N'
         DO I=1,3
            ID(I) = 'N'
         ENDDO
         DO I=1,3
            IF (DABS(TE(I,I) - ONE) < EPS1) THEN
               ID(I) = 'Y'
            ELSE
               ID(I) = 'N'
            ENDIF
         ENDDO
         IF ((ID(1) == 'Y') .AND. (ID(2) == 'Y') .AND. (ID(3) == 'Y')) THEN
            TE_IDENT = 'Y'
         ENDIF

      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Calculate remainder of TE for elements other than ROD

! Calculate V13, vector from G.P.-1 to G.P.-3. For BAR, BEAM, BUDH, USER1 the V13 vector is the v vector = XEB(ELGP+1,i)

begn: IF (TYPE /= 'ROD     ') THEN

         IF ((TYPE == 'BAR     ') .OR. (TYPE == 'BEAM    ') .OR. (TYPE == 'USER1   ')) THEN
            ROWNUM = ELGP + 1
         ELSE
            ROWNUM = 3
         ENDIF
         DO I=1,3
            V13(I) = XEB(ROWNUM,I) - XEB(1,I)
         ENDDO

! Calculate VX x V13 and unit vector in elem z dir. (Col. 3 of TE). If MAGZ is equal to zero, then vector from G.P. 1
! to G.P. 3 is parallel to vector from G.P.-1 to G.P.-2 so write error and quit.

         CALL CROSS ( VX, V13, VZ )
         MAGZ = DSQRT(VZ(1)*VZ(1) + VZ(2)*VZ(2) + VZ(3)*VZ(3))
         IF (MAGZ <=  EPS1) THEN
            IF ((TYPE == 'BAR     ')  .OR. (TYPE == 'BEAM    ')) THEN
               WRITE(ERR,1905) TRIM(TYPE), EID
               WRITE(F06,1905) TRIM(TYPE), EID
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               FATAL_ERR = FATAL_ERR + 1
               RETURN
            ELSE
               WRITE(ERR,1906) TYPE, EID
               WRITE(F06,1906) TYPE, EID
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               FATAL_ERR = FATAL_ERR + 1
               RETURN
            ENDIF
         ENDIF
         DO I=1,3
            TE(3,I) = VZ(I)/MAGZ
         ENDDO

         CALL CROSS ( VZ, VX, VY )                      ! Calc unit vector in Ye dir. (from VZ (cross) VX): If MAGY = 0 quit
         MAGY = DSQRT(VY(1)*VY(1) + VY(2)*VY(2) + VY(3)*VY(3))
         IF(MAGY <= EPS1) THEN
            WRITE(ERR,1912) EID,TYPE
            WRITE(F06,1912) EID,TYPE
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
            RETURN
         ENDIF
         DO I=1,3
            TE(2,I) = VY(I)/MAGY
         ENDDO
                                                       ! Now set TE_IDENT to be 'Y' if TE is an identity matrix.
         TE_IDENT = 'N'
         DO I=1,3
            ID(I) = 'N'
         ENDDO
         DO I=1,3
            IF (DABS(TE(I,I) - ONE) < EPS1) THEN
               ID(I) = 'Y'
            ELSE
               ID(I) = 'N'
            ENDIF
         ENDDO
         IF ((ID(1) == 'Y') .AND. (ID(2) == 'Y') .AND. (ID(3) == 'Y')) THEN
            TE_IDENT = 'Y'
         ENDIF

      ENDIF begn


! ----------------------------------------------------------------------------------------------------------------------------------
! Use TE to get array of elem coords in local system.

      XEL(1,1) = ZERO
      XEL(1,2) = ZERO
      XEL(1,3) = ZERO

      DO I=2,ELGP
         DO J=1,3
            XEL(I,J) = ZERO
            DO K=1,3
               XEL(I,J) = XEL(I,J) + (XEB(I,K) - XEB(1,K))*TE(J,K)
            ENDDO
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************
 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1904 FORMAT(' *ERROR  1904: ',A,I8,' HAS LENGTH BETWEEN ITS ELEMENT ENDS (INCL EFFECTS OF OFFSETS) = ',1ES9.2,'. TOO SMALL')

 1905 FORMAT(' *ERROR  1905: V VECTOR ON ',A,' ELEMENT ',I8,' IS PARALLEL TO VECTOR FROM END A TO END B')

 1906 FORMAT(' *ERROR  1906: ',A,I8,' HAS INTERNAL GRID 3 (FOR V VECTOR) TOO CLOSE TO LINE BETWEEN INTERNAL GRIDS 1 AND 2')

 1912 FORMAT(' *ERROR  1912: CANNOT CALCULATE VECTOR IN ELEMENT Y DIRECTION FOR ELEMENT ',I8,' TYPE ',A)

 1938 FORMAT(' *ERROR  1938: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT CALCULATE PROPER VECTOR IN ELEMENT LOCAL COORDINATE DIRECTION ',A,' FOR ',A,' ELEMENT ',I8,'.' &
                    ,/,14X,' THE VECTOR COMPONENTS CALCULATED WERE ',3(1ES14.6))

 1944 FORMAT(' *ERROR  1944: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE VX VECTOR FOR ',A,' ELEMENT ',I8,' WAS LEFT UNSORTED. IT MUST BE SORTED TO DETERMINE VY, VZ')

 1954 FORMAT(' *ERROR  1954: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' DIMENSION OF ARRAYS OFFDIS, OFFSET ARE ONLY ',I8,' BUT MUST BE ',I8,' FOR ELEM TYPE ',A)

 1959 FORMAT(' *ERROR  1959: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,15X,A,I8,' HAS LENGTH (INCL EFFECTS OF OFFSETS) = ',1ES9.2,'. SHOULD BE ZERO')










! **********************************************************************************************************************************

      END SUBROUTINE ELMGM1


      SUBROUTINE ELMGM1_BUSH ( INT_ELEM_ID, WRITE_WARN )

! Calculates and checks some elem geometry for ROD, BAR, BEAM, triangles and provides a transformation matrix (TE) to transform the
! element stiffness matrix in the element system to the basic coordinate system. Calculates grid point coords in local coord system.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MELGP, MOFFSET, NCORD, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE PARAMS, ONLY                :  EPSIL
      USE MODEL_STUF, ONLY            :  BGRID, BUSH_CID, BUSH_OCID, BUSH_VVEC, CAN_ELEM_TYPE_OFFSET, CORD, EID, ELEM_LEN_12,      &
                                         ELEM_LEN_AB, ELGP, EOFF, GRID, NUM_EMG_FATAL_ERRS, OFFDIS, OFFDIS_L, OFFDIS_O, OFFDIS_B,  &
                                         OFFDIS_G, OFFDIS_GA_GB, RCORD, TE, TE_GA_GB, TE_IDENT, TYPE, XEB, XEL

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE VECTOR_GEOMETRY, ONLY       :  CROSS, GEN_T0L
      USE SORTING, ONLY               :  CALC_VEC_SORT_ORDER
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      IMPLICIT NONE

      CHARACTER( 1*BYTE)              :: ID(3)              ! Used in deciding whether TE_IDENT = 'Y' or 'N'
      CHARACTER( 5*BYTE)              :: SORT_ORDER         ! Order in which the VX(i) have been sorted in subr CALC_VEC_SORT_ORDER
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMGM1_BUSH'
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN         ! If 'Y" write warning messages, otherwise do not
      CHARACTER(1*BYTE)               :: CORD_FND           ! = 'Y' if coord sys ID on CONM2 defined, 'N' otherwise
      CHARACTER(1*BYTE)               :: DO_IT              ! = 'Y' execute code that follows (see where DO_IT is initialized)

      INTEGER(LONG), INTENT(IN)       :: INT_ELEM_ID        ! Internal element ID for which
      INTEGER(LONG)                   :: ACID_G             ! Actual coordinate system ID
      INTEGER(LONG)                   :: I,J,K              ! DO loop indices
      INTEGER(LONG)                   :: I3_IN(3)           ! Integer array used in sorting VX.

      INTEGER(LONG)                   :: I3_OUT(3)          ! Integer array giving order of VX comps. If VX is in the order with
!                                                             comp 2 smallest then comp 3 then comp 1 then I3_OUT is 2, 3, 1

      INTEGER(LONG)                   :: ICID               ! Internal coord sys no. corresponding to an actual coord sys no.
      INTEGER(LONG)                   :: ROWNUM             ! A row number in an array


      REAL(DOUBLE)                    :: DX1(3)             ! Array used in intermediate calc's
      REAL(DOUBLE)                    :: DX2(3)             ! Array used in intermediate calc's
      REAL(DOUBLE)                    :: EPS1               ! A small number to compare to real zero
      REAL(DOUBLE)                    :: LX(3)              ! Distances
      REAL(DOUBLE)                    :: MAGY               ! Magnitude of vector VY
      REAL(DOUBLE)                    :: MAGZ               ! Magnitude of vector VZ
      REAL(DOUBLE)                    :: PHID, THETAD       ! Outputs from subr GEN_T0L
      REAL(DOUBLE)                    :: TET_GA_GB(3,3)     ! Transpose of TE_GA_GB
      REAL(DOUBLE)                    :: T0G(3,3)           ! Matrix to transform offsets from global to basic  coords
      REAL(DOUBLE)                    :: TG0(3,3)           ! Matrix to transform offsets from basic  to global coords
      REAL(DOUBLE)                    :: T0I(3,3)           ! Matrix to transform offsets from BUSH OCID to basic coords
      REAL(DOUBLE)                    :: TET(3,3)           ! Transpose of TE: UEL = TE*UEB
      REAL(DOUBLE)                    :: VX(3)              ! A vector in the elem x dir
      REAL(DOUBLE)                    :: VY(3)              ! A vector in the elem y dir
      REAL(DOUBLE)                    :: VZ(3)              ! A vector in the elem z dir
      REAL(DOUBLE)                    :: V13(3)             ! A vector from grid 1 to grid 3 (for BAR, BEAM or USER1 it is V vector)
      REAL(DOUBLE)                    :: XIB(2,3)           ! Coords at ends of BUSH elem: XIB(1,J) should = XIB(2,J) for 0 len elem



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Initialize

      DO I=1,MELGP
        DO J=1,3
           XEL(I,J) = ZERO
        ENDDO
      ENDDO

      DO I=1,3
        DO J=1,3
           TE(I,J) = ZERO
        ENDDO
      ENDDO

      DO I=1,3
         VX(I) = ZERO
      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------
! Make sure the number of elem grids is not larger than OFFDIS is dimensioned for

      IF ((CAN_ELEM_TYPE_OFFSET == 'Y') .AND. (ELGP > MOFFSET)) THEN
         WRITE(ERR,1954) SUBR_NAME, MOFFSET, ELGP, TYPE
         WRITE(F06,1954) SUBR_NAME, MOFFSET, ELGP, TYPE
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Calculate the TE coord transformation matrix

      IF (ELEM_LEN_12 < .0001) THEN
         IF (BUSH_CID < 0) THEN                               ! Elem len < .0001 so CID must be > 0
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1801) TRIM(TYPE), EID, '0.0001'
            WRITE(F06,1801) TRIM(TYPE), EID, '0.0001'
         ENDIF
      ENDIF

bcid: IF (BUSH_CID > 0) THEN                               ! Get transformation of BUSH_CID to basic and put it into TE

         CORD_FND = 'N'
         ICID = 0
         DO J=1,NCORD
            IF (BUSH_CID == CORD(J,2)) THEN
               CORD_FND = 'Y'
               ICID = J
               EXIT
            ENDIF
         ENDDO

         IF (CORD_FND == 'Y') THEN                         ! NOTE: This BUSH TE transforms a vector in BUSH_CID to basic
            TET(1,1) = RCORD(ICID, 4)   ;   TET(1,2) = RCORD(ICID, 5)   ;   TET(1,3) = RCORD(ICID, 6)
            TET(2,1) = RCORD(ICID, 7)   ;   TET(2,2) = RCORD(ICID, 8)   ;   TET(2,3) = RCORD(ICID, 9)
            TET(3,1) = RCORD(ICID,10)   ;   TET(3,2) = RCORD(ICID,11)   ;   TET(3,3) = RCORD(ICID,12)
            DO I=1,3
               DO J=1,3
                  TE(I,J) = TET(J,I)
               ENDDO
            ENDDO
            write(f06,*)
         ELSE                                              ! Error, could not find BUSH_CID coord system in RCORD
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1822) 'COORD SYSTEM ', BUSH_CID, TYPE, EID
            WRITE(F06,1822) 'COORD SYSTEM ', BUSH_CID, TYPE, EID
         ENDIF

      ELSE IF (BUSH_CID == 0) THEN                         ! BUSH_CID is basic so TE is the identity matrix

         DO I=1,3
            DO J=1,3
               TE(I,J) = ZERO
            ENDDO
            TE(I,I) = ONE
         ENDDO

      ELSE IF (BUSH_CID < 0) THEN                          ! This means no BUSH_CID was found so we look for the v vector

         IF (BUSH_VVEC /= 0) THEN                          ! A v-vector was specified for this bUSH element

            VX(1) = XEB(2,1) - XEB(1,1)                    ! When no BUSH_CID the x axis is along line between the 2 grids
            VX(2) = XEB(2,2) - XEB(1,2)
            VX(3) = XEB(2,3) - XEB(1,3)
            DO I=1,3
               TE(1,I) = VX(I)/ELEM_LEN_12                 ! For BUSH use length bet 2 grids since length bet elem ends is zero
            ENDDO
            DO I=1,3
               V13(I) = XEB(ELGP+1,I) - XEB(1,I)
            ENDDO
            CALL CROSS ( VX, V13, VZ )
            MAGZ = DSQRT(VZ(1)*VZ(1) + VZ(2)*VZ(2) + VZ(3)*VZ(3))
            IF (MAGZ <=  EPS1) THEN
               WRITE(ERR,1906) TYPE, EID
               WRITE(F06,1906) TYPE, EID
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               FATAL_ERR = FATAL_ERR + 1
               RETURN
            ENDIF
            DO I=1,3
               TE(3,I) = VZ(I)/MAGZ
            ENDDO

            CALL CROSS ( VZ, VX, VY )                      ! Calc unit vector in Ye dir.from VZ (cross) VX: If MAGY = 0 quit
            MAGY = DSQRT(VY(1)*VY(1) + VY(2)*VY(2) + VY(3)*VY(3))
            IF(MAGY <= EPS1) THEN
               WRITE(ERR,1912) EID,TYPE
               WRITE(F06,1912) EID,TYPE
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               FATAL_ERR = FATAL_ERR + 1
               RETURN
            ENDIF
            DO I=1,3
               TE(2,I) = VY(I)/MAGY
            ENDDO

         ELSE IF (BUSH_VVEC == 0) THEN                     ! No v-vec or CID specified for BUSH elem so elem x axis is GA->GB
!                                                                  --
            VX(1) = XEB(2,1) - XEB(1,1)
            VX(2) = XEB(2,2) - XEB(1,2)
            VX(3) = XEB(2,3) - XEB(1,3)

            DO I=1,3
               TE(1,I) = VX(I)/ELEM_LEN_12
            ENDDO

            DO I=1,3
               I3_IN(I)   = I
               I3_OUT(I)  = I3_IN(I)
            ENDDO
            CALL CALC_VEC_SORT_ORDER (VX,SORT_ORDER,I3_OUT)! Use this rather than SORT_INT1_REAL1 - didn't work for vec 10.,0.,0.
            IF (SORT_ORDER == '     ') THEN                ! Subr CALC_VEC_SORT_ORDER did not find a sort order
               FATAL_ERR = FATAL_ERR + 1
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               WRITE(ERR,1944) SUBR_NAME, TYPE, EID
               WRITE(F06,1944) SUBR_NAME, TYPE, EID
               RETURN
            ENDIF
                                                           ! See notes on "Some Basic Vector Operations In IDL" on web site:
                                                           ! http://fermi.jhuapl.edu/s1r/idl/s1rlib/vectors/v_basic.html
            VY(I3_OUT(1)) =  ZERO                          !  (a) Component of VY in direction of min VX is set to zero
            VY(I3_OUT(2)) =  VX(I3_OUT(3))                 !  (b) Other 2 VY(i) are corresponding VX(i) switched with one x(-1)
            VY(I3_OUT(3)) = -VX(I3_OUT(2))
            MAGY  = DSQRT(VY(1)*VY(1) + VY(2)*VY(2) + VY(3)*VY(3))


            IF (DABS(MAGY) < EPS1) THEN
               FATAL_ERR = FATAL_ERR + 1
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               WRITE(ERR,1938) SUBR_NAME, 'Y', TYPE, EID, (VY(I),I=1,3)
               WRITE(F06,1938) SUBR_NAME, 'Y', TYPE, EID, (VY(I),I=1,3)
               RETURN
            ENDIF

            DO I=1,3
               TE(2,I) = VY(I)/MAGY
            ENDDO

            CALL CROSS ( VX, VY, VZ )

            MAGZ  = DSQRT(VZ(1)*VZ(1) + VZ(2)*VZ(2) + VZ(3)*VZ(3))

            IF (DABS(MAGZ) < EPS1) THEN
               FATAL_ERR = FATAL_ERR + 1
               NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
               WRITE(ERR,1938) SUBR_NAME, 'Z', TYPE, EID, (VZ(I),I=1,3)
               WRITE(F06,1938) SUBR_NAME, 'Z', TYPE, EID, (VZ(I),I=1,3)
               RETURN
            ENDIF

            DO I=1,3
               TE(3,I) = VZ(I)/MAGZ
            ENDDO

         ENDIF

      ENDIF bcid

      TE_IDENT = 'N'                                       ! Now set TE_IDENT to be 'Y' if TE is an identity matrix.
      DO I=1,3
         ID(I) = 'N'
      ENDDO
      DO I=1,3
         IF (DABS(TE(I,I) - ONE) < EPS1) THEN
            ID(I) = 'Y'
         ELSE
            ID(I) = 'N'
         ENDIF
      ENDDO
      IF ((ID(1) == 'Y') .AND. (ID(2) == 'Y') .AND. (ID(3) == 'Y')) THEN
         TE_IDENT = 'Y'
      ENDIF


! ----------------------------------------------------------------------------------------------------------------------------------
! Use TE to get array of elem coords in local system.

      XEL(1,1) = ZERO
      XEL(1,2) = ZERO
      XEL(1,3) = ZERO

      DO I=2,ELGP
         DO J=1,3
            XEL(I,J) = ZERO
            DO K=1,3
               XEL(I,J) = XEL(I,J) + (XEB(I,K) - XEB(1,K))*TE(J,K)
            ENDDO
         ENDDO
      ENDDO

      DO I=1,3
         DO J=1,3
            TET(I,J) = TE(J,I)
         ENDDO
      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------
! Use TE to get array of elem coords in local system.

      XEL(1,1) = ZERO
      XEL(1,2) = ZERO
      XEL(1,3) = ZERO

      DO I=2,ELGP
         DO J=1,3
            XEL(I,J) = ZERO
            DO K=1,3
               XEL(I,J) = XEL(I,J) + (XEB(I,K) - XEB(1,K))*TE(J,K)
            ENDDO
         ENDDO
      ENDDO

! ----------------------------------------------------------------------------------------------------------------------------------
! Calculate a coord transformation matrix, TE_GA_GB, that will transform a vector whose x axis is along the GA-GB axis to basic:
! This will be used for only the x axis when the BUSH offset is along the line GA-GB in order to specify the offsets whenever
! BUSH_OCID is -1 (default) or blank. The other y and z aves of TE_GA_GB don't matter except thet they be normal to x. So use the
! procedure outlined in "Some Basic Vector Operations In IDL" (below) to find the y and z axes of TE_GA_GB
! This needs to be done if GA and GB are not coincident but also this coord transformation is used in calculating element forces
! from nodal loads (i.e. in subr OFP3_ELFE_1D)

      IF (ELEM_LEN_12 > .0001D0) THEN

         VX(1) = XEB(2,1) - XEB(1,1)                          ! Unit vector in element X direction
         VX(2) = XEB(2,2) - XEB(1,2)
         VX(3) = XEB(2,3) - XEB(1,3)

         DO I=1,3
            TE_GA_GB(1,I) = VX(I)/ELEM_LEN_12                 ! Row 1 of TE_GA_GB
         ENDDO                                                ! -----------------

         DO I=1,3
            I3_IN(I)   = I
            I3_OUT(I)  = I3_IN(I)
         ENDDO
         CALL CALC_VEC_SORT_ORDER ( VX, SORT_ORDER, I3_OUT)   ! Use this rather than SORT_INT1_REAL1 - didn't work for vec 10., 0., 0.
         IF (SORT_ORDER == '     ') THEN                      ! Subr CALC_VEC_SORT_ORDER did not find a sort order so fatal error
            FATAL_ERR = FATAL_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            WRITE(ERR,1944) SUBR_NAME, TYPE, EID
            WRITE(F06,1944) SUBR_NAME, TYPE, EID
            RETURN
         ENDIF
                                                              ! See notes on "Some Basic Vector Operations In IDL" on web site:
                                                              ! http://fermi.jhuapl.edu/s1r/idl/s1rlib/vectors/v_basic.html
         VY(I3_OUT(1)) =  ZERO                                !  (a) Component of VY in direction of min VX is set to zero
         VY(I3_OUT(2)) =  VX(I3_OUT(3))                       !  (b) Other 2 VY(i) are corresponding VX(i) switched with one x(-1)
         VY(I3_OUT(3)) = -VX(I3_OUT(2))
         MAGY  = DSQRT(VY(1)*VY(1) + VY(2)*VY(2) + VY(3)*VY(3))


         IF (DABS(MAGY) < EPS1) THEN
            FATAL_ERR = FATAL_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            WRITE(ERR,1938) SUBR_NAME, 'Y', TYPE, EID, (VY(I),I=1,3)
            WRITE(F06,1938) SUBR_NAME, 'Y', TYPE, EID, (VY(I),I=1,3)
            RETURN
         ENDIF

         DO I=1,3
            TE_GA_GB(2,I) = VY(I)/MAGY                        ! Row 2 of TE_GA_GB
         ENDDO                                                ! -----------------

         CALL CROSS ( VX, VY, VZ )

         MAGZ  = DSQRT(VZ(1)*VZ(1) + VZ(2)*VZ(2) + VZ(3)*VZ(3))

         IF (DABS(MAGZ) < EPS1) THEN
            FATAL_ERR = FATAL_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            WRITE(ERR,1938) SUBR_NAME, 'Z', TYPE, EID, (VZ(I),I=1,3)
            WRITE(F06,1938) SUBR_NAME, 'Z', TYPE, EID, (VZ(I),I=1,3)
            RETURN
         ENDIF

         DO I=1,3
            TE_GA_GB(3,I) = VZ(I)/MAGZ                        ! Row 3 of TE_GA_GB
         ENDDO                                                ! -----------------

         DO I=1,3                                             ! Transpose of TE_GA_GB
            DO J=1,3
               TET_GA_GB(I,J) = TE_GA_GB(J,I)
            ENDDO
         ENDDO


      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
! Offsets for BUSH are specified in a unique coord system. It may or may not be basic or global. So we first transform the input
! offset values to basic and then transform them to global. This way, when we process the offsets in subr ELEM_TRANSFORM_LBG,
! we can treat the final offset values for the BUSH the same as we do for BAR, BEAM or ROD.
! The coord system for BUSH offsets is BUSH_OCID which can be:
!  (1) -1 means offset lies on the line between GA and GB, or
!  (2)  a positive number indicating a coord system number

      DO I=1,ELGP
         DO J=1,3                                          ! Initialize
            OFFDIS_B(I,J) = ZERO
            OFFDIS_G(I,J) = ZERO
         ENDDO
      ENDDO

      IF (EOFF(INT_ELEM_ID) == 'Y') THEN

         IF (BUSH_OCID >= 0) THEN
            IF (BUSH_OCID > 0) THEN                        ! BUSH offsets are defined in coord system BUSH_OCID
               CORD_FND = 'N'
               ICID = 0
               DO J=1,NCORD
                  IF (BUSH_OCID == CORD(J,2)) THEN
                     CORD_FND = 'Y'
                     ICID = J
                     EXIT
                  ENDIF
               ENDDO
               IF (CORD_FND == 'Y') THEN
                  CALL GEN_T0L ( BGRID(1), ICID, THETAD, PHID, T0I )
                  DO J=1,3
                     OFFDIS_B(1,J) = T0I(J,1)*OFFDIS_O(1,1) + T0I(J,2)*OFFDIS_O(1,2) + T0I(J,3)*OFFDIS_O(1,3)
                  ENDDO
                  write(f06,*)
               ELSE
                  NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1822) 'COORD SYSTEM ', BUSH_OCID, TYPE, EID
                  WRITE(F06,1822) 'COORD SYSTEM ', BUSH_OCID, TYPE, EID
               ENDIF

            ELSE IF (BUSH_OCID == 0) THEN                  ! Offset is in basic coords
               DO J=1,3
                  OFFDIS_B(1,J) = OFFDIS_O(1,J)
               ENDDO

            ENDIF

         ELSE IF (BUSH_OCID < 0) THEN                      ! Offsets are along the line between GA and GB

            DX1(1) = OFFDIS_O(1,1)                         ! Transform GA offsets (at 1st grid point relative to GA-GB axis) to basic
            DX1(2) = OFFDIS_O(1,2)
            DX1(3) = OFFDIS_O(1,3)
            CALL MATMULT_FFF (TET_GA_GB, DX1, 3, 3, 1, DX2 )
            OFFDIS_B(1,1) = DX2(1)
            OFFDIS_B(1,2) = DX2(2)
            OFFDIS_B(1,3) = DX2(3)

            DX1(1) = OFFDIS_O(2,1)                         ! Transform GA offsets (at 2nd grid point relative to GA-GB axis) to basic
            DX1(2) = OFFDIS_O(2,2)
            DX1(3) = OFFDIS_O(2,3)
            CALL MATMULT_FFF (TET_GA_GB, DX1, 3, 3, 1, DX2 )
            OFFDIS_B(2,1) = DX2(1)
            OFFDIS_B(2,2) = DX2(2)
            OFFDIS_B(2,3) = DX2(3)

         ENDIF

      ENDIF

! Now that we have the offsets at grid 1 in basic coords we can calc offset values for grid 2 and then transform offsets to global.

      IF (EOFF(INT_ELEM_ID) == 'Y') THEN

         IF (BUSH_OCID /= -1) THEN                         ! We need to calc OFFDIS_B for grid 2. If OCID = -1 then that OFFDIS
                                                           ! was already calculated in subr ELMDAT1 since it lay along line GA-GB
            DO I=1,3
               XEB(3,I)      = XEB(1,I) + OFFDIS_B(1,I)    ! Put BUSH basic coords in XEB (has extra row for more than the 2 grids)
               OFFDIS_B(2,I) = XEB(3,I) - XEB(2,I)         ! Offset from GB is location of BUSH minus location of GB
            ENDDO

         ENDIF

         DO I=1,2
            ACID_G = GRID(BGRID(I),3)                      ! Get global coord sys for this grid
            IF (ACID_G /= 0) THEN                          ! Global is not basic so need to transform offset from basic to global
               ICID = 0
               DO J=1,NCORD
                  IF (ACID_G == CORD(J,2)) THEN
                     ICID = J
                     EXIT
                  ENDIF
               ENDDO
               CALL GEN_T0L ( BGRID(I), ICID, THETAD, PHID, T0G )
               DO J=1,3                                    ! We want transpose of T0G since we are transforming from basic to global
                  DO K=1,3
                     TG0(J,K) = T0G(K,J)
                  ENDDO
               ENDDO
               DO J=1,3
                  OFFDIS_G(I,J) = TG0(J,1)*OFFDIS_B(I,1) + TG0(J,2)*OFFDIS_B(I,2) + TG0(J,3)*OFFDIS_B(I,3)
               ENDDO
            ELSE                                           ! Global was basic so no transformation of coords needed
               DO J=1,3
                  OFFDIS_G(I,J) = OFFDIS_B(I,J)
               ENDDO
            ENDIF
         ENDDO

         DO I=1,ELGP                                       ! Put global offsets in OFFDIS since these are the ones needed when
            DO J=1,3                                       ! offsets are applied in subr ELEM_TRANSFORM_LBG
               OFFDIS(I,J) = OFFDIS_G(I,J)
            ENDDO
         ENDDO

      ELSE                                                 ! There are no offsets so set OFFDIS_B to zero for both grids

         DO I=1,2
            DO J=1,3
               OFFDIS_B(I,J) = ZERO
            ENDDO
         ENDDO

      ENDIF

! Calculate offsets in local, CID, axes

      DO I=1,3
         DX1(I) = OFFDIS_B(1,I)
      ENDDO
      CALL MATMULT_FFF ( TE, DX1, 3, 3, 1, DX2 )

      DO I=1,3
         OFFDIS_L(1,I) = DX2(I)
         DX1(I) = OFFDIS_B(2,I)
      ENDDO
      CALL MATMULT_FFF ( TE, DX1, 3, 3, 1, DX2 )

      DO I=1,3
         OFFDIS_L(2,I) = DX2(I)
      ENDDO

! Calculate offsets in axes along and normal to BUSH axes (i.e. axes that have x along line between GA-GB)
! This will be needed in subr OFP3_ELFE_1D when BUSH element forces are calculated from node forces

      DO I=1,3
         DX1(I) = OFFDIS_B(1,I)
      ENDDO
      CALL MATMULT_FFF ( TE_GA_GB, DX1, 3, 3, 1, DX2 )

      DO I=1,3
         OFFDIS_GA_GB(1,I) = DX2(I)
         DX1(I) = OFFDIS_B(2,I)
      ENDDO
      CALL MATMULT_FFF ( TE_GA_GB, DX1, 3, 3, 1, DX2 )

      DO I=1,3
         OFFDIS_GA_GB(2,I) = DX2(I)
      ENDDO






! ----------------------------------------------------------------------------------------------------------------------------------
      LX(1) = ( XEB(2,1) + OFFDIS_B(2,1) ) - ( XEB(1,1) + OFFDIS_B(1,1) )
      LX(2) = ( XEB(2,2) + OFFDIS_B(2,2) ) - ( XEB(1,2) + OFFDIS_B(1,2) )
      LX(3) = ( XEB(2,3) + OFFDIS_B(2,3) ) - ( XEB(1,3) + OFFDIS_B(1,3) )
      ELEM_LEN_AB = DSQRT( LX(1)*LX(1) + LX(2)*LX(2) + LX(3)*LX(3) )
                                                           ! If ELEM_LEN_AB is not equal to zero then write error and return.
      IF (ELEM_LEN_AB > .0001D0) THEN
         WRITE(ERR,1959) SUBR_NAME, TYPE, EID, ELEM_LEN_AB
         WRITE(F06,1959) SUBR_NAME, TYPE, EID, ELEM_LEN_AB
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         RETURN
      ENDIF

! ----------------------------------------------------------------------------------------------------------------------------------
      IF (DEBUG(110) > 0) THEN
         CALL DEBUG_ELMGM1_FOR_BUSH
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1801 FORMAT(' *ERROR  1801 ',A,'element ',I8,' has length (between grids) less than ',A,' so it must have a CID specified')

 1822 FORMAT(' *ERROR  1822: ',A,I8,' ON ',A,I8,' IS UNDEFINED')

 1906 FORMAT(' *ERROR  1906: ',A,I8,' HAS INTERNAL GRID 3 (FOR V VECTOR) TOO CLOSE TO LINE BETWEEN INTERNAL GRIDS 1 AND 2')

 1912 FORMAT(' *ERROR  1912: CANNOT CALCULATE VECTOR IN ELEMENT Y DIRECTION FOR ELEMENT ',I8,' TYPE ',A)

 1938 FORMAT(' *ERROR  1938: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' CANNOT CALCULATE PROPER VECTOR IN ELEMENT LOCAL COORDINATE DIRECTION ',A,' FOR ',A,' ELEMENT ',I8,'.' &
                    ,/,14X,' THE VECTOR COMPONENTS CALCULATED WERE ',3(1ES14.6))

 1944 FORMAT(' *ERROR  1944: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' THE VX VECTOR FOR ',A,' ELEMENT ',I8,' WAS LEFT UNSORTED. IT MUST BE SORTED TO DETERMINE VY, VZ')

 1954 FORMAT(' *ERROR  1954: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' DIMENSION OF ARRAYS OFFDIS, OFFSET ARE ONLY ',I8,' BUT MUST BE ',I8,' FOR ELEM TYPE ',A)

 1959 FORMAT(' *ERROR  1959: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,15X,A,I8,' HAS LENGTH (INCL EFFECTS OF OFFSETS) = ',1ES9.2,'. SHOULD BE ZERO')





! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEBUG_ELMGM1_FOR_BUSH

      IMPLICIT NONE

! **********************************************************************************************************************************

      WRITE(F06,*)
      WRITE(F06,98720)
      WRITE(F06,'(A,I8)') ' BUSH element number ', EID
      WRITE(F06,'(A)')    ' ============================'
      WRITE(F06,*)

      IF (BUSH_CID >= 0) THEN
         WRITE(F06,'(A,I8)') '    The element coordinate system will be BUSH_CID coord system   ',bush_cid
         WRITE(F06,*)
      ELSE
         WRITE(F06,'(A)') '    The element coordinate system will be determined from the 2 grids and the specified VVEC'
         WRITE(F06,*)
      ENDIF

      IF (EOFF(INT_ELEM_ID) == 'Y') THEN

         WRITE(F06,*) '   OFFDIS_O array of offsets based on input on CBUSH and ELEM_LEN_12 (before any coord transformations)'
         WRITE(F06,'(A,3(1ES14.6))') '                                     End A  = ', (OFFDIS_O(1,j),j=1,3)
         WRITE(F06,'(A,3(1ES14.6))') '                                     End B  = ', (OFFDIS_O(2,j),j=1,3)
         WRITE(F06,*)

         WRITE(F06,*) '   OFFDIS_B array of offsets in basic coords:'
         WRITE(F06,'(A,3(1ES14.6))') '                                     End A  = ', (OFFDIS_B(1,j),j=1,3)
         WRITE(F06,'(A,3(1ES14.6))') '                                     End B  = ', (OFFDIS_B(2,j),j=1,3)
         WRITE(F06,*)

         WRITE(F06,*) '   OFFDIS_G array of offsets transformed to global coords:'
         WRITE(F06,'(A,3(1ES14.6))') '                                     END A  = ', (OFFDIS_G(1,J),J=1,3)
         WRITE(F06,'(A,3(1ES14.6))') '                                     End B  = ', (OFFDIS_G(2,j),j=1,3)
         WRITE(F06,*)

         WRITE(F06,*) '   OFFDIS   array of offsets that will be used as global offsets (must equal OFFDIS_G above)'
         WRITE(F06,'(A,3(1ES14.6))') '                                     End A  = ', (OFFDIS(1,j),j=1,3)
         WRITE(F06,'(A,3(1ES14.6))') '                                     End B  = ', (OFFDIS(2,j),j=1,3)
         WRITE(F06,*)

      ENDIF

      WRITE(F06,'(A,3(1ES14.6))') '    ELEM_LEN_AB                             = ',ELEM_LEN_AB
      WRITE(F06,'(A,3(1ES14.6))') '    ELEM_LEN_12                             = ',ELEM_LEN_12

      WRITE(F06,*)

      WRITE(F06,98799)
      WRITE(F06,*)

! **********************************************************************************************************************************
98720 FORMAT(' __________________________________________________________________________________________________________________',&
             '_________________'                                                                                               ,//,&
             ' ::::::::::::::::::::::::::::::::::::::::::START DEBUG(110) OUTPUT FROM SUBROUTINE ELMGM1::::::::::::::::::::::::::',&
             ':::::::::::::::::',/)

98799 FORMAT(' :::::::::::::::::::::::::::::::::::::::::::END DEBUG(110) OUTPUT FROM SUBROUTINE ELMGM1:::::::::::::::::::::::::::',&
             ':::::::::::::::::'                                                                                                ,/,&
             ' __________________________________________________________________________________________________________________',&
             '_________________',/)

! **********************************************************************************************************************************

      END SUBROUTINE DEBUG_ELMGM1_FOR_BUSH

      END SUBROUTINE ELMGM1_BUSH


      SUBROUTINE ELMGM2 ( WRITE_WARN )

! Calcs and checks elem geometry for quad elems and provides a transformation matrix ( TE ) to transfer the elem stiffness matrix
! in the elem system to the basic coordinate system. Calculates grid point coords in local coord system.
! To define the elem coordinate system, a mean plane is defined which lies midway between the grid points (HBAR is mean dist).
! The elem z direction is in the direction of the cross product of the diagonals (V13 x V24). Initially, the x axis is along
! side 1-2 of the elem projection onto the mean plane. For elems thet are not rectangular, the x,y axes are rotated such that x
! splits the angle between the diagonals.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG, ERR, F06, WRT_BUG, WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, MEFE, MEWE, MELGP, FATAL_ERR, WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, TWO
      USE PARAMS, ONLY                :  EPSIL, QUADAXIS, SUPWARN, QUAD4TYP
      USE MODEL_STUF, ONLY            :  AGRID, BMEANT, EID, ELGP, EMG_IFE, EMG_IWE, EMG_RWE, ERR_SUB_NAM, NUM_EMG_FATAL_ERRS,     &
                                         HBAR, MXWARP, QUAD_DELTA, QUAD_GAMMA, QUAD_THETA, TE, TE_IDENT, TYPE, WARP_WARN, XEB, XEL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE VECTOR_GEOMETRY, ONLY       :  CROSS, PLANE_COORD_TRANS_21
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMGM2'
      CHARACTER(30*BYTE)              :: NAME(15)          ! Names for BUG output purposes
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not

      INTEGER(LONG)                   :: DIAG_GRID1        ! Used for error output purposes
      INTEGER(LONG)                   :: DIAG_GRID2        ! Used for error output purposes
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: QUAD_GEOM_ERR = 0 ! Local error count
      INTEGER(LONG)                   :: ID(3)             ! ID(i) is set to 1 if the i-th diagonal of TE is 1.0
      INTEGER(LONG)                   :: IPNT              ! An internal grid (1,2,3 or 4) of an elem
      INTEGER(LONG)                   :: SIDE_GRID1        ! Used for error output purposes
      INTEGER(LONG)                   :: SIDE_GRID2        ! Used for error output purposes


      REAL(DOUBLE)                    :: V12B(3)           ! Vector from G.P. 1 to G.P. 2 in basic coords
      REAL(DOUBLE)                    :: V13B(3)           ! Vector from G.P. 1 to G.P. 3 in basic coords (a diagonal)
      REAL(DOUBLE)                    :: V24B(3)           ! Vector from G.P. 2 to G.P. 4 in basic coords (a diagonal)
      REAL(DOUBLE)                    :: V13BM             ! Mag of V13B
      REAL(DOUBLE)                    :: V24BM             ! Mag of V24B
      REAL(DOUBLE)                    :: DHBAR             ! DABS(HBAR)
      REAL(DOUBLE)                    :: EPS1              ! A small number to compare to real zero
      REAL(DOUBLE)                    :: EPS4              ! A small number to compare to real zero
      REAL(DOUBLE)                    :: IVEC(3)           ! A vector in the elem x dir
      REAL(DOUBLE)                    :: JVEC(3)           ! A vector in the elem y dir
      REAL(DOUBLE)                    :: KVEC(3)           ! A vector in the elem z dir
      REAL(DOUBLE)                    :: L12               ! Length of side 1-2 of the elem in the mean plane
      REAL(DOUBLE)                    :: L23               ! Length of side 2-3 of the elem in the mean plane
      REAL(DOUBLE)                    :: L34               ! Length of side 3-4 of the elem in the mean plane
      REAL(DOUBLE)                    :: L41               ! Length of side 4-1 of the elem in the mean plane
      REAL(DOUBLE)                    :: MAGI              ! Magnitude of vector IVEC
      REAL(DOUBLE)                    :: MAGJ              ! Magnitude of vector JVEC
      REAL(DOUBLE)                    :: MAGK              ! Magnitude of vector KVEC
      REAL(DOUBLE)                    :: CT_QD(3,3)        ! Coord transf matrix which will rotate a vector thru an angle QUAD_DELTA
      REAL(DOUBLE)                    :: TE_12(3,3)        ! TE matrix for this elem if local x is parallel to side 1-2
      REAL(DOUBLE)                    :: TE_SD(3,3)        ! TE matrix for this elem if local x splits angle between the two diags
      REAL(DOUBLE)                    :: X12               ! (X1 - X2) in elem mean plane in local elem coords
      REAL(DOUBLE)                    :: X13               ! (X1 - X3) in elem mean plane in local elem coords
      REAL(DOUBLE)                    :: X14               ! (X1 - X4) in elem mean plane in local elem coords
      REAL(DOUBLE)                    :: X23               ! (X2 - X3) in elem mean plane in local elem coords
      REAL(DOUBLE)                    :: X24               ! (X2 - X4) in elem mean plane in local elem coords
      REAL(DOUBLE)                    :: X34               ! (X3 - X4) in elem mean plane in local elem coords
      REAL(DOUBLE)                    :: X3                ! -X13
      REAL(DOUBLE)                    :: X4                ! -X14
      REAL(DOUBLE)                    :: Y3                !
      REAL(DOUBLE)                    :: Y4                !
      REAL(DOUBLE)                    :: Y34               ! (Y3 - Y4) in elem mean plane in local elem coords
      REAL(DOUBLE)                    :: VAR(15)           ! Variables for BUG output purposes

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)
      EPS4 = EPSIL(4)

! Initialize XEL to zero

      DO I=1,MELGP
         DO J=1,3
            XEL(I,J) = ZERO
         ENDDO
      ENDDO

! **********************************************************************************************************************************
! Calculate elem z direction from cross products of diagonals

! Generate vectors from G.P 1 to G.P 3 and from G.P. 2 to G.P. 4 (diagonals)

      DO I=1,3
         V13B(I) = XEB(3,I) - XEB(1,I)
         V24B(I) = XEB(4,I) - XEB(2,I)
      ENDDO

      CALL CROSS ( V13B, V24B, KVEC )
      MAGK = DSQRT(KVEC(1)*KVEC(1) + KVEC(2)*KVEC(2) + KVEC(3)*KVEC(3))

! If MAGK = 0 then diagonals are parallel so write error and quit

      IF (MAGK <=  EPS1) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1911) TYPE, EID
            WRITE(F06,1911) TYPE, EID
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1911
            ENDIF
         ENDIF
         RETURN
      ENDIF

! Unit vector in elem local z direction is 3rd row of TE

      DO I=1,3
         KVEC(I)     = KVEC(I)/MAGK
         TE_12(3,I) = KVEC(I)
      ENDDO

! **********************************************************************************************************************************
! Calc initial elem x dir along side 1-2 of the elem projection onto the mean plane.

      DO I=1,3
         V12B(I) = XEB(2,I) - XEB(1,I)
      ENDDO

! HBAR is one half of the projection of V12B in z direction

      HBAR = HALF*((V12B(1)*KVEC(1) + V12B(2)*KVEC(2) + V12B(3)*KVEC(3)))

! Now calculate initial x direction along side 1-2 of the elem projection onto the mean plane.

      DO I=1,3
         IVEC(I) = V12B(I) - TWO*HBAR*KVEC(I)
      ENDDO

! If initial MAGI = 0 then write error and quit.
      MAGI  = DSQRT(IVEC(1)*IVEC(1) + IVEC(2)*IVEC(2) + IVEC(3)*IVEC(3))
      IF (MAGI <= EPS1) THEN
         SIDE_GRID1 = 1
         SIDE_GRID2 = 2
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
            WRITE(F06,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1908
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = SIDE_GRID1
               EMG_IFE(NUM_EMG_FATAL_ERRS,3) = SIDE_GRID2
            ENDIF
         ENDIF
         RETURN
      ENDIF

! Unit vector in initial elem x direction

      DO I=1,3
         IVEC(I)     = IVEC(I)/MAGI                        ! Unit vector along side 1-2 in the mean plane (NOT from G.P. 1-2)
         TE_12(1,I) = IVEC(I)
      ENDDO

! **********************************************************************************************************************************
! Calculate unit vector in initial elem. y direction (from KVEC x IVEC):

      CALL CROSS ( KVEC, IVEC, JVEC )
      MAGJ = DSQRT(JVEC(1)*JVEC(1) + JVEC(2)*JVEC(2) + JVEC(3)*JVEC(3))

! If MAGJ =0 then write error and quit.

      IF (MAGJ <= EPS1) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1912) TYPE, EID
            WRITE(F06,1912) TYPE, EID
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1912
            ENDIF
         ENDIF
         RETURN
      ENDIF

      DO I=1,3
         JVEC(I)    = JVEC(I)/MAGJ
         TE_12(2,I) = JVEC(I)
      ENDDO

! **********************************************************************************************************************************
! Perform some geometry checks on the quad element

! First, calculate XEL coords of grids in local element coord system (relative to node 1) with local x along side 1-2.

      XEL(1,1) = ZERO
      XEL(1,2) = ZERO
      XEL(1,3) = ZERO

      IF ((TYPE == 'QUAD8   ') .OR.                                                                                                &
         ((TYPE == 'QUAD4   ') .AND. ((QUAD4TYP == 'MITC4 ') .OR. (QUAD4TYP == 'MITC4+')))) THEN

                                                           ! The z coordinate of grid points in the
                                                           ! XEL element coordinate system can be non-zero if it's warped.
         DO I=2,ELGP
            DO J=1,3
               XEL(I,J) = ZERO
               DO K=1,3
                  XEL(I,J) = XEL(I,J) + (XEB(I,K) - XEB(1,K))*TE_12(J,K)
               ENDDO
            ENDDO
         ENDDO

      ELSE

         XEL(2,3) = ZERO                                   ! All z coords in mean plane are zero by definition of mean plane
         XEL(3,3) = ZERO                                   ! even though using TE to calc them may not yield zero's
         XEL(4,3) = ZERO

         DO I=2,ELGP
            DO J=1,2                                       ! Only calc x znd y coords. z coords in mean plane are zero by definition
               XEL(I,J) = ZERO
               DO K=1,3
                  XEL(I,J) = XEL(I,J) + (XEB(I,K) - XEB(1,K))*TE_12(J,K)
               ENDDO
            ENDDO
         ENDDO

      ENDIF

      QUAD_GEOM_ERR = 0
      CALL QUAD_GEOM_CHECK
      IF (QUAD_GEOM_ERR > 0) THEN
         RETURN
      ENDIF

! **********************************************************************************************************************************
! Now TE_12 is for elem coord system with x along projected side 1-2. We need to rotate x-y (about z) to get x in a
! direction which splits the angle between the two diagonals. QUAD_THETA is the angle between side 1-2 and diagonal from
! point 1 to point 3. QUAD_GAMMA is the angle between side 1-2 and the diagonal from point 2 to point 4. The rotation
! about z is thru an angle of QUAD_DELTA = (QUAD_THETA - QUAD_GAMMA)/2.

! Find QUAD_THETA from the dot product of vector along side 1-2 and  diagonal from point 1 to point 3 (in the mean plane)
! Find QUAD_GAMMA from the dot product of vector along side 1-2 and  diagonal from point 2 to point 4 (in the mean plane).
! Use ABS to get the acute angle

      QUAD_THETA = DACOS(( V13B(1)*IVEC(1) + V13B(2)*IVEC(2) + V13B(3)*IVEC(3))/V13BM)
      QUAD_GAMMA = DACOS((-V24B(1)*IVEC(1) - V24B(2)*IVEC(2) - V24B(3)*IVEC(3))/V24BM)
      QUAD_DELTA = (QUAD_THETA - QUAD_GAMMA)/TWO

      CALL PLANE_COORD_TRANS_21 ( QUAD_DELTA, CT_QD, SUBR_NAME )
      CALL MATMULT_FFF ( CT_QD, TE_12, 3, 3, 3, TE_SD )

! Select how the local x axis is to be

      IF (QUADAXIS == 'SPLITD') THEN

         DO I=1,3
            DO J=1,3
               TE(I,J) = TE_SD(I,J)
            ENDDO
         ENDDO

      ELSE

         DO I=1,3
            DO J=1,3
               TE(I,J) = TE_12(I,J)
            ENDDO
         ENDDO

      ENDIF

! Now TE is final transformation from basic to elem coordinates. That is, UEL = TE*UEB

      IF ((DEBUG(6) == 1) .AND. (WRT_BUG(0) == 1)) THEN

         WRITE(BUG,*) '  Coord transformation matrix that rotates a vector through angle QUAD_DELTA'
         WRITE(BUG,*) '  --------------------------------------------------------------------------'
         DO I=1,3
            WRITE(BUG,90003) (CT_QD(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         WRITE(BUG,*) '  TE matrix if local x is along side 1-2'
         WRITE(BUG,*) '  --------------------------------------'
         DO I=1,3
            WRITE(BUG,90003) (TE_12(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         WRITE(BUG,*) '  TE matrix if local x splits the angle between the diagonals'
         WRITE(BUG,*) '  -----------------------------------------------------------'
         DO I=1,3
            WRITE(BUG,90003) (TE_SD(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         IF (QUADAXIS == 'SPLITD') THEN
            WRITE(BUG,*) '  Final TE matrix - local x splits angle between diagonals:'
            WRITE(BUG,*) '  --------------------------------------------------------'
         ELSE
            WRITE(BUG,*) '  Final TE matrix - local x parallel to side 1-2:'
            WRITE(BUG,*) '  ----------------------------------------------'
         ENDIF
         DO I=1,3
            WRITE(BUG,90003) (TE(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)
         CALL CHECK_TE_MATRIX ( TE, 'TE' )
      ENDIF


! **********************************************************************************************************************************
! Now set TE_IDENT to be 'Y' if TE is an identity matrix. TE will be an identity matrix if the diagonal terms are unity.

      TE_IDENT = 'N'
      DO I=1,3
         IF (DABS(TE(I,I) - ONE) < EPS1) THEN
            ID(I) = 1
         ELSE
            ID(I) = 0
         ENDIF
      ENDDO
      IF ((ID(1) == 1) .AND. (ID(2) == 1) .AND. (ID(3) == 1)) THEN
         TE_IDENT = 'Y'
      ENDIF

! **********************************************************************************************************************************
! Calculate XEL coords of grids in local element coord system (relative to node 1).

      XEL(1,1) = ZERO
      XEL(1,2) = ZERO
      XEL(1,3) = ZERO

      IF ((TYPE == 'QUAD8   ') .OR.                                                                                                &
         ((TYPE == 'QUAD4   ') .AND. ((QUAD4TYP == 'MITC4 ') .OR. (QUAD4TYP == 'MITC4+')))) THEN

                                                           ! The z coordinate of grid points in the
                                                           ! XEL element coordinate system can be non-zero if it's warped.
         DO I=2,ELGP
            DO J=1,3
               XEL(I,J) = ZERO
               DO K=1,3
                  XEL(I,J) = XEL(I,J) + (XEB(I,K) - XEB(1,K))*TE(J,K)
               ENDDO
            ENDDO
         ENDDO

      ELSE

         XEL(2,3) = ZERO                                   ! All z coords in mean plane are zero by definition of mean plane
         XEL(3,3) = ZERO                                   ! even though using TE to calc them may not yield zero's
         XEL(4,3) = ZERO

         DO I=2,ELGP
            DO J=1,2                                       ! Only calc x znd y coords. z coords in mean plane are zero by definition
               XEL(I,J) = ZERO
               DO K=1,3
                  XEL(I,J) = XEL(I,J) + (XEB(I,K) - XEB(1,K))*TE(J,K)
               ENDDO
            ENDDO
         ENDDO

      ENDIF


      IF ((DEBUG(6) == 1) .AND. (WRT_BUG(0) == 1)) THEN
         WRITE(BUG,*) '  Grid coords in mean plane - using final coord system, TE:'
         WRITE(BUG,*) '  --------------------------------------------------------'
         DO I=1,4
            WRITE(BUG,90003) (XEL(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)
      ENDIF

! **********************************************************************************************************************************
! If HBAR is nonzero, we need to calculate transformation from mean plane to the grid points. This is used only for
! the membrane element. BMEANT is the transpose of the B matrix for the QDMEM1 elem.

      IF (DABS(HBAR) > MXWARP) THEN
         CALL CALC_BMEANT
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1908 FORMAT(' *ERROR  1908: ',A,' ELEMENT ',I8,' HAS LENGTH = ZERO ON THE SIDE THAT HAS INTERNAL GRIDS ',I2,' AND ',I2)

 1911 FORMAT(' *ERROR  1911: ',A,' ELEMENT ',I8,' HAS ITS 2 DIAGONALS PARALLEL.')

 1912 FORMAT(' *ERROR  1912: CANNOT CALCULATE VECTOR IN ELEMENT Y DIRECTION FOR ',A,' ELEMENT ',I8)

90003 FORMAT(3(1ES14.6))




! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE QUAD_GEOM_CHECK

! Checks QUAD geometry. The coords of the 4 points have local x from side 1 to side 2

      IMPLICIT NONE

      INTEGER(LONG)                   :: CW_ERR            ! Error indicator for CW/CCW check

      REAL(DOUBLE)                    :: K13VEC(3)         ! A vector resulting from V12 cross V23
      REAL(DOUBLE)                    :: K24VEC(3)         ! A vector resulting from V23 cross V34
      REAL(DOUBLE)                    :: K31VEC(3)         ! A vector resulting from V34 cross V41
      REAL(DOUBLE)                    :: K42VEC(3)         ! A vector resulting from V41 cross V12
      REAL(DOUBLE)                    :: V12L(3)           ! Vector from G.P. 1 to G.P. 2 in basic coords
      REAL(DOUBLE)                    :: V23L(3)           ! Vector from G.P. 2 to G.P. 3 in basic coords
      REAL(DOUBLE)                    :: V34L(3)           ! Vector from G.P. 3 to G.P. 4 in basic coords
      REAL(DOUBLE)                    :: V41L(3)           ! Vector from G.P. 4 to G.P. 1 in basic coords

! **********************************************************************************************************************************
! Variables used in checking geometry

      X12 = -(V12B(1)*IVEC(1) + V12B(2)*IVEC(2) + V12B(3)*IVEC(3))
      X13 = -(V13B(1)*IVEC(1) + V13B(2)*IVEC(2) + V13B(3)*IVEC(3))
      X24 = -(V24B(1)*IVEC(1) + V24B(2)*IVEC(2) + V24B(3)*IVEC(3))
      X14 =  X12 + X24
      X23 =  X13 - X12
      X34 =  X14 - X13
      Y3  =  (V13B(1)*JVEC(1) + V13B(2)*JVEC(2) + V13B(3)*JVEC(3))
      Y4  =  (V24B(1)*JVEC(1) + V24B(2)*JVEC(2) + V24B(3)*JVEC(3))
      Y34 =  Y3 - Y4
      L12 =  DABS(X12)
      L23 =  DSQRT(X23*X23 + Y3*Y3)
      L34 =  DSQRT(X34*X34 + Y34*Y34)
      L41 =  DSQRT(X14*X14 + Y4*Y4)

      X3  = -X13
      X4  = -X14

      IF ((DEBUG(6) == 1) .AND. (WRT_BUG(0) == 1)) THEN
         WRITE(BUG,*) '  Variables used in checking quad geometry:'
         WRITE(BUG,*) '  ----------------------------------------'
         NAME( 1) = ' X12 = -V12B(t)*IVEC       =  '  ;  VAR( 1) = X12
         NAME( 2) = ' X13 = -V13B(t)*IVEC       =  '  ;  VAR( 2) = X13
         NAME( 3) = ' X24 = -V24B(t)*IVEC       =  '  ;  VAR( 3) = X24
         NAME( 4) = ' X14 =  X12 + X24          =  '  ;  VAR( 4) = X14
         NAME( 5) = ' X23 =  X13 - X12          =  '  ;  VAR( 5) = X23
         NAME( 6) = ' X34 =  X14 - X13          =  '  ;  VAR( 6) = X34
         NAME( 7) = ' X3  = -X13                =  '  ;  VAR( 7) = X3
         NAME( 8) = ' Y3  =  V13B(t)*JVEC       =  '  ;  VAR( 8) = Y3
         NAME( 9) = ' X4  = -X14                =  '  ;  VAR( 9) = X4
         NAME(10) = ' Y4  =  V24B(t)*JVEC       =  '  ;  VAR(10) = Y4
         NAME(11) = ' Y34 =  Y3 - Y4            =  '  ;  VAR(11) = Y34
         NAME(12) = ' L12 = Side 1-2 length     =  '  ;  VAR(12) = L12
         NAME(13) = ' L23 = Side 2-3 length     =  '  ;  VAR(13) = L23
         NAME(14) = ' L34 = Side 3-4 length     =  '  ;  VAR(14) = L34
         NAME(15) = ' L41 = Side 4-1 length     =  '  ;  VAR(15) = L41
         DO I=1,15
            WRITE(BUG,90001) NAME(I),VAR(I)
         ENDDO
         WRITE(BUG,*)
      ENDIF

! Make sure that all side lengths are finite

      IF (DABS(L12) < EPS1) THEN
         QUAD_GEOM_ERR = QUAD_GEOM_ERR + 1
         SIDE_GRID1 = 1
         SIDE_GRID2 = 2
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
            WRITE(F06,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1908
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = SIDE_GRID1
               EMG_IFE(NUM_EMG_FATAL_ERRS,3) = SIDE_GRID2
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(L23) < EPS1) THEN
         QUAD_GEOM_ERR = QUAD_GEOM_ERR + 1
         SIDE_GRID1 = 2
         SIDE_GRID2 = 3
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
            WRITE(F06,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1908
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = SIDE_GRID1
               EMG_IFE(NUM_EMG_FATAL_ERRS,3) = SIDE_GRID2
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(L34) < EPS1) THEN
         QUAD_GEOM_ERR = QUAD_GEOM_ERR + 1
         SIDE_GRID1 = 3
         SIDE_GRID2 = 4
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
            WRITE(F06,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1908
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = SIDE_GRID1
               EMG_IFE(NUM_EMG_FATAL_ERRS,3) = SIDE_GRID2
            ENDIF
         ENDIF
      ENDIF

      IF (DABS(L41) < EPS1) THEN
         QUAD_GEOM_ERR = QUAD_GEOM_ERR + 1
         SIDE_GRID1 = 4
         SIDE_GRID2 = 1
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
            WRITE(F06,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1908
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = SIDE_GRID1
               EMG_IFE(NUM_EMG_FATAL_ERRS,3) = SIDE_GRID2
            ENDIF
         ENDIF
      ENDIF

! Check that element is numered clockwise or counter clockwise. This can be skipped if DEBUG(194) = 1 or 3

      IF ((DEBUG(194) == 1) .OR. (DEBUG(194) == 3)) THEN

         DO I=1,3
            V12L(I) = XEL(2,I) - XEL(1,I)
            V23L(I) = XEL(3,I) - XEL(2,I)
            V34L(I) = XEL(4,I) - XEL(3,I)
            V41L(I) = XEL(1,I) - XEL(4,I)
         ENDDO

         CALL CROSS ( V12L, V23L, K13VEC )
         CALL CROSS ( V23L, V34L, K24VEC )
         CALL CROSS ( V34L, V41L, K31VEC )
         CALL CROSS ( V41L, V12L, K42VEC )

         CW_ERR = 0
         IF      ((K13VEC(3) > ZERO) .AND. (K24VEC(3) > ZERO) .AND. (K31VEC(3) > ZERO) .AND. (K42VEC(3) > ZERO)) THEN
            CONTINUE
         ELSE IF ((K13VEC(3) < ZERO) .AND. (K24VEC(3) < ZERO) .AND. (K31VEC(3) < ZERO) .AND. (K42VEC(3) < ZERO)) THEN
            CONTINUE
         ELSE
            CW_ERR = CW_ERR + 1
         ENDIF

         IF (CW_ERR > 0) THEN
            QUAD_GEOM_ERR = QUAD_GEOM_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
            IF (WRT_ERR > 0) THEN
               WRITE(ERR,1910) TYPE, EID
               WRITE(F06,1910) TYPE, EID
            ELSE
               IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
                  ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
                  EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1910
               ENDIF
            ENDIF
         ENDIF

      ENDIF

! Calculate and check diagonal lengths

      DO I=1,3
         V13B(I) = XEB(3,I) - XEB(1,I)
         V24B(I) = XEB(4,I) - XEB(2,I)
      ENDDO
      V13BM = DSQRT(V13B(1)*V13B(1) + V13B(2)*V13B(2) + V13B(3)*V13B(3))
      V24BM = DSQRT(V24B(1)*V24B(1) + V24B(2)*V24B(2) + V24B(3)*V24B(3))
      IF ((DEBUG(6) == 1) .AND. (WRT_BUG(0) == 1)) THEN
         WRITE(BUG,*) '  V13BM, V24BM diagonal lengths (in mean plane):'
         WRITE(BUG,*) '  ----------------------------------------------'
         WRITE(BUG,90003) V13BM,V24BM
         WRITE(BUG,*)
      ENDIF

      IF (V13BM < EPS1) THEN
         QUAD_GEOM_ERR = QUAD_GEOM_ERR + 1
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         DIAG_GRID1 = 1
         DIAG_GRID2 = 3
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1909) TYPE, EID, DIAG_GRID1, DIAG_GRID2
            WRITE(F06,1909) TYPE, EID, DIAG_GRID1, DIAG_GRID2
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1909
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = DIAG_GRID1
               EMG_IFE(NUM_EMG_FATAL_ERRS,3) = DIAG_GRID2
            ENDIF
         ENDIF
      ENDIF

      IF (V24BM < EPS1) THEN
         QUAD_GEOM_ERR = QUAD_GEOM_ERR + 1
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         DIAG_GRID1 = 2
         DIAG_GRID2 = 4
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1909) TYPE, EID, DIAG_GRID1, DIAG_GRID2
            WRITE(F06,1909) TYPE, EID, DIAG_GRID1, DIAG_GRID2
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1909
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = DIAG_GRID1
               EMG_IFE(NUM_EMG_FATAL_ERRS,3) = DIAG_GRID2
            ENDIF
         ENDIF
      ENDIF

! Print warning if all points not in a plane by more than WARP_WARN. Use average diagonal length times EPS4 as measure.

      WARP_WARN = EPS4*(V13BM + V24BM)/TWO
      DHBAR = DABS(HBAR)
      IF (DHBAR > WARP_WARN) THEN
         WARN_ERR = WARN_ERR + 1
         IF ((WRT_ERR > 0) .AND. (WRITE_WARN == 'Y')) THEN
            WRITE(ERR,1915) TYPE, EID, DHBAR, WARP_WARN
            IF (SUPWARN == 'N') THEN
               WRITE(F06,1915) TYPE, EID, DHBAR, WARP_WARN
            ENDIF
         ELSE
            IF (WARN_ERR <= MEWE) THEN
               EMG_IWE(WARN_ERR,1) = 1915
               EMG_RWE(WARN_ERR,1) = DHBAR
               EMG_RWE(WARN_ERR,2) = WARP_WARN
            ENDIF
         ENDIF
      ENDIF

! Check interior angles. This can be skipped if DEBUG(194) = 2 or 3

      IF ((DEBUG(194) == 2) .OR. (DEBUG(194) == 3)) THEN

         IPNT = 0
         IF      (Y4  < ZERO) THEN
            IPNT = 1
         ELSE IF (Y3  < ZERO) THEN
            IPNT = 2
         ELSE IF (X14 < (X12 + X23*Y4/Y3)) THEN
            IPNT = 3
         ELSE IF (X13 > (Y3/Y4)*X14) THEN
            IPNT = 4
         ENDIF
         IF (IPNT > 0) THEN
            QUAD_GEOM_ERR = QUAD_GEOM_ERR + 1
            NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
            FATAL_ERR = FATAL_ERR + 1
            IF (WRT_ERR > 0) THEN
               WRITE(ERR,1914) TYPE, EID, IPNT
               WRITE(F06,1914) TYPE, EID, IPNT
            ELSE
               IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
                  ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
                  EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1914
                  EMG_IFE(NUM_EMG_FATAL_ERRS,2) = IPNT
               ENDIF
            ENDIF
         ENDIF

      ENDIF

      RETURN

! **********************************************************************************************************************************
 1908 FORMAT(' *ERROR  1908: ',A,' ELEMENT ',I8,' HAS LENGTH = ZERO ON THE SIDE THAT HAS GRIDS ',I2,' AND ',I2)

 1909 FORMAT(' *ERROR  1909: ',A,' ELEMENT ',I8,' HAS LENGTH = ZERO FOR THE DIAGONAL BETWEEN INTERNAL GRIDS ',I2,' AND ',I2)

 1910 FORMAT(' *ERROR  1910: ',A,' ELEMENT ',I8,' IS NOT NUMBERED IN A CLOCKWISE OR COUNTER CLOCKWISE MANNER')

 1914 FORMAT(' *ERROR  1914: ',A,' ELEMENT ',I8,' HAS INTERIOR ANGLE >= 180 DEG FOR INTERNAL GRID ',I2)

 1915 FORMAT(' *WARNING    : ',A,' ELEMENT ',I8,' IS NOT PLANAR. GRIDS ARE OUT OF PLANE BY +/-',1ES9.2,'. MAX SUGGESTED IS ',1ES9.2)

90001 FORMAT(3X,A,(1ES19.11))

90003 FORMAT(3(1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE QUAD_GEOM_CHECK

! ##################################################################################################################################

      SUBROUTINE CALC_BMEANT

      IMPLICIT NONE

      REAL(DOUBLE)                    :: BMEAN(12,8)       ! Matrix that transforms the 8 nodal forces for an assumed flat element
!                                                            to the 12 forces for the actual warped element
      REAL(DOUBLE)                    :: DELTA1            ! sin(theta2 - quad_gamma) = sin(TH1 - GAM)
      REAL(DOUBLE)                    :: DELTA2            ! sin(theta2 + quad_gamma) = sin(TH1 + GAM)
      REAL(DOUBLE)                    :: SIN_TH1           ! sin(theta1), TH1 = interior angle of quad at node 1
      REAL(DOUBLE)                    :: COS_TH1           ! cos(theta1)
      REAL(DOUBLE)                    :: SIN_TH2           ! sin(theta2). TH2 = interior angle of quad at node 2
      REAL(DOUBLE)                    :: COS_TH2           ! cos(theta2)
      REAL(DOUBLE)                    :: SIN_GAM           ! sin(quad_gamma ). GAM = angle of side 3-4 from side 1-2
      REAL(DOUBLE)                    :: COS_GAM           ! cos(quad_gamma )
      REAL(DOUBLE)                    :: CTN_TH1           ! cot(theta1)
      REAL(DOUBLE)                    :: CTN_TH2           ! cot(theta2)
      REAL(DOUBLE)                    :: BLK(3,2)          ! One grid-by-node block of BMEAN
      INTEGER(LONG)                   :: IG, JN            ! Grid (row block) and flat node (col block) indices

! **********************************************************************************************************************************
      DO I=1,12
         DO J=1,8
            BMEAN(I,J) = ZERO
         ENDDO
      ENDDO

      SIN_TH1 =  Y4/L41
      COS_TH1 = -X14/L41

      SIN_TH2 =  Y3/L23
      COS_TH2 =  X23/L23

      SIN_GAM = (Y4 - Y3)/L34
      COS_GAM = (X3 - X4)/L34

      CTN_TH1 =  COS_TH1/SIN_TH1
      CTN_TH2 =  COS_TH2/SIN_TH2

      DELTA1 = SIN_TH2*COS_GAM - COS_TH2*SIN_GAM
      DELTA2 = SIN_TH1*COS_GAM + COS_TH1*SIN_GAM

      IF ((DEBUG(6) == 1) .AND. (WRT_BUG(0) == 1)) THEN
         NAME( 1) = 'SIN(THETA1)                     =  '  ;  VAR( 1) = SIN_TH1
         NAME( 2) = 'COS(THETA1)                     =  '  ;  VAR( 2) = COS_TH1
         NAME( 3) = 'SIN(THETA2)                     =  '  ;  VAR( 3) = SIN_TH2
         NAME( 4) = 'COS(THETA2)                     =  '  ;  VAR( 4) = COS_TH2
         NAME( 5) = 'SIN(QUAD_GAMMA)                 =  '  ;  VAR( 5) = SIN_GAM
         NAME( 6) = 'COS(QUAD_GAMMA)                 =  '  ;  VAR( 6) = COS_GAM
         NAME( 7) = 'CTN(THETA1)                     =  '  ;  VAR( 7) = CTN_TH1
         NAME( 8) = 'CTN(THETA2)                     =  '  ;  VAR( 8) = CTN_TH2
         NAME( 9) = 'DELTA1 = SIN(THETA2-QUAD_GAMMA) =  '  ;  VAR( 9) = DELTA1
         NAME(10) = 'DELTA2 = SIN(THETA1+QUAD_GAMMA) =  '  ;  VAR(10) = DELTA2
         WRITE(BUG,*) '  Variables used in creating BMEAN matrix:'
         WRITE(BUG,*) '  ---------------------------------------'
         DO I=1,10
            WRITE(BUG,90001) NAME(I),VAR(I)
         ENDDO
         WRITE(BUG,*)
      ENDIF

      BMEAN( 1,1) =  ONE
      BMEAN( 2,2) =  ONE
      BMEAN( 4,3) =  ONE
      BMEAN( 5,4) =  ONE
      BMEAN( 7,5) =  ONE
      BMEAN( 8,6) =  ONE
      BMEAN(10,7) =  ONE
      BMEAN(11,8) =  ONE
      BMEAN( 3,1) =  HBAR/L12
      BMEAN( 3,2) = -HBAR*(CTN_TH1/L12 - ONE/(L41*SIN_TH1))
      BMEAN( 3,3) = -BMEAN( 3,1)
      BMEAN( 3,4) = -HBAR*CTN_TH2/L12
      BMEAN( 3,7) = -HBAR*SIN_GAM/(L41*DELTA2)
      BMEAN( 3,8) = -HBAR*COS_GAM/(L41*DELTA2)

      BMEAN( 6,1) = -BMEAN( 3,1)
      BMEAN( 6,2) =  HBAR*CTN_TH1/L12
      BMEAN( 6,3) =  BMEAN( 3,1)
      BMEAN( 6,4) =  HBAR*(CTN_TH2/L12 - ONE/(L23*SIN_TH2))
      BMEAN( 6,5) =  HBAR*SIN_GAM/(L23*DELTA1)
      BMEAN( 6,6) =  HBAR*COS_GAM/(L23*DELTA1)

      BMEAN( 9,4) =  HBAR/(L23*SIN_TH2)
      BMEAN( 9,5) = -HBAR*(SIN_GAM/L23 + SIN_TH2/L34)/DELTA1
      BMEAN( 9,6) = -HBAR*(COS_GAM/L23 + COS_TH2/L34)/DELTA1
      BMEAN( 9,7) =  HBAR*SIN_TH1/(L34*DELTA2)
      BMEAN( 9,8) = -HBAR*COS_TH1/(L34*DELTA2)

      BMEAN(12,2) = -HBAR/(L41*SIN_TH1)
      BMEAN(12,5) =  HBAR*SIN_TH2/(L34*DELTA1)
      BMEAN(12,6) =  HBAR*COS_TH2/(L34*DELTA1)
      BMEAN(12,7) = -HBAR*(SIN_TH1/L34 - SIN_GAM/L41)/DELTA2
      BMEAN(12,8) =  HBAR*(COS_TH1/L34 + COS_GAM/L41)/DELTA2

! The terms above are in the axes with x along side 1-2 (X12...Y4 from QUAD_GEOM_CHECK, TE_12). With QUADAXIS = 'SPLITD' the
! membrane matrices BMEANT multiplies are in the final axes TE = CT_QD*TE_12, so rotate BMEAN to them: each 3x2 block becomes
! CT_QD*block*CT_QD(1:2,1:2)'. Without this a warped quad with QUAD_DELTA /= 0 is not rigid: a rigid rotation about an
! in-plane axis gives grid forces in proportion to HBAR.

      IF (QUADAXIS == 'SPLITD') THEN
         DO IG=0,3
            DO JN=0,3
               BLK = BMEAN(3*IG+1:3*IG+3,2*JN+1:2*JN+2)
               BMEAN(3*IG+1:3*IG+3,2*JN+1:2*JN+2) = MATMUL(CT_QD, MATMUL(BLK, TRANSPOSE(CT_QD(1:2,1:2))))
            ENDDO
         ENDDO
      ENDIF

      DO I=1,8
         DO J=1,12
            BMEANT(I,J) = BMEAN(J,I)
         ENDDO
      ENDDO

      RETURN

! **********************************************************************************************************************************
90001 FORMAT(3X,A,(1ES14.6))

91827 FORMAT(6(1ES15.6))

! **********************************************************************************************************************************

      END SUBROUTINE CALC_BMEANT

! **********************************************************************************************************************************

      END SUBROUTINE ELMGM2


      SUBROUTINE ELMGM3 ( WRITE_WARN )

! Calculates and checks elem geometry for 3D elems and provides a transformation matrix ( TE ) to transfer the elem stiffness matrix
! in the elem system to the basic coordinate system. Calculates grid point coords in local coord system.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, BUG, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MEFE, MELGP
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, HALF, ONE, TWO
      USE PARAMS, ONLY                :  EPSIL, HEXAXIS
      USE MODEL_STUF, ONLY            :  EID, ELGP, EMG_IFE, ERR_SUB_NAM, HEXA_DELTA, HEXA_GAMMA, HEXA_THETA,                      &
                                         NUM_EMG_FATAL_ERRS, TE, TE_IDENT, TYPE, XEB, XEL
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE DATE_TIME_UTILS, ONLY       :  OURTIM
      USE VECTOR_GEOMETRY, ONLY       :  CROSS, PLANE_COORD_TRANS_21
      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELMGM3'
      CHARACTER(LEN=*), INTENT(IN)    :: WRITE_WARN        ! If 'Y" write warning messages, otherwise do not

      INTEGER(LONG)                   :: SIDE_GRID1        ! Used for error output purposes
      INTEGER(LONG)                   :: SIDE_GRID2        ! Used for error output purposes
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: ID(3)             ! ID(i) is set to 1 if the i-th diagonal of TE is 1.0


      REAL(DOUBLE)                    :: EPS1              ! A small number to compare to real zero
      REAL(DOUBLE)                    :: IVEC(3)           ! A vector in the elem x dir
      REAL(DOUBLE)                    :: JVEC(3)           ! A vector in the elem y dir
      REAL(DOUBLE)                    :: HEXA_HBAR         ! Warp of HEXA element (only used in calc initial x direction along
!                                                            side 1-2 of the elem projection onto the mean plane)
      REAL(DOUBLE)                    :: KVEC(3)           ! A vector in the elem z dir
      REAL(DOUBLE)                    :: MAGI              ! Magnitude of vector IVEC
      REAL(DOUBLE)                    :: MAGJ              ! Magnitude of vector JVEC
      REAL(DOUBLE)                    :: MAGK              ! Magnitude of vector KVEC
      REAL(DOUBLE)                    :: CT_QD(3,3)        ! Coord transf matrix which will rotate a vector in a coord sys that has
!                                                            x axis along side 1-2 thru angle HEXA_DELTA to get a vector in the
!                                                            element local coord sys
      REAL(DOUBLE)                    :: TE_12(3,3)        ! TE matrix for this elem if local x is parallel to side 1-2
      REAL(DOUBLE)                    :: TE_SD(3,3)        ! TE matrix for this elem if local x splits angle between the two diags
      REAL(DOUBLE)                    :: XEBM(4,3)         ! Basic coordinates of 4 points that lie midway between the 4 grids
!                                                            1-4 and the 4 grids 5-8 of the element. Used to construct a mean plane.
!                                                            These will be referred to as points AM, BM, CM, DM
      REAL(DOUBLE)                    :: V12B(3)           ! Vector from G.P. 1 to G.P. 2 in basic coords
      REAL(DOUBLE)                    :: V13B(3)           ! Vector from G.P. 1 to G.P. 3 in basic coords (a diagonal)
      REAL(DOUBLE)                    :: V24B(3)           ! Vector from G.P. 2 to G.P. 4 in basic coords (a diagonal)
      REAL(DOUBLE)                    :: V13BM             ! Mag of V13B
      REAL(DOUBLE)                    :: V24BM             ! Mag of V24B

      INTRINSIC                       :: DABS



! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

! Initialize XEL to zero

      DO I=1,MELGP
         DO J=1,3
            XEL(I,J) = ZERO
         ENDDO
      ENDDO

! Calc coords of a mean plane between grids 1-4 and grids 5-8. These will be referred to as points AM, BM, CM, DM:

      DO J=1,3
         XEBM(1,J) = (XEB(1,J) + XEB(5,J))/TWO
         XEBM(2,J) = (XEB(2,J) + XEB(6,J))/TWO
         XEBM(3,J) = (XEB(3,J) + XEB(7,J))/TWO
         XEBM(4,J) = (XEB(4,J) + XEB(8,J))/TWO
      ENDDO

! **********************************************************************************************************************************
! Calculate elem z direction from cross products of diagonals

! Generate vectors from G.P AM to G.P CM and from G.P. BM to G.P. DM (diagonals in the mean (M) plane)

      DO I=1,3
         V13B(I) = XEBM(3,I) - XEBM(1,I)
         V24B(I) = XEBM(4,I) - XEBM(2,I)
      ENDDO

      V13BM = DSQRT(V13B(1)*V13B(1) + V13B(2)*V13B(2) + V13B(3)*V13B(3))
      V24BM = DSQRT(V24B(1)*V24B(1) + V24B(2)*V24B(2) + V24B(3)*V24B(3))

      CALL CROSS ( V13B, V24B, KVEC )
      MAGK = DSQRT(KVEC(1)*KVEC(1) + KVEC(2)*KVEC(2) + KVEC(3)*KVEC(3))

! If MAGK = 0 then diagonals are parallel so write error and quit

      IF (MAGK <=  EPS1) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1911) TYPE, EID
            WRITE(F06,1911) TYPE, EID
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1911
            ENDIF
         ENDIF
         RETURN
      ENDIF

! Unit vector in elem local z direction is 3rd row of TE

      DO I=1,3
         KVEC(I)    = KVEC(I)/MAGK
         TE_12(3,I) = KVEC(I)
      ENDDO

! **********************************************************************************************************************************
! Calc initial elem x dir along side 1-2 of the element. First, get vector from pt 1 to 2:

      DO I=1,3
         V12B(I) = XEB(2,I) - XEB(1,I)
      ENDDO

! HEXA_HBAR is one half of the projection of V12B in z direction

      HEXA_HBAR = HALF*((V12B(1)*KVEC(1) + V12B(2)*KVEC(2) + V12B(3)*KVEC(3)))

! Now calculate initial x direction along side 1-2 of the elem projection onto the mean plane.

      DO I=1,3
         IVEC(I) = V12B(I) - TWO*HEXA_HBAR*KVEC(I)
      ENDDO

! If initial MAGI = 0 then write error and quit.

      MAGI  = DSQRT(IVEC(1)*IVEC(1) + IVEC(2)*IVEC(2) + IVEC(3)*IVEC(3))
      IF (MAGI <= EPS1) THEN
         SIDE_GRID1 = 1
         SIDE_GRID2 = 2
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
            WRITE(F06,1908) TYPE, EID, SIDE_GRID1, SIDE_GRID2
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1908
               EMG_IFE(NUM_EMG_FATAL_ERRS,2) = SIDE_GRID1
               EMG_IFE(NUM_EMG_FATAL_ERRS,3) = SIDE_GRID2
            ENDIF
         ENDIF
         RETURN
      ENDIF

! Unit vector in initial elem x direction

      DO I=1,3
         IVEC(I)     = IVEC(I)/MAGI                        ! Unit vector along side 1-2 in the mean plane (NOT from G.P. 1-2)
         TE_12(1,I) = IVEC(I)
      ENDDO

! **********************************************************************************************************************************
! Calculate unit vector in initial elem. y dir. (from VZ (cross) VX):

      CALL CROSS ( KVEC, IVEC, JVEC )
      MAGJ = DSQRT(JVEC(1)*JVEC(1) + JVEC(2)*JVEC(2) + JVEC(3)*JVEC(3))

! If MAGJ =0 then write error and quit.

      IF (MAGJ <= EPS1) THEN
         NUM_EMG_FATAL_ERRS = NUM_EMG_FATAL_ERRS + 1
         FATAL_ERR = FATAL_ERR + 1
         IF (WRT_ERR > 0) THEN
            WRITE(ERR,1912) TYPE, EID
            WRITE(F06,1912) TYPE, EID
         ELSE
            IF (NUM_EMG_FATAL_ERRS <= MEFE) THEN
               ERR_SUB_NAM(NUM_EMG_FATAL_ERRS) = SUBR_NAME
               EMG_IFE(NUM_EMG_FATAL_ERRS,1) = 1912
            ENDIF
         ENDIF
         RETURN
      ENDIF

      DO I=1,3
         JVEC(I)     = JVEC(I)/MAGJ
         TE_12(2,I) = JVEC(I)
      ENDDO

! **********************************************************************************************************************************
! Now TE_12 is for elem coord system with x along projected side 1-2. We need to rotate x-y (about z) to get x in a
! direction which splits the angle between the two diagonals. HEXA_THETA is the angle between side 1-2 and diagonal from
! point 1 to point 3. HEXA_GAMMA is the angle between side 1-2 and the diagonal from point 2 to point 4. The rotation
! about z is thru an angle of HEXA_DELTA = (HEXA_THETA - HEXA_GAMMA)/2.

! Find HEXA_THETA from the dot product of vector along side 1-2 and  diagonal from point 1 to point 3 (in the mean plane)
! Find HEXA_GAMMA from the dot product of vector along side 1-2 and  diagonal from point 2 to point 4 (in the mean plane).
! Use ABS to get the acute angle

      HEXA_THETA = DACOS(( V13B(1)*IVEC(1) + V13B(2)*IVEC(2) + V13B(3)*IVEC(3))/V13BM)
      HEXA_GAMMA = DACOS((-V24B(1)*IVEC(1) - V24B(2)*IVEC(2) - V24B(3)*IVEC(3))/V24BM)
      HEXA_DELTA = (HEXA_THETA - HEXA_GAMMA)/TWO

      CALL PLANE_COORD_TRANS_21 ( HEXA_DELTA, CT_QD, SUBR_NAME )
      CALL MATMULT_FFF ( CT_QD, TE_12, 3, 3, 3, TE_SD )

! Select how the local x axis is to be

      IF (HEXAXIS == 'SPLITD') THEN

         DO I=1,3
            DO J=1,3
               TE(I,J) = TE_SD(I,J)
            ENDDO
         ENDDO

      ELSE

         DO I=1,3
            DO J=1,3
               TE(I,J) = TE_12(I,J)
            ENDDO
         ENDDO

      ENDIF

! Now TE is final transformation from basic to elem coordinates. That is, UEL = TE*UEB

      IF ((DEBUG(6) == 2) .AND. (WRT_BUG(0) == 1)) THEN

         WRITE(BUG,*) '  Coord transformation matrix that rotates a vector through angle HEXA_DELTA'
         WRITE(BUG,*) '  --------------------------------------------------------------------------'
         DO I=1,3
            WRITE(BUG,90003) (CT_QD(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         WRITE(BUG,*) '  TE matrix if local x is along side 1-2'
         WRITE(BUG,*) '  --------------------------------------'
         DO I=1,3
            WRITE(BUG,90003) (TE_12(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         WRITE(BUG,*) '  TE matrix if local x splits the angle between the diagonals'
         WRITE(BUG,*) '  -----------------------------------------------------------'
         DO I=1,3
            WRITE(BUG,90003) (TE_SD(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)

         IF (HEXAXIS == 'SPLITD') THEN
            WRITE(BUG,*) '  Final TE matrix - local x splits angle between diagonals:'
            WRITE(BUG,*) '  --------------------------------------------------------'
         ELSE
            WRITE(BUG,*) '  Final TE matrix - local x parallel to side 1-2:'
            WRITE(BUG,*) '  ----------------------------------------------'
         ENDIF
         DO I=1,3
            WRITE(BUG,90003) (TE(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)
         CALL CHECK_TE_MATRIX ( TE, 'TE' )
      ENDIF

! **********************************************************************************************************************************
! Now set TE_IDENT to be 'Y' if TE is an identity matrix. TE will be an identity matrix if the diagonal terms are unity.

      TE_IDENT = 'N'
      DO I=1,3
         IF (DABS(TE(I,I) - ONE) < EPS1) THEN
            ID(I) = 1
         ELSE
            ID(I) = 0
         ENDIF
      ENDDO
      IF ((ID(1) == 1) .AND. (ID(2) == 1) .AND. (ID(3) == 1)) THEN
         TE_IDENT = 'Y'
      ENDIF

! **********************************************************************************************************************************
! Calculate XEL coords of grids in local element coord system (relative to node 1).

      XEL(1,1) = ZERO
      XEL(1,2) = ZERO
      XEL(1,3) = ZERO

      DO I=2,ELGP
         DO J=1,3
            XEL(I,J) = ZERO
            DO K=1,3
               XEL(I,J) = XEL(I,J) + (XEB(I,K) - XEB(1,K))*TE(J,K)
            ENDDO
         ENDDO
      ENDDO

      IF ((DEBUG(6) == 1) .AND. (WRT_BUG(0) == 1)) THEN
         WRITE(BUG,*) '  Grid coords in mean plane - using final coord system, TE:'
         WRITE(BUG,*) '  --------------------------------------------------------'
         DO I=1,4
            WRITE(BUG,90003) (XEL(I,J),J=1,3)
         ENDDO
         WRITE(BUG,*)
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1908 FORMAT(' *ERROR  1908: ',A,' ELEMENT ',I8,' HAS LENGTH = ZERO ON THE SIDE THAT HAS INTERNAL GRIDS ',I2,' AND ',I2)

 1911 FORMAT(' *ERROR  1911: ELEMENT ',I8,' TYPE ',A,' HAS DIAG FROM INTERNAL GRIDS ',I2,' TO ',I2,                                &
                           ' PARALLEL TO DIAG FROM INTERNAL GRIDS ',I2,' TO ',I2)

 1912 FORMAT(' *ERROR  1912: CANNOT CALCULATE VECTOR IN ELEMENT Y DIRECTION FOR ELEMENT ',I8,' TYPE ',A)

 1913 FORMAT(' *ERROR  1913: BAD GEOMETRY FOR ELEMENT ',I8,' TYPE ',A,'. ELEMENT EITHER NOT NUMBERED CLOCKWISE/COUNTER CLOCKWISE.' &
                          ,' OR AN INTERNAL ANGLE IS > 180 DEG')

 1939 FORMAT(' *ERROR  1939: INTERNAL GRIDS 1 AND 2 ON ELEMENT ',I8,' TYPE ',A,' ARE COINCIDENT')

90003 FORMAT(3(1ES14.6))

! **********************************************************************************************************************************

      END SUBROUTINE ELMGM3


      SUBROUTINE CHECK_TE_MATRIX ( TE_IN, NAME_IN )

! Checks TE_IN(t)*TE_IN to see if the element transformation matrix times its transpose is an identity matrix

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE IOUNT1, ONLY                :  BUG

      USE FULL_MATRIX_ALGEBRA, ONLY   :  MATMULT_FFF_T

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: NAME_IN           ! Name for output purposes

      INTEGER(LONG)                   :: I,J               ! DO loop indices

      REAL(DOUBLE), INTENT(IN)        :: TE_IN(3,3)        ! Input TE matrix
      REAL(DOUBLE)                    :: DUM1(3,3)         ! Dummy matrix in a matrix multiply
      REAL(DOUBLE)                    :: DUM2(3,3)         ! Dummy matrix in a matrix multiply

! **********************************************************************************************************************************
      DO I=1,3
         DO J=1,3
            DUM1(I,J) = TE_IN(I,J)
         ENDDO
      ENDDO

      CALL MATMULT_FFF_T ( DUM1, TE_IN, 3, 3, 3, DUM2 )
      WRITE(BUG,1) NAME_IN, NAME_IN
      DO I=1,3
         WRITE(BUG,2) (DUM2(I,J),J=1,3)
      ENDDO
      WRITE(BUG,*)

      RETURN

! **********************************************************************************************************************************
    1 FORMAT('  Check on ',A,'(t)*',A,' to see if it is the identity matrix:',/,                                                   &
             '  -----------------------------------------------------------------')

    2 FORMAT(3(1ES14.6))

      END SUBROUTINE CHECK_TE_MATRIX

   END MODULE ELEMENT_GEOMETRY_PREPARATION
