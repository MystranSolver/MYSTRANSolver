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

   MODULE ELEMENT_MODEL_INDEXING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ELESORT, ELEM_PROP_MATL_IIDS, ELSAVE

   CONTAINS

      SUBROUTINE ELESORT

! Performs 3 functions:

!   (1) Generates ESORT1 array (actual elem ID's) and ESORT2 array (internal elem no.) for all elements.

!   (2) Sorts ESORT1, ESORT2, EPNT, ETYPE, EOFF so that elem ID's in ESORT1 are in numerically increasing order.

!   (3) Sorts RIGID_ELEM_IDS so that the actual ID's of rigid elements are in numerically increasing order

!   (4) Checks for redundant element ID's - elastic as well as rigid elements

! At the beginning of the subroutine, ESORT1(I) is set to element ID's in the order in which they were encountered
! in the Bulk Data Deck and ESORT2(I) = I. EPNT(I) gives the location in EDAT where connection data starts for element
! ESORT1(I) and ETYPE(I) is the element type for element ESORT1(I).

! This subr then sorts ESORT1, ESORT2, EPNT and ETYPE (together) so that the actual element numbers in ESORT1 are in
! numerical order. Then ESORT2(I) is the position, in the Bulk Data Deck (BDD), where actual elem ESORT1(I) was located.

! For example:

!    EID's in    |    At Beginning of Subr  |   At End of Subroutine:
!     Order      |        (after (1))       |
!    as input    |                          |
!    in BDD    I | ESORT1 ESORT2 EPNT ETYPE | ESORT1 ESORT2 EPNT ETYPE
!   -------------|--------------------------|-------------------------
!       31     1 |     31      1    1    B1 |     11      2    9    E1
!       11     2 |     11      2    9    E1 |     21      4   19    T2
!       41     3 |     41      3   15    R1 |     31      1    1    B1
!       21     4 |     21      4   19    T2 |     41      3   15    R1
!       61     5 |     61      5   24    Q2 |     51      6   30    B1
!       51     6 |     51      6   30    B1 |     61      5   24    Q2

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, ELESORT_RUN, NELE, NRIGEL
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EDAT, EOFF, EPNT, ESORT1, ESORT2, ETYPE, RIGID_ELEM_IDS
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG

      USE SORTING, ONLY               :  SORT_INT1, SORT_INT3_CHAR2
      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELESORT'

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: IERROR            ! Error count




! **********************************************************************************************************************************
! (1) Generate ESORT1 and ESORT2. Initially, set ESORT1(I) = elem ID's in order read in Bulk Data and ESORT2(I) = I.

      DO I=1,NELE
         ESORT1(I) = EDAT(EPNT(I))
         ESORT2(I) = I
      ENDDO

      IF (DEBUG(7) == 1) THEN                              ! Debug output before sorting ESORT1, ESORT2, EPNT, ETYPE
         WRITE(F06,1001)
         WRITE(F06,1011)
         DO I=1,NELE
            WRITE(F06,1021) I,ESORT1(I),ESORT2(I),EPNT(I),ETYPE(I)
         ENDDO
         WRITE(F06,*)
      ENDIF

! (2) Sort ESORT1, ESORT2, EPNT, ETYPE so that ESORT1 has actual element ID's in numerically increasing order

      IF (NELE > 1) THEN
         CALL SORT_INT3_CHAR2 ( SUBR_NAME, 'ESORT1, ESORT2, EPNT, ETYPE, EOFF', NELE, ESORT1, ESORT2, EPNT, ETYPE, EOFF )
      ENDIF

      IF (DEBUG(7) == 1) THEN                              ! Debug output after sorting ESORT1, ESORT2, EPNT, ETYPE
         WRITE(F06,1002)
         WRITE(F06,1011)
         DO I=1,NELE
            WRITE(F06,1021) I,ESORT1(I),ESORT2(I),EPNT(I),ETYPE(I)
         ENDDO
         WRITE(F06,*)
      ENDIF

! (3) Sort RIGID_ELEM_IDS

      IF (NRIGEL > 1) THEN
         CALL SORT_INT1 ( SUBR_NAME, 'RIGID_ELEM_IDS', NRIGEL, RIGID_ELEM_IDS )
      ENDIF

! (4) Check for duplicate element numbers

      IERROR = 0

      DO I=1,NELE-1                                        ! Elastic elements
         IF (ESORT1(I) == ESORT1(I+1)) THEN
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1400) ESORT1(I+1)
            WRITE(F06,1400) ESORT1(I+1)
         ENDIF
      ENDDO

      DO I=1,NRIGEL-1                                      ! Rigid elements
         IF (RIGID_ELEM_IDS(I) == RIGID_ELEM_IDS(I+1)) THEN
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1400) RIGID_ELEM_IDS(I+1)
            WRITE(F06,1400) RIGID_ELEM_IDS(I+1)
         ENDIF
      ENDDO

      IF (IERROR > 0) THEN
         WRITE(ERR,1408) IERROR
         WRITE(F06,1408) IERROR
         CALL OUTA_HERE ( 'Y' )                                    ! Duplicate elem numbers, so quit
      ENDIF

! Set ELESORT_RUN so subrs which need to know if ELESORT subr has run, will know it has run

      ELESORT_RUN = 'Y'



      RETURN

! **********************************************************************************************************************************
 1001 FORMAT('             Element arrays before sorting in subroutine ELESORT:')

 1002 FORMAT('             Element arrays after sorting in subroutine ELESORT:')

 1011 FORMAT('            I   ESORT1(I)   ESORT2(I)     EPNT(I)       ETYPE(I)')

 1021 FORMAT(1X,4I12,10X,A)

 1400 FORMAT(' *ERROR  1400: ELEMENT NUMBER ',I8,' IS A DUPLICATE ELEMENT NUMBER.')

 1408 FORMAT(' PROCESSING TERMINATED DUE TO ABOVE ',I8,' ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE ELESORT


      SUBROUTINE ELEM_PROP_MATL_IIDS

! Fix up routine for element property and material ID's. When the Bulk Data File was read, element and material property ID's
! were actual ID values. WE need to convert these to consequetive integer ID's so that the element property and material arrays
! can be accessed sequentially. This subr performs that function.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, IN4FIL_NUM, NUM_IN4_FILES
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, DEDAT_Q4_SHELL_KEY, DEDAT_T3_SHELL_KEY, DEDAT_Q8_SHELL_KEY, FATAL_ERR,      &
                                         MPCOMP0, MPCOMP_PLIES, NCMASS, NELE, NMATL, NPBAR, NPBEAM,                                &
                                         NPBUSH, NPCOMP, NPELAS, NPMASS, NPROD, npshear, NPSHEL, NPSOLID, NPUSER1, NPUSERIN
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  CMASS, ETYPE, EPNT, EDAT, PELAS, PROD, PBAR, PBEAM, PBUSH, PCOMP, PMASS, PSHEAR,          &
                                         PSHEL, PSOLID, PUSER1, PUSERIN, MATL

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME        = 'ELEM_PROP_MATL_IIDS'
      CHARACTER( 6*BYTE)              :: CMASS_TYPE        ! 'CMASS1', 'CMASS2' 'CMASS3' or 'CMASS4'
      CHARACTER( 1*BYTE)              :: FOUND             ! Used to indicate if a prop or mat'l ID was found
      CHARACTER( 1*BYTE)              :: FOUND_PCOMP       ! Used to indicate if a PCOMP prop ID was found
      CHARACTER( 1*BYTE)              :: FOUND_PSHEL       ! Used to indicate if a PSHELL prop ID was found
      CHARACTER( 8*BYTE)              :: NAME = 'MATERIAL' ! Used for output error message
      CHARACTER( 8*BYTE)              :: PROPERTY_NAME     ! Name of an element property card

      INTEGER(LONG)                   :: EPNTK             ! Value from array EPNT at the row for this internal elem ID. It is the
!                                                            row number in array EDAT where data begins for this element.

      INTEGER(LONG)                   :: I,J,K             ! DO loop indices
      INTEGER(LONG)                   :: IERROR            ! Cumulative error count
      INTEGER(LONG)                   :: IN4_ID            !
      INTEGER(LONG)                   :: MATERIAL_ID       ! Material ID
      INTEGER(LONG)                   :: PCOMP_INDEX       ! Index into PCOMP array
      INTEGER(LONG)                   :: PCOMP_PLIES       ! Number of plies in PCOMP array
      INTEGER(LONG)                   :: PROPERTY_ID       ! Property ID




! **********************************************************************************************************************************
      PROPERTY_NAME = '        '
      IERROR = 0

! Process PID'S on element connection cards to internal values

      DO I = 1,NELE
         EPNTK = EPNT(I)

         FOUND = 'N'

         IF      (ETYPE(I) == 'BAR     ') THEN             ! Process property ID's for BAR elements
            PROPERTY_NAME(1:4) = 'PBAR'
            PROPERTY_ID = EDAT(EPNTK+1)
            DO J = 1,NPBAR
               IF (PBAR(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

         ELSE IF (ETYPE(I) == 'BEAM    ') THEN             ! Process property ID's for BEAM elements
            PROPERTY_NAME(1:4) = 'PBEAM'
            PROPERTY_ID = EDAT(EPNTK+1)
            DO J = 1,NPBEAM
               IF (PBEAM(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

         ELSE IF (ETYPE(I)(1:4) == 'BUSH') THEN            ! Process property ID's for BUSH elements
            PROPERTY_NAME(1:5) = 'PBUSH'
            PROPERTY_ID = EDAT(EPNTK+1)
            DO J = 1,NPBUSH
               IF (PBUSH(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

         ELSE IF (ETYPE(I)(1:4) == 'ELAS') THEN            ! Process property ID's for ELAS elements
            PROPERTY_NAME(1:5) = 'PELAS'
            PROPERTY_ID = EDAT(EPNTK+1)
            DO J = 1,NPELAS
               IF (PELAS(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
                                                           ! Process property ID's for solid elements
         ELSE IF ((ETYPE(I) == 'HEXA8   ') .OR. (ETYPE(I) == 'HEXA20  ') .OR.                                                      &
                  (ETYPE(I) == 'PENTA6  ') .OR. (ETYPE(I) == 'PENTA15 ') .OR.                                                      &
                  (ETYPE(I) == 'TETRA4  ') .OR. (ETYPE(I) == 'TETRA10 ')) THEN
            PROPERTY_NAME(1:6) = 'PSOLID'
            PROPERTY_ID = EDAT(EPNTK+1)
            DO J = 1,NPSOLID
               IF (PSOLID(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF
                                                           ! Process property ID's for shell elements
         ELSE IF ((ETYPE(I)(1:5) == 'QUAD4') .OR. (ETYPE(I)(1:5) == 'TRIA3') .OR. (ETYPE(I)(1:5) == 'QUAD8')) THEN
            PROPERTY_ID = EDAT(EPNTK+1)
            FOUND_PSHEL = 'N'
            DO J = 1,NPSHEL                                ! Search PSHEL to see if we find ID
               IF (PSHEL(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND_PSHEL = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND_PSHEL == 'Y') THEN
               PROPERTY_NAME(1:6) = 'PSHELL'
               IF      (ETYPE(I)(1:5) == 'QUAD4') THEN     ! Set flag to indicate, in EDAT, that this element refers to a PSHELL
                  EDAT(EPNTK+DEDAT_Q4_SHELL_KEY) = 1       !     flag for QUAD's using PSHELL is 1
               ELSE IF (ETYPE(I)(1:5) == 'TRIA3') THEN
                  EDAT(EPNTK+DEDAT_T3_SHELL_KEY) = 1       !     flag for TRIA's using PSHELL is 1
               ELSE IF (ETYPE(I)(1:5) == 'QUAD8') THEN
                  EDAT(EPNTK+DEDAT_Q8_SHELL_KEY) = 1       !     flag for QUAD8's using PSHELL is 1
               ENDIF
            ENDIF

            FOUND_PCOMP = 'N'
            DO J = 1,NPCOMP                                ! Search PCOMP to see if we find ID
               IF (PCOMP(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND_PCOMP = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND_PCOMP == 'Y') THEN
               PROPERTY_NAME(1:5) = 'PCOMP'
               IF      (ETYPE(I)(1:5) == 'QUAD4') THEN     ! Set flag to indicate, in EDAT, that this element refers to a PSHELL
                  EDAT(EPNTK+DEDAT_Q4_SHELL_KEY) = 2       !     flag for QUAD's using PCOMP is 2
               ELSE IF (ETYPE(I)(1:5) == 'TRIA3') THEN
                  EDAT(EPNTK+DEDAT_T3_SHELL_KEY) = 2       !     flag for TRIA's using PCOMP  is 2
               ELSE IF (ETYPE(I)(1:5) == 'QUAD8') THEN
                  EDAT(EPNTK+DEDAT_Q8_SHELL_KEY) = 2       !     flag for QUAD8's using PCOMP  is 2
               ENDIF
            ENDIF

            IF ((FOUND_PSHEL == 'N') .AND. (FOUND_PCOMP == 'N')) THEN
               WRITE(ERR,1404) 'PSHELL/PCOMP', PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) 'PSHELL/PCOMP', PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

            IF ((FOUND_PSHEL == 'Y') .AND. (FOUND_PCOMP == 'Y')) THEN
               WRITE(ERR,1406) PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1406) PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

         ELSE IF (ETYPE(I) == 'SHEAR   ') THEN             ! Process property ID's for SHEAR elements
            PROPERTY_NAME(1:6) = 'PSHEAR'
            PROPERTY_ID = EDAT(EPNTK+1)
            DO J = 1,NPSHEAR
               IF (PSHEAR(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

         ELSE IF (ETYPE(I) == 'ROD     ') THEN             ! Process property ID's for ROD elements
            PROPERTY_NAME(1:4) = 'PROD'
            PROPERTY_ID = EDAT(EPNTK+1)
            DO J = 1,NPROD
               IF (PROD(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

         ELSE IF (ETYPE(I) == 'USER1   ') THEN             ! Process property ID's for USER1 elements
            PROPERTY_NAME(1:6) = 'PUSER1'
            PROPERTY_ID = EDAT(EPNTK+1)
            DO J = 1,NPUSER1
               IF (PUSER1(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

         ELSE IF (ETYPE(I) == 'USERIN  ') THEN             ! Process property ID's for USERIN elements
            PROPERTY_NAME(1:7) = 'PUSERIN'
            PROPERTY_ID = EDAT(EPNTK+1)
            DO J = 1,NPUSERIN
               IF (PUSERIN(J,1) == PROPERTY_ID) THEN
                  EDAT(EPNTK+1) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, ETYPE(I), EDAT(EPNTK)
               IERROR = IERROR + 1
               FATAL_ERR = FATAL_ERR + 1
            ENDIF

         ELSE IF (ETYPE(I) == 'PLOTEL  ') THEN
            CONTINUE

         ELSE

            WRITE(ERR,1403) SUBR_NAME,ETYPE(I)
            WRITE(F06,1403) SUBR_NAME,ETYPE(I)
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )                         ! Coding error (elem type not valid), so quit

         ENDIF

      ENDDO

! **********************************************************************************************************************************
! Process MID'S on property cards to internal pointers

      PROPERTY_NAME(1:4) = 'PBAR'                          ! Process material ID's on PBAR
      DO I = 1,NPBAR
         FOUND = 'N'
         MATERIAL_ID = PBAR(I,2)
         DO J = 1,NMATL
            IF (MATL(J,1) == MATERIAL_ID) THEN
               PBAR(I,2) = J
               FOUND = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (FOUND == 'N') THEN
            WRITE(ERR,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PBAR(I,1)
            WRITE(F06,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PBAR(I,1)
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO

      PROPERTY_NAME(1:4) = 'PBEAM'                         ! Process material ID's on PBEAM
      DO I = 1,NPBEAM
         FOUND = 'N'
         MATERIAL_ID = PBEAM(I,2)
         DO J = 1,NMATL
            IF (MATL(J,1) == MATERIAL_ID) THEN
               PBEAM(I,2) = J
               FOUND = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (FOUND == 'N') THEN
            WRITE(ERR,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PBEAM(I,1)
            WRITE(F06,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PBEAM(I,1)
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO

      PROPERTY_NAME(1:4) = 'PSHEAR'                          ! Process material ID's on PSHEAR
      DO I = 1,NPSHEAR
         FOUND = 'N'
         MATERIAL_ID = PSHEAR(I,2)
         DO J = 1,NMATL
            IF (MATL(J,1) == MATERIAL_ID) THEN
               PSHEAR(I,2) = J
               FOUND = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (FOUND == 'N') THEN
            WRITE(ERR,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PSHEAR(I,1)
            WRITE(F06,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PSHEAR(I,1)
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO

      PROPERTY_NAME(1:4) = 'PROD'                          ! Process material ID's on PROD
      DO I = 1,NPROD
         FOUND = 'N'
         MATERIAL_ID = PROD(I,2)
         DO J = 1,NMATL
            IF (MATL(J,1) == MATERIAL_ID) THEN
               PROD(I,2) = J
               FOUND = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (FOUND == 'N') THEN
            WRITE(ERR,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PROD(I,1)
            WRITE(F06,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PROD(I,1)
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO

     PROPERTY_NAME(1:6) = 'PSHELL'                         ! Process material ID's on PSHELL
      DO I = 1,NPSHEL
         DO K = 2,5
            FOUND = 'N'
            MATERIAL_ID = PSHEL(I,K)
            IF (MATERIAL_ID == 0) CYCLE
            DO J = 1,NMATL
               IF (MATL(J,1) == MATERIAL_ID) THEN
                  PSHEL(I,K) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PSHEL(I,1)
               WRITE(F06,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PSHEL(I,1)
               FATAL_ERR = FATAL_ERR + 1
               IERROR = IERROR + 1
            ELSE
               CYCLE
            ENDIF
         ENDDO
      ENDDO

      PROPERTY_NAME(1:5) = 'PCOMP'                          ! Process material ID's on PCOMP
      DO I = 1,NPCOMP
         PCOMP_PLIES = PCOMP(I,5)
         DO K = 1,PCOMP_PLIES
            PCOMP_INDEX = MPCOMP0 + MPCOMP_PLIES*(K-1) + 1
            FOUND = 'N'
            MATERIAL_ID = PCOMP(I,PCOMP_INDEX)
            DO J = 1,NMATL
               IF (MATL(J,1) == MATERIAL_ID) THEN
                  PCOMP(I,PCOMP_INDEX) = J
                  FOUND = 'Y'
                  EXIT
               ENDIF
            ENDDO
            IF (FOUND == 'N') THEN
               WRITE(ERR,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PCOMP(I,1)
               WRITE(F06,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PCOMP(I,1)
               FATAL_ERR = FATAL_ERR + 1
               IERROR = IERROR + 1
            ELSE
               CYCLE
            ENDIF
         ENDDO
      ENDDO

      PROPERTY_NAME(1:6) = 'PSOLID'                        ! Process material ID's on PSOLID
      DO I = 1,NPSOLID
         FOUND = 'N'
         MATERIAL_ID = PSOLID(I,2)
            DO J = 1,NMATL
            IF (MATL(J,1) == MATERIAL_ID) THEN
               PSOLID(I,2) = J
               FOUND = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (FOUND == 'N') THEN
            WRITE(ERR,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PSOLID(I,1)
            WRITE(F06,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PSOLID(I,1)
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO

      PROPERTY_NAME(1:6) = 'PUSER1'                        ! Process material ID's on PUSER1
      DO I = 1,NPUSER1
         FOUND = 'N'
         MATERIAL_ID = PUSER1(I,2)
         DO J = 1,NMATL
            IF (MATL(J,1) == MATERIAL_ID) THEN
               PUSER1(I,2) = J
               FOUND = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (FOUND == 'N') THEN
            WRITE(ERR,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PUSER1(I,1)
            WRITE(F06,1401) NAME, MATERIAL_ID, PROPERTY_NAME, PUSER1(I,1)
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO

      PROPERTY_NAME(1:7) = 'PUSERIN'                       ! Process material ID's on PUSERIN (actually IN4 ID)
      DO I = 1,NPUSERIN
         FOUND = 'N'
         IN4_ID = PUSERIN(I,2)
         DO J = 1,NUM_IN4_FILES
            IF (IN4FIL_NUM(J) == IN4_ID) THEN
               PUSERIN(I,2) = J
               FOUND = 'Y'
               EXIT
            ENDIF
         ENDDO
         IF (FOUND == 'N') THEN
            WRITE(ERR,1401) 'IN4FIL_NUM', IN4_ID, PROPERTY_NAME, PUSERIN(I,1)
            WRITE(F06,1401) 'IN4FIL_NUM', IN4_ID, PROPERTY_NAME, PUSERIN(I,1)
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF
      ENDDO

! **********************************************************************************************************************************
! Check that all PMASS ID's have been defined

      PROPERTY_NAME(1:5) = 'PMASS'
      DO I=1,NCMASS

         PROPERTY_ID = CMASS(I,3)

         IF      (CMASS(I,2) == 1) THEN
            CMASS_TYPE = 'CMASS1'
         ELSE IF (CMASS(I,2) == 2) THEN
            CMASS_TYPE = 'CMASS2'
         ELSE IF (CMASS(I,2) == 3) THEN
            CMASS_TYPE = 'CMASS3'
         ELSE IF (CMASS(I,2) == 4) THEN
            CMASS_TYPE = 'CMASS4'
         ENDIF

         FOUND = 'N'

j_do:    DO J=1,NPMASS
            IF (PMASS(J,1) == PROPERTY_ID) THEN
               FOUND = 'Y'
               EXIT j_do
            ENDIF
         ENDDO j_do

         IF (FOUND == 'N') THEN
            WRITE(ERR,1404) PROPERTY_NAME, PROPERTY_ID, CMASS_TYPE, CMASS(I,1)
            WRITE(F06,1404) PROPERTY_NAME, PROPERTY_ID, CMASS_TYPE, CMASS(I,1)
            IERROR = IERROR + 1
            FATAL_ERR = FATAL_ERR + 1
         ENDIF

      ENDDO

! **********************************************************************************************************************************
! Quit if IERROR > 0

      IF (IERROR > 0) THEN
         CALL OUTA_HERE ( 'Y' )
      ENDIF



      RETURN

! **********************************************************************************************************************************
 1401 FORMAT(' *ERROR  1401: ',A,I8,' ON ',A,I8,' UNDEFINED')

 1403 FORMAT(' *ERROR  1403: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ELEMENT TYPE ',A,' NOT VALID')

 1404 FORMAT(' *ERROR  1404: ',A,I8,' ON C',A,I8,' UNDEFINED')

 1406 FORMAT(' *ERROR  1406: BOTH A PSHELL AND A PCOMP WITH ID ',I8,' WERE FOUND FOR C',A,I8,'. THEY CANNOT HAVE SAME ID.')


! **********************************************************************************************************************************

      END SUBROUTINE ELEM_PROP_MATL_IIDS


      SUBROUTINE ELSAVE

! Saves element data to file LINK1G.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR,     F06,     L1G
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, DATA_NAM_LEN, MMATL, MPBAR, MPBEAM, MPBUSH, MPELAS, MPROD, MPSHEL,          &
                                         MPSOLID, MPUSER1,MPUSERIN, MRMATLC, MRPBAR, MRPBEAM, MRPBUSH, MRPELAS, MRPROD, MPSHEAR,   &
                                         MRPSHEAR, MRPSHEL, MRPUSER1, NBAROFF, NBUSHOFF, NEDAT, NELE, NMATANGLE, NMATL, MPCOMP0,   &
                                         MRPCOMP0, MPCOMP_PLIES, MRPCOMP_PLIES, MUSERIN_MAT_NAMES, NPBAR, NPBEAM, NPBUSH, NPCOMP,  &
                                         NPELAS, NPLATEOFF, NPLATETHICK, NPROD, NPSHEAR, NPSHEL, NPSOLID, NPUSER1, NPUSERIN, NVVEC
      USE PARAMS, ONLY                :  CBMIN3, CBMIN4, IORQ1M, IORQ1S, IORQ1B, IORQ2B, IORQ2T
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  BAROFF, BUSHOFF, EDAT, EOFF, EPNT, ESORT1, ESORT2, ETYPE, MATANGLE, MATL, RMATL,PBAR,     &
                                         RPBAR, PBEAM, RPBEAM, PBUSH, RPBUSH, PCOMP, RPCOMP, PELAS, RPELAS, PROD, RPROD, PSHEAR,   &
                                         RPSHEAR, PSHEL, RPSHEL, PSOLID, PUSER1, RPUSER1, PUSERIN, PLATEOFF, PLATETHICK,           &
                                         USERIN_MAT_NAMES, VVEC
      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELSAVE'
      CHARACTER(LEN=DATA_NAM_LEN)     :: DATA_SET_NAME     ! A data set name for output purposes

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: PCOMP_PLIES       ! Number of plies in 1 PCOMP entry incl sym plies not explicitly defined




! **********************************************************************************************************************************
! Write element data

      DATA_SET_NAME = 'ETYPE, EPNT, ESORT1, ESORT2, EOFF'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NELE
      DO I = 1,NELE
         WRITE(L1G) ETYPE(I)
         WRITE(L1G) EPNT(I)
         WRITE(L1G) ESORT1(I)
         WRITE(L1G) ESORT2(I)
         WRITE(L1G) EOFF(I)
      ENDDO

      DATA_SET_NAME = 'EDAT'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NEDAT
      DO I = 1,NEDAT
         WRITE(L1G) EDAT(I)
      ENDDO

! Write element parameters

      DATA_SET_NAME = 'ELEM PARAMETERS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) IORQ1M
      WRITE(L1G) IORQ1S
      WRITE(L1G) IORQ1B
      WRITE(L1G) IORQ2B
      WRITE(L1G) IORQ2T
      WRITE(L1G) CBMIN3
      WRITE(L1G) CBMIN4

! Write BAR, BEAM v vectors

      DATA_SET_NAME = 'V VECTORS IN GLOBAL COORDS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NVVEC
      DO I=1,NVVEC
         WRITE(L1G) (VVEC(I,J),J=1,3)
      ENDDO

! Write BAR, BEAM offsets

      DATA_SET_NAME = 'BAR, BEAM OFFSETS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NBAROFF
      DO I = 1,NBAROFF
         DO J = 1,6
            WRITE(L1G) BAROFF(I,J)
         ENDDO
      ENDDO

! Write BUSH offsets

      DATA_SET_NAME = 'BUSH OFFSETS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NBUSHOFF
      DO I = 1,NBUSHOFF
         DO J = 1,6
            WRITE(L1G) BUSHOFF(I,J)
         ENDDO
      ENDDO

! Write plate offsets

      DATA_SET_NAME = 'PLATE OFFSETS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPLATEOFF
      DO I = 1,NPLATEOFF
         WRITE(L1G) PLATEOFF(I)
      ENDDO

! Write plate thicknesses

      DATA_SET_NAME = 'PLATE THICKNESSES FROM CONNECTION ENTRIES'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPLATETHICK
      DO I = 1,NPLATETHICK
         WRITE(L1G) PLATETHICK(I)
      ENDDO

! Write property data

      DATA_SET_NAME = 'PBAR, RPBAR'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPBAR
      WRITE(L1G) MPBAR
      WRITE(L1G) MRPBAR
      DO I = 1,NPBAR
         DO J=1,MPBAR
            WRITE(L1G) PBAR(I,J)
         ENDDO
         DO J = 1,MRPBAR
            WRITE(L1G) RPBAR(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PBEAM, RPBEAM'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPBEAM
      WRITE(L1G) MPBEAM
      WRITE(L1G) MRPBEAM
      DO I = 1,NPBEAM
         DO J=1,MPBEAM
            WRITE(L1G) PBEAM(I,J)
         ENDDO
         DO J = 1,MRPBEAM
            WRITE(L1G) RPBEAM(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PBUSH, RPBUSH'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPBUSH
      WRITE(L1G) MPBUSH
      WRITE(L1G) MRPBUSH
      DO I = 1,NPBUSH
         DO J=1,MPBUSH
            WRITE(L1G) PBUSH(I,J)
         ENDDO
         DO J = 1,MRPBUSH
            WRITE(L1G) RPBUSH(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PROD, RPROD'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPROD
      WRITE(L1G) MPROD
      WRITE(L1G) MRPROD
      DO I = 1,NPROD
         DO J=1,MPROD
            WRITE(L1G) PROD(I,J)
         ENDDO
         DO J = 1,MRPROD
            WRITE(L1G) RPROD(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PELAS, RPELAS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPELAS
      WRITE(L1G) MPELAS
      WRITE(L1G) MRPELAS
      DO I = 1,NPELAS
         DO J=1,MPELAS
            WRITE(L1G) PELAS(I,J)
         ENDDO
         DO J = 1,MRPELAS
            WRITE(L1G) RPELAS(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PSHEAR, RPSHEAR'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPSHEAR
      WRITE(L1G) MPSHEAR
      WRITE(L1G) MRPSHEAR
      DO I = 1,NPSHEAR
         DO J = 1,MPSHEAR
            WRITE(L1G) PSHEAR(I,J)
         ENDDO
         DO J = 1,MRPSHEAR
            WRITE(L1G) RPSHEAR(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PSHEL, RPSHEL'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPSHEL
      WRITE(L1G) MPSHEL
      WRITE(L1G) MRPSHEL
      DO I = 1,NPSHEL
         DO J = 1,MPSHEL
            WRITE(L1G) PSHEL(I,J)
         ENDDO
         DO J = 1,MRPSHEL
            WRITE(L1G) RPSHEL(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PCOMP, RPCOMP'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPCOMP
      WRITE(L1G) MPCOMP0
      WRITE(L1G) MPCOMP_PLIES
      WRITE(L1G) MRPCOMP0
      WRITE(L1G) MRPCOMP_PLIES
      DO I = 1,NPCOMP
         PCOMP_PLIES = PCOMP(I,5)
         WRITE(L1G) PCOMP_PLIES
         DO J = 1,MPCOMP0+MPCOMP_PLIES*PCOMP_PLIES
            WRITE(L1G) PCOMP(I,J)
         ENDDO
         DO J = 1,MRPCOMP0+MRPCOMP_PLIES*PCOMP_PLIES
            WRITE(L1G) RPCOMP(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PSOLID'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPSOLID
      WRITE(L1G) MPSOLID
      DO I = 1,NPSOLID
         DO J = 1,MPSOLID
            WRITE(L1G) PSOLID(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PUSER1, RPUSER1'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPUSER1
      WRITE(L1G) MPUSER1
      WRITE(L1G) MRPUSER1
      DO I = 1,NPUSER1
         DO J=1,MPUSER1
            WRITE(L1G) PUSER1(I,J)
         ENDDO
         DO J = 1,MRPUSER1
            WRITE(L1G) RPUSER1(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PUSERIN'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPUSERIN
      WRITE(L1G) MPUSERIN
      DO I = 1,NPUSERIN
         DO J=1,MPUSERIN
            WRITE(L1G) PUSERIN(I,J)
         ENDDO
      ENDDO

! Write USERIN_MAT_NAMES

      DATA_SET_NAME = 'USERIN_MAT_NAMES'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPUSERIN                                  ! NOTE: there are as many USERIN_MAT_NAMES as PUSERIN entries
      WRITE(L1G) MUSERIN_MAT_NAMES
      DO I = 1,NPUSERIN
         DO J=1,MUSERIN_MAT_NAMES
            WRITE(L1G) USERIN_MAT_NAMES(I,J)
         ENDDO
      ENDDO

! Write material data

      DATA_SET_NAME = 'MATL, RMATL'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NMATL
      WRITE(L1G) MMATL
      WRITE(L1G) MRMATLC
      DO I = 1,NMATL
         DO J=1,MMATL
            WRITE(L1G) MATL(I,J)
         ENDDO
         DO J = 1,MRMATLC
            WRITE(L1G) RMATL(I,J)
         ENDDO
      ENDDO

! Write material property angles

      DATA_SET_NAME = 'MATERIAL PROPERTY ANGLES'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NMATANGLE
      DO I = 1,NMATANGLE
         WRITE(L1G) MATANGLE(I)
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE ELSAVE

   END MODULE ELEMENT_MODEL_INDEXING
