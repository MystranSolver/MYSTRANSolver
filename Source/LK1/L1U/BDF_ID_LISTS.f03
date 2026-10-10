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

! THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
! LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN
! NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
! LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
! OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
! THE SOFTWARE.
! _______________________________________________________________________________________________________

! End MIT license text.

      MODULE BDF_ID_LISTS

! The ID list of a Bulk Data entry (SPC1, ASET1/OMIT1, USET1, PARTVEC1, SPOINT) runs from a field of the parent line through the
! fields 2-9 of its continuations. An item of the list is an ID or a range "ID1 THRU ID2" (ID2 >= ID1); items may be mixed in any
! order and a range may run from one line to the next, so "1 THRU 2 21" is the IDs 1, 2 and 21 and "2 THRU 5 7" is 2, 3, 4, 5
! and 7. Blank fields are skipped. Anything else in a list field (a "THRU" that does not stand between two IDs, a word or a real
! number) stops the run with an error, so that no ID of the list is lost without a message (before, the readers took a range
! only as the whole entry "ID1 THRU ID2" and dropped the IDs after it).

      IMPLICIT NONE

      PRIVATE

      PUBLIC :: READ_ID_LIST

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE READ_ID_LIST ( CARD, LARGE_FLD_INP, FIRST_FLD, COUNT_ONLY, NRANGE, RANGES, NBAD )

! CARD is the parent line on input; the continuations are read here (on output CARD is the line after the last continuation, as
! after subr NEXTC). FIRST_FLD is the field of the parent line where the list begins.
! COUNT_ONLY = .TRUE. in the first pass over the Bulk Data (subr LOADB0): the continuations are read with NEXTC0/NEXTC20 and no
! message is written; the items are the same as in the second pass, so a count made from them (SPOINT) is the one used there.
! Otherwise errors are written and counted in FATAL_ERR, and subr CRDERR is called for each line (with the field errors of the
! parent line that the caller found before calling, e.g. in its components field).
! RANGES(1:2,1:NRANGE) are the items in their order, a single ID as (ID, ID). NBAD is the number of fields in error.

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE IOUNT1, ONLY                :  ERR, F06
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, FATAL_ERR, IERRFL, JCARD_LEN, JF
      USE BDF_CARD_CONTINUATIONS, ONLY:  MKJCARD, NEXTC, NEXTC0, NEXTC2, NEXTC20
      USE BDF_FIELD_VALIDATION, ONLY  :  CRDERR, I4FLD
      USE BDF_SET_SYNTAX, ONLY        :  TOKCHK

      CHARACTER(LEN=*), INTENT(INOUT) :: CARD              ! The parent line of the entry
      CHARACTER(LEN=*), INTENT(IN)    :: LARGE_FLD_INP     ! If 'Y', CARD is large field format
      INTEGER(LONG)   , INTENT(IN)    :: FIRST_FLD         ! The field of the parent line where the list begins
      LOGICAL         , INTENT(IN)    :: COUNT_ONLY        ! First pass: no messages, NEXTC0/NEXTC20
      INTEGER(LONG)   , INTENT(OUT)   :: NRANGE            ! Number of items
      INTEGER(LONG)   , ALLOCATABLE, INTENT(OUT) :: RANGES(:,:)! The items (ID1, ID2)
      INTEGER(LONG)   , INTENT(OUT)   :: NBAD              ! Number of fields in error

      CHARACTER(LEN=BD_ENTRY_LEN)     :: CHILD             ! A continuation read by NEXTC2/NEXTC20
      CHARACTER(LEN=JCARD_LEN)        :: JCARD(10)         ! The 10 fields of a line
      CHARACTER(LEN=JCARD_LEN)        :: ENTRY_NAME        ! Field 1 of the parent line
      CHARACTER(LEN=JCARD_LEN)        :: ENTRY_ID          ! Field 2 of the parent line
      CHARACTER(LEN=JCARD_LEN)        :: FLD               ! A field, left adjusted
      CHARACTER(LEN=8)                :: TOKTYP            ! Output of TOKCHK
      INTEGER(LONG), ALLOCATABLE      :: GROW(:,:)         ! RANGES when it is enlarged
      INTEGER(LONG)                   :: ICONT             ! 1 if the next line is a continuation
      INTEGER(LONG)                   :: IERR              ! Output of NEXTC
      INTEGER(LONG)                   :: ID                ! An ID read from a field
      INTEGER(LONG)                   :: J                 ! A field
      INTEGER(LONG)                   :: J1                ! The first list field of a line
      INTEGER(LONG)                   :: LINE              ! The line of the entry (1 = parent)
      LOGICAL                         :: PENDING_THRU      ! A "THRU" waits for its ending ID
      LOGICAL                         :: LAST_SINGLE       ! The last item is one ID that a "THRU" may extend
      LOGICAL                         :: OK                ! The field holds an ID

      NRANGE = 0
      NBAD   = 0
      ALLOCATE ( RANGES(2,16) )
      PENDING_THRU = .FALSE.
      LAST_SINGLE  = .FALSE.

      CALL MKJCARD ( 'READ_ID_LIST', CARD, JCARD )
      ENTRY_NAME = JCARD(1)
      ENTRY_ID   = JCARD(2)
      LINE = 1
      J1   = FIRST_FLD

      DO

         DO J=J1,9
            FLD = ADJUSTL(JCARD(J))
            IF (FLD == ' ') CYCLE
            CALL TOKCHK ( FLD(1:8), TOKTYP )

            IF (TOKTYP == 'THRU    ') THEN                 ! "THRU": the last item must be one ID, and an ID must follow
               IF (LAST_SINGLE .AND. (.NOT.PENDING_THRU)) THEN
                  PENDING_THRU = .TRUE.
               ELSE
                  NBAD = NBAD + 1
                  IF (.NOT.COUNT_ONLY) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1214) TRIM(ENTRY_NAME), TRIM(ENTRY_ID), JF(J), LINE
                     WRITE(F06,1214) TRIM(ENTRY_NAME), TRIM(ENTRY_ID), JF(J), LINE
                  ENDIF
                  PENDING_THRU = .FALSE.
                  LAST_SINGLE  = .FALSE.
               ENDIF
               CYCLE
            ENDIF

            OK = .FALSE.                                   ! An ID is an integer to both TOKCHK and I4FLD (the same in both passes)
            IF (TOKTYP == 'INTEGER ') THEN
               CALL I4FLD ( JCARD(J), JF(J), ID )
               OK = (IERRFL(J) == 'N')
            ELSE IF (.NOT.COUNT_ONLY) THEN
               CALL I4FLD ( JCARD(J), JF(J), ID )          ! Its message, if it is not an integer to I4FLD either
               IF (IERRFL(J) == 'N') THEN
                  FATAL_ERR = FATAL_ERR + 1
                  WRITE(ERR,1215) TRIM(ENTRY_NAME), TRIM(ENTRY_ID), JF(J), LINE, TRIM(FLD)
                  WRITE(F06,1215) TRIM(ENTRY_NAME), TRIM(ENTRY_ID), JF(J), LINE, TRIM(FLD)
               ENDIF
            ENDIF

            IF (.NOT.OK) THEN
               NBAD = NBAD + 1
               IF (PENDING_THRU) NRANGE = NRANGE - 1       ! The range is dropped with its start (the run stops anyway)
               PENDING_THRU = .FALSE.
               LAST_SINGLE  = .FALSE.
            ELSE IF (PENDING_THRU) THEN
               PENDING_THRU = .FALSE.
               LAST_SINGLE  = .FALSE.
               IF (ID < RANGES(1,NRANGE)) THEN
                  NBAD = NBAD + 1
                  IF (.NOT.COUNT_ONLY) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,1128) TRIM(ENTRY_NAME), TRIM(ENTRY_ID), RANGES(1,NRANGE), ID
                     WRITE(F06,1128) TRIM(ENTRY_NAME), TRIM(ENTRY_ID), RANGES(1,NRANGE), ID
                  ENDIF
                  NRANGE = NRANGE - 1
               ELSE
                  RANGES(2,NRANGE) = ID
               ENDIF
            ELSE
               IF (NRANGE == SIZE(RANGES,2)) THEN
                  ALLOCATE ( GROW(2,2*NRANGE) )
                  GROW(:,1:NRANGE) = RANGES(:,1:NRANGE)
                  CALL MOVE_ALLOC ( GROW, RANGES )
               ENDIF
               NRANGE = NRANGE + 1
               RANGES(1,NRANGE) = ID
               RANGES(2,NRANGE) = ID
               LAST_SINGLE = .TRUE.
            ENDIF
         ENDDO

         IF (.NOT.COUNT_ONLY) CALL CRDERR ( CARD )

         IF (COUNT_ONLY) THEN                              ! The next line, if it is a continuation
            IF (LARGE_FLD_INP == 'N') THEN
               CALL NEXTC0  ( CARD, ICONT, IERR )
            ELSE
               CALL NEXTC20 ( CARD, ICONT, IERR, CHILD )
               CARD = CHILD
            ENDIF
         ELSE
            IF (LARGE_FLD_INP == 'N') THEN
               CALL NEXTC  ( CARD, ICONT, IERR )
            ELSE
               CALL NEXTC2 ( CARD, ICONT, IERR, CHILD )
               CARD = CHILD
            ENDIF
         ENDIF
         IF (ICONT /= 1) EXIT
         CALL MKJCARD ( 'READ_ID_LIST', CARD, JCARD )
         LINE = LINE + 1
         J1   = 2

      ENDDO

      IF (PENDING_THRU) THEN                               ! The list ends with "THRU"
         NBAD   = NBAD + 1
         NRANGE = NRANGE - 1
         IF (.NOT.COUNT_ONLY) THEN
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,1216) TRIM(ENTRY_NAME), TRIM(ENTRY_ID)
            WRITE(F06,1216) TRIM(ENTRY_NAME), TRIM(ENTRY_ID)
         ENDIF
      ENDIF

      RETURN

! **********************************************************************************************************************************
 1128 FORMAT(' *ERROR  1128: ON ',A,' ',A,' THE IDs MUST BE IN INCREASING ORDER FOR THRU OPTION (',I0,' THRU ',I0,')')

 1214 FORMAT(' *ERROR  1214: ON ',A,' ',A,' "THRU" IN FIELD ',I0,' OF LINE ',I0,' DOES NOT FOLLOW AN ID. AN ID LIST HOLDS IDs',  &
             ' AND RANGES "ID1 THRU ID2"')

 1215 FORMAT(' *ERROR  1215: ON ',A,' ',A,' FIELD ',I0,' OF LINE ',I0,' IS NOT AN ID: "',A,'"')

 1216 FORMAT(' *ERROR  1216: ON ',A,' ',A,' THE ID LIST ENDS WITH "THRU"; A RANGE IS "ID1 THRU ID2"')

! **********************************************************************************************************************************

      END SUBROUTINE READ_ID_LIST

      END MODULE BDF_ID_LISTS
