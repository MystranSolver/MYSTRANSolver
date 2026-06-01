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

      SUBROUTINE EC_PARTN ( CARD1, IERR )

! EC_PARTN reads in the Exec Control entry PARTN (partition a matrix). The NASTRAN form of the entry is:

!     PARTN A,CP,RP/A11,A21,A12,A22/SYM/TYPE/F11,F21,F12,F22

! The MYSTRAN form will be the same except only the data up to the 1st slash (/), if it exists, is processed:

!     PARTN A,CP,RP

! Matrix (A) will be partitioned using the column (CP) and row (RP) partitioning vectors.  The partition will be what is output as
! the OUTPUT4 matrix. PARTN must follow an OUTPUT4 requesting output of the matrix A

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG
      USE SCONTR, ONLY                :  EC_ENTRY_LEN
      USE IOUNT1, ONLY                :  ERR, F06, MOU4, OU4, OU4_ELM_OTM, OU4_GRD_OTM
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE OUTPUT4_MATRICES, ONLY      :  NUM_OU4_REQUESTS, NUM_PARTN_REQUESTS, OU4_PART_VEC_NAMES, OU4_PART_MAT_NAMES,             &
                                         ACT_OU4_MYSTRAN_NAMES, ACT_OU4_OUTPUT_NAMES

      USE EC_PARTN_USE_IFs                                 ! Added 2019/07/14
      USE TO_UPPER_Interface

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: CARD1             ! Card read in LOADE and shifted to begin in col 1
      CHARACTER(LEN=LEN(CARD1))       :: CARD2             ! CARD1 truncated at $ (trailing comment) if there is one
      CHARACTER(LEN=EC_ENTRY_LEN)     :: DATA_80(3)        ! Temp slot for holding data until lead/trail blanks stripped
      CHARACTER(16*BYTE)              :: DATA_16(3)        ! Matrix name read from OUTPUT4 entry
      CHARACTER( 1*BYTE)              :: FOUND             ! 'Y' if we found something we were looking for

      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Error indicator. If CHAR not found, IERR set to 1
      INTEGER(LONG)                   :: DATA_BEG          ! Column where data begins (after OUTPUT4)
      INTEGER(LONG)                   :: DATA_END          ! Column where data ends (after ist slash or at $)
      INTEGER(LONG)                   :: COMMA_COL(3)      ! Column where comma is found in CARD2
      INTEGER(LONG)                   :: I,J               ! DO loop index
      INTEGER(LONG)                   :: JBEG              ! Beg col in data




! **********************************************************************************************************************************
      IF ((DEBUG(197) == 2) .OR. (DEBUG(197) == 3)) THEN
         WRITE(F06,'(A)') '******************************************************************************************************'
      ENDIF

! Initialize

      IERR =  0

      DO I=1,3
         DATA_16(I)(1:) = ' '
         DATA_80(I)(1:) = ' '
      ENDDO

      CARD2(1:) = ' '

! Find where data begins after PARTN. NOTE: input CARD1 was shifted in subr LOADE so that the "PARTN" began in col 1. Thus the
! data will begin after col 6 (after "PARTN")

      DATA_BEG = 1
      DO I=6,EC_ENTRY_LEN
         IF (CARD1(I:I) == ' ') THEN
            CYCLE
         ELSE
            DATA_BEG = I
            EXIT
         ENDIF
      ENDDO

! Find where data we want ends (no data will be kept past a "/" or a "$" character)

      DATA_END = EC_ENTRY_LEN
      DO I=6,EC_ENTRY_LEN
         IF ((CARD1(I:I) == "/") .OR. (CARD2(I:I) == "$")) THEN
            DATA_END = I - 1
            EXIT
         ENDIF
      ENDDO

! The data we want is in CARD1 from DATA_BEG thru DATA_END

      CARD2(1:) = CARD1(DATA_BEG:DATA_END)

! Find commas (there should be no 3rd comma - this is checked below)

      COMMA_COL(1) =                INDEX ( CARD2             , "," )
      COMMA_COL(2) = COMMA_COL(1) + INDEX ( CARD2(COMMA_COL(1)+1:), "," )
      COMMA_COL(3) = COMMA_COL(2) + INDEX ( CARD2(COMMA_COL(2)+1:), "," )

! Make sure there are 2 commas

      IF ((COMMA_COL(1) == 0) .OR. (COMMA_COL(2) == COMMA_COL(1))) THEN
         IERR = IERR + 1
         WRITE(ERR,1031) '2'
         WRITE(F06,1031) '2'
      ENDIF

! Make sure there are not more than 2 commas

      IF (COMMA_COL(3) /= COMMA_COL(2)) THEN
         IERR = IERR + 1
         WRITE(ERR,1032) '2'
         WRITE(F06,1032) '2'
      ENDIF

! Process data from CARD2 if there were no errors

nerr: IF (IERR == 0) THEN

         NUM_PARTN_REQUESTS = NUM_PARTN_REQUESTS + 1
                                                           ! DATA_80(1-3) is the data from CARD2 bet commas from data beg to //:
         DATA_80(1)(1:) = CARD2 ( 1              : COMMA_COL(1)-1 ) ! This is the name of the matrix to be partitioned
         DATA_80(2)(1:) = CARD2 ( COMMA_COL(1)+1 : COMMA_COL(2)-1 ) ! This is the name of the col partitioning vector
         DATA_80(3)(1:) = CARD2 ( COMMA_COL(2)+1 : DATA_END       ) ! This is the name of the row partitioning vector

         DO I=1,3                                          ! Get 16 char matrix names by stripping leading/trailing blanks
            IF (DATA_80(I)(1:) /= ' ') THEN
               DO J=1,EC_ENTRY_LEN
                  IF (DATA_80(I)(J:J) == ' ') THEN
                     CYCLE
                  ELSE
                     JBEG = J
                     EXIT
                  ENDIF
               ENDDO
               DO J=EC_ENTRY_LEN,JBEG,-1
                  IF (DATA_80(I)(J:J) == ' ') THEN
                     CYCLE
                  ELSE
                     EXIT
                  ENDIF
               ENDDO
               DATA_16(I)(1:) = DATA_80(I)(JBEG:JBEG+15)
            ELSE
               DATA_16(I)(1:) = ' '
            ENDIF
         ENDDO

         FOUND = 'N'
         DO I=1,NUM_OU4_REQUESTS                           ! Set names of the matrix to be partitioned, the partitions and the vecs
            ! DATA_16(1) was upper-cased when reading the card.
            ! Here, upper-case the allowed names so that comparison is case-insensitive
            ! and thus allows KRRcb, MRRcb, KLR(t).
            IF (DATA_16(1) == TO_UPPER(ACT_OU4_MYSTRAN_NAMES(I))) THEN
               FOUND = 'Y'
               OU4_PART_MAT_NAMES(I,1) = DATA_16(1)
               OU4_PART_VEC_NAMES(I,1) = DATA_16(2)
               OU4_PART_VEC_NAMES(I,2) = DATA_16(3)
            ENDIF
         ENDDO

         IF (FOUND == 'N') THEN                            ! Matrix to be partitioned is not an OU4 matrix or an OU4 request
            IERR = IERR + 1
            WRITE(ERR,1040) DATA_16(1)
            WRITE(F06,1040) DATA_16(1)
         ENDIF

      ENDIF nerr

      IF ((DEBUG(197) == 2) .OR. (DEBUG(197) == 3)) THEN

         WRITE(F06,'(A)') 'Debug output from subr EC_PARTN'
         WRITE(F06,'(A)') '---------------------------------'
         WRITE(F06,*)
         WRITE(F06,99881) 'CARD1', CARD1
         WRITE(F06,*)
         WRITE(F06,99881) 'CARD2', CARD2
         WRITE(F06,*)

         WRITE(F06,99882) DATA_BEG
         WRITE(F06,*)

         WRITE(F06,99883) (COMMA_COL(I),I=1,2)
         WRITE(F06,*)

         IF (COMMA_COL(3) /= COMMA_COL(2)) THEN
            WRITE(F06,99884) COMMA_COL(3)
            WRITE(F06,*)
         ENDIF

         DO I=1,NUM_OU4_REQUESTS

            IF (ACT_OU4_MYSTRAN_NAMES(I)(1:16) == DATA_80(1)(1:16)) THEN
               WRITE(F06,99901) I, 'Matrix to be partitioned: OU4_PART_MAT_NAMES(I,1)', OU4_PART_MAT_NAMES(I,1)
               WRITE(F06,*)
               WRITE(F06,99901) I, 'Partitioning vec 1 name : OU4_PART_VEC_NAMES(I,1)', OU4_PART_VEC_NAMES(I,1)
               WRITE(F06,99901) I, 'Partitioning vec 2 name : OU4_PART_VEC_NAMES(I,2)', OU4_PART_VEC_NAMES(I,2)
            ENDIF

         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************
 1025 FORMAT(' *ERROR  1025: THE PARTN ENTRY FORMAT REQUIRES ONE SLASH (/) AFTER THE SPECIFICATION OF THE PARTITIONING VECTORS')

 1031 FORMAT(' *ERROR  1031: EXEC CONTROL ENTRY OUTPUT4 MUST HAVE ',A,' COMMAS (AS IN PARTN A,CP,RP).')

 1032 FORMAT(' *ERROR  1032: EXEC CONTROL ENTRY OUTPUT4 CANNOT HAVE MORE THAN ',A,' COMMAS.')

 1039 FORMAT(' *ERROR  1039: A SLASH (/) MUST FOLLOW THE 3 PIECES OF DATA (MATRIX NAME, 2 PARTITION VECTORS) SEPARATED BY 2 COMMAS')

 1040 FORMAT(' *ERROR  1040: PARTN REQUEST HAS THE "',A,'" MATRIX REQUESTED TO BE PARTITIONED. HOWEVER, EITHER THIS IS NOT A'      &
                    ,/,14X,' VALID OUTPUT4 MATRIX NAME OR THE OUTPUT4 REQUEST FOR IT WAS NOT FOUND BEFORE THE PARTN REQUEST')

99881 FORMAT(A,':'/,A)

99882 FORMAT('Data begins in OUTPUT4 entry in col ',I3)

99883 FORMAT('COMMA_COL(1:2) = ',3(I4,','),I4,' (cols where commas exist)')

99884 FORMAT('There are more than 4 commas. A 5th one was found in col ',I4)

99901 FORMAT(' I =',I3,2X,A,2X,'"',A,'"')

! **********************************************************************************************************************************

      END SUBROUTINE EC_PARTN
