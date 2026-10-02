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

   MODULE SORTING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: SORT_INT1, SORT_INT1_REAL1, SORT_INT1_REAL3, SORT_INT2, SORT_INT2_REAL1, SORT_INT3_CHAR2, SORT_REAL1_INT1, SORT_GRID_RGRID, SORT_TDOF, CALC_VEC_SORT_ORDER

   CONTAINS

      SUBROUTINE SORTLEN ( NLEN, JCT )

! Calculates shell sort length parameter, JCT

      USE PENTIUM_II_KIND, ONLY       :  LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR
      USE SCONTR, ONLY                :  BLNK_SUB_NAM
      USE CONSTANTS_1, ONLY           :  TWO
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORTLEN'

      INTEGER(LONG), INTENT(IN)       :: NLEN              ! Length of the array that will be sorted in the calling procedure
      INTEGER(LONG), INTENT(OUT)      :: JCT               ! Sort parameter to be used by calling procedure
      INTEGER(LONG)                   :: MAX_JCT           ! Max practical value of JCT to use in sort by the calling procedure.
!                                                            Values of JCT > MAX_JCT will not cause any error, but will introduce
!                                                            inefficiency into the sort (a DO loop will run excessively).


      INTRINSIC DLOG



! **********************************************************************************************************************************
! Initialize outputs

      JCT = 0

! MAX_JCT is the largest integer value of JCT that is practical to use in the shell sort routine in the calling
! routine. Larger values cause unnecessary DO looping in that routine.

      MAX_JCT  = FLOOR(  (DLOG(DBLE(NLEN)+1.D0)) / (DLOG(TWO))  )

! Calculate shell sort parameter JCT based on array size (NLEN)

      IF (NLEN <= 5) THEN
         JCT = 1
      ELSE IF ((NLEN >       5) .AND. (NLEN <=      13)) THEN ! Add      8
         JCT = 2
      ELSE IF ((NLEN >      13) .AND. (NLEN <=      29)) THEN ! Add     16
         JCT = 3
      ELSE IF ((NLEN >      29) .AND. (NLEN <=      61)) THEN ! Add     32
         JCT = 4
      ELSE IF ((NLEN >      61) .AND. (NLEN <=     125)) THEN ! Add     64
         JCT = 5
      ELSE IF ((NLEN >     125) .AND. (NLEN <=     253)) THEN ! Add    128
         JCT = 6
      ELSE IF ((NLEN >     253) .AND. (NLEN <=     509)) THEN ! Add    256
         JCT = 7
      ELSE IF ((NLEN >     509) .AND. (NLEN <=    1021)) THEN ! Add    512
         JCT = 8
      ELSE IF ((NLEN >    1022) .AND. (NLEN <=    2045)) THEN ! Add   1024
         JCT = 9
      ELSE IF ((NLEN >    2045) .AND. (NLEN <=    4093)) THEN ! Add   2048
         JCT = 10
      ELSE IF ((NLEN >    4093) .AND. (NLEN <=    8189)) THEN ! Add   4096
         JCT = 11
      ELSE IF ((NLEN >    8189) .AND. (NLEN <=   16381)) THEN ! Add   8192
         JCT = 12
      ELSE IF ((NLEN >   16381) .AND. (NLEN <=   32765)) THEN ! Add  16384
         JCT = 13
      ELSE IF ((NLEN >   32765) .AND. (NLEN <=   65533)) THEN ! Add  32768
         JCT = 14
      ELSE IF ((NLEN >   65533) .AND. (NLEN <=  131069)) THEN ! Add  65536
         JCT = 15
      ELSE IF ((NLEN >  131069) .AND. (NLEN <=  262141)) THEN ! Add 131072
         JCT = 16
      ELSE IF ((NLEN >  262141) .AND. (NLEN <=  524285)) THEN ! Add 262144
         JCT = 17
      ELSE IF ((NLEN >  524285) .AND. (NLEN <= 1048573)) THEN ! Add 524288
         JCT = 18
      ELSE
         JCT = MAX_JCT
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE SORTLEN


      SUBROUTINE SORT_INT1 ( CALLING_SUBR, MESSAG, NSIZE, IARRAY )

! Performs shell sort on integer array IARRAY (of size NSIZE) to put it into numerically increasing order

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_INT1'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to be written out if this subr fails to sort
      CHARACTER( 3*BYTE)              :: SORTED            ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(IN)       :: NSIZE             ! No. rows in arrays IARRAY, RARRAY
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY(NSIZE)     ! Integer array
      INTEGER(LONG)                   :: I,K,M             ! DO loop indices
      INTEGER(LONG)                   :: IDUM              ! Dummy values in IARRAY used when switching IARRAY rows during sort.
      INTEGER(LONG)                   :: IFLIP             ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT               ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM              ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                 ! An array index
      INTEGER(LONG)                   :: SORTPK            ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM          ! How many times the sort has to be performed in order for the data




! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (IARRAY(M) > IARRAY(N)) THEN
                     IDUM      = IARRAY(M)
                     IARRAY(M) = IARRAY(N)
                     IARRAY(N) = IDUM
                     IFLIP = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO

! Make sure array is sorted. If not, repeat algorithm SORT_MAX times

         SORTED = 'YES'                                    ! Make sure array is sorted
         DO I=1,NSIZE-1
            IF (IARRAY(I) > IARRAY(I+1)) THEN
               SORTED = 'NO '
               EXIT
            ENDIF
         ENDDO

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')

! **********************************************************************************************************************************

      END SUBROUTINE SORT_INT1


      SUBROUTINE SORT_INT1_REAL1 ( CALLING_SUBR, MESSAG, NSIZE, IARRAY, RARRAY )

! Performs shell sort on integer array IARRAY to put it into numerically increasing order and sorts real array RARRAY
! along with IARRAY. Both arrays ore of size NSIZE

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_INT1_REAL1'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to be written out if this subr fails to sort
      CHARACTER( 3*BYTE)              :: SORTED            ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(IN)       :: NSIZE             ! No. rows in arrays IARRAY, RARRAY
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY(NSIZE)     ! Array of integer values
      INTEGER(LONG)                   :: I,K,M             ! DO loop indices
      INTEGER(LONG)                   :: IDUM              ! Dummy values in IARRAY used when switching IARRAY rows during sort.
      INTEGER(LONG)                   :: IFLIP             ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT               ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM              ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                 ! An array index
      INTEGER(LONG)                   :: SORTPK            ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM          ! How many times the sort has to be performed in order for the data


      REAL(DOUBLE),  INTENT(INOUT)    :: RARRAY(NSIZE)     ! Array of real values
      REAL(DOUBLE)                    :: RDUM              ! Dummy values in RARRAY used when switching RARRAY rows during the sort.



! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (IARRAY(M) > IARRAY(N)) THEN
                     IDUM      = IARRAY(M)
                     RDUM      = RARRAY(M)
                     IARRAY(M) = IARRAY(N)
                     RARRAY(M) = RARRAY(N)
                     IARRAY(N) = IDUM
                     RARRAY(N) = RDUM
                     IFLIP = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         SORTED = 'YES'                                    ! Make sure array is sorted.
chk_sort:DO I=1,NSIZE-1
            IF (IARRAY(I) > IARRAY(I+1)) THEN
            SORTED = 'NO '
            EXIT chk_sort
            ENDIF
         ENDDO chk_sort

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                              ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')

! **********************************************************************************************************************************

      END SUBROUTINE SORT_INT1_REAL1


      SUBROUTINE SORT_INT1_REAL3 ( CALLING_SUBR, MESSAG, NSIZE, IARRAY, RARRAY )

! Performs shell sort on integer array IARRAY to put it into numerically increasing order and sorts real array RARRAY (with 3 cols)
! along with IARRAY. Both arrays have NSIZE rows

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_INT1_REAL3'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to be written out if this subr fails to sort
      CHARACTER( 3*BYTE)              :: SORTED            ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(IN)       :: NSIZE             ! No. rows in arrays IARRAY, RARRAY
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY(NSIZE)     ! Array of integer values
      INTEGER(LONG)                   :: I,K,M             ! DO loop indices
      INTEGER(LONG)                   :: IDUM              ! Dummy values in IARRAY used when switching IARRAY rows during sort.
      INTEGER(LONG)                   :: IFLIP             ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT               ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM              ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                 ! An array index
      INTEGER(LONG)                   :: SORTPK            ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM          ! How many times the sort has to be performed in order for the data


      REAL(DOUBLE),  INTENT(INOUT)    :: RARRAY(NSIZE,3)   ! Array of real values
      REAL(DOUBLE)                    :: RDUM              ! Dummy values in RARRAY used when switching RARRAY rows during the sort.



! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (IARRAY(M) > IARRAY(N)) THEN
                     IDUM      = IARRAY(M)
                     IARRAY(M) = IARRAY(N)
                     IARRAY(N) = IDUM
                     DO I=1,3
                        RDUM        = RARRAY(M,I)
                        RARRAY(M,I) = RARRAY(N,I)
                        RARRAY(N,I) = RDUM
                     ENDDO
                     IFLIP = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         SORTED = 'YES'                                    ! Make sure array is sorted.
chk_sort:DO I=1,NSIZE-1
            IF (IARRAY(I) > IARRAY(I+1)) THEN
            SORTED = 'NO '
            EXIT chk_sort
            ENDIF
         ENDDO chk_sort

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                              ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')

! **********************************************************************************************************************************

      END SUBROUTINE SORT_INT1_REAL3


      SUBROUTINE SORT_INT2 ( CALLING_SUBR, MESSAG, NSIZE, IARRAY1, IARRAY2 )

! Performs shell sort on integer array IARRAY1 to put it into numerically increasing order and sorts integer array
! IARRAY2 along with IARRAY1. Both arrays are of size NSIZE

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_INT2'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to be written out if this subr fails to sort
      CHARACTER( 3*BYTE)              :: SORTED            ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(IN)       :: NSIZE             ! No. rows in arrays IARRAY, RARRAY
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY1(NSIZE)    ! Integer array
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY2(NSIZE)    ! Integer array
      INTEGER(LONG)                   :: I,K,M             ! DO loop indices
      INTEGER(LONG)                   :: IDUM1,IDUM2       ! Dummy values in IARRAY used when switching IARRAY rows during sort.
      INTEGER(LONG)                   :: IFLIP             ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT               ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM              ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                 ! An array index
      INTEGER(LONG)                   :: SORTPK            ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM          ! How many times the sort has to be performed in order for the data




! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (IARRAY1(M) > IARRAY1(N)) THEN
                     IDUM1      = IARRAY1(M)
                     IDUM2      = IARRAY2(M)
                     IARRAY1(M) = IARRAY1(N)
                     IARRAY2(M) = IARRAY2(N)
                     IARRAY1(N) = IDUM1
                     IARRAY2(N) = IDUM2
                     IFLIP = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         SORTED = 'YES'                                    ! Make sure array is sorted.
chk_sort:DO I=1,NSIZE-1
            IF (IARRAY1(I) > IARRAY1(I+1)) THEN
            SORTED = 'NO '
            EXIT chk_sort
            ENDIF
         ENDDO chk_sort

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                              ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')

! **********************************************************************************************************************************

      END SUBROUTINE SORT_INT2


      SUBROUTINE SORT_INT2_REAL1 ( CALLING_SUBR, MESSAG, NSIZE, IARRAY1, IARRAY2, RARRAY )

! Performs shell sort on integer array IARRAY1 to put it into numerically increasing order and sorts integer array IARRAY2 and real
! array RARRAY along with IARRAY1. All arrays are of size NSIZE

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_INT2_REAL1'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to be written out if this subr fails to sort
      CHARACTER( 3*BYTE)              :: SORTED            ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(IN)       :: NSIZE             ! No. rows in arrays IARRAY, RARRAY
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY1(NSIZE)    ! Array of integer values
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY2(NSIZE)    ! Array of integer values
      INTEGER(LONG)                   :: I,K,M             ! DO loop indices
      INTEGER(LONG)                   :: IDUM1,IDUM2       ! Dummy values in IARRAY used when switching IARRAY rows during sort.
      INTEGER(LONG)                   :: IFLIP             ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT               ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM              ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                 ! An array index
      INTEGER(LONG)                   :: SORTPK            ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM          ! How many times the sort has to be performed in order for the data


      REAL(DOUBLE),  INTENT(INOUT)    :: RARRAY(NSIZE)     ! Array of real values
      REAL(DOUBLE)                    :: RDUM              ! Dummy values in RARRAY used when switching RARRAY rows during the sort



! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (IARRAY1(M) > IARRAY1(N)) THEN
                     IDUM1       = IARRAY1(M)
                     IDUM2       = IARRAY2(M)
                     RDUM        = RARRAY(M)
                     IARRAY1(M)  = IARRAY1(N)
                     IARRAY2(M)  = IARRAY2(N)
                     RARRAY(M)   = RARRAY(N)
                     IARRAY1(N)  = IDUM1
                     IARRAY2(N)  = IDUM2
                     RARRAY(N)   = RDUM
                     IFLIP       = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         SORTED = 'YES'                                    ! Make sure array is sorted.
chk_sort:DO I=1,NSIZE-1
            IF (IARRAY1(I) > IARRAY1(I+1)) THEN
            SORTED = 'NO '
            EXIT chk_sort
            ENDIF
         ENDDO chk_sort

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                              ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')

! **********************************************************************************************************************************

      END SUBROUTINE SORT_INT2_REAL1


      SUBROUTINE SORT_INT3 ( CALLING_SUBR, MESSAG, NSIZE, IARRAY1, IARRAY2, IARRAY3 )

! Performs shell sort on integer array IARRAY1 to put it into numerically increasing order and sorts integer arrays
! IARRAY2 and IARRAY3 along with IARRAY1. All 3 arrays are of size NSIZE

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_INT3'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to be written out if this subr fails to sort
      CHARACTER( 3*BYTE)              :: SORTED            ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(IN)       :: NSIZE             ! No. rows in arrays IARRAY, RARRAY
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY1(NSIZE)    ! Integer array
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY2(NSIZE)    ! Integer array
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY3(NSIZE)    ! Integer array
      INTEGER(LONG)                   :: I,K,M             ! DO loop indices
      INTEGER(LONG)                   :: IDUM1,IDUM2,IDUM3 ! Dummy values in IARRAY used when switching IARRAY rows during sort.
      INTEGER(LONG)                   :: IFLIP             ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT               ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM              ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                 ! An array index
      INTEGER(LONG)                   :: SORTPK            ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM          ! How many times the sort has to be performed in order for the data




! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (IARRAY1(M) > IARRAY1(N)) THEN
                     IDUM1      = IARRAY1(M)
                     IDUM2      = IARRAY2(M)
                     IDUM3      = IARRAY3(M)
                     IARRAY1(M) = IARRAY1(N)
                     IARRAY2(M) = IARRAY2(N)
                     IARRAY3(M) = IARRAY3(N)
                     IARRAY1(N) = IDUM1
                     IARRAY2(N) = IDUM2
                     IARRAY3(N) = IDUM3
                     IFLIP = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         SORTED = 'YES'                                    ! Make sure array is sorted.
chk_sort:DO I=1,NSIZE-1
            IF (IARRAY1(I) > IARRAY1(I+1)) THEN
            SORTED = 'NO '
            EXIT chk_sort
            ENDIF
         ENDDO chk_sort

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                              ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')

! **********************************************************************************************************************************

      END SUBROUTINE SORT_INT3


      SUBROUTINE SORT_INT3_CHAR2 ( CALLING_SUBR, MESSAG, NSIZE, IARRAY1, IARRAY2, IARRAY3, CARRAY1, CARRAY2 )

! Performs shell sort on integer array IARRAY1 to put it into numerically increasing order and sorts integer arrays
! IARRAY2, IARRAY3 and character arrays CARRAY1, CARRAY2 along with IARRAY1. All arrays are of size NSIZE

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_INT3_CHAR2'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr
      INTEGER(LONG), INTENT(IN)       :: NSIZE             ! No. rows in arrays IARRAY, RARRAY
      CHARACTER(LEN=*), INTENT(INOUT) :: CARRAY1(NSIZE)    ! Character array
      CHARACTER(LEN=*), INTENT(INOUT) :: CARRAY2(NSIZE)    ! Character array
      CHARACTER(LEN=LEN(CARRAY1))     :: CDUM1             ! Dummy values in CARRAY1 used when switching CARRAY1 rows during sort.
      CHARACTER(LEN=LEN(CARRAY2))     :: CDUM2             ! Dummy values in CARRAY2 used when switching CARRAY2 rows during sort.
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to be written out if this subr fails to sort
      CHARACTER(3*BYTE)               :: SORTED            ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(INOUT)    :: IARRAY1(NSIZE)    ! Integer array
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY2(NSIZE)    ! Integer array
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY3(NSIZE)    ! Integer array
      INTEGER(LONG)                   :: I,K,M             ! DO loop indices
      INTEGER(LONG)                   :: IDUM1,IDUM2,IDUM3 ! Dummy values in IARRAY used when switching IARRAY rows during sort.
      INTEGER(LONG)                   :: IFLIP             ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT               ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM              ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                 ! An array index
      INTEGER(LONG)                   :: SORTPK            ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM          ! How many times the sort has to be performed in order for the data




! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (IARRAY1(M) > IARRAY1(N)) THEN
                     IDUM1      = IARRAY1(M)
                     IDUM2      = IARRAY2(M)
                     IDUM3      = IARRAY3(M)
                     CDUM1      = CARRAY1(M)
                     CDUM2      = CARRAY2(M)
                     IARRAY1(M) = IARRAY1(N)
                     IARRAY2(M) = IARRAY2(N)
                     IARRAY3(M) = IARRAY3(N)
                     CARRAY1(M)  = CARRAY1(N)
                     IARRAY1(N) = IDUM1
                     IARRAY2(N) = IDUM2
                     IARRAY3(N) = IDUM3
                     CARRAY1(N) = CDUM1
                     CARRAY2(N) = CDUM2
                     IFLIP = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                 CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         SORTED = 'YES'                                    ! Make sure array is sorted.
chk_sort:DO I=1,NSIZE-1
            IF (IARRAY1(I) > IARRAY1(I+1)) THEN
            SORTED = 'NO '
            EXIT chk_sort
            ENDIF
         ENDDO chk_sort

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                              ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')

! **********************************************************************************************************************************

      END SUBROUTINE SORT_INT3_CHAR2


      SUBROUTINE SORT_REAL1_INT1 ( CALLING_SUBR, MESSAG, NSIZE, RARRAY, IARRAY )

! Performs shell sort on real array RARRAY to put it into numerically increasing order and sorts integer array IARRAY
! along with RARRAY. Both arrays are of size NSIZE

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_REAL1_INT1'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to be written out if this subr fails to sort
      CHARACTER( 3*BYTE)              :: SORTED            ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(IN)       :: NSIZE             ! No. rows in arrays IARRAY, RARRAY
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY(NSIZE)     ! Integer array
      INTEGER(LONG)                   :: I,K,M             ! DO loop indices
      INTEGER(LONG)                   :: IDUM              ! Dummy values in IARRAY used when switching IARRAY rows during sort.
      INTEGER(LONG)                   :: IFLIP             ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT               ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM              ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                 ! An array index
      INTEGER(LONG)                   :: SORTPK            ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM          ! How many times the sort has to be performed in order for the data


      REAL(DOUBLE),  INTENT(INOUT)    :: RARRAY(NSIZE)     ! Array of real values
      REAL(DOUBLE)                    :: RDUM              ! Dummy values in RARRAY used when switching RARRAY rows during sort.



! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (RARRAY(M) > RARRAY(N)) THEN
                     RDUM      = RARRAY(M)
                     IDUM      = IARRAY(M)
                     RARRAY(M) = RARRAY(N)
                     IARRAY(M) = IARRAY(N)
                     RARRAY(N) = RDUM
                     IARRAY(N) = IDUM
                     IFLIP = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         SORTED = 'YES'                                    ! Make sure array is sorted.
chk_sort:DO I=1,NSIZE-1
            IF (RARRAY(I) > RARRAY(I+1)) THEN
            SORTED = 'NO '
            EXIT chk_sort
            ENDIF
         ENDDO chk_sort

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                              ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')

! **********************************************************************************************************************************

      END SUBROUTINE SORT_REAL1_INT1


      SUBROUTINE SORT_GRID_RGRID ( CALLING_SUBR, MESSAG, NSIZE, IARRAY, RARRAY )

! Performs shell sort on integer IARRAY (GRID) and real RARRAY (RGRID) to put them into an order where the GRID ID
! (in column 1 of IARRAY) is in numerically increasing order

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MGRID, MRGRID
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_GRID_RGRID'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR        ! Subr that called this subr
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG              ! Message to be written out if this subr fails to sort
      CHARACTER( 3*BYTE)              :: SORTED              ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(IN)       :: NSIZE               ! No. rows in arrays IARRAY, RARRAY
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY(NSIZE,MGRID) ! Array GRID
      INTEGER(LONG)                   :: I,K,M               ! DO loop indices
      INTEGER(LONG)                   :: IDUM1               ! Dummy values in IARRAY used when switching IARRAY rows during sort
      INTEGER(LONG)                   :: IFLIP               ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT                 ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM                ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                   ! An array index
      INTEGER(LONG)                   :: SORTPK              ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM            ! How many times the sort has to be performed in order for the data
!                                                              to be in sort order. SORT_MAX is max value


      REAL(DOUBLE),  INTENT(INOUT)    :: RARRAY(NSIZE,MRGRID)! Array RGRID
      REAL(DOUBLE)                    :: RDUM1               ! Dummy values in RARRAY used when switching RARRAY rows during sort



! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (IARRAY(M,1) > IARRAY(N,1)) THEN
                     DO I=1,MGRID
                        IDUM1       = IARRAY(M,I)
                        IARRAY(M,I) = IARRAY(N,I)
                        IARRAY(N,I) = IDUM1
                     ENDDO
                     DO I=1,MRGRID
                        RDUM1       = RARRAY(M,I)
                        RARRAY(M,I) = RARRAY(N,I)
                        RARRAY(N,I) = RDUM1
                     ENDDO
                     IFLIP = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO

         SORTED = 'YES'                                    ! Make sure array is sorted
chk_sort:DO I=1,NSIZE-1
            IF (IARRAY(I,1) > IARRAY(I+1,1)) THEN
            SORTED = 'NO '
            EXIT chk_sort
            ENDIF
         ENDDO chk_sort

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')

! **********************************************************************************************************************************

      END SUBROUTINE SORT_GRID_RGRID


      SUBROUTINE SORT_TDOF ( CALLING_SUBR, MESSAG, NSIZE, IARRAY, ICOL )

! Performs shell sort on 2D integer array IARRAY1 such that column ICOL is in numerically increasing order. IARRAY
! has MTDOF columns (which is the number of columns in array TDOF)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, MTDOF
      USE PARAMS, ONLY                :  SORT_MAX
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'SORT_TDOF'
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this subr
      CHARACTER(LEN=*), INTENT(IN)    :: MESSAG            ! Message to be written out if this subr fails to sort
      CHARACTER( 3*BYTE)              :: SORTED              ! = 'YES' if array is sorted

      INTEGER(LONG), INTENT(IN)       :: NSIZE               ! No. rows in arrays IARRAY, RARRAY
      INTEGER(LONG), INTENT(INOUT)    :: IARRAY(NSIZE,MTDOF) ! Integer array
      INTEGER(LONG)                   :: I,J,K,M             ! DO loop indices
      INTEGER(LONG), INTENT(IN)       :: ICOL                ! Col ICOL will be the col in numerical order after sort
      INTEGER(LONG)                   :: IDUM(MTDOF)         ! Dummy values in IARRAY used when switching IARRAY rows during sort.
      INTEGER(LONG)                   :: IFLIP               ! Indicates whether two values have been switched in sort order.
      INTEGER(LONG)                   :: JCT                 ! Shell sort parameter returned from subroutine SORTLEN.
      INTEGER(LONG)                   :: MAXM                ! NSIZE - SORTPK
      INTEGER(LONG)                   :: N                   ! An array index
      INTEGER(LONG)                   :: SORTPK              ! Intermediate variable used in setting a DO loop range.
      INTEGER(LONG)                   :: SORT_NUM            ! How many times the sort has to be performed in order for the data




! **********************************************************************************************************************************
! Call SORTLEN to calculate the shell sort parameter JCT

      SORT_NUM = 1

      CALL SORTLEN ( NSIZE, JCT )

outer:DO                                                      ! Run sort until array is sorted or SORT_MAX is exceeded

         DO K = JCT,1,-1                                      ! Do the sort
            SORTPK = 2**K - 1
            MAXM  = NSIZE - SORTPK
            IFLIP = 0
            DO
               DO M = 1,MAXM
                  N = M + SORTPK
                  IF (IARRAY(M,ICOL) > IARRAY(N,ICOL)) THEN
                     DO J = 1,MTDOF
                        IDUM(J)     = IARRAY(M,J)
                        IARRAY(M,J) = IARRAY(N,J)
                        IARRAY(N,J) = IDUM(J)
                     ENDDO
                     IFLIP = 1
                  ENDIF
               ENDDO
               IF (IFLIP == 1) THEN
                  IFLIP = 0
                  CYCLE
               ELSE
                  EXIT
               ENDIF
            ENDDO
         ENDDO


         SORTED = 'YES'                                    ! Make sure array is sorted.
chk_sort:DO I=1,NSIZE-1
            IF (IARRAY(I,ICOL) > IARRAY(I+1,ICOL)) THEN
            SORTED = 'NO '
            EXIT chk_sort
            ENDIF
         ENDDO chk_sort

         IF (SORTED == 'YES') THEN                         ! If array is sorted, exit outer loop
            EXIT outer
         ELSE                                              ! If array is not sorted, sort again if SORT_MAX num times not exceeded
            SORT_NUM = SORT_NUM + 1
            IF (SORT_NUM <= SORT_MAX) THEN
               CYCLE outer
            ELSE
               WRITE(ERR, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               WRITE(F06, 914) SUBR_NAME, CALLING_SUBR, SORT_MAX, MESSAG
               FATAL_ERR = FATAL_ERR + 1
               CALL OUTA_HERE ( 'Y' )                      ! Can't get array sorted after trying SORT_MAX times, so quit
            ENDIF
         ENDIF

      ENDDO outer



      RETURN

! **********************************************************************************************************************************
  914 FORMAT(' *ERROR   914: SUBROUTINE ',A,', CALLED BY SUBR ',A                                                                  &
                    ,/,14X,' HAS MADE ',I8,' UNSUCCESSFUL ATTEMPTS TO SORT ARRAY(S) ',A                                            &
                    ,/,14X,' THE MAX NUMBER OF SORT ATTEMPTS CAN BE INCREASED WITH BULK DATA PARAM SORT_MAX')


! **********************************************************************************************************************************

      END SUBROUTINE SORT_TDOF


      SUBROUTINE CALC_VEC_SORT_ORDER ( VEC, SORT_ORDER, SORT_INDICES )

! Determines the order of the 3 components of a vector if they were arranged from lowest value to largest value. The values are not
! actually sorted.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_VEC_SORT_ORDER'
      CHARACTER( 5*BYTE), INTENT(OUT) :: SORT_ORDER        ! Order in which the VX(i) have been sorted. If none of the tests below
!                                                            are satisfied, SORT_ORDER is returned as null

      INTEGER(LONG), INTENT(OUT)      :: SORT_INDICES(3)   ! Indices of VEC in the order from lowest value component to highest


      REAL(DOUBLE), INTENT(IN)        :: VEC(3)            ! A 3 component vector



! **********************************************************************************************************************************
      SORT_ORDER = '     '

      IF ((DABS(VEC(1)) <= DABS(VEC(2))) .AND. (DABS(VEC(1)) <= DABS(VEC(3))) .AND. (DABS(VEC(2)) <= DABS(VEC(3)))) THEN
         SORT_ORDER      = '1 2 3'
         SORT_INDICES(1) = 1
         SORT_INDICES(2) = 2
         SORT_INDICES(3) = 3
      ENDIF

      IF ((DABS(VEC(1)) <= DABS(VEC(2))) .AND. (DABS(VEC(1)) <= DABS(VEC(3))) .AND. (DABS(VEC(3)) <= DABS(VEC(2)))) THEN
         SORT_ORDER      = '1 3 2'
         SORT_INDICES(1) = 1
         SORT_INDICES(2) = 3
         SORT_INDICES(3) = 2
      ENDIF

      IF ((DABS(VEC(2)) <= DABS(VEC(3))) .AND. (DABS(VEC(2)) <= DABS(VEC(1))) .AND. (DABS(VEC(1)) <= DABS(VEC(3)))) THEN
         SORT_ORDER      = '2 1 3'
         SORT_INDICES(1) = 2
         SORT_INDICES(2) = 1
         SORT_INDICES(3) = 3
      ENDIF

      IF ((DABS(VEC(2)) <= DABS(VEC(3))) .AND. (DABS(VEC(2)) <= DABS(VEC(1))) .AND. (DABS(VEC(3)) <= DABS(VEC(1)))) THEN
         SORT_ORDER      = '2 3 1'
         SORT_INDICES(1) = 2
         SORT_INDICES(2) = 3
         SORT_INDICES(3) = 1
      ENDIF

      IF ((DABS(VEC(3)) <= DABS(VEC(1))) .AND. (DABS(VEC(3)) <= DABS(VEC(2))) .AND. (DABS(VEC(1)) <= DABS(VEC(2)))) THEN
         SORT_ORDER      = '3 1 2'
         SORT_INDICES(1) = 3
         SORT_INDICES(2) = 1
         SORT_INDICES(3) = 2
      ENDIF

      IF ((DABS(VEC(3)) <= DABS(VEC(1))) .AND. (DABS(VEC(3)) <= DABS(VEC(2))) .AND. (DABS(VEC(2)) <= DABS(VEC(1)))) THEN
         SORT_ORDER      = '3 2 1'
         SORT_INDICES(1) = 3
         SORT_INDICES(2) = 2
         SORT_INDICES(3) = 1
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE CALC_VEC_SORT_ORDER


   END MODULE SORTING
