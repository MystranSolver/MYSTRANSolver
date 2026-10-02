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

   MODULE DOF_ARRAY_INDEXING

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ARRAY_SIZE_ERROR_1, CALC_TDOF_ROW_START, GET_ARRAY_ROW_NUM, GET_GRID_NUM_COMPS, GET_I_MAT_FROM_I2_MAT, GET_UG_123_IN_GRD_ORD

   CONTAINS

      SUBROUTINE ARRAY_SIZE_ERROR_1 ( INP_SUBR_NAME, NTERM_VAL, MATIN_NAME )

! Print error and quit when a subr tries to exceed allocated number of terms when storing/retrieving terms in an array

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ARRAY_SIZE_ERROR_1'
      CHARACTER(LEN=*), INTENT(IN)    :: INP_SUBR_NAME     ! Subroutine in which the error was detected
      CHARACTER(LEN=*), INTENT(IN)    :: MATIN_NAME        ! Name of matrix (for output message purposes)

      INTEGER(LONG), INTENT(IN)       :: NTERM_VAL         ! Size of the array that was exceeded




! **********************************************************************************************************************************
      WRITE(ERR,937) INP_SUBR_NAME, NTERM_VAL, MATIN_NAME
      WRITE(F06,937) INP_SUBR_NAME, NTERM_VAL, MATIN_NAME
      FATAL_ERR = FATAL_ERR + 1
      CALL OUTA_HERE ( 'Y' )                               ! Coding error (attempt to exceed allocated array size), so quit



      RETURN

! **********************************************************************************************************************************
  937 FORMAT(' *ERROR   937: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' ATTEMPT TO EXCEED ALLOCATED ARRAY SIZE = ',I12                                                        &
                    ,/,14X,' WHILE STORING OR RETRIEVING TERMS IN ARRAY(S) = ',A)

! **********************************************************************************************************************************

      END SUBROUTINE ARRAY_SIZE_ERROR_1


      SUBROUTINE CALC_TDOF_ROW_START ( PRTDEB )

! Calculates the row number in array TDOF where a grid (or SPOINT) DOF data is to begin

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, NGRID
      USE IOUNT1, ONLY                :  ERR, F06
      USE TIMDAT, ONLY                :  TSEC
      USE DOF_TABLES, ONLY            :  TDOF_ROW_START
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MODEL_STUF, ONLY            :  GRID_ID

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'CALC_TDOF_ROW_START'
      CHARACTER(LEN=*), INTENT(IN)    :: PRTDEB            ! If 'Y' then print debug info if DEBUG(183) also > 0

      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: NUM_COMPS         ! Number of displ components (1 for SPOINT, 6 for physical grid)




! **********************************************************************************************************************************
      IF ((DEBUG(183) > 0) .AND. (PRTDEB == 'Y')) THEN
         WRITE(F06,88677)
      ENDIF

      TDOF_ROW_START(1) = 1
      DO I=2,NGRID
         CALL GET_GRID_NUM_COMPS ( I-1, NUM_COMPS, SUBR_NAME )
         TDOF_ROW_START(I) = TDOF_ROW_START(I-1) + NUM_COMPS
      ENDDO

      IF ((DEBUG(183) > 0) .AND. (PRTDEB == 'Y')) THEN
         DO I=1,NGRID
            WRITE(F06,88678) I, GRID_ID(I), NUM_COMPS, TDOF_ROW_START(I)
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************
88677 FORMAT(47x,'GRID_ID(I)      NUM_COMPS   TDOF_ROW_START')

88678 FORMAT(41X,3I14)

! **********************************************************************************************************************************

      END SUBROUTINE CALC_TDOF_ROW_START


      SUBROUTINE GET_ARRAY_ROW_NUM ( ARRAY_NAME, CALLING_SUBR, ASIZE, ARRAY, EXT_ID, ROW_NUM )

! Searches integer array ARRAY to find an external (actual) ID (EXT_ID) in order to find the row number where it exists.
! If EXT_ID is not found, ROW_NUM is set to -1 to indicate an error. ARRAY must be sorted into numerical order

! The algorithm uses HI and LO to bound the range in ARRAY where EXT_ID may be found. Initially, HI is set to the size
! of ARRAY and LO is set at 0 and the first estimate of the location of EXT_ID is near the middle of the range at:

!                              N = (HI + LO +1)/2                 (1)

! The algorithm then iterates until EXT_ID = ARRAY(N) by modifying HI and LO as follows:
!   (a) If EXT_ID < ARRAY(N) then HI is lowered   to N and a new N is calculated from (1) and the procedure repeated
!   (b) If EXT_ID > ARRAY(N) then LO is increased to N and a new N is calculated from (1) and the procedure repeated

!                          TMP_N = (TMP_HI + TMP_LO + 1.D0)/2.D0  (2)


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, f06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ONE, TWO

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_ARRAY_ROW_NUM'
      CHARACTER(LEN=*), INTENT(IN)    :: ARRAY_NAME        ! Name of array to be searched
      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Name of subr that called this one

      INTEGER(LONG), INTENT(IN)       :: ASIZE             ! Size of ARRAY
      INTEGER(LONG), INTENT(IN)       :: ARRAY(ASIZE)      ! Array to search
      INTEGER(LONG), INTENT(IN)       :: EXT_ID            ! External (actual) ID to find in ARRAY
      INTEGER(LONG), INTENT(OUT)      :: ROW_NUM           ! Internal ID (row in ARRAY) where EXT_ID exists
      INTEGER(LONG)                   :: HI, LO            ! Used to bound the range of N where EXT_ID is expected to be found
      INTEGER(LONG)                   :: LAST              ! Previous value of N in the search
      INTEGER(LONG)                   :: N                 ! When the search is completed, N is the ROW_NUM we ara looking for


      INTEGER(LONG)                   :: TMP_N             ! Real value of (DBL_HI + DBL_LO + 1.D0)/2.D0
      INTEGER(LONG)                   :: TMP_HI            ! Real value of HI
      INTEGER(LONG)                   :: TMP_LO            ! Real value of LO



! **********************************************************************************************************************************


! Check to see if our binary search will have an integer overflow
! Pick a more appropriate version if so, or raise an error
     IF ((ASIZE+1) > (HUGE(TMP_N)-1)/2) THEN
        ! use the double version I guess? Just error for now.
        FATAL_ERR = FATAL_ERR + 1
        WRITE(ERR,9003) CALLING_SUBR, ARRAY_NAME
        WRITE(F06,9003) CALLING_SUBR, ARRAY_NAME
        CALL OUTA_HERE ( 'Y' )
     ENDIF

! Initialize outputs

      ROW_NUM = 0

! Calc outputs

      HI     = ASIZE
      LO     = 0
      TMP_HI = HI
      TMP_LO = LO
      LAST   = 0

      DO
         N = (TMP_HI + TMP_LO + 1)/2
      !    N     = FLOOR(DBL_N)
         IF (N == LAST) THEN
            ROW_NUM = -1
            RETURN
         ENDIF
         LAST = N
         IF      (EXT_ID <  ARRAY(N)) THEN
           HI     = N
           TMP_HI = HI
           CYCLE
         ELSE IF (EXT_ID >  ARRAY(N)) THEN
           LO     = N
           TMP_LO = LO
           CYCLE
         ELSE IF (EXT_ID == ARRAY(N)) THEN
           EXIT
         ENDIF
      ENDDO

      ROW_NUM = N



      RETURN


 9003 FORMAT(' *ERROR   9003: ERROR IN SUBROUTINE ',A                                                                   &
      ,/,14X,' INPUT ARRAY ',A,' HAS TOO MANY ELEMENTS AND WILL OVERFLOW INTEGER(LONG) INDEX, MAX SUPPORTED ELEMS: 1073741823')


! **********************************************************************************************************************************

      END SUBROUTINE GET_ARRAY_ROW_NUM




      SUBROUTINE GET_GRID_NUM_COMPS ( IGRID, NUM_COMPS, CALLING_SUBR )

! Gets the number of components for a "grid" number from array GRID by testing GRID(IGRID,6)
! If GRID(I,6) is 6 then this is a physical grid with 6 comps of displ and if GRID(I,6) is 1 then this is an SPOINT.

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE MODEL_STUF, ONLY            :  GRID

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this one

      INTEGER(LONG), INTENT(IN)       :: IGRID             ! An internal grid number
      INTEGER(LONG), INTENT(OUT)      :: NUM_COMPS         ! 6 if GRID_NUM is an physical grid, 1 if an SPOINT

! **********************************************************************************************************************************

      NUM_COMPS = GRID(IGRID, 6)

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE GET_GRID_NUM_COMPS


      SUBROUTINE GET_I2_MAT_FROM_I_MAT ( MAT_NAME, NROWS, NTERMS, I_MAT, I2_MAT )

! This subr does the inverse of subr GET_I_MAT_FROM_I2_MAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_I2_MAT_FROM_I_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME          ! Matrix name

      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in MAT
      INTEGER(LONG), INTENT(IN)       :: NTERMS            ! Number of matrix terms that should be in MAT
      INTEGER(LONG), INTENT(IN)       :: I_MAT(NROWS+1)    ! Row indicators for terms in matrix MAT
      INTEGER(LONG), INTENT(OUT)      :: I2_MAT(NTERMS)    ! Row numbers for terms in matrix MAT
      INTEGER(LONG)                   :: I,J,K             ! DO loop indices or counters
      INTEGER(LONG)                   :: NUM_IN_ROW_I      ! Number of nonzero terms in row I




! **********************************************************************************************************************************
! Initialize

      DO I=1,NTERMS
         I2_MAT(I) = 0
      ENDDO

! Get I2_MAT

      K = 0
      IF (NTERMS > 0) THEN
         DO I=1,NROWS
            NUM_IN_ROW_I = I_MAT(I+1) - I_MAT(I)
            DO J=1,NUM_IN_ROW_I
               K = K + 1
               IF (K > NTERMS)  CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, K, MAT_NAME )
               I2_MAT(K) = I
            ENDDO
         ENDDO
      ENDIF



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE GET_I2_MAT_FROM_I_MAT



      SUBROUTINE GET_I_MAT_FROM_I2_MAT ( MAT_NAME, NROWS, NTERMS, I2_MAT, I_MAT )

! I_MAT is the sparse compressed row format for sparse array MAT. I2_MAT has a row number for each of the terms in sparse MAT.
! I2_MAT is generally used when some sparse matrices are written to a file since it then has a row, col, value for all terms in the
! sparse MAT. I_MAT is usually the way the row numbers are stored in memory for sparse matrices.
! This subr creates I_MAT from a given I2_MAT

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR
      USE TIMDAT, ONLY                :  TSEC

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_I_MAT_FROM_I2_MAT'
      CHARACTER(LEN=*), INTENT(IN)    :: MAT_NAME          ! Matrix name

      INTEGER(LONG), INTENT(IN)       :: NROWS             ! Number of rows in MAT
      INTEGER(LONG), INTENT(IN)       :: NTERMS            ! Number of matrix terms that should be in MAT
      INTEGER(LONG), INTENT(IN)       :: I2_MAT(NTERMS)    ! Row numbers for terms in matrix MAT
      INTEGER(LONG), INTENT(OUT)      :: I_MAT(NROWS+1)    ! Row numbers for terms in matrix MAT
      INTEGER(LONG)                   :: I,K               ! DO loop indices or counters
      INTEGER(LONG)                   :: IROW              ! Integer row value read from FILNAM
      INTEGER(LONG)                   :: IROW_OLD          ! Previous value of IROW
      INTEGER(LONG)                   :: KTERM             ! Count of number of nonzero terms read from FILNAM
      INTEGER(LONG)                   :: MAT_ERR    = 0    ! Error count




! **********************************************************************************************************************************
      IF (NTERMS > 0) THEN

         KTERM    = 0
         IROW_OLD = 0
         I_MAT(1) = 1

k_do1:   DO K = 1,NTERMS
            IROW = I2_MAT(K)
            IF (IROW > IROW_OLD) THEN
               DO I=IROW_OLD+1,IROW
                  I_MAT(I+1) = I_MAT(I)
               ENDDO
               IROW_OLD = IROW
            ELSE IF (IROW < IROW_OLD) THEN
               MAT_ERR   = MAT_ERR + 1
               FATAL_ERR = FATAL_ERR + 1
               WRITE(ERR,926) SUBR_NAME, MAT_NAME, IROW_OLD, IROW
               WRITE(F06,926) SUBR_NAME, MAT_NAME, IROW_OLD, IROW
               CYCLE k_do1
            ENDIF
            I_MAT(IROW+1) = I_MAT(IROW+1) + 1
            KTERM         = KTERM + 1
         ENDDO k_do1

         IF (MAT_ERR /= 0) THEN
            WRITE(ERR,9996) SUBR_NAME,MAT_ERR
            WRITE(F06,9996) SUBR_NAME,MAT_ERR
            CALL OUTA_HERE ( 'Y' )                         ! Quit due to above errors in k_do1 loop
         ENDIF

         IF (IROW < NROWS) THEN                            ! Fill out remainder of I_MAT, if needed
            DO I=IROW+1,NROWS
               I_MAT(I+1) = I_MAT(I)
            ENDDO
         ENDIF

      ELSE                                                 ! MAT is null so set I_MAT

         DO I=1,NROWS+1
            I_MAT(I) = 1
         ENDDO

      ENDIF



      RETURN

! **********************************************************************************************************************************
  926 FORMAT(' *ERROR   926: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' MATRIX ',A,' IS STORED INCORRECTLY. THE MATRIX MUST BE STORED BY ROWS IN NUMERICAL ORDER '            &
                    ,/,14X,' HOWEVER, ROW ',I12,' IS STORED BEFORE ROW ',I12)

 9996 FORMAT(/,' PROCESSING ABORTED IN SUBROUTINE ',A,' DUE TO ABOVE ',I8,' ERRORS')

! **********************************************************************************************************************************

      END SUBROUTINE GET_I_MAT_FROM_I2_MAT



      SUBROUTINE GET_UG_123_IN_GRD_ORD ( IERR )

! Gets the 3 translation (T1, T2, T3) displ values fro the G-set displ vector (UG_COL) for 1 subcase and puts those values into an
! NGRID x 3 array where T1 is col1, T2 is col 2 and T3 is col3

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, FATAL_ERR, NDOFG, NGRID
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  GRID_ID
      USE DOF_TABLES, ONLY            :  TDOFI
      USE COL_VECS, ONLY              :  UG_COL
      USE MISC_MATRICES, ONLY         :  UG_T123_MAT

      USE FILE_LIFECYCLE, ONLY   :  OUTA_HERE
      USE SORTING, ONLY               :  SORT_INT1_REAL3

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'GET_UG_123_IN_GRD_ORD'


      INTEGER(LONG), INTENT(OUT)      :: IERR              ! Local error indicator
      INTEGER(LONG)                   :: GRDS_GLOBAL(NGRID)!
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: IGRID             ! Count of grids (1 to NGRID)
      INTEGER(LONG)                   :: IDOFG             ! Count of G-set DOF's  (1 to NDOFG)
      INTEGER(LONG)                   :: NUM_COMPS         ! 6 if GRID_NUM is an physical grid, 1 if an SPOINT


! **********************************************************************************************************************************
      IERR = 0

! Put the grid numbers that are in TDOFI (which are in global order) into array GRDS_GLOBAL

      IGRID = 0
      IDOFG = 0
      DO I = 1,NGRID
         CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
         DO J = 1,NUM_COMPS
            IDOFG = IDOFG + 1
            IF (IDOFG > NDOFG) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, IDOFG, 'TDOFI' )
            IF (J == 1) THEN
               IGRID = IGRID + 1
               IF (IGRID > NGRID) CALL ARRAY_SIZE_ERROR_1 ( SUBR_NAME, IGRID, 'GRDS_GLOBAL' )
               GRDS_GLOBAL(IGRID) = TDOFI(IDOFG,1)
            ENDIF
         ENDDO
      ENDDO
! Make sure IGRID = NGRID and IDOFG = NDOFG. Otherwise something is wrong in the code above. That code was checked to make sure
! IGRID did not attempt to exceed NGRID and IDOFG did not attempt to exceed NDOFG. Now we make sure they are equal.

      IF (IGRID /= NGRID) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,966) SUBR_NAME, 'IGRID', IGRID, 'NGRID', NGRID
         WRITE(F06,966) SUBR_NAME, 'IGRID', IGRID, 'NGRID', NGRID
         CALL OUTA_HERE ( 'Y' )
      ENDIF

      IF (IDOFG /= NDOFG) THEN
         FATAL_ERR = FATAL_ERR + 1
         WRITE(ERR,966) SUBR_NAME, 'IDOFG', IDOFG, 'NDOFG', NDOFG
         WRITE(F06,966) SUBR_NAME, 'IDOFG', IDOFG, 'NDOFG', NDOFG
         CALL OUTA_HERE ( 'Y' )
      ENDIF

! Put the data from UG_COL into an array that has 3 cols (for each of the T1, T2, T3 translation freedoms) and rows equal to NDOFG.

      IDOFG = 0
      DO I = 1,NGRID
         CALL GET_GRID_NUM_COMPS ( I, NUM_COMPS, SUBR_NAME )
j_do2:   DO J = 1,NUM_COMPS
            IDOFG = IDOFG + 1
            IF (J <= 3) THEN                               ! We only want the 3 translations (or 1 if SPOINT)
               UG_T123_MAT(I,J) = UG_COL(IDOFG)
            ELSE
               CYCLE j_do2
            ENDIF
         ENDDO j_do2
      ENDDO
! Sort GRDS_GLOBAL and UG_T123_MAT so that GRDS_GLOBAL is in numerical grid order

      CALL SORT_INT1_REAL3 ( SUBR_NAME, 'GRDS_GLOBAL and UG_T123_MAT', NGRID, GRDS_GLOBAL, UG_T123_MAT )


      RETURN

! **********************************************************************************************************************************
  966 FORMAT(' *ERROR   966: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' VARIABLE ',A,' = ',I8,' BUT SHOULD BE THE SAME AS VARIABLE ',A,' = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE GET_UG_123_IN_GRD_ORD


   END MODULE DOF_ARRAY_INDEXING
