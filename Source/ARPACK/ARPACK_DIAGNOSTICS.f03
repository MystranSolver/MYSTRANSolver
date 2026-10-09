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

   MODULE ARPACK_DIAGNOSTICS

   IMPLICIT NONE

   PRIVATE

   PUBLIC :: ARPACK_INFO_MSG

   CONTAINS

      SUBROUTINE ARPACK_INFO_MSG ( SUBNAME, INFO, IPARAM, LWORKL, NEV, NCV )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  PROG_NAME, FATAL_ERR, NDOFL, WARN_ERR
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06
      USE PARAMS, ONLY                :  DARPACK, SUPWARN
      USE MODEL_STUF, ONLY            :  EIG_N2

      IMPLICIT NONE

      CHARACTER(  6*BYTE), INTENT(IN):: SUBNAME            ! Name of subr w/ msg
      CHARACTER(  5*BYTE)            :: ERROR_NUM          ! Character value of an error number
      CHARACTER(113*BYTE)            :: MSG(21,4)          ! Message written

      INTEGER, INTENT(IN)            :: INFO               ! Error msg number from ARPACK routine to do eigenvalue/vec calc
      INTEGER, INTENT(IN)            :: IPARAM(11)         ! Integer array from subr DSBAND (in module ARPACK_LANCZOS_1)
      INTEGER, INTENT(IN)            :: LWORKL             ! Length of array WORKL in subr dsband
      INTEGER, INTENT(IN)            :: NEV                !
      INTEGER, INTENT(IN)            :: NCV                !
      INTEGER                        :: I,J                ! DO loop indices

! **********************************************************************************************************************************
      DO I=1,21
         DO J=1,2
            MSG(I,J)(1:) = ' '
         ENDDO
      ENDDO

      MSG( 1,1) = 'NORMAL EXIT'
      MSG( 2,1) = 'MAX NUMBER OF ITERATIONS TAKEN. ALL POSSIBLE EIGENVALUES HAVE BEEN FOUND. THE NUMBER OF CONVERGED VALUES IS:'
      MSG( 2,2) = 'THE MAX NUMBER OF ITERATIONS CAN BE INCREASED THROUGH USE OF BULK DATA PARAM MXITERL ENTRY'
      MSG( 3,1) = 'NO SHIFTS COULD BE APPLIED DURING A CYCLE OF THE IMPLICITLY RESTARTED ARNOLDI ITERATION. ONE POSSIBILITY IS TO'
      MSG( 3,2) = 'INCREASE THE SIZE OF NCV RELATIVE TO NEV. MYSTRAN CALCS A VALUE TO USE FOR NCV = MIN(EIG_NCVFACL*NEV, N) WHERE:'
      MSG( 3,3) = 'NEV = NUMBER OF EIGENVALUES REQUESTED AND N IS THE ROW OR COLUMN SIZE OF THE MATRIX (NUMBER OF A-SET DOF''s)'
      MSG( 3,4) = 'TO INCREASE NCV, EIG_NCVFACL CAN BE CHANGED FROM ITS DEFAULT VALUE ON THE BULK DATA EIGRL ENTRY.'
      MSG( 4,1) = 'MATRIX SIZE N INPUT TO ARPACK SUBR DSBAND, MUST BE > 0.'
      MSG( 5,1) = 'NEV (THE NUMBER OF EIGENVALUES TO BE CALCULATED), INPUT TO ARPACK SUBR DSBAND, MUST BE > 0.'
      MSG( 6,1) = 'NCV, INPUT TO ARPACK SUBR DSBAND, MUST BE > NEV AND LESS THAN OR EQUAL TO NDOFA.'
      MSG( 7,1) = 'MXITER (THE MAX NUM OF ARNOLDI UPDATE ITERATIONS ALLOWED), INPUT TO ARPACK SUBR DSBAND, MUST BE > 0'
      MSG( 8,1) = '"WHICH", INPUT TO ARPACK SUBR DSBAND, MUST BE ONE OF LM, SM, LA, SA OR BE.'
      MSG( 9,1) = 'BMAT, INPUT TO ARPACK SUBR DSBAND, MUST BE ONE OF I OR G'
      MSG(10,1) = 'LWORKL (LEN OF WORK ARRAY WORKL), INPUT TO ARPACK SUBR DSBAND, IS NOT SUFFICIENT. MUST BE AT LEAST NCV*(NCV+8)'
      MSG(11,1) = 'THE ALGORITHM HAS FAILED TO FIND ALL THE EIGENVALUES IN A TOTAL OF 30*N ITERATIONS (N = MATRIX COL OR ROW SIZE)'
      MSG(12,1) = 'THE INITIAL RESIDUAL VECTOR FOR THE ARNOLDI (LANCZOS) PROCESS IS ZERO.'
      MSG(13,1) = 'IPARAM(7), INPUT TO ARPACK SUBR DSBAND, MUST BE 1,2,3,4,5.'
      MSG(14,1) = 'IPARAM(7) AND BMAT, INPUT TO ARPACK SUBR DSBAND, ARE INCOMPATABLE.'
      MSG(15,1) = 'IPARAM(1), INPUT TO ARPACK SUBR DSBAND, MUST BE EQUAL TO 0 OR 1.'
      MSG(16,1) = 'NEV AND WHICH = BE, INPUT TO ARPACK SUBR DSBAND, ARE INCOMPATABLE.'
      MSG(17,1) = 'HOWMNY, INPUT TO ARPACK SUBR DSBAND, MUST BE ONE OF A OR S IF RVEC = .TRUE.'
      MSG(18,1) = 'ARPACK SUBR DSEUPD DID NOT FIND ANY EIGENVALUES TO SUFFICIENT ACCURACY.'
      MSG(19,1) = 'HOWMNY, INPUT TO ARPACK SUBR DSBAND, = S NOT YET IMPLEMENTED.'
      MSG(20,1) = 'ARPACK SUBR DSEUPD GOT A DIFFERENT COUNT OF THE NUMBER OF CONVERGED RITZ VALUES THAN DSAUPD GOT. THIS INDICATES'
      MSG(20,2) = 'THE USER PROBABLY MADE AN ERROR IN PASSING DATA FROM DSAUPD TO DSEUPD OR THE DATA WAS MODIFIED BEFORE DSEUPD.'
      MSG(21,1) = 'ARPACK SUBR DSAITR COULD NOT BUILD AN ARNOLDI FACTORIZATION. IPARAM(5) RETURNS THE SIZE OF THE CURRENT ARNOLDI'
      MSG(21,2) = 'FACTORIZATION. THE USER IS ADVISED TO CHECK THAT ENOUGH WORKSPACE AND ARRAY STORAGE HAS BEEN ALLOCATED.'


      IF ((SUBNAME == 'DSAUPD') .OR. (SUBNAME == 'dsaupd')) THEN

         IF      (INFO ==     1) THEN                      ! MSG( 2)

            ERROR_NUM = ' 7101'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,702) ERROR_NUM, SUBNAME, INFO, MSG( 2,1)
            WRITE(ERR,1104) IPARAM(5)
            WRITE(ERR,1103) MSG( 2,2)
            WRITE(ERR,*)
            WRITE(F06,702) ERROR_NUM, SUBNAME, INFO, MSG( 2,1)
            WRITE(F06,1104) IPARAM(5)
            WRITE(F06,1103) MSG( 2,2)
            WRITE(F06,*)

         ELSE IF (INFO ==     3) THEN                      ! MSG( 3)

            ERROR_NUM = ' 7103'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,702) ERROR_NUM, SUBNAME, INFO, MSG( 3,1)
            WRITE(ERR,1103) MSG( 3,2)
            WRITE(ERR,1103) MSG( 3,3)
            WRITE(ERR,1103) MSG( 3,4)
            WRITE(ERR,*)
            WRITE(F06,702) ERROR_NUM, SUBNAME, INFO, MSG( 3,1)
            WRITE(F06,1103) MSG( 3,2)
            WRITE(F06,1103) MSG( 3,3)
            WRITE(F06,1103) MSG( 3,4)
            WRITE(F06,*)

         ELSE IF (INFO ==    -1) THEN                      ! MSG( 4), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7104'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 4,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 4,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -2) THEN                      ! MSG( 5), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7105'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 5,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 5,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -3) THEN                      ! MSG( 6), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7106'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 6,1)
            WRITE(ERR,71006) NCV, NEV, NDOFL
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 6,1)
            WRITE(F06,71006) NCV, NEV, NDOFL
            WRITE(F06,*)

         ELSE IF (INFO ==    -4) THEN                      ! MSG( 7), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7107'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 7,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 7,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -5) THEN                      ! MSG( 8), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7108'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 8,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 8,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -6) THEN                      ! MSG( 9), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7109'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 9,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 9,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -7) THEN                      ! MSG(10), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7110'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(10,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(10,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -8) THEN                      ! MSG(11)

            ERROR_NUM = ' 7111'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,702) ERROR_NUM, SUBNAME, INFO, MSG(11,1)
            WRITE(ERR,*)
            WRITE(F06,702) ERROR_NUM, SUBNAME, INFO, MSG(11,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -9) THEN                      ! MSG(12)

            ERROR_NUM = ' 7112'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,702) ERROR_NUM, SUBNAME, INFO, MSG(12,1)
            WRITE(ERR,*)
            WRITE(F06,702) ERROR_NUM, SUBNAME, INFO, MSG(12,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -10) THEN                      ! MSG(13), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7113'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(13,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(13,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -11) THEN                      ! MSG(14), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7114'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(14,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(14,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -12) THEN                      ! MSG(15), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7115'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(15,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(15,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -13) THEN                      ! MSG(16), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7116'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(16,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(16,1)
            WRITE(F06,*)

         ELSE IF (INFO == -9999) THEN                      ! MSG(21)

            ERROR_NUM = ' 7121'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,702) ERROR_NUM, SUBNAME, INFO, MSG(21,1)
            WRITE(ERR,1103) MSG(21,2)
            WRITE(ERR,1105) IPARAM(5)
            WRITE(ERR,*)
            WRITE(F06,702) ERROR_NUM, SUBNAME, INFO, MSG(21,1)
            WRITE(F06,1103) MSG(21,2)
            WRITE(F06,1105) IPARAM(5)
            WRITE(F06,*)

         ELSE

            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,9991) INFO
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9991) INFO
            ENDIF

         ENDIF

      ELSE IF ((SUBNAME == 'DSEUPD') .OR. (SUBNAME == 'dseupd')) THEN

         IF      (INFO ==    -1) THEN                      ! MSG( 4), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7204'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 4,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 4,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -2) THEN                      ! MSG( 5), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7205'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 5,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 5,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -3) THEN                      ! MSG( 6), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7206'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 6,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 6,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -5) THEN                      ! MSG( 8), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7208'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 8,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 8,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -6) THEN                      ! MSG( 9), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7209'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG( 9,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG( 9,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -7) THEN                      ! MSG(10), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7210'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(10,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(10,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -8) THEN                      ! MSG(11)

            ERROR_NUM = ' 7211'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,702) ERROR_NUM, SUBNAME, INFO, MSG(11,1)
            WRITE(ERR,*)
            WRITE(F06,702) ERROR_NUM, SUBNAME, INFO, MSG(11,1)
            WRITE(F06,*)

         ELSE IF (INFO ==    -9) THEN                      ! MSG(12)

            ERROR_NUM = ' 7212'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,702) ERROR_NUM, SUBNAME, INFO, MSG(12,1)
            WRITE(ERR,*)
            WRITE(F06,702) ERROR_NUM, SUBNAME, INFO, MSG(12,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -10) THEN                      ! MSG(13), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7213'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(13,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(13,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -11) THEN                      ! MSG(14), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7214'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(14,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(14,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -12) THEN                      ! MSG(16), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7216'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(16,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(16,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -14) THEN                      ! MSG(18)

            ERROR_NUM = ' 7218'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,702) ERROR_NUM, SUBNAME, INFO, MSG(18,1)
            WRITE(ERR,*)
            WRITE(F06,702) ERROR_NUM, SUBNAME, INFO, MSG(18,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -15) THEN                      ! MSG(17), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7217'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(17,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(17,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -16) THEN                      ! MSG(19), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7219'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(19,1)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(19,1)
            WRITE(F06,*)

         ELSE IF (INFO ==   -17) THEN                      ! MSG(20), coding error (wrong input to subr DSBAND)

            ERROR_NUM = ' 7220'
            FATAL_ERR = FATAL_ERR + 1
            WRITE(ERR,710) ERROR_NUM, SUBNAME, INFO, MSG(20,1)
            WRITE(ERR,1103) MSG(20,2)
            WRITE(ERR,*)
            WRITE(F06,710) ERROR_NUM, SUBNAME, INFO, MSG(20,1)
            WRITE(F06,1103) MSG(20,2)
            WRITE(F06,*)

         ELSE

            WARN_ERR = WARN_ERR + 1
            WRITE(ERR,9991) INFO
            IF (SUPWARN == 'N') THEN
               WRITE(F06,9991) INFO
            ENDIF

         ENDIF

      ELSE

         WARN_ERR = WARN_ERR + 1
         WRITE(ERR,9992) SUBNAME
         IF (SUPWARN == 'N') THEN
            WRITE(F06,9992) SUBNAME
         ENDIF

      ENDIF

! **********************************************************************************************************************************
  100 FORMAT(' *INFORMATION: THE FOLLOWING MESSAGE WAS FOUND BY ARPACK SUBR ',A,' CALLED BY ARPACK DRIVER DSBAND IN ',A,           &
                           ' SUBR EIG_LANCZOS.'                                                                                    &
                    ,/,14X,' IT RELATES TO AN ERROR, LISTED BELOW, FOR ARPACK "INFO" = ',I6,' IN THE LANCZOS EIGENVALUE METHOD:',/)

  710 FORMAT(' *ERROR ',A,': PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' FOR LANCZOS EIGENVALUE EXTRACTION WITH ARPACK "INFO" = ',I6,':'                                       &
                    ,/,15X,A)

  702 FORMAT(' *ERROR ',A,': ERROR IN SUBROUTINE ',A,' FOR LANCZOS EIGENVALUE EXTRACTION WITH ARPACK "INFO" = ',I6,/,15X,A)

 1103 FORMAT(15X,A)

 1104 FORMAT(15X,I8)

 1105 FORMAT(15X,'IPARAM(5) = ',I8)

 1106 FORMAT(15X,'ARPACK "MODE" = ',I3)

 9991 FORMAT(' *WARNING    : INVALID INPUT FOR VARIABLE INFO = ',I6,' IN ARPACK SUBR DSBAND_INFO_MSG.'                             &
                    ,/,14X,' CANNOT WRITE ERROR MESSAGES FROM ARPACK ROUTINES FOR LANCZOS EIGENVALUE EXTRACTION',/)

 9992 FORMAT(' *WARNING    : INVALID INPUT FOR VARIABLE SUBNAME = ',A,' IN ARPACK SUBR DSBAND_INFO_MSG.'                           &
                    ,/,14X,' CANNOT WRITE ERROR MESSAGES FROM ARPACK ROUTINES FOR LANCZOS EIGENVALUE EXTRACTION',/)

71006 FORMAT( 14X,' NCV = ',I8,',  NEV = ',I8,',  NDOFL = ',I8)

! **********************************************************************************************************************************

      END SUBROUTINE ARPACK_INFO_MSG

   END MODULE ARPACK_DIAGNOSTICS
