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

! MYSTRAN stands for "MY STRuctural ANalysis and is a general purpose finite element analysis program
!                     -- ---        --
      PROGRAM MYSTRAN

! This is the main program. The major functions are to call the primary subroutines that comprise MYTRAN. The primary subroutines
! are LINK0, LINK1, LINK2, LINK3, LINK4, LINK5, LINK6, LINK9:

!  LINK 0:
!     Read input data deck and check for errors and possible restart
!     Process grid and coordinate system input data
!     Process Case Control output requests
!     Forms degree of freedom (DOF) tables
!     Process concentrated mass input data
!     Calculates rigid body mass properties (Grid Point Weight Generator)
!     Process temperature and pressure load input data intp arrays needed for element load calculations

!  LINK 1:
!     Process MPC's and rigid elements into sparse array RMG
!     Process all applied forces (grid forces and moments, gravity, pressure, thermal, centrifugal, scalar) into sparse array PG)
!     Formulate the G-set sparse stiffness and mass arrays KGG and MGG
!     Formulate the G-set sparse differential stiffness array KGGD

!  LINK 2:
!     Reduce the G-set stiffness, mass, load and constraint matrices to the L-set

!  LINK 3 (for statics problems only):
!     Solve for the L-set displacements

!  LINK 4 (for eigenvalue problems only):
!     Solve for the L-set eigenvalues and eigenvectors

!  LINK 5:
!     Build the A-set displacements back up to the G-set through use of the constraint matrices

!  LINK 6 (for Craig-Bampton substructure model generation only):
!     Builds a Craig-Bampton model from the input physical model. It is a normal eigenvalue run with SUPORT Bulk Data that defines
!     a boundary where the model will connect with other models. All of the matrices needed to couple this model to ones connected
!     at its boundary are calculated along with many output transformation matrices used by Craig-Bampton analysts.

!  LINK 9:
!     Use the G-set displacements from LINK5 to solve for the outputs requested in Case Control.


      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE IOUNT1, ONLY                :  BUG, ERR, F06, IN0, IN1, L1A, NEU, SC1

      USE IOUNT1, ONLY                :  F06FIL, IN0FIL, INFILE, NEUFIL

      USE IOUNT1, ONLY                :  IN0_MSG, IN1_MSG, L1A_MSG

      USE IOUNT1, ONLY                :  BUGSTAT, BUGSTAT_OLD, ERRSTAT, ERRSTAT_OLD, OP2STAT

      USE IOUNT1, ONLY                :  LEN_INPUT_FNAME, LEN_RESTART_FNAME,LINK1A, BUGOUT, RESTART_FILNAM

      USE SCONTR, ONLY                :  COMM, FATAL_ERR, LINKNO_START, LSETLN, LSETS, LSUB, NDOFL, NSETS, NSUB, NTSUB,            &
                                         PROG_NAME, RESTART, SETLEN, SOL_NAME, WARN_ERR

      USE TIMDAT, ONLY                :  HOUR, MINUTE, SEC, SFRAC, YEAR, MONTH, DAY, STIME,                                        &
                                         START_HOUR, START_MINUTE, START_SEC, START_SFRAC, START_YEAR, START_MONTH, START_DAY

      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP, NL_ITER_NUM, NL_MAXITER, NL_NORM, NL_NUM_LOAD_STEPS

      USE PARAMS, ONLY                :  EPSIL, RELINK3, SUPWARN
      USE CONSTANTS_1, ONLY           :  ZERO
      USE DEBUG_PARAMETERS, ONLY      :  DEBUG
      USE MYSTRAN_Version, ONLY       :  MYSTRAN_VER_NUM, MYSTRAN_VER_MONTH, MYSTRAN_VER_DAY, MYSTRAN_VER_YEAR, MYSTRAN_AUTHOR
      USE OUTPUT4_MATRICES, ONLY      :  ACT_OU4_MYSTRAN_NAMES, ACT_OU4_OUTPUT_NAMES, ALLOW_OU4_MYSTRAN_NAMES,                     &
                                         ALLOW_OU4_OUTPUT_NAMES, HAS_OU4_MAT_BEEN_PROCESSED, NUM_OU4_REQUESTS,                     &
                                         NUM_OU4_VALID_NAMES, SUBR_WHEN_TO_WRITE_OU4_MATS
      USE COL_VECS, ONLY              :  UG_COL

      USE DATE_TIME_UTILS, ONLY       :  OURDAT, OURTIM, TIME_INIT
      USE INI_FILE, ONLY              :  READ_INI
      USE STARTUP_FILES, ONLY         :  MYSTRAN_FILES, PROCESS_INCLUDE_FILES, READ_CL, READ_INPUT_FILE_NAME
      USE FILE_LIFECYCLE, ONLY        :  CLOSE_LIJFILES, CLOSE_OUTFILES, FILE_CLOSE, FILE_INQUIRE, FILE_OPEN, WRITE_FILNAM
      USE INPUT_FILE_MECHANICS, ONLY  :  IS_THIS_A_RESTART
      USE INPUT_STAGE_DISPATCH, ONLY  :  LOADE0
      USE TEMP_FILE_READERS, ONLY     :  READ_L1A
      USE LINK0_LINK1_MOD, ONLY       :  LINK0, LINK1
      USE LINK2_MOD, ONLY             :  LINK2
      USE LINK3_MOD, ONLY             :  LINK3
      USE LINK4_MOD, ONLY             :  LINK4
      USE LINK6_MOD, ONLY             :  LINK6
      USE RIGID_BODY_STORAGE_LIFECYCLE, ONLY:  DEALLOCATE_RBGLOBAL
      USE LINK5_MOD, ONLY             :  LINK5
      USE FILE_LIFECYCLE, ONLY        :  OUTA_HERE
      USE LINK9_MOD, ONLY             :  LINK9
      USE RESTART_FILE_IO, ONLY       :  RESTART_DATA_FOR_L3
      USE NONLINEAR_PARAMETER_LIFECYCLE, ONLY:  DEALLOCATE_NL_PARAMS
      USE VECTOR_METRICS, ONLY        :  VECTOR_NORM
      USE STARTUP_SPECIALS, ONLY      :  PRINT_BUILD_INFO, SET_BLAS_THREADS

      IMPLICIT NONE


      CHARACTER( 1*BYTE)              :: INI_EXIST         ! 'Y' if file MYSTRAN.INI exists or 'N' otherwise


      INTEGER(LONG)                   :: I,J                ! DO loop indices
      INTEGER(LONG)                   :: I1                 ! Length of file name of LINK1A without the extension
      INTEGER(LONG)                   :: ITER_MAX           ! Naximum iterations before quiting
      INTEGER(LONG)                   :: NUM_INCL_FILES     ! Number of INCLUDE statements found in the Bulk data file
      INTEGER(LONG)                   :: NUM_LOAD_STEPS     ! Number of steps to divide the load into (1 unless SOL is nonlinear)
      INTEGER(LONG)                   :: OUNT(2)            ! File units to write messages to

      REAL(DOUBLE)                    :: CPU_SECS           ! CPU time for job
      REAL(DOUBLE)                    :: EPS1               ! Small number with which to compare to 0
      REAL(DOUBLE) , ALLOCATABLE      :: PERCENT_CHANGE(:,:)! Used to test convergence of UG_NORM
      REAL(DOUBLE)                    :: TIME_END           ! Time job is ended
      REAL(DOUBLE)                    :: TIME_START         ! Time job is started
      REAL(DOUBLE) , ALLOCATABLE      :: UG_NORM(:,:)       ! Norm of UG_COL vector
      CHARACTER(LEN=LEN(INFILE))      :: ARG1               ! First argument. Used to check for --notice
      INTEGER(LONG)                   :: NC_ARG1            ! Number of characters in ARG1.

! **********************************************************************************************************************************
      EPS1 = EPSIL(1)

      CALL CPU_TIME ( TIME_START )

! Default units for writing errors the screen (until LINK1A is read)

      OUNT(1) = SC1
      OUNT(2) = SC1

! Open SC1

      OPEN(SC1)

! Cap threads used by the linked BLAS provider (e.g. OpenBLAS). MYSTRAN
! drives BLAS from its own serial loops, so per-call thread fan-out is
! pure overhead. This is a no-op unless a recognized threaded BLAS
! provider was detected at configure time.

      CALL SET_BLAS_THREADS ( 1_LONG )

! Set time initializing parameters

      CALL TIME_INIT

! Write logo and copyright notice to screen. Then write MYSTRAN start time/date

      WRITE(SC1,117) PROG_NAME, MYSTRAN_VER_NUM, MYSTRAN_VER_MONTH, MYSTRAN_VER_DAY, MYSTRAN_VER_YEAR, MYSTRAN_AUTHOR

! Check for --notice and --licenses

      CALL READ_CL (ARG1, NC_ARG1)
      IF (NC_ARG1 == 9 .AND. ARG1(1:9) == "--version") THEN
         CALL PRINT_BUILD_INFO (.FALSE.)
         STOP
      END IF
      IF (NC_ARG1 == 10 .AND. ARG1(1:10) == "--licenses") THEN
         CALL PRINT_BUILD_INFO (.TRUE.)
         STOP
      END IF

! Read data in initialization file, MYSTRAN.INI, if it exists.

      CALL READ_INI ( INI_EXIST )

      CALL OURTIM
      CALL OURDAT
      WRITE(SC1,150) MONTH, DAY, YEAR, HOUR, MINUTE, SEC, SFRAC
      START_HOUR   = HOUR
      START_MINUTE = MINUTE
      START_SEC    = SEC
      START_SFRAC  = SFRAC
      START_YEAR   = YEAR
      START_MONTH  = MONTH
      START_DAY    = DAY

      STIME = (10**8)*MONTH+(10**6)*DAY+(10**4)*HOUR+(10**2)*MINUTE+SEC

! Read input data filename (from command line) and calc LEN_INPUT_FNAME.

      CALL READ_INPUT_FILE_NAME ( INI_EXIST )



      IF (RESTART == 'N') THEN
         LEN_RESTART_FNAME = LEN_INPUT_FNAME
      ENDIF

! Write file name to screen

      CALL WRITE_FILNAM ( INFILE, SC1, 4 )

! If MYSTRAN is to start in LINK 0 then do the following. It can be started in a later LINK if MYSTRAN was aborted just as a new
! LINK was beginning and then warm restarted in the next LINK if the INI file is modified to have LINK_START = that LINK number

      IF (LINKNO_START == 0) THEN
                                                           ! Open input data file (FEM model) for reading
         CALL FILE_OPEN ( IN1, INFILE, OUNT, 'OLD', IN1_MSG, 'NEITHER', 'FORMATTED', 'READ', 'REWIND', 'N', 'N' )

         CALL IS_THIS_A_RESTART                            ! Only check if RESTART entry is in E.C. (need for subr MYSTRAN_FILES)
         REWIND (IN1)
         CALL MYSTRAN_FILES ( START_MONTH, START_DAY, START_YEAR, START_HOUR, START_MINUTE, START_SEC, START_SFRAC )
                                                           ! Process INCLUDE entries in whole DAT file here to create the IN0 file
         CALL PROCESS_INCLUDE_FILES ( NUM_INCL_FILES )
         IF (NUM_INCL_FILES > 0) THEN
            CALL FILE_CLOSE ( IN1, INFILE, 'KEEP' )
            CALL FILE_CLOSE ( IN0, INFILE, 'KEEP' )
            CALL FILE_OPEN ( IN1, IN0FIL, OUNT, 'OLD', IN0_MSG, 'NEITHER', 'FORMATTED', 'READ', 'REWIND', 'N', 'N' )
            INFILE(1:) = IN0FIL(1:)
         ELSE
            CALL FILE_CLOSE ( IN0, INFILE, 'DELETE' )
            REWIND (IN1)
         ENDIF

         CALL LOADE0                                       ! Now perform the other functions of old LOADE0 (F06, etc now open)

         IF (RESTART == 'N') THEN                          ! Normal execution (not a RESTART)


            BUGSTAT_OLD = BUGSTAT                          ! Default value from module IOUNT1
            ERRSTAT_OLD = ERRSTAT                          ! Default value from module IOUNT1

         ELSE                                              ! This is a RESTART

            I1 = LEN_RESTART_FNAME
            LINK1A(1:I1)  = RESTART_FILNAM(1:I1)
            LINK1A(I1+1:) = 'L1A'
            CALL READ_L1A ( 'KEEP' )

            BUGSTAT_OLD = BUGSTAT                          ! Old value from LINK1A
            ERRSTAT_OLD = ERRSTAT                          ! Old value from LINK1A

            NSUB  = 0                                      ! Initialize items read from L1A that could be changed in the RESTART
            NTSUB = 0

         ENDIF

         LOAD_ISTEP = 1                                    ! Want 1 in LINK0 - also when the steps: DO loop runs for I=1 below
         CALL LINK0                                        ! Call LINK's to process input and solve problem.

      ENDIF

      IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
         ALLOCATE(UG_NORM(NL_NUM_LOAD_STEPS,NL_MAXITER))
         ALLOCATE(PERCENT_CHANGE(NL_NUM_LOAD_STEPS,NL_MAXITER))
         DO I=1,NL_NUM_LOAD_STEPS
            DO J=1,NL_MAXITER
               UG_NORM(I,J) = ZERO
               PERCENT_CHANGE(I,J) = ZERO
            ENDDO
         ENDDO
      ENDIF

      IF (RESTART == 'N') THEN

         IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
            NUM_LOAD_STEPS = NL_NUM_LOAD_STEPS
            ITER_MAX       = NL_MAXITER
         ELSE IF ( SOL_NAME(1:8) == 'BUCKLING') THEN
            NUM_LOAD_STEPS = 2                             ! 1st "load" step is linear statics, 2nd is buckling eigens
            ITER_MAX       = 0
         ELSE
            NUM_LOAD_STEPS = 1
            ITER_MAX       = 0
         ENDIF

steps:   DO I=1,NUM_LOAD_STEPS

            NL_ITER_NUM = 0
            LOAD_ISTEP = I

iters:      DO

               NL_ITER_NUM = NL_ITER_NUM + 1

               CALL LINK1

               CALL LINK2
               IF (SOL_NAME(1:8) == 'BUCKLING') THEN
                  IF (I == 1) THEN                         ! This is the linear static portion of BUCKLING
                     IF (NDOFL > 0) THEN
                        CALL LINK3
                     ELSE                                  ! This is eigens for BUCKLING
                        COMM(3) = 'C'                      ! --- Need to do this, otherwise LINK5 will fail
                     ENDIF
                  ELSE
                     CALL LINK4
                  ENDIF
               ELSE IF ((SOL_NAME(1:7) == 'STATICS') .OR. (SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
                  IF (NDOFL > 0) THEN
                     CALL LINK3
                  ELSE
                     COMM(3) = 'C'                         ! Need to do this, otherwise LINK5 will fail
                  ENDIF
               ELSE IF ((SOL_NAME(1:5) == 'MODES') .OR. (SOL_NAME(1:12) == 'GEN CB MODEL')) THEN
                  CALL LINK4
               ENDIF

               IF (SOL_NAME(1:12) == 'GEN CB MODEL') THEN
                  CALL LINK6
               ENDIF
               CALL DEALLOCATE_RBGLOBAL ( 'G ' )

               CALL LINK5
               IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
                  CALL CALC_CONVERGENCE
                  IF (DABS(PERCENT_CHANGE(LOAD_ISTEP,NL_ITER_NUM)) < EPSIL(3)) THEN
                     EXIT iters
                  ELSE
                     IF (NL_ITER_NUM < ITER_MAX) THEN
                        CYCLE iters
                     ELSE
                        WRITE(ERR,30) ITER_MAX, EPSIL(3)
                        WRITE(F06,30) ITER_MAX, EPSIL(3)
                        FATAL_ERR = FATAL_ERR + 1
                        CALL OUTA_HERE ( 'Y' )
                     ENDIF
                  ENDIF
               ENDIF

               IF (SOL_NAME(1:8) == 'NLSTATIC') THEN
                  IF (NL_ITER_NUM > ITER_MAX) THEN
                     FATAL_ERR = FATAL_ERR + 1
                     WRITE(ERR,20) ITER_MAX
                     WRITE(F06,20) ITER_MAX
                     CALL OUTA_HERE ( 'Y' )
                  ENDIF
               ELSE
                  EXIT iters
               ENDIF

            ENDDO iters

            IF ((SOL_NAME(1:8) == 'DIFFEREN') .OR. (SOL_NAME(1:8) == 'NLSTATIC')) THEN
               WRITE(F06,1001) LOAD_ISTEP
               WRITE(F06,1002)
               DO J=1,NL_ITER_NUM
                  WRITE(F06,1003) J, UG_NORM(I,J), PERCENT_CHANGE(I,J)
               ENDDO
            ENDIF

            CALL LINK9 ( I )

         ENDDO steps

      ELSE                                                 ! This is a restart

         IF (RELINK3 == 'Y') THEN
            IF (SOL_NAME(1:7) == 'STATICS') THEN
               CALL RESTART_DATA_FOR_L3
               CALL LINK3
               CALL LINK9 ( 1 )
               CALL LINK5
            ENDIF
         ENDIF

      ENDIF

      CALL DEALLOCATE_NL_PARAMS

! If any OUTPUT4 matrices have been requested for output, make sure they were all processed

      DO I=1,NUM_OU4_REQUESTS
         DO J=1,NUM_OU4_VALID_NAMES
            IF (ACT_OU4_MYSTRAN_NAMES(I) == ALLOW_OU4_MYSTRAN_NAMES(J)) THEN
               IF (HAS_OU4_MAT_BEEN_PROCESSED(I,1) /= 'Y') THEN
                  IF (ACT_OU4_MYSTRAN_NAMES(I) == ACT_OU4_OUTPUT_NAMES(I)) THEN
                     WRITE(ERR,10) SUBR_WHEN_TO_WRITE_OU4_MATS(J), ALLOW_OU4_MYSTRAN_NAMES(J)
                     WRITE(F06,10) SUBR_WHEN_TO_WRITE_OU4_MATS(J), ALLOW_OU4_MYSTRAN_NAMES(J)
                  ELSE
                     WRITE(ERR,10) SUBR_WHEN_TO_WRITE_OU4_MATS(J), ALLOW_OU4_OUTPUT_NAMES(J)
                     WRITE(F06,10) SUBR_WHEN_TO_WRITE_OU4_MATS(J), ALLOW_OU4_OUTPUT_NAMES(J)
                  ENDIF
               ENDIF
            ENDIF
         ENDDO
      ENDDO

! Set close status for output files

      IF (BUGOUT == 'Y') THEN
         BUGSTAT = 'KEEP'
      ELSE
         IF (BUGSTAT_OLD == 'KEEP    ') THEN
            BUGSTAT = 'KEEP'
         ELSE
            BUGSTAT = 'DELETE'
         ENDIF
      ENDIF

      IF ((FATAL_ERR > 0) .OR. (WARN_ERR > 0)) THEN
         ERRSTAT = 'KEEP'
      ELSE
         IF (ERRSTAT_OLD == 'KEEP    ') THEN
            ERRSTAT = 'KEEP'
         ELSE
            ERRSTAT = 'DELETE'
         ENDIF
      ENDIF


!xx   IF (WRITE_NEU) THEN
!xx      CALL FILE_CLOSE ( NEU, NEUFIL, 'KEEP' )
!xx   ENDIF

! Write MYSTRAN END to BUG, ERR, F06 and then close those files

      WRITE(F06,*)
      CALL CPU_TIME ( TIME_END )
      CPU_SECS = TIME_END - TIME_START
      WRITE(F06,200) CPU_SECS

      CALL OURTIM
      CALL OURDAT
      WRITE(BUG,152) MONTH,DAY,YEAR,HOUR,MINUTE,SEC,SFRAC
      WRITE(ERR,152) MONTH,DAY,YEAR,HOUR,MINUTE,SEC,SFRAC
      WRITE(F06,152) MONTH,DAY,YEAR,HOUR,MINUTE,SEC,SFRAC
      WRITE(F06,*)
      WRITE(F06,*) "                                         * * * END OF JOB * * *"

      IF (( DEBUG(193) == 100) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of MAIN' )
      ENDIF

      CALL CLOSE_OUTFILES ( BUGSTAT, ERRSTAT, OP2STAT )

! Close LIJ files

      CALL CLOSE_LIJFILES ( 'FILE_STAT' )

! Write MYSTRAN END and output file name to screen

      WRITE(SC1,151) MONTH,DAY,YEAR,HOUR,MINUTE,SEC,SFRAC
      CALL WRITE_FILNAM ( F06FIL, SC1, 4 )

      IF (WARN_ERR > 0) THEN
   !xx   WRITE(SC1, * )
         IF (SUPWARN == 'N') THEN
            WRITE(SC1,156) WARN_ERR
         ELSE
            WRITE(SC1,157) WARN_ERR
         ENDIF
      ENDIF

      WRITE(SC1,155) CPU_SECS

      STOP

! **********************************************************************************************************************************
   10 FORMAT(' *ERROR    10: PROGRAMMING ERROR IN SUBROUTINE ',A                                                                   &
                    ,/,14X,' REQUESTED OUTPUT4 MATRIX ',A,' HAS NOT BEEN PROCESSED.',/)

   20 FORMAT(' *ERROR    20: MAXIMUM ITERATIONS, NL_MAXITER = ',I8,' EXCEEDED')

   30 FORMAT(' *ERROR    30: MAXIMUM ITERATIONS = ',I8,' HAVE BEEN TAKEN WITHOUT CONVERGING ON THE NONLINEAR DISPL.'               &
                    ,/,14X,' BULK DATA PARAM ITER_MAX CAN BE USED TO INCREASE THE NUMBER OF ITERATIONS OR PARAM EPSIL 3 = ',1ES14.6&
                    ,/,14X,' (THE ITERATION % CHANGE CRITERIA) CAN BE CHANGED')

   99 FORMAT(/,' PROCESSING TERMINATED DUE TO ',I8,' INPUT ERRORS. CHECK OUTPUT FILE FOR ERROR MESSAGES')

  111 FORMAT(/' Initialization file ',A,' does not exist. internal defaults will be used',/)

  117 FORMAT(1X,A,' Version',5(1X,A))

  150 FORMAT(/,' >> MYSTRAN BEGIN  : ',I2,'/',I2,'/',I4,' at ',I2,':',I2,':',I2,'.',I3,'. The input file is:')

  151 FORMAT(/,' >> MYSTRAN END    : ',I2,'/',I2,'/',I4,' at ',I2,':',I2,':',I2,'.',I3,'. The output file is:')

  152 FORMAT(/,' >> MYSTRAN END    : ',I2,'/',I2,'/',I4,' at ',I2,':',I2,':',I2,'.',I3)

  155 FORMAT(' MYSTRAN terminated normally. Total CPU time = ',1ES9.2,' seconds',/,'                    ========')

  156 FORMAT(' Check F06 output file for ',I8,' warning messages')

  157 FORMAT(' Check ERR output file for ',I8,' warning messages')

  200 FORMAT('    Total CPU time = ',1ES9.2,' seconds')

 1001 FORMAT(40X,'Convergence data for load step ',I4,':',/,                                                                       &
             40X,'-----------------------------------')

 1002 FORMAT(35X,'Iter Num   Norm of UG vec   % chge since last')

 1003 FORMAT(32X,I9,1ES18.6,7X,E9.2,'%')


! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE CALC_CONVERGENCE

      USE SCONTR, ONLY                :  NDOFG
      USE CONSTANTS_1, ONLY           :  ZERO, ONE_HUNDRED

      IMPLICIT NONE

      INTEGER(LONG)                   :: IERR

! **********************************************************************************************************************************
      PERCENT_CHANGE(LOAD_ISTEP,1) = ONE_HUNDRED
      CALL VECTOR_NORM ( UG_COL, NDOFG, NL_NORM, UG_NORM(LOAD_ISTEP,NL_ITER_NUM), IERR )
      IF (NL_ITER_NUM > 1) THEN
!zzzz    IF (DABS(UG_NORM(LOAD_ISTEP,NL_ITER_NUM-1)) > EPS1) THEN
            PERCENT_CHANGE(LOAD_ISTEP,NL_ITER_NUM) = 1.0D2*(UG_NORM(LOAD_ISTEP,NL_ITER_NUM) - UG_NORM(LOAD_ISTEP,NL_ITER_NUM-1))   &
                                                    /UG_NORM(LOAD_ISTEP,NL_ITER_NUM-1)
!zzzz    ELSE
!zzzz       IF (DABS(UG_NORM(LOAD_ISTEP,NL_ITER_NUM) - UG_NORM(LOAD_ISTEP,NL_ITER_NUM-1)) > EPS1) THEN
!zzzz          PERCENT_CHANGE(LOAD_ISTEP,NL_ITER_NUM) = ONE_HUNDRED                ! Indeterminant but set it to a large number
!zzzz       ELSE
!zzzz          PERCENT_CHANGE(LOAD_ISTEP,NL_ITER_NUM) = ZERO
!zzzz       ENDIF
!zzzz    ENDIF
      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE CALC_CONVERGENCE

      END PROGRAM MYSTRAN
