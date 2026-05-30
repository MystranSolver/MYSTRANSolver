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

      SUBROUTINE LINK2

! Reduces G-set matrices to L-set. Performs the following major functions:

!   1) Calls routines to reduce KGG, MGG, PG to N-set (needs to solve for GMN = RMM(-1)*RMN first)
!   2) Calls routines to reduce KNN, MNN, PN to F-set
!   3) Calls routines to reduce KFF, MFF, PF to A-set
!   4) Calls routines to reduce KAA, MAA, PA to L-set

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE

      USE IOUNT1, ONLY                :  WRT_BUG, ERR, F06, L1A, ERRSTAT, SC1
      USE IOUNT1, ONLY                :  L2G,     L2H    , L2I   ,  L2O    , L2P   ,  L2Q
      USE IOUNT1, ONLY                :  LINK2G,  LINK2H , LINK2I , LINK2O , LINK2P , LINK2Q
      USE IOUNT1, ONLY                :  L2G_MSG, L2H_MSG, L2I_MSG, L2O_MSG, L2P_MSG, L2Q_MSG

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, COMM, FATAL_ERR, LINKNO, MBUG, NSUB,                                        &
                                         NDOFG, NDOFM, NDOFN, NDOFS, NDOFSA, NDOFF, NDOFO, NDOFA, NDOFL, NDOFR,                    &
                                         NTERM_KGG , NTERM_KNN , NTERM_KNM , NTERM_KMM ,                                           &
                                         NTERM_KFF , NTERM_KFS , NTERM_KSS ,                                                       &
                                         NTERM_KAA , NTERM_KAO , NTERM_KOO ,                                                       &
                                         NTERM_KLL , NTERM_KRL , NTERM_KRR ,                                                       &
                                         NTERM_KGGD, NTERM_KNND, NTERM_KNMD, NTERM_KMMD,                                           &
                                         NTERM_KFFD, NTERM_KFSD, NTERM_KSSD,                                                       &
                                         NTERM_KAAD, NTERM_KAOD, NTERM_KOOD,                                                       &
                                         NTERM_KLLD, NTERM_KRLD, NTERM_KRRD,                                                       &
                                                     NTERM_MAA , NTERM_MLL ,                                                       &
                                                     NTERM_PA  , NTERM_PL  ,                                                       &
                                         NTERM_RMG , SOL_NAME  , WARN_ERR

      USE CONSTANTS_1, ONLY           :  ONE
      USE RIGID_BODY_DISP_MATS, ONLY  :  RBGLOBAL_GSET
      USE PARAMS, ONLY                :  CUSERIN, CUSERIN_XSET, EQCHK_OUTPUT, EQCHK_NORM, PRTSTIFF, PRTSTIFD, PRTMASS, PRTFOR,     &
                                         SUPINFO, SUPWARN
      USE NONLINEAR_PARAMS, ONLY      :  LOAD_ISTEP
      USE DEBUG_PARAMETERS

      USE SPARSE_MATRICES, ONLY       :  I_KGG, J_KGG, KGG, I_KGGD, J_KGGD, KGGD,                                                  &
                                         I_KAA, J_KAA, KAA, I_MAA, J_MAA, MAA, I_PA, J_PA , PA,                                    &
                                         I_KLL, J_KLL, KLL, I_MLL, J_MLL, MLL, I_PL, J_PL , PL
      USE FULL_MATRICES, ONLY         :  KGG_FULL
      USE SPARSE_MATRICES, ONLY       :  SYM_KGG
      USE OUTPUT4_MATRICES, ONLY      :  NUM_OU4_REQUESTS
      USE LINK2_USE_IFs
      USE LINK_MESSAGE_Interface


      IMPLICIT NONE

      CHARACTER, PARAMETER            :: CR13 = CHAR(13)   ! This causes a carriage return simulating the "+" action in a FORMAT
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK2'
      CHARACTER(  1*BYTE)             :: CLOSE_IT          ! Input to subr READ_MATRIX_i. 'Y'/'N' whether to close a file or not
      CHARACTER( 8*BYTE)              :: CLOSE_STAT        ! What to do with file when it is closed

      INTEGER(LONG)                   :: NROWS             ! Value of DOF size to pass to subr WRITE_USERIN_BD_CARDS
      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG), PARAMETER        :: P_LINKNO  = 1     ! Prior LINK no's that should have run before this LINK can execute

      REAL(DOUBLE)                    :: KGG_DIAG(NDOFG)   ! Diagonal of KGG
      REAL(DOUBLE)                    :: KGG_MAX_DIAG      ! Max diag term from KGG
      REAL(DOUBLE)                    :: KGGD_DIAG(NDOFG)  ! Diagonal of KGGD
      REAL(DOUBLE)                    :: KGGD_MAX_DIAG     ! Max diag term from KGGD

! **********************************************************************************************************************************
      LINKNO = 2

! Set time initializing parameters

      CALL TIME_INIT

! Initialize WRT_BUG

      DO I=0,MBUG-1
         WRT_BUG(I) = 0
      ENDDO

! Get date and time, write to screen

      CALL OURDAT
      CALL OURTIM
      WRITE(SC1,152) LINKNO

! Write info to text files

      WRITE(F06,150) LINKNO
      WRITE(ERR,150) LINKNO

! Read LINK1A file

      CALL READ_L1A ( 'KEEP' )
! Check COMM for successful completion of prior LINKs

      IF (COMM(P_LINKNO) /= 'C') THEN
         WRITE(ERR,9998) P_LINKNO,P_LINKNO,LINKNO
         WRITE(F06,9998) P_LINKNO,P_LINKNO,LINKNO
         FATAL_ERR = FATAL_ERR + 1
         CALL OUTA_HERE ( 'Y' )                            ! Prior LINK's didn't complete, so quit
      ENDIF

! Deallocate the sparse matrices that will be generated in LINK2 (just to make sure none get missed later)

      CALL DEALLOCATE_LINK2_ARRAYS ( 1 )

! Write G-set stiff matrix diagonal and stats, if requested.
! NOTE: call this subr even if PRTSTIFFD(1) = 0 since we need KGG_DIAG, KGG_MAX_DIAG for the equilibrium check
! Do equilibrium check on the G-set stiffness matrix, if requested

      IF (DEBUG(57) > 0) THEN
         CALL ALLOCATE_FULL_MAT ( 'KGG_FULL', NDOFG, NDOFG, SUBR_NAME )
         CALL SPARSE_CRS_TO_FULL ( 'KGG', NTERM_KGG, NDOFG, NDOFG, 'Y', I_KGG, J_KGG, KGG, KGG_FULL )
         WRITE(F06,*) ' G-Set Stiffness Matrix, KGG'
         DO I=1,NDOFG
            IF (DEBUG(57) == 2) THEN
               WRITE(F06,*) '  ROW',I
            ENDIF
            WRITE(F06,3006) (KGG_FULL(I,J),J=1,NDOFG)
            WRITE(F06,*)
         ENDDO
 3006    FORMAT(6(1ES15.6))
      ENDIF
      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         CALL GET_MATRIX_DIAG_STATS ( 'KGGD', 'G ', NDOFG, NTERM_KGGD, I_KGGD, J_KGGD, KGGD, PRTSTIFD(1), KGGD_DIAG , KGGD_MAX_DIAG)
      ELSE
         CALL GET_MATRIX_DIAG_STATS ( 'KGG ', 'G ', NDOFG, NTERM_KGG , I_KGG , J_KGG , KGG , PRTSTIFD(1), KGG_DIAG  , KGG_MAX_DIAG )
         IF (EQCHK_OUTPUT(1) > 0) THEN
            CALL LINK_MESSAGE('EQUILIBRIUM CHECK ON KGG                  ')
            CALL STIFF_MAT_EQUIL_CHK ( EQCHK_OUTPUT(1), 'G ', SYM_KGG, NDOFG, NTERM_KGG, I_KGG, J_KGG, KGG, KGG_DIAG, KGG_MAX_DIAG,&
                                       RBGLOBAL_GSET )
         ENDIF
      ENDIF

! If NDOFR > 0 make sure A-set equil check is done

      IF (NDOFR > 0) THEN
         IF (EQCHK_OUTPUT(4) == 0) THEN
            EQCHK_OUTPUT(4) = 2
         ENDIF
      ENDIF

! Reduce G-set to N and M-sets.  After REDUCE_N_FS has run, deallocate array GRID_ID that has been
! allocated since LINK1 (needed in subr REDUCE_N_FS if AUTOSPC moves any F-set DOF's to the SA set due to small diagonal terms)

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         NTERM_KNND = 0
         NTERM_KNMD = 0
         NTERM_KMMD = 0
      ELSE
         NTERM_KNN  = 0
         NTERM_KNM  = 0
         NTERM_KMM  = 0
      ENDIF
      CALL LINK_MESSAGE('REDUCE G-SET TO N, M-SETS                   ')
      CALL REDUCE_G_NM

      WRITE(ERR,2001) NDOFM
      WRITE(ERR,2002) NDOFN
      IF (SUPINFO == 'N') THEN
         WRITE(F06,2001) NDOFM
         WRITE(F06,2002) NDOFN
      ENDIF

! Reduce N-set to F and S-sets.

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         NTERM_KFFD = 0
         NTERM_KFSD = 0
         NTERM_KSSD = 0
      ELSE
         NTERM_KFF  = 0
         NTERM_KFS  = 0
         NTERM_KSS  = 0
      ENDIF
      CALL LINK_MESSAGE('REDUCE N-SET TO F, S-SETS                   ')
      CALL REDUCE_N_FS
      WRITE(ERR,2003) NDOFS
      WRITE(ERR,2013) NDOFSA
      WRITE(ERR,2004) NDOFF
      IF (SUPINFO == 'N') THEN
         WRITE(F06,2003) NDOFS
         WRITE(F06,2013) NDOFSA
         WRITE(F06,2004) NDOFF
      ENDIF

! Reduce F-set to A and O-sets

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         NTERM_KAAD = 0
         NTERM_KAOD = 0
         NTERM_KOOD = 0
      ELSE
         NTERM_KAA  = 0
         NTERM_KAO  = 0
         NTERM_KOO  = 0
      ENDIF
      CALL LINK_MESSAGE('REDUCE F-SET TO A, O-SETS                   ')
      CALL REDUCE_F_AO
      WRITE(ERR,2005) NDOFO
      WRITE(ERR,2006) NDOFA
      IF (SUPINFO == 'N') THEN
         WRITE(F06,2005) NDOFO
         WRITE(F06,2006) NDOFA
      ENDIF

! Write A-set arrays to files if this is a CB run

      IF (SOL_NAME == 'GEN CB MODEL') THEN

         CALL LINK_MESSAGE('WRITE KAA STIFFNESS ARRAYS TO FILE')
         IF (NTERM_KAA > 0) THEN
            CALL WRITE_MATRIX_1 ( LINK2O, L2O, 'Y', 'KEEP', L2O_MSG, 'KAA', NTERM_KAA, NDOFA, I_KAA, J_KAA, KAA )
         ENDIF

         CALL LINK_MESSAGE('WRITE MAA MASS ARRAYS TO FILE')
         IF (NTERM_MAA > 0) THEN
            CALL WRITE_MATRIX_1 ( LINK2P, L2P, 'Y', 'KEEP', L2P_MSG, 'MAA', NTERM_MAA, NDOFA, I_MAA, J_MAA, MAA )
         ENDIF

         CALL LINK_MESSAGE('WRITE PA  LOAD ARRAYS TO FILE')
         IF (NTERM_PA > 0) THEN
            CALL WRITE_MATRIX_1 ( LINK2Q, L2Q, 'Y', 'KEEP', L2Q_MSG, 'PA' , NTERM_PA , NDOFA, I_PA , J_PA , PA  )
         ENDIF

      ENDIF

! Reduce A-set to L and R-sets

      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN
         NTERM_KLLD = 0
         NTERM_KRLD = 0
         NTERM_KRRD = 0
      ELSE
         NTERM_KLL  = 0
         NTERM_KRL  = 0
         NTERM_KRR  = 0
      ENDIF
      CALL LINK_MESSAGE('REDUCE A-SET TO L, R-SETS                   ')
      CALL REDUCE_A_LR
      CALL DEALLOCATE_RBGLOBAL ( 'G ' )
      WRITE(ERR,2007) NDOFR
      WRITE(ERR,2008) NDOFL
      IF (SUPINFO == 'N') THEN
         WRITE(F06,2007) NDOFR
         WRITE(F06,2008) NDOFL
      ENDIF



      IF ((SOL_NAME(1:8) == 'BUCKLING') .AND. (LOAD_ISTEP == 2)) THEN

         CONTINUE

      ELSE

                                                           ! Write L-set arrays to files
         CALL LINK_MESSAGE('WRITE L SET ARRAYS TO FILE')

         IF (NTERM_KLL > 0) THEN
            CLOSE_IT   = 'Y'
            CALL WRITE_MATRIX_1 ( LINK2G, L2G, CLOSE_IT, 'KEEP', L2G_MSG, 'KLL', NTERM_KLL, NDOFL, I_KLL, J_KLL, KLL )
         ENDIF

         IF (NTERM_MLL > 0) THEN
            CLOSE_IT   = 'Y'
            CALL WRITE_MATRIX_1 ( LINK2I, L2I, CLOSE_IT, 'KEEP', L2I_MSG, 'MLL', NTERM_MLL, NDOFL, I_MLL, J_MLL, MLL )
         ENDIF

         IF (NTERM_PL > 0) THEN
            CLOSE_IT   = 'Y'
            CLOSE_STAT = 'KEEP'
            CALL WRITE_MATRIX_1 ( LINK2H, L2H, CLOSE_IT, 'KEEP', L2H_MSG, 'PL ', NTERM_PL , NDOFL, I_PL , J_PL , PL  )
         ENDIF

         WRITE(ERR,146) 'KLL STIFFNESS MATRIX IS ', NTERM_KLL
         WRITE(ERR,146) 'MLL MASS MATRIX IS      ', NTERM_MLL
         WRITE(ERR,146) 'PL LOAD MATRIX IS       ', NTERM_PL
         IF (SUPINFO == 'N') THEN
            WRITE(F06,146) 'KLL STIFFNESS MATRIX IS ', NTERM_KLL
            WRITE(F06,146) 'MLL MASS MATRIX IS      ', NTERM_MLL
            WRITE(F06,146) 'PL LOAD MATRIX IS       ', NTERM_PL
         ENDIF

         IF (DEBUG(206) > 0) THEN
            IF (NTERM_KLL > 0) THEN
               CALL WRITE_MATRIX_MARKET_SPARSE ( 'KLL', NDOFL, NDOFL, 'Y', NTERM_KLL, I_KLL, J_KLL, KLL )
            ENDIF
            IF (NTERM_MLL > 0) THEN
               CALL WRITE_MATRIX_MARKET_SPARSE ( 'MLL', NDOFL, NDOFL, 'Y', NTERM_MLL, I_MLL, J_MLL, MLL )
            ENDIF
            IF (NTERM_PL > 0) THEN
               CALL WRITE_MATRIX_MARKET_SPARSE ( 'PL', NDOFL, NSUB, 'N', NTERM_PL, I_PL, J_PL, PL )
            ENDIF
         ENDIF

         IF (NUM_OU4_REQUESTS > 0) THEN                    ! Call OUTPUT4 processor to process output requests for OUTPUT4 matrices

            CALL LINK_MESSAGE('WRITE OUTPUT4 NATRICES      ')
            WRITE(F06,*)
            CALL OUTPUT4_PROC ( SUBR_NAME )

!xx         WRITE(SC1, * ) '     DEALLOCATE SOME ARRAYS'
      !xx   WRITE(SC1, * )                                 ! Advance 1 line for screen messages
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating KGG', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KGG' )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating MGG', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MGG' )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating PG ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'PG'  )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating KAA', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KAA' )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating MAA', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MAA' )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating PA ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'PA'  )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating KRL', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KRL' )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating KRR', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'KRR' )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating MRL', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MRL' )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating MRR', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'MRR' )
!xx         WRITE(SC1,12345,ADVANCE='NO') '       Deallocating PL ', CR13  ;   CALL DEALLOCATE_SPARSE_MAT ( 'PL'  )
!xx
         ENDIF

         IF (CUSERIN == 'Y') THEN                          ! Write out B.D. CUSERIN card images for use in INPUT4 for this substr

            NROWS = 0

            IF      (CUSERIN_XSET(1:2) == 'A ') THEN
               NROWS = NDOFA
            ELSE IF (CUSERIN_XSET(1:2) == 'G ') THEN
               NROWS = NDOFG
            ELSE IF (CUSERIN_XSET(1:2) == 'L ') THEN
               NROWS = NDOFL
            ELSE
               IF (SOL_NAME(1:12) /= 'GEN CB MODEL') THEN
                  WARN_ERR = WARN_ERR + 1
                  WRITE(ERR,2888) CUSERIN_XSET
                  IF (SUPWARN == 'N') THEN
                     WRITE(F06,2888) CUSERIN_XSET
                  ENDIF
               ENDIF
            ENDIF

            IF (NROWS > 0) THEN                               ! We don't want to call WRITE_USERIN_BD_CARDS at this time unless the
!                                                            CUSERIN request is for one of the displ sets above
               CALL WRITE_USERIN_BD_CARDS ( NROWS, CUSERIN_XSET )

            ENDIF

         ENDIF

      ENDIF

! Deallocate the sparse matrices not needed any more

      CALL DEALLOCATE_LINK2_ARRAYS ( 2 )

! Process is now complete so set COMM(LINKNO). Also set COMM(3) based on NDOFL since LINK3 will be skipped in main if NDOFL = 0

      COMM(LINKNO) = 'C'
      IF (NDOFL == 0) THEN
         COMM(3) = 'C'
      ENDIF

! Write data to L1A

      CALL WRITE_L1A ( 'KEEP', 'Y' )

! Check allocation status of allocatable arrays, if requested

      IF (DEBUG(100) > 0) THEN
         CALL CHK_ARRAY_ALLOC_STAT
         IF (DEBUG(100) > 1) THEN
            CALL WRITE_ALLOC_MEM_TABLE ( 'at the end of '//SUBR_NAME )
         ENDIF
      ENDIF

! Write LINK2 end to F06

      CALL OURTIM
      WRITE(F06,151) LINKNO

! Close files

      IF (( DEBUG(193) == 2) .OR. (DEBUG(193) == 999)) THEN
         CALL FILE_INQUIRE ( 'near end of LINK2' )
      ENDIF
      ERRSTAT = 'KEEP'

! Write LINK2 end to screen
      WRITE(SC1,153) LINKNO

! **********************************************************************************************************************************
  146 FORMAT(' *INFORMATION: NUMBER OF NONZERO TERMS IN THE ',A,'                = ',I12,/)

  150 FORMAT(/,' >> LINK',I3,' BEGIN',/)

  151 FORMAT(/,' >> LINK',I3,' END',/)

  152 FORMAT(/,' >> LINK',I3,' BEGIN')

  153 FORMAT(  ' >> LINK',I3,' END')

 2001 FORMAT(' *INFORMATION: NUMBER OF  M SET DEGREES OF FREEDOM (NDOFM)                            = ',I12)

 2002 FORMAT(' *INFORMATION: NUMBER OF  N SET DEGREES OF FREEDOM (NDOFN)                            = ',I12)

 2003 FORMAT(' *INFORMATION: NUMBER OF  S SET DEGREES OF FREEDOM (NDOFS)                            = ',I12)

 2013 FORMAT(' *INFORMATION: NUMBER OF SA SET DEGREES OF FREEDOM (NDOFSA)                           = ',I12)

 2004 FORMAT(' *INFORMATION: NUMBER OF  F SET DEGREES OF FREEDOM (NDOFF)                            = ',I12)

 2005 FORMAT(' *INFORMATION: NUMBER OF  O SET DEGREES OF FREEDOM (NDOFO)                            = ',I12)

 2006 FORMAT(' *INFORMATION: NUMBER OF  A SET DEGREES OF FREEDOM (NDOFA)                            = ',I12)

 2007 FORMAT(' *INFORMATION: NUMBER OF  R SET DEGREES OF FREEDOM (NDOFR)                            = ',I12)

 2008 FORMAT(' *INFORMATION: NUMBER OF  L SET DEGREES OF FREEDOM (NDOFL)                            = ',I12,/)

 2888 FORMAT(' *WARNING    : CANNOT OUTPUT CUSERIN BD CARDS SINCE THE DOF SET REQUESTED ("',A,'") IS NOT ONE PROGRAMMED')

 9998 FORMAT(' *ERROR  9998: COMM ',I3,' INDICATES UNSUCCESSFUL LINK ',I2,' COMPLETION.'                                           &
                    ,/,14X,' FATAL ERROR - CANNOT START LINK ',I2)

12345 FORMAT(A,10X,A)

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE DEALLOCATE_LINK2_ARRAYS ( WHICH )

      IMPLICIT NONE

      INTEGER(LONG), INTENT(IN)       :: WHICH             ! Which set of matrices to deallocate

! **********************************************************************************************************************************
      CALL DEALLOCATE_SPARSE_MAT ( 'CG_LTM'  )
      CALL DEALLOCATE_SPARSE_MAT ( 'DLR'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'DLRt'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'HMN'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'IF_LTM'  )
      CALL DEALLOCATE_SPARSE_MAT ( 'IRR'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KAA'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KAAD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KAO'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KAOD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KFF'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KFFD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KFS'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KFSD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KFSDe'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'KFSe'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KLLs'    )
      IF (SOL_NAME(1:8) /= 'BUCKLING') THEN
         CALL DEALLOCATE_SPARSE_MAT ( 'KLLD'    )
         CALL DEALLOCATE_SPARSE_MAT ( 'KLLDn'   )
      ENDIF
      CALL DEALLOCATE_SPARSE_MAT ( 'KMM'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KMMD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KMN'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KMND'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KMSM'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KMSMn'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'KMSMs'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'KNM'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KNMD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KNN'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KNND'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KOO'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KOOD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KOODs'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'KOOs'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KRL'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KRLD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KRR'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KRRcb'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'KRRcbn'  )
      CALL DEALLOCATE_SPARSE_MAT ( 'KRRcbs'  )
      CALL DEALLOCATE_SPARSE_MAT ( 'KRRD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KSF'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KSFD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KSS'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'KSSD'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KSSDe'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'KSSe'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'KXX  '   )
      CALL DEALLOCATE_SPARSE_MAT ( 'LMN'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'LTM'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MAA'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MAO'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MFF'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MFS'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MLLn'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'MLLs'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'MLR'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MMM'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MMN'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MNM'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MNN'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MOO'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MPF0 '   )
      CALL DEALLOCATE_SPARSE_MAT ( 'MRL'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MRN  '   )
      CALL DEALLOCATE_SPARSE_MAT ( 'MRR'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MRRcb'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'MRRcbn'  )
      CALL DEALLOCATE_SPARSE_MAT ( 'MSF'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MSS'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'MXX  '   )
      CALL DEALLOCATE_SPARSE_MAT ( 'MXXn'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'PA'      )
      CALL DEALLOCATE_SPARSE_MAT ( 'PF'      )
      CALL DEALLOCATE_SPARSE_MAT ( 'PF_TMP'  )
      CALL DEALLOCATE_SPARSE_MAT ( 'PFYS'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'PFYS1'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'PHIXA'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'PHIXG'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'PHIZL'   )
      CALL DEALLOCATE_SPARSE_MAT ( 'PHIZL1'  )
      CALL DEALLOCATE_SPARSE_MAT ( 'PHIZL1t' )
      CALL DEALLOCATE_SPARSE_MAT ( 'PHIZL2'  )
      CALL DEALLOCATE_SPARSE_MAT ( 'PM'      )
      CALL DEALLOCATE_SPARSE_MAT ( 'PN'      )
      CALL DEALLOCATE_SPARSE_MAT ( 'PO'      )
      CALL DEALLOCATE_SPARSE_MAT ( 'PR'      )
      CALL DEALLOCATE_SPARSE_MAT ( 'PS'      )
      CALL DEALLOCATE_SPARSE_MAT ( 'QSYS'    )
      CALL DEALLOCATE_SPARSE_MAT ( 'RMM'     )
      CALL DEALLOCATE_SPARSE_MAT ( 'RMN'     )

      IF (WHICH  == 1) THEN
         IF (SOL_NAME(1:8) /= 'BUCKLING') THEN
            CALL DEALLOCATE_SPARSE_MAT ( 'KLL'     )
            CALL DEALLOCATE_SPARSE_MAT ( 'MLL'     )
         ENDIF
         CALL DEALLOCATE_SPARSE_MAT ( 'PL'      )
      ELSE IF (WHICH == 2) THEN
         CALL DEALLOCATE_SPARSE_MAT ( 'GMN'     )
         CALL DEALLOCATE_SPARSE_MAT ( 'GMNt'    )
         CALL DEALLOCATE_SPARSE_MAT ( 'GOA'     )
         CALL DEALLOCATE_SPARSE_MAT ( 'GOAt'    )
         CALL DEALLOCATE_SPARSE_MAT ( 'KGG'     )
         CALL DEALLOCATE_SPARSE_MAT ( 'KGGD'    )
         CALL DEALLOCATE_SPARSE_MAT ( 'MGG'     )
         CALL DEALLOCATE_SPARSE_MAT ( 'PG'      )
         CALL DEALLOCATE_SPARSE_MAT ( 'RMG'     )
      ENDIF

! **********************************************************************************************************************************

      END SUBROUTINE DEALLOCATE_LINK2_ARRAYS

      END SUBROUTINE LINK2
