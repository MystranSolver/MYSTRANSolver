! ##################################################################################################################################
! Begin MIT license text.
! _______________________________________________________________________________________________________

! Copyright 2022 Dr William R Case, Jr (mystransolver@gmail.com)
! Copyright 2026 Bruno Borber (pre-factored matrix support for adaptive Lanczos)

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

      MODULE DSBAND_PREFAC_MOD
      IMPLICIT NONE
      PRIVATE
      PUBLIC :: DSBAND_PREFAC
      CONTAINS
c
c     DSBAND_PREFAC: Version of DSBAND that assumes the matrix has already been factored.
c
c     For BANDED solver: RFAC must contain the factored form (output from dgbtrf/dpbtrf)
c                        IWORK must contain the pivot indices from factorization
c     For SPARSE solver: SLU_FACTORS must be set up by prior call to SYM_MAT_DECOMP_SUPRLU
c
c     This routine is used by EIG_LANCZOS_ARPACK_ADAPTIVE to avoid re-factoring the same
c     matrix when iterating with increasing NEV values.
c
      subroutine dsband_prefac( rvec, howmny, select, d, z, ldz, sigma,
     &           n, lda, rfac, kl, ku, which, bmat, nev,
     &           tol, resid, ncv, v, ldv, iparam, workd, workl,
     &           lworkl, iwork, info, info_lapack, dtbsv_msg, piters)

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06, SC1
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, SOL_NAME,
     &                                   NTERM_KLLDn, NTERM_MLLn,
     &                                   NTERM_KMSMn, NTERM_ALL
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  EIG_MSGLVL, EIG_LAP_MAT_TYPE
      USE SuperLU_STUF, ONLY          :  SLU_FACTORS, SLU_INFO
      USE PARAMS, ONLY                :  SOLLIB
      USE SPARSE_MATRICES, ONLY       :  I_KLLDn, J_KLLDn, KLLDn,
     &                                   I_MLLn , J_MLLn , MLLn,
     &                                   I_KMSMn, J_KMSMn, KMSMn

      USE ARPACK_LANCZOS_EIG, ONLY    :  dsaupd, dseupd, cr13_a

      USE SPARSE_FULL_MULTIPLICATION, ONLY: MATMULT_SFF
      USE ARPACK_DIAGNOSTICS, ONLY : ARPACK_INFO_MSG
      USE SUPERLU_ADAPTERS, ONLY     :  FBS_SUPRLU

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'DSBAND_PREFAC'

c
c     %------------------%
c     | Scalar Arguments |
c     %------------------%
c
      character        which*2, bmat, howmny
      character*1      dtbsv_msg
      character*80     pagetitle
      integer          n, lda, kl, ku, nev, ncv, ldv,
     &                 ldz, lworkl, info
      integer          info_lapack, dsaupd_loop_count
      real(double)
     &                 tol, sigma
      logical          rvec
c
c     %-----------------%
c     | Array Arguments |
c     %-----------------%
c
      integer          iparam(*), iwork(*)
      logical          select(*)
      real(double)  :: workd1(n), workd2(n)
      real(double)
     &                 d(*), resid(*), v(ldv,*), z(ldz,*),
     &                 rfac(lda,*), workd(*), workl(*)
      integer nev_user
      real(double) eig_old(n)

      integer iter, iter_old

      integer piters

c
c     %--------------%
c     | Local Arrays |
c     %--------------%
c
      integer          ipntr(14)
      REAL(DOUBLE)     KLLDn_DIAG(N), MLLn_DIAG(N), KMSMn_DIAG(N)
c
c     %---------------%
c     | Local Scalars |
c     %---------------%
c
      integer          ido, i, j, type, imid, itop, ibot, ierr
c
c     %------------%
c     | Parameters |
c     %------------%
c
      REAL(DOUBLE)
     &                  one, zero
      parameter        (one = 1.0, zero = 0.0)
c
c     %-----------------------%
c     | Executable Statements |
c     %-----------------------%
c

      ierr = 0
! **********************************************************************************************************************************

      do j=1,n
         eig_old(j) = 0.d0
      enddo

      nev_user = nev

c     %----------------------------------------------------------------%
c     | Set type of the problem to be solved.                          |
c     |                                                                |
c     | type = 3 --> Solving generalized problem in regular mode.      |
c     | type = 4 --> Solving generalized problem in shift-invert mode. |
c     %----------------------------------------------------------------%
c
      if      ( iparam(7) .eq. 1 ) then
         type = 1                                                            ! ARPACK "mode" = 1
      else if ( iparam(7) .eq. 2 ) then
         type = 3                                                            ! ARPACK "mode" = 2
      else if ( iparam(7) .eq. 3 .and. bmat .eq. 'G') then
         type = 4                                                            ! ARPACK "mode" = 3
      else                                                                   ! Err msg given in BD_EIGRL
         go to 9000
      end if
c
c     %------------------------%
c     | Initialize the reverse |
c     | communication flag.    |
c     %------------------------%
c
      ido   = 0
c
c     %----------------%
c     | Exact shift is |
c     | used.          |
c     %----------------%
c
      iparam(1) = 1
c
c     %-----------------------------------%
c     | Both matrices A and M are stored  |
c     | between rows itop and ibot.  Imid |
c     | is the index of the row that      |
c     | stores the diagonal elements.     |
c     %-----------------------------------%
c
      if      (eig_lap_mat_type(1:3) == 'DGB') then
         itop = kl + 1
         ibot = 2*kl + ku + 1
         imid = kl + ku + 1
      else if (eig_lap_mat_type(1:3) == 'DPB') then
         itop = 1
         ibot = ku + 1
         imid = ibot
      endif

c     %-------------------------------------------------------------------%
c     | FACTORIZATION IS SKIPPED IN THIS VERSION.                         |
c     | The caller must ensure that:                                      |
c     |   - For BANDED: RFAC contains factored matrix, IWORK has pivots   |
c     |   - For SPARSE: SLU_FACTORS is set up from SYM_MAT_DECOMP_SUPRLU  |
c     %-------------------------------------------------------------------%

      IF (EIG_MSGLVL > 0) THEN
         IF (SOL_NAME(1:8) == 'BUCKLING') THEN
            DO I=1,N
               IF (I_KLLDn(I) == I_KLLDn(I+1)) THEN
                  KLLDn_DIAG(I) = ZERO
               ELSE
                  KLLDn_DIAG(I) = KLLDn(I_KLLDn(I))
               ENDIF
               IF (I_KMSMn(I) == I_KMSMn(I+1)) THEN
                  KMSMn_DIAG(I) = ZERO
               ELSE
                  KMSMn_DIAG(I) = KMSMn(I_KMSMn(I))
               ENDIF
            ENDDO
         ELSE
            DO I=1,N
               IF (I_MLLn(I) == I_MLLn(I+1)) THEN
                  MLLn_DIAG(I) = ZERO
               ELSE
                  MLLn_DIAG(I) = MLLn(I_MLLn(I))
               ENDIF
               IF (I_KMSMn(I) == I_KMSMn(I+1)) THEN
                  KMSMn_DIAG(I) = ZERO
               ELSE
                  KMSMn_DIAG(I) = KMSMn(I_KMSMn(I))
               ENDIF
            ENDDO
         ENDIF
      ENDIF
c
c     %--------------------------------------------%
c     |  M A I N   L O O P (reverse communication) |
c     %--------------------------------------------%
c
      iter_old          = 0
      dsaupd_loop_count = 0
      write(sc1, * )
  90  continue

      IF (EIG_MSGLVL > 0) THEN
         WRITE(F06,99990)
         WRITE(F06,98710) dsaupd_loop_count, ido, type,
     &                    eig_lap_mat_type(1:3)
      ENDIF
c
      iter = 0
      call dsaupd ( ido, bmat, n, which, nev, tol, resid, ncv,
     &              v, ldv, iparam, ipntr, workd, workl, lworkl,
     &              info, nev_user, sigma, eig_old, iter, piters )
c
      if (iter > iter_old) then
         dsaupd_loop_count = 1
         iter_old = iter
      else
         dsaupd_loop_count = dsaupd_loop_count + 1
      endif
      write(sc1,12345,advance='no') iter+1,dsaupd_loop_count,ido,cr13_a

! **********************************************************************
      if (ido .eq. -1) then
c
         if ( type .eq. 3 ) then
c
c           %-------------------------------------%
c           | Perform  y <--- OP*x = inv[MB]*AB*x |
c           | MB = [K - sigma*M] and AB = M       |
c           | to force the starting vector into   |
c           | the range of OP.                    |
c           %-------------------------------------%
c
            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn', n, n, nterm_klldn, 'N',
     &                             i_klldn, j_klldn, klldn,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn', n, n, nterm_mlln, 'N',
     &                             i_mlln, j_mlln, mlln,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo

            IF (EIG_MSGLVL > 0) CALL ARP_DEB_PREFAC(1,N,IDO,IPNTR,
     &         KLLDn_DIAG,MLLn_DIAG,KMSMn_DIAG,WORKD1,WORKD2,RFAC,IMID)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)

            IF(SOLLIB(1:6) == 'SPARSE') THEN

               SLU_INFO = 0
               call FBS_SUPRLU ( SUBR_NAME, 'KMSMn', n,
     &                        NTERM_KMSMn, I_KMSMn, J_KMSMn, KMSMn,
     &                        0, workd(ipntr(2)), SLU_INFO )

            ELSE

               if      (eig_lap_mat_type(1:3) == 'DGB') then

                  call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                          iwork, workd(ipntr(2)), n, ierr )

               else if (eig_lap_mat_type(1:3) == 'DPB') then

                  call dpbtrs ( 'U', n, ku, 1, rfac, ku+1,
     &                          workd(ipntr(2)), n, ierr )

               endif
               IF (EIG_MSGLVL > 0) CALL ARP_DEB_PREFAC(2,N,IDO,IPNTR,
     &         KLLDn_DIAG,MLLn_DIAG,KMSMn_DIAG,WORKD1,WORKD2,RFAC,IMID)

               if (ierr .ne. 0) then
                  info_lapack = ierr
                  go to 9000
               end if

            ENDIF

c
         else if ( type .eq. 4 ) then
c
c           %-----------------------------------------%
c           | Perform y <-- OP*x                      |
c           |           = inv[A-SIGMA*M]*M            |
c           | to force the starting vector into the   |
c           | range of OP.                            |
c           %-----------------------------------------%
c
            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn' , n, n, nterm_klldn, 'N',
     &                             i_klldn, j_klldn, klldn,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn'  , n, n, nterm_mlln , 'N',
     &                             i_mlln, j_mlln, mlln,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB_PREFAC(1,N,IDO,IPNTR,
     &         KLLDn_DIAG,MLLn_DIAG,KMSMn_DIAG,WORKD1,WORKD2,RFAC,IMID)

            IF(SOLLIB(1:6) == 'SPARSE') THEN

               SLU_INFO = 0
               call FBS_SUPRLU ( SUBR_NAME, 'KMSMn', n,
     &                        NTERM_KMSMn, I_KMSMn, J_KMSMn, KMSMn,
     &                        0, workd(ipntr(2)), SLU_INFO )

            ELSE

               if      (eig_lap_mat_type(1:3) == 'DGB') then

                  call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                         iwork, workd(ipntr(2)), n, ierr )

               else if (eig_lap_mat_type(1:3) == 'DPB') then

                  call dpbtrs ( 'U', n, ku, 1, rfac, ku+1,
     &                           workd(ipntr(2)), n, ierr )
               endif
               IF (EIG_MSGLVL > 0) CALL ARP_DEB_PREFAC(2,N,IDO,IPNTR,
     &         KLLDn_DIAG,MLLn_DIAG,KMSMn_DIAG,WORKD1,WORKD2,RFAC,IMID)

               if (ierr .ne. 0) then
                  info_lapack = ierr
                  go to 9000
               end if

            ENDIF


         endif
c
! **********************************************************************
      else if (ido .eq. 1) then
c
         if ( type .eq. 3 ) then
c
c           %----------------------------------------------%
c           | Perform  y <--- OP*x = inv[MB]*AB*x          |
c           |                      = inv[K - sigma*M]*M*x  |
c           %----------------------------------------------%
c
            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn' , n, n, nterm_klldn, 'N',
     &                             i_klldn, j_klldn, klldn,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn'  , n, n, nterm_mlln , 'N',
     &                             i_mlln , j_mlln, mlln,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB_PREFAC(1,N,IDO,IPNTR,
     &         KLLDn_DIAG,MLLn_DIAG,KMSMn_DIAG,WORKD1,WORKD2,RFAC,IMID)
            call dcopy(n, workd(ipntr(2)), 1, workd(ipntr(1)), 1)

            IF(SOLLIB(1:6) == 'SPARSE') THEN

               SLU_INFO = 0
               call FBS_SUPRLU ( SUBR_NAME, 'KMSMn', n,
     &                        NTERM_KMSMn, I_KMSMn, J_KMSMn, KMSMn,
     &                        0, workd(ipntr(2)), SLU_INFO )

            ELSE

               if      (eig_lap_mat_type(1:3) == 'DGB') then

                  call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                          iwork, workd(ipntr(2)), n, ierr )

               else if (eig_lap_mat_type(1:3) == 'DPB') then

                  call dpbtrs ( 'U', n, ku, 1, rfac, ku+1,
     &                          workd(ipntr(2)), n, ierr )

               endif
               IF (EIG_MSGLVL > 0) CALL ARP_DEB_PREFAC(2,N,IDO,IPNTR,
     &         KLLDn_DIAG,MLLn_DIAG,KMSMn_DIAG,WORKD1,WORKD2,RFAC,IMID)

               if (ierr .ne. 0) then
                  info_lapack = ierr
                  go to 9000
               end if

            ENDIF
c
         else if ( type .eq. 4 ) then
c
c           %----------------------------------------%
c           | Perform y <-- inv[AB-sigma*MB]*(MB*x). |
c           |               inv[K - sigma*M]*m*x     |
c           | (MB*x) has been computed and stored    |
c           | in workd(ipntr(3)).                    |
c           %----------------------------------------%
c

            call dcopy(n, workd(ipntr(3)), 1, workd(ipntr(2)), 1)

            IF(SOLLIB(1:6) == 'SPARSE') THEN

               SLU_INFO = 0
               call FBS_SUPRLU ( SUBR_NAME, 'KMSMn', n,
     &                        NTERM_KMSMn, I_KMSMn, J_KMSMn, KMSMn,
     &                        0, workd(ipntr(2)), SLU_INFO )

            ELSE

               if      (eig_lap_mat_type(1:3) == 'DGB') then

                  call dgbtrs ('Notranspose', n, kl, ku, 1, rfac, lda,
     &                          iwork, workd(ipntr(2)), n, ierr )

               else if (eig_lap_mat_type(1:3) == 'DPB') then

                  call dpbtrs ( 'U', n, ku, 1, rfac, ku+1,
     &                          workd(ipntr(2)), n, ierr )

               endif
               IF (EIG_MSGLVL > 0) CALL ARP_DEB_PREFAC(2,N,IDO,IPNTR,
     &         KLLDn_DIAG,MLLn_DIAG,KMSMn_DIAG,WORKD1,WORKD2,RFAC,IMID)

               if (ierr .ne. 0) then
                  info_lapack = ierr
                  go to 9000
               end if

            ENDIF

c
         end if
c
! **********************************************************************
      else if (ido .eq. 2) then
c
c        %----------------------------------%
c        |        Perform y <-- B*x         |
c        %----------------------------------%
c
         if (type == 3) then

            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo
            call matmult_sff ( 'KMSMn', n, n, nterm_kmsmn, 'N',
     &                          i_kmsmn, j_kmsmn, kmsmn,
     &                         'workd1', n, 1, workd1, 'N',
     &                         'workd2', one,  workd2 )

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB_PREFAC(3,N,IDO,IPNTR,
     &         KLLDn_DIAG,MLLn_DIAG,KMSMn_DIAG,WORKD1,WORKD2,RFAC,IMID)

         else if (type == 4) then

            do j=1,n
               workd1(j) = workd(ipntr(1)+j-1)
            enddo

            if (sol_name(1:8) == 'BUCKLING') then
               call matmult_sff ( 'KLLDn' , n, n, nterm_klldn, 'N',
     &                             i_klldn, j_klldn, klldn,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2',-one,  workd2 )
            else
               call matmult_sff ( 'MLLn'  , n, n, nterm_mlln , 'N',
     &                             i_mlln , j_mlln, mlln,
     &                            'workd1', n, 1, workd1, 'N',
     &                            'workd2', one,  workd2 )
            endif

            do j=1,n
               workd(ipntr(2)+j-1) = workd2(j)
            enddo
            IF (EIG_MSGLVL > 0) CALL ARP_DEB_PREFAC(1,N,IDO,IPNTR,
     &         KLLDn_DIAG,MLLn_DIAG,KMSMn_DIAG,WORKD1,WORKD2,RFAC,IMID)

         endif

c
      else
c
c        %--------------------------------------%
c        | Either we have convergence, or error |
c        %--------------------------------------%
c
         if ( info .lt. 0) then
c
            call arpack_info_msg ('dsaupd',info,iparam,lworkl,nev,ncv)
            go to 9000
c
         else
c
            if ((info == 1) .or. (info == 3)) then
               call arpack_info_msg('dsaupd',info,iparam,lworkl,nev,ncv)
            endif
c
            if (iparam(5) .gt. 0) then
c
               call dseupd ( rvec, 'A', select, d, z, ldz, sigma,
     &                  bmat, n, which, nev, tol, resid, ncv, v, ldv,
     &                  iparam, ipntr, workd, workl, lworkl, info )
c
               if ( info .ne. 0) then
c
              call arpack_info_msg ('dseupd',info,iparam,lworkl,nev,ncv)
                  go to 9000
c
               end if
c
            end if
c
         end if
c
         go to 9000
c
      end if
c
c     %----------------------------------------%
c     | L O O P  B A C K to call DSAUPD again. |
c     %----------------------------------------%
c
      go to 90
c
 9000 continue
c
      Write(f06,*) ! blank line after eigen iteration results
! **********************************************************************************************************************************
12345 format(5X,'Iteration',i4,' Rev comm loop',i4,' with IDO =',i3,a)

 4907 FORMAT(/,22X,A
     &      ,/,7X,'1',12X,'2',12X,'3',12X,'4',12X,'5',12X,'6',12X,
     &        '7',12X,'8',12X,'9',12X,'10')

 4908 FORMAT(10(1X,1ES12.5))

98710 FORMAT(' dsaupd loop count = ',I4,' ido = ',i4,', "type" = ',I3,
     &       ', using ',a,' LAPACK matrices',/)

99990 FORMAT('**********************************************************
     &***************************')

      RETURN

! **********************************************************************************************************************************

! ##################################################################################################################################

      CONTAINS

! ##################################################################################################################################

      SUBROUTINE ARP_DEB_PREFAC ( WHICH_DEB, N_DEB, IDO_DEB, IPNTR_DEB,
     &                    KLLDn_DIAG_DEB, MLLn_DIAG_DEB, KMSMn_DIAG_DEB,
     &                    WORKD1_DEB, WORKD2_DEB, RFAC_DEB, IMID_DEB )

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE SCONTR, ONLY                :  PROG_NAME, FATAL_ERR, WARN_ERR
      USE IOUNT1, ONLY                :  WRT_BUG, WRT_ERR, ERR, F06

      IMPLICIT NONE

      INTEGER, INTENT(IN)            :: IDO_DEB
      INTEGER, INTENT(IN)            :: IPNTR_DEB(14)
      INTEGER, INTENT(IN)            :: N_DEB
      INTEGER, INTENT(IN)            :: WHICH_DEB
      INTEGER, INTENT(IN)            :: IMID_DEB
      INTEGER                        :: II, JJ
      REAL(DOUBLE), INTENT(IN)       :: KLLDn_DIAG_DEB(N_DEB)
      REAL(DOUBLE), INTENT(IN)       :: MLLn_DIAG_DEB(N_DEB)
      REAL(DOUBLE), INTENT(IN)       :: KMSMn_DIAG_DEB(N_DEB)
      REAL(DOUBLE), INTENT(IN)       :: WORKD1_DEB(N_DEB)
      REAL(DOUBLE), INTENT(IN)       :: WORKD2_DEB(N_DEB)
      REAL(DOUBLE), INTENT(IN)       :: RFAC_DEB(LDA,*)

! **********************************************************************************************************************************
      IF (WHICH_DEB == 1) THEN

         WRITE(F06,98700) IDO_DEB, IPNTR_DEB(1), IPNTR_DEB(2)
         IF (EIG_MSGLVL > 1) THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,98711)
               DO JJ=1,N_DEB
                  WRITE(F06,98713) JJ, KLLDn_DIAG_DEB(JJ),
     &                             WORKD1_DEB(JJ), WORKD2_DEB(JJ)
               ENDDO
            ELSE
               WRITE(F06,98712)
               DO JJ=1,N_DEB
                  WRITE(F06,98713) JJ, MLLn_DIAG_DEB(JJ),
     &                             WORKD1_DEB(JJ), WORKD2_DEB(JJ)
               ENDDO
            ENDIF
            WRITE(F06,*)
         ENDIF

      ELSE IF (WHICH_DEB == 2) THEN

         WRITE(F06,98700) IDO_DEB, IPNTR_DEB(1), IPNTR_DEB(2)
         IF (EIG_MSGLVL > 1) THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,98721)
            ELSE
               WRITE(F06,98722)
            ENDIF
            DO JJ=1,N_DEB
               WRITE(F06,98723) JJ, RFAC_DEB(IMID_DEB,JJ),
     &                          WORKD1_DEB(JJ), WORKD2_DEB(JJ)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ELSE IF (WHICH_DEB == 3) THEN

         WRITE(F06,98700) IDO_DEB, IPNTR_DEB(1), IPNTR_DEB(2)
         IF (EIG_MSGLVL > 1) THEN
            IF (SOL_NAME(1:8) == 'BUCKLING') THEN
               WRITE(F06,98731)
            ELSE
               WRITE(F06,98732)
            ENDIF
            DO JJ=1,N_DEB
               WRITE(F06,98733) JJ, KMSMn_DIAG_DEB(JJ),
     &                          WORKD1_DEB(JJ), WORKD2_DEB(JJ)
            ENDDO
            WRITE(F06,*)
         ENDIF

      ENDIF

! **********************************************************************************************************************************
98700 FORMAT(' ido, ipntr(1), ipntr(2) = ',I9,3x,2i3,/)

98711 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF                KLL                  WORKD1(J)      WOR
     &KD2(J)',/,' ---------        -------------------        ----------
     &--   ------------')

98712 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF                MLL                  WORKD1(J)      WOR
     &KD2(J)',/,' ---------        -------------------        ----------
     &--   ------------')

98713 FORMAT(I8,12X,1ES14.6,9X,1ES14.6,1X,1ES14.6)

98721 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF  tridiag factor of [KLL-sig*MLL]    WORKD1(J)      WOR
     &KD2(J)',/,' ---------  -------------------------------  ----------
     &--   ------------')

98722 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF  tridiag factor of [KLL-sig*KLL]    WORKD1(J)      WOR
     &KD2(J)',/,' ---------  -------------------------------  ----------
     &--   ------------')

98723 FORMAT(I8,12X,1ES14.6,9X,1ES14.6,1X,1ES14.6)

98731 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF          [KLL-sigma*MLL]            WORKD1(J)      WOR
     &KD2(J)',/,' ---------        -------------------        ----------
     &--   ------------')

98732 FORMAT('                  Diagonal of matrix:',/,
     &' L-set DOF          [KLL-sigma*KLL]            WORKD1(J)      WOR
     &KD2(J)',/,' ---------        -------------------        ----------
     &--   ------------')

98733 FORMAT(I8,12X,1ES14.6,9X,1ES14.6,1X,1ES14.6)

! **********************************************************************************************************************************

      END SUBROUTINE ARP_DEB_PREFAC

      end subroutine dsband_prefac

      END MODULE DSBAND_PREFAC_MOD
