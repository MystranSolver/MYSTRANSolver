! ##################################################################################################################################

      MODULE LAPACK_STD_EIG_1
! --- lapack_surgery begin --- !

      USE PENTIUM_II_KIND, ONLY          :  BYTE, LONG, DOUBLE
      USE LAPACK_STD_EIG_1_HELPER, ONLY  :  DSYEV_HELPER  => DSYEV,                                             &
     &                                      DSYTRD_HELPER => DSYTRD,                                            &
     &                                      DORGTR_HELPER => DORGTR

      CONTAINS

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DSYEV( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      CHARACTER          JOBZ, UPLO
      INTEGER            INFO, LDA, LWORK, N
      REAL(DOUBLE)       A( LDA, * ), W( * ), WORK( * )
      CALL DSYEV_HELPER( JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO )
      RETURN
      END SUBROUTINE DSYEV
! --- lapack_peeloff end --- !

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DSYTRD( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
      REAL(DOUBLE)       A( LDA, * ), D( * ), E( * ), TAU( * ),
     $                   WORK( * )
      CALL DSYTRD_HELPER( UPLO, N, A, LDA, D, E, TAU, WORK, LWORK,
     $                    INFO )
      RETURN
      END SUBROUTINE DSYTRD
! --- lapack_peeloff end --- !

! ##################################################################################################################################

! --- lapack_peeloff begin --- !
      SUBROUTINE DORGTR( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
      CHARACTER          UPLO
      INTEGER            INFO, LDA, LWORK, N
      REAL(DOUBLE)       A( LDA, * ), TAU( * ), WORK( LWORK )
      CALL DORGTR_HELPER( UPLO, N, A, LDA, TAU, WORK, LWORK, INFO )
      RETURN
      END SUBROUTINE DORGTR
! --- lapack_peeloff end --- !

! --- lapack_surgery end --- !
      END MODULE LAPACK_STD_EIG_1
