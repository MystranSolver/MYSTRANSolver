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

      SUBROUTINE LINK1_RESTART_DATA

! Reads data from files LINK1B, LINK1G, LINK1K, LINK1Q, LINK1Y (created in LINK1) needed in LINK1 restart

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  ERR, F06   , FILE_NAM_MAXLEN,                                                        &
                                         L1B    , L1G    , L1K    , L1Q   , L1Y     ,                                              &
                                         LINK1B , LINK1G , LINK1K , LINK1Q, LINK1Y  ,                                              &
                                         L1B_MSG, L1G_MSG, L1K_MSG, L1Q_MSG, L1Y_MSG,                                              &
                                         L1BSTAT, L1GSTAT, L1KSTAT, L1YSTAT

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, DATA_NAM_LEN, MCORD, MRCORD, MGRID, MRGRID, NBAROFF, NCORD,                 &
                                         NCONM2, NEDAT, NELE, NGRID, NMATANGLE, NMATL, NPBAR, NPBEAM, NPDAT, NPELAS,NPROD, NPSHEL, &
                                         NPSOLID, NPUSER1, NSEQ, NSUB, NTCARD, NTDAT, NTSUB, NVVEC, MCONM2, MMATL, MPBAR, MPBEAM,  &
                                         MPELAS, MPLOAD4_3D_DATA, MPROD, MPSHEL, MPSOLID, MPUSER1, MRCONM2, MRMATLC, MRPBAR,       &
                                         MRPBEAM, MRPELAS, MRPROD, MRPSHEL, MRPUSER1, NPLATEOFF, NPLOAD4_3D, NPUSERIN
      USE SCONTR, ONLY                :  MPBEAM_STATIONS  

      USE SCONTR, ONLY                :  NPLATETHICK, NPBUSH, MPBUSH, MRPBUSH, NPSHEAR, MPSHEAR, MRPSHEAR, NPCOMP, MPCOMP0,        &
                                         MPCOMP_PLIES, MRPCOMP0, MRPCOMP_PLIES, MPUSERIN, MUSERIN_MAT_NAMES

      USE TIMDAT, ONLY                :  TSEC
      USE PARAMS, ONLY                :  CBMIN3, CBMIN4, IORQ1M, IORQ1S, IORQ1B, IORQ2B, IORQ2T
      USE MODEL_STUF, ONLY            :  CORD, RCORD
      USE MODEL_STUF, ONLY            :  CONM2, RCONM2
      USE MODEL_STUF, ONLY            :  BAROFF, EDAT, EOFF, EPNT, ESORT1, ESORT2, ETYPE, PLATEOFF, VVEC
      USE MODEL_STUF, ONLY            :  GRID, RGRID, GRID_ID, GRID_SEQ, INV_GRID_SEQ

      USE MODEL_STUF, ONLY            :  MATL, RMATL, PBAR, RPBAR, PBEAM, RPBEAM, PELAS, RPELAS, PPNT, PDATA,                      &
                                         PLOAD4_3D_DATA, PTYPE, PROD, RPROD, PSHEL, PSOLID, RPSHEL, PUSER1, RPUSER1, MATANGLE

      USE MODEL_STUF, ONLY            :  GTEMP, TDATA, TPNT
      USE MODEL_STUF, ONLY            :  PLATETHICK, PBUSH, RPBUSH, PSHEAR, RPSHEAR, PCOMP, RPCOMP, PUSERIN, USERIN_MAT_NAMES
      USE MODEL_STUF, ONLY            :  PBEAM_NSTATIONS, PBEAM_XL, PBEAM_RPROPS		
      USE LINK1_RESTART_DATA_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'LINK1_RESTART_DATA'
      CHARACTER(FILE_NAM_MAXLEN*BYTE) :: FILNAM            ! Name of a file that is to be opened for reading
      CHARACTER(132*BYTE)             :: MESSAG            ! Char message for file name
      CHARACTER(LEN=DATA_NAM_LEN)     :: NAME_Is           ! Name of data actually read from file
      CHARACTER(LEN=DATA_NAM_LEN)     :: NAME_ShouldBe     ! Name of data that should be read from file

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: INT2              ! An integer value read from a file from a file
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when opening or reading a file
      INTEGER(LONG)                   :: OUNT(2)           ! File units to write messages to. Input to subr UNFORMATTED_OPEN
      INTEGER(LONG)                   :: PCOMP_PLIES       ! Number of plies in 1 PCOMP entry incl sym plies not explicitly defined
      INTEGER(LONG)                   :: REC_NO            ! Record number of a record read from a file
      INTEGER(LONG)                   :: UNT               ! Unit number of a file to be read




! **********************************************************************************************************************************
! Make units for writing errors the error file and output file

      OUNT(1) = ERR
      OUNT(2) = F06
!-----------------------------------------------------------------------------------------------------------------------------------
! Open L1B and read data

      FILNAM = LINK1B
      UNT    = L1B
      MESSAG = L1B_MSG

      CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

! Read GRID, RGRID data

      NAME_ShouldBe = 'GRID, RGRID'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1B, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NGRID)     CALL DATA_SET_SIZE_ERROR ( LINK1B, NAME_Is, 'NGRID', NGRID, INT2 )

      DO I=1,NGRID
         DO J=1,MGRID
            READ(UNT,IOSTAT=IOCHK) GRID(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRGRID
            READ(UNT,IOSTAT=IOCHK) RGRID(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO
      DO I=1,NGRID
         GRID_ID(I) = GRID(I,1)
      ENDDO

! Read CORD, RCORD data

      NAME_ShouldBe = 'COORDINATE SYSTEM DATA'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1B, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NCORD)     CALL DATA_SET_SIZE_ERROR ( LINK1B, NAME_Is, 'NCORD', NCORD, INT2 )

      DO I=1,NCORD
         DO J=1,MCORD
            READ(UNT,IOSTAT=IOCHK) CORD(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRCORD
            READ(UNT,IOSTAT=IOCHK) RCORD(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read GRID_SEQ, INV_GRID_SEQ data

      NAME_ShouldBe = 'GRID_SEQ, INV_GRID_SEQ'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, FILNAM, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NGRID)     CALL DATA_SET_SIZE_ERROR ( FILNAM, NAME_Is, 'NGRID', NGRID, INT2 )

      DO I=1,NGRID
         READ(UNT,IOSTAT=IOCHK) GRID_SEQ(I), INV_GRID_SEQ(I)                   ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

      CALL FILE_CLOSE ( L1B, LINK1B, L1BSTAT )

!-----------------------------------------------------------------------------------------------------------------------------------
! Open L1G

      FILNAM = LINK1G
      UNT    = L1G
      MESSAG = L1G_MSG

      CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

! Read ETYPE, EPNT, ESORT1 ESORT,2, EOFF

      NAME_ShouldBe = 'ETYPE, EPNT, ESORT1, ESORT2, EOFF'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NELE) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NELE', NELE, INT2 )
      DO I = 1,NELE

         READ(UNT,IOSTAT=IOCHK) ETYPE(I)                                       ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) EPNT(I)                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) ESORT1(I)                                      ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) ESORT2(I)                                      ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

         READ(UNT,IOSTAT=IOCHK) EOFF(I)                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      ENDDO

! Read EDAT

      NAME_ShouldBe = 'EDAT'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NEDAT) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NEDAT', NEDAT, INT2 )
      DO I = 1,NEDAT
         READ(UNT,IOSTAT=IOCHK) EDAT(I)                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read element PARAMETERS

      NAME_ShouldBe = 'ELEM PARAMETERS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) IORQ1M                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) IORQ1S                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) IORQ1B                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) IORQ2B                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) IORQ2T                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) CBMIN3                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

      READ(UNT,IOSTAT=IOCHK) CBMIN4                                            ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )

! Read BAR, BEAM v vectors

      NAME_ShouldBe = 'V VECTORS IN GLOBAL COORDS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NVVEC) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NVVEC', NVVEC, INT2 )
      DO I=1,NVVEC
         READ(UNT,IOSTAT=IOCHK) (VVEC(I,J),J=1,3)                              ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read BAR, BEAM offsets

      NAME_ShouldBe = 'BAR, BEAM OFFSETS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NBAROFF) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NBAROFF', NBAROFF, INT2 )
      DO I = 1,NBAROFF
         DO J = 1,6
            READ(UNT,IOSTAT=IOCHK) BAROFF(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
          ENDDO
      ENDDO

! Read plate offsets

      NAME_ShouldBe = 'PLATE OFFSETS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPLATEOFF) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPLATEOFF', NPLATEOFF, INT2 )
      DO I = 1,NPLATEOFF
         READ(UNT,IOSTAT=IOCHK) PLATEOFF(I)                                    ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read plate thicknesses from connection entries

      NAME_ShouldBe = 'PLATE THICKNESSES FROM CONNECTION ENTRIES'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPLATETHICK) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPLATETHICK', NPLATETHICK, INT2 )
      DO I = 1,NPLATETHICK
         READ(UNT,IOSTAT=IOCHK) PLATETHICK(I)                                  ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

! Read PBAR, RPBAR

      NAME_ShouldBe = 'PBAR, RPBAR'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPBAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPBAR', NPBAR, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPBAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPBAR', MPBAR, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPBAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPBAR', MRPBAR, INT2 )

      DO I=1,NPBAR
         DO J=1,MPBAR
            READ(UNT,IOSTAT=IOCHK) PBAR(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPBAR
            READ(UNT,IOSTAT=IOCHK) RPBAR(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PBEAM, RPBEAM

      NAME_ShouldBe = 'PBEAM, RPBEAM'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPBEAM) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPBEAM', NPBEAM, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPBEAM) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPBEAM', MPBEAM, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPBEAM) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPBEAM', MRPBEAM, INT2 )

      DO I=1,NPBEAM
         DO J=1,MPBEAM
            READ(UNT,IOSTAT=IOCHK) PBEAM(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         READ(UNT,IOSTAT=IOCHK) PBEAM_NSTATIONS(I)                             ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         DO J=1,MPBEAM_STATIONS
            READ(UNT,IOSTAT=IOCHK) PBEAM_XL(I,J)                               ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MPBEAM_STATIONS
            READ(UNT,IOSTAT=IOCHK) PBEAM_RPROPS(I,J,1)                         ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            READ(UNT,IOSTAT=IOCHK) PBEAM_RPROPS(I,J,2)                         ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            READ(UNT,IOSTAT=IOCHK) PBEAM_RPROPS(I,J,3)                         ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            READ(UNT,IOSTAT=IOCHK) PBEAM_RPROPS(I,J,4)                         ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            READ(UNT,IOSTAT=IOCHK) PBEAM_RPROPS(I,J,5)                         ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            READ(UNT,IOSTAT=IOCHK) PBEAM_RPROPS(I,J,6)                         ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPBEAM
            READ(UNT,IOSTAT=IOCHK) RPBEAM(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PBUSH, RPBUSH

      NAME_ShouldBe = 'PBUSH, RPBUSH'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPBUSH) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPBUSH', NPBUSH, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPBUSH) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPBUSH', MPBUSH, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPBUSH) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPBUSH', MRPBUSH, INT2 )

      DO I=1,NPBUSH
         DO J=1,MPBUSH
            READ(UNT,IOSTAT=IOCHK) PBUSH(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPBUSH
            READ(UNT,IOSTAT=IOCHK) RPBUSH(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PROD, RPROD

      NAME_ShouldBe = 'PROD, RPROD'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPROD) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPROD', NPROD, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPROD) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPROD', MPROD, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPROD) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPROD', MRPROD, INT2 )

      DO I = 1,NPROD
         DO J=1,MPROD
            READ(UNT,IOSTAT=IOCHK) PROD(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPROD
            READ(UNT,IOSTAT=IOCHK) RPROD(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PELAS, RPELAS

      NAME_ShouldBe = 'PELAS, RPELAS'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPELAS) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPELAS', NPELAS, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPELAS) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPELAS', MPELAS, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPELAS) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPELAS', MRPELAS, INT2 )

      DO I = 1,NPELAS
         DO J=1,MPELAS
            READ(UNT,IOSTAT=IOCHK) PELAS(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPELAS
            READ(UNT,IOSTAT=IOCHK) RPELAS(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PSHEAR, RPSHEAR

      NAME_ShouldBe = 'PSHEAR, RPSHEAR'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPSHEAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPSHEAR', NPSHEAR, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPSHEAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPSHEAR', MPSHEAR, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPSHEAR) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPSHEAR', MRPSHEAR, INT2 )

      DO I = 1,NPSHEAR
         DO J=1,MPSHEAR
            READ(UNT,IOSTAT=IOCHK) PSHEAR(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPSHEAR
            READ(UNT,IOSTAT=IOCHK) RPSHEAR(I,J)                                ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PSHEL, RPSHEL

      NAME_ShouldBe = 'PSHEL, RPSHEL'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPSHEL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPSHEL', NPSHEL, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPSHEL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPSHEL', MPSHEL, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPSHEL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPSHEL', MRPSHEL, INT2 )

      DO I = 1,NPSHEL
         DO J=1,MPSHEL
            READ(UNT,IOSTAT=IOCHK) PSHEL(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPSHEL
            READ(UNT,IOSTAT=IOCHK) RPSHEL(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PCOMP, RPCOMP

      NAME_ShouldBe = 'PCOMP, RPCOMP'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPCOMP) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPCOMP', NPCOMP, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPCOMP0) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPCOMP0', MPCOMP0, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPCOMP_PLIES) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPCOMP_PLIES', MPCOMP_PLIES, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPCOMP0) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPCOMP0', MRPCOMP0, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPCOMP_PLIES) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPCOMP_PLIES', MRPCOMP_PLIES, INT2 )

      DO I = 1,NPCOMP

         READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         PCOMP_PLIES = INT2

         DO J=1,MPCOMP0+MPCOMP_PLIES*PCOMP_PLIES
            READ(UNT,IOSTAT=IOCHK) PCOMP(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO

         DO J=1,MRPCOMP0+MRPCOMP_PLIES*PCOMP_PLIES
            READ(UNT,IOSTAT=IOCHK) RPCOMP(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO

      ENDDO

! Read PSOLID

      NAME_ShouldBe = 'PSOLID'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPSOLID) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPSOLID', NPSOLID, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPSOLID) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPSOLID', MPSOLID, INT2 )

      DO I = 1,NPSOLID
         DO J=1,MPSOLID
            READ(UNT,IOSTAT=IOCHK) PSOLID(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PUSER1, RPUSER1

      NAME_ShouldBe = 'PUSER1, RPUSER1'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPUSER1) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPUSER1', NPUSER1, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPUSER1) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPUSER1', MPUSER1, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRPUSER1) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRPUSER1', MRPUSER1, INT2 )

      DO I = 1,NPUSER1
         DO J=1,MPUSER1
            READ(UNT,IOSTAT=IOCHK) PUSER1(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRPUSER1
            READ(UNT,IOSTAT=IOCHK) RPUSER1(I,J)                                ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read PUSERIN

      NAME_ShouldBe = 'PUSERIN'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPUSERIN) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPUSERIN', NPUSERIN, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MPUSERIN) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MPUSERIN', MPUSERIN, INT2 )

      DO I = 1,NPUSERIN
         DO J=1,MPUSERIN
            READ(UNT,IOSTAT=IOCHK) PUSERIN(I,J)                                ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read USERIN_MAT_NAMES

      NAME_ShouldBe = 'USERIN_MAT_NAMES'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NPUSERIN) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NPUSERIN', NPUSERIN, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MUSERIN_MAT_NAMES) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MUSERIN_MAT_NAMES', MUSERIN_MAT_NAMES, INT2 )

      DO I = 1,NPUSERIN
         DO J=1,MUSERIN_MAT_NAMES
            READ(UNT,IOSTAT=IOCHK) USERIN_MAT_NAMES(I,J)                       ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read material data from L1G

      NAME_ShouldBe = 'MATL, RMATL'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NMATL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NMATL', NMATL, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MMATL) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MMATL', MMATL, INT2 )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= MRMATLC) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'MRMATLC', MRMATLC, INT2 )

      DO I = 1,NMATL
         DO J=1,MMATL
            READ(UNT,IOSTAT=IOCHK) MATL(I,J)                                   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRMATLC
            READ(UNT,IOSTAT=IOCHK) RMATL(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO

! Read material property angles

      NAME_ShouldBe = 'MATERIAL PROPERTY ANGLES'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1G, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NMATANGLE) CALL DATA_SET_SIZE_ERROR ( LINK1G, NAME_Is, 'NMATANGLE', NMATANGLE, INT2 )
      DO I = 1,NMATANGLE
         READ(UNT,IOSTAT=IOCHK) MATANGLE(I)                                    ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      ENDDO

      CALL FILE_CLOSE ( L1G, LINK1G, L1GSTAT )

!-----------------------------------------------------------------------------------------------------------------------------------
! Open L1K and read data

      FILNAM = LINK1K
      UNT    = L1K
      MESSAG = L1K_MSG

      IF (NTCARD > 0) THEN

         CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

! Read TPNT

         NAME_ShouldBe = 'TPNT'
         REC_NO = 0

         READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1K, NAME_Is )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NELE) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NELE', NELE, INT2 )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NTSUB) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NTSUB', NTSUB, INT2 )

         DO I=1,NELE
            DO J=1,NTSUB
               READ(UNT,IOSTAT=IOCHK) TPNT(I,J)                                ; REC_NO = REC_NO + 1
               CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            ENDDO
         ENDDO

! Read TDATA

         NAME_ShouldBe = 'TDATA'
         REC_NO = 0

         READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1K, NAME_Is )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NTDAT) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NTDAT', NTDAT, INT2 )

         DO I=1,NTDAT
            READ(UNT,IOSTAT=IOCHK) TDATA(I)                                    ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO

! Read GTEMP

         NAME_ShouldBe = 'GTEMP'
         REC_NO = 0

         READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1K, NAME_Is )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NGRID) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NGRID', NGRID, INT2 )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NTSUB) CALL DATA_SET_SIZE_ERROR ( LINK1K, NAME_Is, 'NTSUB', NTSUB, INT2 )

         DO I=1,NGRID
            DO J=1,NTSUB
               READ(UNT,IOSTAT=IOCHK) GTEMP(I,J)                               ; REC_NO = REC_NO + 1
               CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            ENDDO
         ENDDO

         CALL FILE_CLOSE ( L1K, LINK1K, L1KSTAT )

      ENDIF

!-----------------------------------------------------------------------------------------------------------------------------------
! Open L1Q and read data

      FILNAM = LINK1Q
      UNT    = L1Q
      MESSAG = L1Q_MSG

      IF (NTCARD > 0) THEN

         CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

! Read PPNT

         NAME_ShouldBe = 'PPNT'
         REC_NO = 0

         READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1Q, NAME_Is )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NELE) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NELE', NELE, INT2 )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NSUB) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NSUB', NSUB, INT2 )

         DO I=1,NELE
            DO J=1,NSUB
               READ(UNT,IOSTAT=IOCHK) PPNT(I,J)                                ; REC_NO = REC_NO + 1
               CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
            ENDDO
         ENDDO

! Read PDATA

         NAME_ShouldBe = 'PDATA'
         REC_NO = 0

         READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1Q, NAME_Is )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NPDAT) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NPDAT', NPDAT, INT2 )

         DO I=1,NPDAT
            READ(UNT,IOSTAT=IOCHK) PDATA(I)                                    ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO

! Read PTYPE

         NAME_ShouldBe = 'PTYPE'
         REC_NO = 0

         READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1Q, NAME_Is )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NPDAT) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NPDAT', NPDAT, INT2 )

         DO I=1,NELE
            READ(UNT,IOSTAT=IOCHK) PTYPE(I)                                    ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO

      ENDIF

! Read PLOAD4_3D_DATA

         NAME_ShouldBe = 'PLOAD4_3D_DATA'
         REC_NO = 0

         READ(UNT,IOSTAT=IOCHK) NAME_Is                                        ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (NAME_Is /= NAME_ShouldBe) CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1Q, NAME_Is )

         READ(UNT,IOSTAT=IOCHK) INT2                                           ; REC_NO = REC_NO + 1
         CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         IF (INT2 /= NPLOAD4_3D) CALL DATA_SET_SIZE_ERROR ( LINK1Q, NAME_Is, 'NPLOAD4_3D', NPLOAD4_3D, INT2 )

         DO I=1,NPLOAD4_3D
            READ(UNT,IOSTAT=IOCHK) (PLOAD4_3D_DATA(I,J),J=1,MPLOAD4_3D_DATA)   ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO

!-----------------------------------------------------------------------------------------------------------------------------------
! Open L1Y and read data

      FILNAM = LINK1Y
      UNT    = L1Y
      MESSAG = L1Y_MSG

      CALL FILE_OPEN ( UNT, FILNAM, OUNT, 'OLD', MESSAG, 'READ_STIME', 'UNFORMATTED', 'READ', 'REWIND', 'Y', 'N' )

! Read CONM2, RCONM2 data

      NAME_ShouldBe = 'CONM2, RCONM2'
      REC_NO = 0

      READ(UNT,IOSTAT=IOCHK) NAME_Is                                           ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (NAME_Is /= NAME_ShouldBe)     CALL DATA_SET_NAME_ERROR ( NAME_ShouldBe, LINK1Y, NAME_Is )

      READ(UNT,IOSTAT=IOCHK) INT2                                              ; REC_NO = REC_NO + 1
      CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
      IF (INT2 /= NCONM2)     CALL DATA_SET_SIZE_ERROR ( LINK1Y, NAME_Is, 'NCONM2', NCONM2, INT2 )

      DO I=1,NCONM2
         DO J=1,MCONM2
            READ(UNT,IOSTAT=IOCHK) CONM2(I,J)                                  ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
         DO J=1,MRCONM2
            READ(UNT,IOSTAT=IOCHK) RCONM2(I,J)                                 ; REC_NO = REC_NO + 1
            CALL READ_CHK ( IOCHK, FILNAM, NAME_ShouldBe, REC_NO, OUNT )
         ENDDO
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE LINK1_RESTART_DATA
