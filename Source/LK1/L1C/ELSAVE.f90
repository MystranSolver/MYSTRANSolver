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

      SUBROUTINE ELSAVE

! Saves element data to file LINK1G.

      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR,     F06,     L1G
      USE SCONTR, ONLY                :  BLNK_SUB_NAM, DATA_NAM_LEN, MMATL, MPBAR, MPBEAM, MPBEAM_STATIONS, MPBUSH, MPELAS, MPROD, MPSHEL,          &
                                         MPSOLID, MPUSER1,MPUSERIN, MRMATLC, MRPBAR, MRPBEAM, MRPBUSH, MRPELAS, MRPROD, MPSHEAR,   &
                                         MRPSHEAR, MRPSHEL, MRPUSER1, NBAROFF, NBUSHOFF, NEDAT, NELE, NMATANGLE, NMATL, MPCOMP0,   &
                                         MRPCOMP0, MPCOMP_PLIES, MRPCOMP_PLIES, MUSERIN_MAT_NAMES, NPBAR, NPBEAM, NPBUSH, NPCOMP,  &
                                         NPELAS, NPLATEOFF, NPLATETHICK, NPROD, NPSHEAR, NPSHEL, NPSOLID, NPUSER1, NPUSERIN, NVVEC
      USE PARAMS, ONLY                :  CBMIN3, CBMIN4, IORQ1M, IORQ1S, IORQ1B, IORQ2B, IORQ2T
      USE TIMDAT, ONLY                :  TSEC
      USE MODEL_STUF, ONLY            :  BAROFF, BUSHOFF, EDAT, EOFF, EPNT, ESORT1, ESORT2, ETYPE, MATANGLE, MATL, RMATL,PBAR,     &
                                         RPBAR, PBEAM, PBEAM_NSTATIONS, PBEAM_XL, PBEAM_RPROPS, RPBEAM, PBUSH, RPBUSH, PCOMP, RPCOMP, PELAS, RPELAS, PROD, RPROD, PSHEAR,   &
                                         RPSHEAR, PSHEL, RPSHEL, PSOLID, PUSER1, RPUSER1, PUSERIN, PLATEOFF, PLATETHICK,           &
                                         USERIN_MAT_NAMES, VVEC
      USE ELSAVE_USE_IFs

      IMPLICIT NONE

      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME = 'ELSAVE'
      CHARACTER(LEN=DATA_NAM_LEN)     :: DATA_SET_NAME     ! A data set name for output purposes

      INTEGER(LONG)                   :: I,J               ! DO loop indices
      INTEGER(LONG)                   :: PCOMP_PLIES       ! Number of plies in 1 PCOMP entry incl sym plies not explicitly defined




! **********************************************************************************************************************************
! Write element data

      DATA_SET_NAME = 'ETYPE, EPNT, ESORT1, ESORT2, EOFF'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NELE
      DO I = 1,NELE
         WRITE(L1G) ETYPE(I)
         WRITE(L1G) EPNT(I)
         WRITE(L1G) ESORT1(I)
         WRITE(L1G) ESORT2(I)
         WRITE(L1G) EOFF(I)
      ENDDO

      DATA_SET_NAME = 'EDAT'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NEDAT
      DO I = 1,NEDAT
         WRITE(L1G) EDAT(I)
      ENDDO

! Write element parameters

      DATA_SET_NAME = 'ELEM PARAMETERS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) IORQ1M
      WRITE(L1G) IORQ1S
      WRITE(L1G) IORQ1B
      WRITE(L1G) IORQ2B
      WRITE(L1G) IORQ2T
      WRITE(L1G) CBMIN3
      WRITE(L1G) CBMIN4

! Write BAR, BEAM v vectors

      DATA_SET_NAME = 'V VECTORS IN GLOBAL COORDS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NVVEC
      DO I=1,NVVEC
         WRITE(L1G) (VVEC(I,J),J=1,3)
      ENDDO

! Write BAR, BEAM offsets

      DATA_SET_NAME = 'BAR, BEAM OFFSETS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NBAROFF
      DO I = 1,NBAROFF
         DO J = 1,6
            WRITE(L1G) BAROFF(I,J)
         ENDDO
      ENDDO

! Write BUSH offsets

      DATA_SET_NAME = 'BUSH OFFSETS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NBUSHOFF
      DO I = 1,NBUSHOFF
         DO J = 1,6
            WRITE(L1G) BUSHOFF(I,J)
         ENDDO
      ENDDO

! Write plate offsets

      DATA_SET_NAME = 'PLATE OFFSETS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPLATEOFF
      DO I = 1,NPLATEOFF
         WRITE(L1G) PLATEOFF(I)
      ENDDO

! Write plate thicknesses

      DATA_SET_NAME = 'PLATE THICKNESSES FROM CONNECTION ENTRIES'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPLATETHICK
      DO I = 1,NPLATETHICK
         WRITE(L1G) PLATETHICK(I)
      ENDDO

! Write property data

      DATA_SET_NAME = 'PBAR, RPBAR'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPBAR
      WRITE(L1G) MPBAR
      WRITE(L1G) MRPBAR
      DO I = 1,NPBAR
         DO J=1,MPBAR
            WRITE(L1G) PBAR(I,J)
         ENDDO
         DO J = 1,MRPBAR
            WRITE(L1G) RPBAR(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PBEAM, RPBEAM'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPBEAM
      WRITE(L1G) MPBEAM
      WRITE(L1G) MRPBEAM
      DO I = 1,NPBEAM
         DO J=1,MPBEAM
            WRITE(L1G) PBEAM(I,J)
         ENDDO
         WRITE(L1G) PBEAM_NSTATIONS(I)
         DO J=1,MPBEAM_STATIONS
            WRITE(L1G) PBEAM_XL(I,J)
         ENDDO
         DO J=1,MPBEAM_STATIONS
            WRITE(L1G) PBEAM_RPROPS(I,J,1)
            WRITE(L1G) PBEAM_RPROPS(I,J,2)
            WRITE(L1G) PBEAM_RPROPS(I,J,3)
            WRITE(L1G) PBEAM_RPROPS(I,J,4)
            WRITE(L1G) PBEAM_RPROPS(I,J,5)
            WRITE(L1G) PBEAM_RPROPS(I,J,6)
         ENDDO
         DO J = 1,MRPBEAM
            WRITE(L1G) RPBEAM(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PBUSH, RPBUSH'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPBUSH
      WRITE(L1G) MPBUSH
      WRITE(L1G) MRPBUSH
      DO I = 1,NPBUSH
         DO J=1,MPBUSH
            WRITE(L1G) PBUSH(I,J)
         ENDDO
         DO J = 1,MRPBUSH
            WRITE(L1G) RPBUSH(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PROD, RPROD'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPROD
      WRITE(L1G) MPROD
      WRITE(L1G) MRPROD
      DO I = 1,NPROD
         DO J=1,MPROD
            WRITE(L1G) PROD(I,J)
         ENDDO
         DO J = 1,MRPROD
            WRITE(L1G) RPROD(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PELAS, RPELAS'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPELAS
      WRITE(L1G) MPELAS
      WRITE(L1G) MRPELAS
      DO I = 1,NPELAS
         DO J=1,MPELAS
            WRITE(L1G) PELAS(I,J)
         ENDDO
         DO J = 1,MRPELAS
            WRITE(L1G) RPELAS(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PSHEAR, RPSHEAR'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPSHEAR
      WRITE(L1G) MPSHEAR
      WRITE(L1G) MRPSHEAR
      DO I = 1,NPSHEAR
         DO J = 1,MPSHEAR
            WRITE(L1G) PSHEAR(I,J)
         ENDDO
         DO J = 1,MRPSHEAR
            WRITE(L1G) RPSHEAR(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PSHEL, RPSHEL'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPSHEL
      WRITE(L1G) MPSHEL
      WRITE(L1G) MRPSHEL
      DO I = 1,NPSHEL
         DO J = 1,MPSHEL
            WRITE(L1G) PSHEL(I,J)
         ENDDO
         DO J = 1,MRPSHEL
            WRITE(L1G) RPSHEL(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PCOMP, RPCOMP'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPCOMP
      WRITE(L1G) MPCOMP0
      WRITE(L1G) MPCOMP_PLIES
      WRITE(L1G) MRPCOMP0
      WRITE(L1G) MRPCOMP_PLIES
      DO I = 1,NPCOMP
         PCOMP_PLIES = PCOMP(I,5)
         WRITE(L1G) PCOMP_PLIES
         DO J = 1,MPCOMP0+MPCOMP_PLIES*PCOMP_PLIES
            WRITE(L1G) PCOMP(I,J)
         ENDDO
         DO J = 1,MRPCOMP0+MRPCOMP_PLIES*PCOMP_PLIES
            WRITE(L1G) RPCOMP(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PSOLID'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPSOLID
      WRITE(L1G) MPSOLID
      DO I = 1,NPSOLID
         DO J = 1,MPSOLID
            WRITE(L1G) PSOLID(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PUSER1, RPUSER1'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPUSER1
      WRITE(L1G) MPUSER1
      WRITE(L1G) MRPUSER1
      DO I = 1,NPUSER1
         DO J=1,MPUSER1
            WRITE(L1G) PUSER1(I,J)
         ENDDO
         DO J = 1,MRPUSER1
            WRITE(L1G) RPUSER1(I,J)
         ENDDO
      ENDDO

      DATA_SET_NAME = 'PUSERIN'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPUSERIN
      WRITE(L1G) MPUSERIN
      DO I = 1,NPUSERIN
         DO J=1,MPUSERIN
            WRITE(L1G) PUSERIN(I,J)
         ENDDO
      ENDDO

! Write USERIN_MAT_NAMES

      DATA_SET_NAME = 'USERIN_MAT_NAMES'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NPUSERIN                                  ! NOTE: there are as many USERIN_MAT_NAMES as PUSERIN entries
      WRITE(L1G) MUSERIN_MAT_NAMES
      DO I = 1,NPUSERIN
         DO J=1,MUSERIN_MAT_NAMES
            WRITE(L1G) USERIN_MAT_NAMES(I,J)
         ENDDO
      ENDDO

! Write material data

      DATA_SET_NAME = 'MATL, RMATL'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NMATL
      WRITE(L1G) MMATL
      WRITE(L1G) MRMATLC
      DO I = 1,NMATL
         DO J=1,MMATL
            WRITE(L1G) MATL(I,J)
         ENDDO
         DO J = 1,MRMATLC
            WRITE(L1G) RMATL(I,J)
         ENDDO
      ENDDO

! Write material property angles

      DATA_SET_NAME = 'MATERIAL PROPERTY ANGLES'
      WRITE(L1G) DATA_SET_NAME
      WRITE(L1G) NMATANGLE
      DO I = 1,NMATANGLE
         WRITE(L1G) MATANGLE(I)
      ENDDO



      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE ELSAVE
