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

      MODULE CC_OUTPUT_DESCRIBERS

      USE PENTIUM_II_KIND, ONLY         :  BYTE, LONG
      USE SCONTR, ONLY                  :  CC_CMD_DESCRIBERS

      IMPLICIT NONE

      SAVE

      ! The following are the default values for the Case Control command
      ! describers.  These are used to check the values in parens () in Case
      ! Control entries to see if they are valid for MYSTRAN.  The ones below
      ! are ones MYSTRAN honors.  If other values are encountered in the user's
      !  Case Control, warning messages are written to the MYSTRAN output file.

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: ACCE_SORT = 'SORT1   '
      LOGICAL                           :: ACCE_F06  = .TRUE.
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: ACCE_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: DISP_SORT = 'SORT1   '
      LOGICAL                           :: DISP_F06  = .TRUE.
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: DISP_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: GPFO_SORT = 'SORT1   '
      LOGICAL                           :: GPFO_F06  = .TRUE.

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: MPCF_SORT = 'SORT1   '
      LOGICAL                           :: MPCF_F06  = .TRUE.
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: MPCF_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: OLOA_SORT = 'SORT1   '
      LOGICAL                           :: OLOA_F06  = .TRUE.
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: OLOA_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: SPCF_SORT = 'SORT1   '
      LOGICAL                           :: SPCF_F06  = .TRUE.
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: SPCF_MAG  = 'MAG     '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: FORC_SORT = 'SORT1   '
      LOGICAL                           :: FORC_F06  = .TRUE.
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: FORC_MAG  = 'MAG     '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: FORC_LOC  = 'CENTER  '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_SORT = 'SORT1   '
      LOGICAL                           :: STRN_F06  = .TRUE.
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_MAG  = 'MAG     '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_OPT  = 'VONMISES'
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_CUR  = 'STRCUR  '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRN_LOC  = 'CENTER  '

      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRE_SORT = 'SORT1   '
      LOGICAL                           :: STRE_F06  = .TRUE.
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRE_MAG  = 'MAG     '
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRE_OPT  = 'VONMISES'
      CHARACTER(LEN(CC_CMD_DESCRIBERS)) :: STRE_LOC  = 'CENTER  '

      END MODULE CC_OUTPUT_DESCRIBERS
