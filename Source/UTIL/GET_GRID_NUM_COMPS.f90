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


      SUBROUTINE GET_GRID_NUM_COMPS ( IGRID, NUM_COMPS, CALLING_SUBR )

! Gets the number of components for a "grid" number from array GRID by testing GRID(IGRID,6)
! If GRID(I,6) is 6 then this is a physical grid with 6 comps of displ and if GRID(I,6) is 1 then this is an SPOINT.

      USE PENTIUM_II_KIND, ONLY       :  LONG
      USE MODEL_STUF, ONLY            :  GRID

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN)    :: CALLING_SUBR      ! Subr that called this one

      INTEGER(LONG), INTENT(IN)       :: IGRID             ! An internal grid number
      INTEGER(LONG), INTENT(OUT)      :: NUM_COMPS         ! 6 if GRID_NUM is an physical grid, 1 if an SPOINT

      INTEGER(LONG)                   :: I
      INTEGER(LONG)                   :: NGRID

      NGRID = SIZE(GRID,1)

! **********************************************************************************************************************************

      ! Most callers pass an internal grid row number, but some report/output paths
      ! can hand us the external grid ID. Accept either form to avoid a bounds crash.
      IF ((IGRID >= 1) .AND. (IGRID <= NGRID)) THEN
         NUM_COMPS = GRID(IGRID, 6)
      ELSE
         NUM_COMPS = 6
         DO I = 1, NGRID
            IF (GRID(I,1) == IGRID) THEN
               NUM_COMPS = GRID(I, 6)
               EXIT
            ENDIF
         ENDDO
      ENDIF

      RETURN

! **********************************************************************************************************************************

      END SUBROUTINE GET_GRID_NUM_COMPS

