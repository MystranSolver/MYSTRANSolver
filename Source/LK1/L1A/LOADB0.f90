!--- cbeam add --- begin!
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

      SUBROUTINE LOADB0

      ! Preliminary reading of the Bulk Data to count several data sizes so
      ! that arrays may be allocated prior to the final reading of the Bulk Data.
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06, IN1                
      USE SCONTR, ONLY                :  BD_ENTRY_LEN, BLNK_SUB_NAM, FATAL_ERR, LCMASS, LDOFG, LELE,                               &
                                         LEDAT, LFORCE, LCONM2, LCORD, LGRAV, LGRID, LGUSERIN, LLOADC, LLOADR,                     &
                                         LMATL, LMPC, LMPCADDC, LMPCADDR, LPBAR, LPBEAM, LPBUSH, LPCOMP, LPCOMP_PLIES, LPDAT,      &
                                         LPELAS, LPLOAD, LPMASS, LPROD, LPSHEL, LPSHEAR, LPSOLID, LPUSER1, LPUSERIN, LRFORCE,      &
                                         LRIGEL, LSEQ, LSLOAD, LSNORM, LSPC, LSPC1, LSPCADDC, LSPCADDR, LSUSERIN, LTDAT,           &
                                         MEDAT_CBAR, MEDAT_CBEAM, MEDAT_CBUSH,                                                     &
                                         MEDAT_CELAS1, MEDAT_CELAS2, MEDAT_CELAS3, MEDAT_CELAS4,                                   &
                                         MEDAT_CQUAD, MEDAT_CQUAD8, MEDAT_CROD, MEDAT_CSHEAR, MEDAT_CTRIA, MEDAT_CUSER1,           &
                                         MEDAT0_CUSERIN, MMPC,                                                                     &
                                         MPDAT_PLOAD2, MPDAT_PLOAD4, MEDAT_PLOTEL, MRBE3, MRSPLINE, MTDAT_TEMPRB, MTDAT_TEMPP1,    &
! add pload1
                                         MPDAT_PLOAD1,                                                                             &
                                         NPBARL, NSPOINT, PROG_NAME
      USE TIMDAT, ONLY                :  TSEC

      USE MODEL_STUF, ONLY            :  GRDSET3, GRDSET7, GRDSET8
      USE PARAMS, ONLY                :  GRIDSEQ
 
      USE LOADB0_USE_IFs

      IMPLICIT NONE
 
      CHARACTER(LEN=LEN(BLNK_SUB_NAM)):: SUBR_NAME   = 'LOADB0'
      CHARACTER( 7*BYTE), PARAMETER   :: END_CARD    = 'ENDDATA'

      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD              ! 16 col field card (either CARD1 if small field or CARD1 + CARD2 if
!                                                            a large field card. This is output from subr FFIELD

      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD1             ! BD card (a small field card or the 1st half of a large field card)
      CHARACTER(LEN=BD_ENTRY_LEN)     :: CARD2             ! 2nd half of a large field card
      CHARACTER( 9*BYTE)              :: DECK_NAME   = 'BULK DATA'
      CHARACTER( 1*BYTE)              :: LARGE_FLD_INP     ! If 'Y', card is in large field format

      INTEGER(LONG)                   :: COMMENT_COL       ! Col on CARD where a comment begins (if one exists)
      INTEGER(LONG)                   :: DELTA_LEDAT       ! Delta number of words to add to LEDAT for an element (for elements
!                                                            where the delta cannot be discerned from the BD card name, such as
!                                                            CHEXA which can have 8 or 20 nodes
      INTEGER(LONG)                   :: DELTA_SLOAD       ! Delta number of words to add to LSLOAD for SLOAD's
      INTEGER(LONG)                   :: DELTA_SPOINT      ! Delta number of words to add to LGRID for SPOINT's
      INTEGER(LONG)                   :: I                 ! DO loop index
      INTEGER(LONG)                   :: ICNT              ! Card count
      INTEGER(LONG)                   :: IERR              ! Error indicatro from subr FFIELD
      INTEGER(LONG)                   :: IOCHK             ! IOSTAT error number when reading a file
      INTEGER(LONG)                   :: ILOAD             ! Number of loads defined on 1 logical B.D. LOAD card
      INTEGER(LONG)                   :: IMPC              ! Number of grid/comp/coeff triplets defined on 1 logical B.D. MPC  card
      INTEGER(LONG)                   :: IMPCADD           ! Number of MPC set ID's defined on 1 B.D. MPCADD card
      INTEGER(LONG)                   :: IRBE3             ! Number of grid/comp/coeff triplets defined on 1 logical B.D. RBE3 card
      INTEGER(LONG)                   :: IRSPLINE          ! Number of fields that can have a grid or comp num on an RSPLINE entry
      INTEGER(LONG)                   :: ISPCADD           ! Number of SPC set ID's defined on 1 B.D. SPCADD card
      INTEGER(LONG)                   :: IPLIES            ! Number of composite layers on 1 B.D. PCOMP card
      INTEGER(LONG)                   :: NG_USERIN         ! Number of grids found on USERIN elems (not incl SPOINT's)
      INTEGER(LONG)                   :: NS_USERIN         ! Number of SPOINT's found on USERIN elems

 

! **********************************************************************************************************************************
      CARD(1:) = ' '

! Set MRBE3, MRSPLINE so that we will make sure that they will be at least 1

      MRBE3    = 1
      MRSPLINE = 1

! Process Bulk Data cards in a loop that runs until either an ENDDATA card is found or when an error or EOF/EOR occurs
    
      ICNT = 0
!xx   REWIND (IN1)
      DO

         ICNT = ICNT + 1
         CALL READ_BDF_LINE(IN1, IOCHK, CARD1)
         CARD(1:)  = CARD1(1:)

         ! Quit if EOF/EOR occurs.
         IF (IOCHK < 0) THEN
            WRITE(ERR,1011) END_CARD
            WRITE(F06,1011) END_CARD
            FATAL_ERR = FATAL_ERR + 1
            CALL OUTA_HERE ( 'Y' )
         ENDIF

         ! Check if error occurs.
         IF (IOCHK > 0) THEN
            WRITE(ERR,1010) DECK_NAME
            WRITE(F06,1010) DECK_NAME
            WRITE(F06,'(A)') CARD1
            FATAL_ERR = FATAL_ERR + 1
            CYCLE
         ENDIF

         ! Remove any comments within the CARD1 by deleting everything
         ! from $ on (after col 1)
         COMMENT_COL = 1
         DO I=2,BD_ENTRY_LEN
            IF (CARD1(I:I) == '$') THEN
               COMMENT_COL = I
               EXIT
            ENDIF
         ENDDO

         IF (COMMENT_COL > 1) THEN
            CARD1(COMMENT_COL:) = ' '
         ENDIF

         ! Determine if the card is large or small format
         LARGE_FLD_INP = 'N'
         DO I=1,8
            IF (CARD1(I:I) == '*') THEN
               LARGE_FLD_INP = 'Y'
            ENDIF
         ENDDO

         ! FFIELD converts free-field card to fixed field and left justifies
         ! data in fields 2-9 and outputs a 10 field, 16 col/field CARD1
          IF ((CARD1(1:1) /= '$')  .AND. (CARD1(1:) /= ' ')) THEN

            IF (LARGE_FLD_INP == 'N') THEN

               CALL FFIELD ( CARD1, IERR )
               CARD(1:) = CARD1(1:)

            ELSE

               ! Read 2nd physical entry for a large field parent B.D. entry
               CALL READ_BDF_LINE(IN1, IOCHK, CARD2)

               IF (IOCHK < 0) THEN
                  WRITE(ERR,1011) END_CARD
                  WRITE(F06,1011) END_CARD
                  FATAL_ERR = FATAL_ERR + 1
                  CALL OUTA_HERE ( 'Y' )
               ENDIF
 
               IF (IOCHK > 0) THEN
                  WRITE(ERR,1010) DECK_NAME
                  WRITE(F06,1010) DECK_NAME
                  WRITE(F06,'(A)') CARD2
                  FATAL_ERR = FATAL_ERR + 1
                  CYCLE
               ENDIF
 
               COMMENT_COL = 1
               DO I=2,BD_ENTRY_LEN
                  IF (CARD2(I:I) == '$') THEN
                     COMMENT_COL = I
                     EXIT
                  ENDIF
               ENDDO

               IF (COMMENT_COL > 1) THEN
                  CARD2(COMMENT_COL:) = ' '
               ENDIF

!              IF (CARD2(1:8) /= CARD1(73:80)) THEN
!                 BACKSPACE(IN1)
!                 CARD2(1:) = ' '
!                 CARD2(1:8) = CARD1(73:80)
!              ENDIF

               IF      (CARD2( 1: 8) == CARD1(73:80)) THEN
                  CONTINUE     !ICONT = 1
               ELSE IF ((CARD2( 1: 8) == '*       ') .AND. (CARD1(73:80) == '        ')) THEN
                  CONTINUE     !ICONT = 1
               ELSE IF ((CARD2( 1: 8) == '        ') .AND. (CARD1(73:80) == '*       ')) THEN
                  CONTINUE     !ICONT = 1
               ELSE
                  BACKSPACE(IN1)
                  CARD2(1:) = ' '
                  CARD2(1:8) = CARD1(73:80)
               ENDIF

               CALL FFIELD2 ( CARD1, CARD2, CARD, IERR )

            ENDIF

            IF (IERR /= 0) THEN
               CYCLE
            ENDIF 

         ENDIF
 
          ! No errors, so process Bulk Data card. No need to check for
          ! imbedded blanks found when FFIELD was run - this will be
          ! checked when LOADB reads the bulk data
          IF      (CARD(1:5) == 'BAROR'   )  THEN
            CALL BD_BAROR0 ( CARD )

         ELSE IF (CARD(1:6) == 'BEAMOR'  )  THEN
            CALL BD_BEAMOR0( CARD )

         ELSE IF ((CARD(1:4) == 'CBAR'    ) .OR. (CARD(1:5) == 'CBEAM'   ))  THEN
            LELE  = LELE + 1
            IF (CARD(1:4) == 'CBAR'    ) THEN
               LEDAT = LEDAT + MEDAT_CBAR
            ELSE
               LEDAT = LEDAT + MEDAT_CBEAM
            ENDIF
            CALL BD_CBAR0 ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:5) == 'CBUSH'   )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CBUSH
            CALL BD_CBUSH0 ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'CELAS1'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CELAS1

         ELSE IF (CARD(1:6) == 'CELAS2'  )  THEN
            LELE   = LELE + 1
            LEDAT  = LEDAT + MEDAT_CELAS2
            LPELAS = LPELAS + 1                            ! CELAS2 has props on conn entry and we create a PELAS for them

         ELSE IF (CARD(1:6) == 'CELAS3'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CELAS3

         ELSE IF (CARD(1:6) == 'CELAS4'  )  THEN
            LELE   = LELE + 1
            LEDAT  = LEDAT + MEDAT_CELAS4
            LPELAS = LPELAS + 1                            ! CELAS4 has props on conn entry and we create a PELAS for them

         ELSE IF (CARD(1:5) == 'CHEXA'   ) THEN
            LELE  = LELE + 1
            CALL BD_CHEXA0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )
            LEDAT = LEDAT + DELTA_LEDAT

         ELSE IF (CARD(1:6) == 'CMASS1'  )  THEN
            LCMASS = LCMASS + 1
  
         ELSE IF (CARD(1:6) == 'CMASS2'  )  THEN
            LCMASS = LCMASS + 1
            LPMASS = LPMASS + 1
  
         ELSE IF (CARD(1:6) == 'CMASS3'  )  THEN
            LCMASS = LCMASS + 1
  
         ELSE IF (CARD(1:6) == 'CMASS4'  )  THEN
            LCMASS = LCMASS + 1
            LPMASS = LPMASS + 1
  
         ELSE IF (CARD(1:5) == 'CONM2'   )  THEN
            LCONM2 = LCONM2 + 1
  
         ELSE IF (CARD(1:6) == 'CONROD'  )  THEN
            LELE  = LELE  + 1
            LPROD = LPROD + 1
            LEDAT = LEDAT + MEDAT_CROD
  
         ELSE IF((CARD(1:6) == 'CORD1C'  ) .OR. (CARD(1:6) == 'CORD1R'  ) .OR. (CARD(1:6) == 'CORD1S'  )) THEN
            LCORD = LCORD + 1
            IF (CARD(41:48) /= '        ') LCORD = LCORD + 1
  
         ELSE IF((CARD(1:6) == 'CORD2C'  ) .OR. (CARD(1:6) == 'CORD2R'  ) .OR. (CARD(1:6) == 'CORD2S'  )) THEN
            LCORD = LCORD + 1
  
         ELSE IF (CARD(1:6) == 'CPENTA'  ) THEN
            LELE  = LELE + 1
            CALL BD_CPENTA0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )
            LEDAT = LEDAT + DELTA_LEDAT

         ELSE IF (CARD(1:6) == 'CQUAD4'  ) THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CQUAD
            CALL BD_CQUAD0 ( CARD, LARGE_FLD_INP )

         ELSE IF (CARD(1:6) == 'CQUAD8'  ) THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CQUAD8
            CALL BD_CQUAD80 ( CARD, LARGE_FLD_INP )
   
         ELSE IF (CARD(1:4) == 'CROD'    )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CROD
  
         ELSE IF (CARD(1:6) == 'CSHEAR'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CSHEAR
  
         ELSE IF (CARD(1:6) == 'CTETRA'  ) THEN
            LELE  = LELE + 1
            CALL BD_CTETRA0 ( CARD, LARGE_FLD_INP, DELTA_LEDAT )
            LEDAT = LEDAT + DELTA_LEDAT

         ELSE IF (CARD(1:6) == 'CTRIA3'  ) THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CTRIA
            CALL BD_CTRIA0 ( CARD, LARGE_FLD_INP )
   
         ELSE IF (CARD(1:6) == 'CUSER1'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_CUSER1

         ELSE IF (CARD(1:7) == 'CUSERIN' )  THEN
            LELE  = LELE + 1
            CALL BD_CUSERIN0 ( CARD, NG_USERIN, NS_USERIN )
            IF (NG_USERIN > LGUSERIN) THEN
               LGUSERIN = NG_USERIN
            ENDIF
            IF (NS_USERIN > LSUSERIN) THEN
               LSUSERIN = NS_USERIN
            ENDIF
                                                           ! LEDAT has "+ 1" term since last record is NUM_BDY_DOF not in MEDAT0
            LEDAT = LEDAT + MEDAT0_CUSERIN + 2*NG_USERIN + NS_USERIN +  1

         ELSE IF (CARD(1:5) == 'DEBUG'   )  THEN
            CALL BD_DEBUG0 ( CARD )
 
         ELSE IF((CARD(1:5) == 'FORCE'   ) .OR. (CARD(1:6) == 'MOMENT'  )) THEN
            LFORCE = LFORCE +1
 
         ELSE IF (CARD(1:4) == 'GRAV'    )  THEN
            LGRAV = LGRAV + 1
 
         ELSE IF (CARD(1:6) == 'GRDSET'  )  THEN
            CALL BD_GRDSET0 ( CARD )
 
         ELSE IF (CARD(1:4) == 'GRID'    )  THEN
            LGRID = LGRID + 1
            LDOFG = LDOFG + 6
 
         ELSE IF (CARD(1:4) == 'LOAD'    )  THEN
            LLOADR = LLOADR + 1
            CALL BD_LOAD0 ( CARD, LARGE_FLD_INP, ILOAD )
            IF (ILOAD > LLOADC) THEN
               LLOADC = ILOAD
            ENDIF

         ELSE IF((CARD(1:4) == 'MAT1'    ) .OR. (CARD(1:4) == 'MAT2'    )  .OR.                                                    &
                 (CARD(1:4) == 'MAT8'    ) .OR. (CARD(1:4) == 'MAT9'    )) THEN
            LMATL = LMATL + 1

         ELSE IF((CARD(1:4) == 'MPC '    ) .OR. (CARD(1:4) == 'MPC*'    )) THEN
            LMPC = LMPC + 1
            CALL BD_MPC0 ( CARD, LARGE_FLD_INP, IMPC )
            IF (IMPC > MMPC) THEN
               MMPC = IMPC
            ENDIF

         ELSE IF (CARD(1:6) == 'MPCADD'  )  THEN
            LMPCADDR = LMPCADDR + 1
            CALL BD_MPCADD0 ( CARD, LARGE_FLD_INP, IMPCADD )
            IF (IMPCADD > LMPCADDC) THEN
               LMPCADDC = IMPCADD
            ENDIF

         ELSE IF (CARD(1:6) == 'PARAM '  )  THEN
            CALL BD_PARAM0 ( CARD )
 
         ELSE IF((CARD(1:5) == 'PBAR '   ) .OR. (CARD(1:5) == 'PBAR*'   ))  THEN
            LPBAR = LPBAR + 1
 
         ELSE IF (CARD(1:5) == 'PBARL'   )  THEN
            LPBAR  = LPBAR  + 1
            NPBARL = NPBARL + 1
 ! --- cbeam_add begin --- !
         ELSE IF (CARD(1:6) == 'PBEAML'  )  THEN
            LPBEAM = LPBEAM + 1
 ! --- cbeam_add end --- !
         ELSE IF (CARD(1:5) == 'PBEAM'   )  THEN
            LPBEAM = LPBEAM + 1
 
         ELSE IF (CARD(1:5) == 'PBUSH'   )  THEN
            LPBUSH = LPBUSH + 1
 
         ELSE IF (CARD(1:5) == 'PCOMP'   )  THEN
            LPCOMP = LPCOMP + 1
            CALL BD_PCOMP0 ( CARD, LARGE_FLD_INP, IPLIES )
            IF (IPLIES > LPCOMP_PLIES) THEN
               LPCOMP_PLIES = IPLIES
            ENDIF

         ELSE IF (CARD(1:6) == 'PCOMP1'  )  THEN
            LPCOMP = LPCOMP + 1
            CALL BD_PCOMP10 ( CARD, LARGE_FLD_INP, IPLIES )
            IF (IPLIES > LPCOMP_PLIES) THEN
               LPCOMP_PLIES = IPLIES
            ENDIF

         ELSE IF (CARD(1:5) == 'PELAS'   )  THEN
            LPELAS = LPELAS + 1
 ! --- cbeam_add begin --- !
         ELSE IF ((CARD(1:6) == 'PLOAD1'  ) .OR. (CARD(1:6) == 'PLOAD2'  ))  THEN
            IF (CARD(1:6) == 'PLOAD1') THEN
               LPDAT  = LPDAT  + MPDAT_PLOAD1               ! Beam/bar PLOAD1 data slots
            ELSE
               LPDAT  = LPDAT  + MPDAT_PLOAD2
            ENDIF
            LPLOAD = LPLOAD + 1
 ! --- cbeam_add end --- !

         ELSE IF (CARD(1:6) == 'PLOAD4'  )  THEN
            LPDAT  = LPDAT  + MPDAT_PLOAD4
            LPLOAD = LPLOAD + 1
 
         ELSE IF (CARD(1:6) == 'PLOTEL'  )  THEN
            LELE  = LELE + 1
            LEDAT = LEDAT + MEDAT_PLOTEL
 
         ELSE IF (CARD(1:5) == 'PMASS'   )  THEN
            LPMASS = LPMASS + 4

         ELSE IF (CARD(1:4) == 'PROD'    )  THEN
            LPROD = LPROD + 1

         ELSE IF (CARD(1:6) == 'PSHEAR'  )  THEN
            LPSHEAR = LPSHEAR + 1
  
         ELSE IF (CARD(1:6) == 'PSHELL'  )  THEN
            LPSHEL = LPSHEL + 1
  
         ELSE IF (CARD(1:6) == 'PSOLID'  )  THEN
            LPSOLID = LPSOLID + 1
  
         ELSE IF (CARD(1:6) == 'PUSER1'  )  THEN
            LPUSER1 = LPUSER1 + 1
  
         ELSE IF (CARD(1:7) == 'PUSERIN' )  THEN
            LPUSERIN = LPUSERIN + 1
  
         ELSE IF (CARD(1:4) == 'RBAR'    ) THEN
            LRIGEL  = LRIGEL + 1

         ELSE IF (CARD(1:4) == 'RBE1'    ) THEN
            LRIGEL  = LRIGEL + 1

         ELSE IF (CARD(1:4) == 'RBE2'    ) THEN
            LRIGEL  = LRIGEL + 1

         ELSE IF (CARD(1:4) == 'RBE3'    )  THEN
            LRIGEL  = LRIGEL + 1
            CALL BD_RBE30 ( CARD, LARGE_FLD_INP, IRBE3 )
            IF (IRBE3 > MRBE3) THEN
               MRBE3 = IRBE3
            ENDIF

         ELSE IF (CARD(1:6) == 'RFORCE'  )  THEN
            LRFORCE = LRFORCE + 1
 
         ELSE IF (CARD(1:7) == 'RSPLINE' )  THEN
            CALL BD_RSPLINE0 ( CARD, LARGE_FLD_INP, IRSPLINE )
            LRIGEL = LRIGEL + 1
            IF (IRSPLINE > MRSPLINE) THEN
               MRSPLINE = IRSPLINE
            ENDIF
 
         ELSE IF (CARD(1:5) == 'SLOAD'   )  THEN
            CALL BD_SLOAD0 ( CARD, DELTA_SLOAD )
            LSLOAD = LSLOAD + DELTA_SLOAD
 
         ELSE IF (CARD(1:5) == 'SEQGP'   )  THEN
            LSEQ = LSEQ + 4                                ! Conservative estimate. There can only be 4 entries per card

         ELSE IF (CARD(1:5) == 'SNORM'   )  THEN
            LSNORM = LSNORM + 1
 
         ELSE IF((CARD(1:4) == 'SPC '    ) .OR. (CARD(1:4) == 'SPC*'    )) THEN
            LSPC = LSPC + 1
 ! --- cbeam_add begin --- !
         ELSE IF (CARD(1:4) == 'SPCD'    )  THEN
            LSPC = LSPC + 1
 ! --- cbeam_add end --- !
         ELSE IF (CARD(1:4) == 'SPC1'    )  THEN
            LSPC1 = LSPC1 + 1

         ELSE IF (CARD(1:6) == 'SPCADD'  )  THEN
            LSPCADDR = LSPCADDR + 1
            CALL BD_SPCADD0 ( CARD, LARGE_FLD_INP, ISPCADD )
            IF (ISPCADD > LSPCADDC) THEN
               LSPCADDC = ISPCADD
            ENDIF

         ELSE IF (CARD(1:6) == 'SPOINT'  )  THEN
            CALL BD_SPOINT0 ( CARD, DELTA_SPOINT )
            NSPOINT = NSPOINT + DELTA_SPOINT               ! DELTA_SPOINT = number of SPOINTS defined on this SPOINT Bulk Data entry
            LGRID   = LGRID   + DELTA_SPOINT               ! Each SPOINT counts as 1 in the number of grids
            LDOFG   = LDOFG   + DELTA_SPOINT

         ELSE IF (CARD(1:6) == 'TEMPRB'  )  THEN
            LTDAT = LTDAT + MTDAT_TEMPRB
 
         ELSE IF (CARD(1:6) == 'TEMPP1'  )  THEN
            LTDAT = LTDAT + MTDAT_TEMPP1

         ELSE IF (CARD(1:7) == 'ENDDATA' )  THEN
            EXIT          
 
         ENDIF
 
      ENDDO

! ! Reset FATAL_ERR to 0. Errors that incremented FATAL_ERR in LOADBO (and routines it calls) will be recorded and 
! ! reported when LOADB runs.
! ! 
! !   FATAL_ERR = 0 

      RETURN

! **********************************************************************************************************************************
  101 FORMAT(A)

 1010 FORMAT(' *ERROR  1010: ERROR READING FOLLOWING ',A,' ENTRY. ENTRY IGNORED')

 1011 FORMAT(' *ERROR  1011: NO ',A10,' ENTRY FOUND BEFORE END OF FILE OR END OF RECORD IN INPUT FILE')

! **********************************************************************************************************************************
 
      END SUBROUTINE LOADB0

