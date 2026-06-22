! ###############################################################################################################################
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

   MODULE BD_PARAM_Interface

   INTERFACE

      SUBROUTINE BD_PARAM ( CARD )

  
      USE PENTIUM_II_KIND, ONLY       :  BYTE, LONG, DOUBLE
      USE IOUNT1, ONLY                :  WRT_ERR, ERR, F06

      USE SCONTR, ONLY                :  BLNK_SUB_NAM, ECHO, FATAL_ERR, IERRFL, JCARD_LEN, JF, MEPSIL, MPBARLU, NUM_USETSTR,       &
                                         WARN_ERR
      USE TIMDAT, ONLY                :  TSEC
      USE CONSTANTS_1, ONLY           :  ZERO, ONE
      USE MACHINE_PARAMS, ONLY        :  MACH_PREC
      USE DOF_TABLES, ONLY            :  TSET_CHR_LEN

      USE PARAMS, ONLY                :  ARP_TOL         , ART_KED         , ART_ROT_KED     , ART_TRAN_KED    ,                   &
                                         ART_MASS        , ART_ROT_MASS    , ART_TRAN_MASS   , AUTOSPC         , AUTOSPC_NSET    , &
                                         AUTOSPC_RAT     , AUTOSPC_INFO    , AUTOSPC_SPCF    , BAILOUT         , CRS_CCS         , &
!---  cbeam add --- begin!
! add cbeam parameter stiffnes modifier
                                         CBMIN3          , CBMIN4          , CBMIN4T         , CBEAMAREA       , CBEAMAREA_PID   , &
                                         CBEAMAREA_VAL   , CBEAMSHR        , CBEAMSHR_PID    , CBEAMSHR_VAL    , &
                                         MCBEAMAREA_PID  , MCBEAMSHR_PID   , NCBEAMAREA_PID  , NCBEAMSHR_PID   , CHKGRDS         , &
!---  cbeam add --- end!
                                         CUSERIN         , CUSERIN_EID     , CUSERIN_IN4     , CUSERIN_PID     , CUSERIN_SPNT_ID , &
                                         CUSERIN_XSET    , CUSERIN_COMPTYP , DARPACK         ,                                     &
                                         DELBAN          , EIGESTL         , EIGNORM2        , ELFORCEN        , EPSERR          , &
                                         EQCHK_REF_GRID  , EQCHK_NORM      , EQCHK_OUTPUT    , EQCHK_TINY      ,                   &
                                         EPSIL           , EMP0_PAUSE      , ESP0_PAUSE      , F06_COL_START   ,                   &
                                         GRDPNT          , GRIDSEQ         , HEXAXIS         ,                                     &
                                         IORQ1M          , IORQ1S          , IORQ1B          , IORQ2B          , IORQ2T          , &
! add lancmeth for next eigen solver choice : DENSE, SUBSPACE, FEAST
                                         ITMAX           , KLLRAT          , KOORAT          , LANCMETH        , MATSPARS        , &
                                         MEMAFAC         , MIN4TRED        , MXALLOCA        , MAXRATIO        ,                   &
                                         MEFMCORD        , MEFMLOC         , MEFMGRID        ,                                     &
                                         MPFOUT          , MXITERI         , MXITERL         , OTMSKIP         , POST            , &
                                         PBARLDEC        , PBARLSHR        , PCOMPEQ         , PCHSPC1         , PCMPTSTM        , &
                                         PRTBASIC        , PRTCONN         , PRTCORD         , PRTDISP         , PRTDLR          , &
                                         PRTDOF          , PRTFOR          , PRTHMN          , PRTGMN          , PRTGOA          , &
                                         PRTCGLTM        , PRTPHIZL        , PRTIFLTM        , PRTKXX          , PRTMXX          , &
                                         PRTOU4          , PRTPHIXA        , PRTMASS         , PRTMASSD        , PRTRMG          , &
                                         PRTSCP          , PRTPSET         , PRTTSET         , PRTUSET         ,                   &
                                         PRTSTIFD        , PRTSTIFF        , PRTUO0          ,                                     &
                                         PRTYS           , PRTQSYS         ,                                                       &
                                         Q4SURFIT        , QUADAXIS        , QUAD4TYP        , RCONDK          , RELINK3         , &
                                         SEQPRT          , SEQQUIT         , SETLKTM         , SETLKTK         , SHRFXFAC        , &
                                         SKIPMGG         , SOLLIB          , SPARSE_FLAVOR   , SPARSTOR        ,                   &
                                         SPC1QUIT        , SORT_MAX        , SPC1SID         , STR_CID                           , &
                                         SUPINFO         , SUPWARN                                                               , &
                                         THRESHK         , THRESHK_LAP     , TINY            ,                                     &
                                         TSTM_DEF        , USR_JCT         , USR_LTERM_KGG   , USR_LTERM_MGG   , WINAMEM         , &
                                         WTMASS
 
      IMPLICIT NONE
 
      CHARACTER(LEN=*), INTENT(IN)    :: CARD              ! A Bulk Data card
      CHARACTER(15*BYTE)              :: PARNAM            ! The name of a parameter
 
      INTEGER(LONG)                   :: LOWER             ! Lower allowable value for an integer parameter
      INTEGER(LONG)                   :: UPPER             ! Upper allowable value for an integer parameter
      END SUBROUTINE BD_PARAM

   END INTERFACE

   END MODULE BD_PARAM_Interface


