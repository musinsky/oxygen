co************************************************************************
      SUBROUTINE BLTREE(NEL,RN,DEL,TPR,MV,BIMP,PDL)
************************************************************************
*     PURPOSE : TREE OF CASCADES (FOR STUDY THE CODE AND DEBUGGING)    *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     RN    - PARAMETER SPECIFYING MINIMUM DISTANCE BETWEEN CENTRES    *
*             OF NUCLEONS IN BOTH NUCLEI (RN=0.4[FM]).                 *
*     DEL   - RADIUS OF STRONG INTERACTION                             *
*     TPR   - PARAMETER SPECIFYING FORMATION TIME OF SECONDARIES       *
*     PDL   - PARAMETER SPECIFYING NUCLEAR TRANSPERANCY; NOT USED.     *
*     BIMP  - INPUT VALUE OF IMPACT PARAMETER; IF BIMP < 0 THEN        *
*             IMPACT PARAMETER IS CALCULATED RANDOMLY BY PROGRAM.      *     
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES:                                                 *
*     NEL      - =0 --> INTERACTION REALIZED                           *
*                =1 --> INTERACTION NOT REALIZED                       *
*     MV       - CURRENT NUMBER OF CASCADE PARTICLES IN ARRAYS         *
*                PM AND IM.                                            *
************************************************************************
***  COMMON BLOCKS  ***                                                *
*      RESCAS:                                                         *
*          AC    - CURRENT ATOMIC NUMBERS OF COLLIDING NUCLEI          *
*          ZC    - THEIR CURRENT CHARGES                               *
*          E     - THEIR CURRENT EXCITATION ENERGY                     *
*          PC    - THIER CURRENT MOMENTUM IN LAB                       *
*          AMC   - THEIR CURRENT ANGULAR MOMENTUN                      *
*      DATIN:                                                          *
*          A     - INITIAL ATOMIC NUMBERS OF COLLIDING NUCLEI          *
*          Z     - INITIAL CHARGES OF COLLIDING NUCLEI                 *
*          WT    - KINETIC ENERGY OF PROJECTILE                        *
*          VPI   - MESON POTENTIAL IN NUCLEUS                          *
*          EPS   - 
*          AR,CR - PARAMETERS FOR WOODS-SAXON DENSITY DISTRIBUTION     *
*          DR    - CUT PARAMETERS FOR DENSITY DISTRIBUTION             *
*          R     - NUCLEAR RADII                                       *
*          TR    - BOUND FERMI ENERGY                                  *
*      LNUC:                                                           *
*          BR    - PARAMETERS OF GAUSSIAN DENSITY DISTRIBUTION         *
*          PIMPAC - IMPACT PARAMETER                                   *
*      CARCIL:                                                         *
*          NK1P  - ARRAY COMPOSING INDECES OF TARGET NUCLEON-PARTNERS  *
*                  IN SUBROUTINE 'TIMEC'.                              *
*          TIME  - ARRAY COMPOSING TIMES NEEDED FOR PROJECTILE-NUCLEONS*
*                  TO PASS.                                            *
*          MK1P  - ARRAY COMPOSING INDICES OF NUCLEON-PARTNERS         *
*                  IN SUBROUTINE 'TIMEB'.                              *
*          MK2P  - ARRAY COMPOSING INDICES OF SECOND NUCLEON-PARTNER   *
*                  FOR PION ABSORPTION IN SUB. 'TIMEB'                 *
*          MK1T  - ARRAY COMPOSING INDICES OF NUCLEON-PARTNERS         *
*                  IN SUBROUTINE 'TIMEA'                               *
*          MK2T  - ARRAY COMPOSING INDICES OF SECOND NUCLEON-PARTNER   *
*                  FOR PION ABSORPTION IN SUB. 'TIMEA'                 *
*          MKAS  - ARRAY COMPOSING INDICES OF CASCADE-PARTNERS         *
*                  IN SUBROUTINE 'TIMED'                               *
*      ACTIVE:                                                         *
*          MPA   - ARRAY OF INDICES INDICATING ABILITY OF PROJECTILE   *
*                  NUCLEONS TO INTERACT WITH TARGET NUCLEONS           *
*          MYP   - ARRAY OF INDICES INDICATING ABILITY OF CASCADE      *
*                  PARTICLES TO INTERACT WITH NUCLEONS OF PROJECTILE   *
*                  NUCLEUS.                                            *
*          MYT   - ARRAY OF INDICES INDICATING ABILITY OF CASCADE      *
*                  PARTICLES TO INTERACT WITH NUCLEONS OF TARGET       *
*                  NUCLEUS.                                            *
*      XYZIZ:                                                          *
*          XYZ   - ARRAY OF COORDINATES OF NUCLEONS IN NUCLEI.         *
*          IZ    - ARRAY OF THEIR CHARGES.                             *
*                                                                      *
*      PMIM:                                                           *
*                   IN SUBROUTINES CALCULATING CHARACTERISTICS OF      *
*                 HADRON-HADRON INTERACTIONS OR DECAY OF RESONANCE     *
*          PM(1,J) -                                                   *
*          PM(2,J) - THESE ARE MOMENTUM COMPONENTS OF PRODUCED PARTICLE,
*          PM(3,J) -                                                   *
*          PM(4,J) - MODULE OF MOMENTUM.                               *
*                                                                      *
*                    IN OTHER SUBROUTINES                              *
*          PM(1,J) - X-COORDINATE OF CASCADE PARTICLE                  *
*          PM(2,J) - Y-COORDINATE OF CASCADE PARTICLE                  *
*          PM(3,J) - Z-COORDINATE OF CASCADE PARTICLE                  *
*          PM(4,J) - SIN OF ZENITH ANGLE CASCADE PARTICLE              *
*          PM(5,J) - COS OF ZENITH ANGLE CASCADE PARTICLE              *
*          PM(6,J) - SIN OF AZIMUTHAL ANGLE CASCADE PARTICLE           *
*          PM(7,J) - COS OF AZIMUTHAL ANGLE CASCADE PARTICLE           *
*          PM(8,J) - KINETIC ENERGY OF CASCADE PARTICLE                *  
*                                 *****                                *
*          PM(9,J) - MASS OF PARTICLE                                  *
*          PM(10,J) - WIDTH OF RESONANCE                               *
*          PM(11,J) - TIME TO THE NEXT INTERACTION OF CASCADE PARTICLE *
*                      IN PROJECTILE NUCLEUS.                          *
*          PM(12,J) - DECAY TIME OF RESONANCE IN PROJECTILE REGION     *
*          PM(13,J) - TIME TO THE NEXT INTERACTION OF CASCADE PARTICLE *
*                      IN TARGET NUCLEUS.                              *
*          PM(14,J) - DECAY TIME OF RESONANCE IN TARGET REGION.        *
*          PM(15,J) - TIME TO THE NEXT CASCADE-CASCADE INTERACTION.    *
*          PM(16,J) - EVOLUTION TIME OF CASCADE PARTICLE IN PROJECTILE *
*                     NUCLEUS.                                         *
*          PM(17,J) - EVOLUTION TIME OF CASCADE PARTICLE IN TARGET     *
*                     NUCLEUS.                                         *
*          PM(18,J) - VALUE INDICATING THE NUMBER OF INTERACTIONS.     *
*          PM(19,J) - FEYNMAN VARIABLE OF CASCADE PARTICLE.            *
************************************************************************  
      parameter (max1=2000)
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/LNUC/BR(2)                                                 
      common/inel/inel              
      COMMON/PMIM/pm(19,max1),IM(5,max1)                                
      COMMON/KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KS1,LST
      COMMON /KSALL/KSALL,KPTAL,KPT0                                    
      COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     1MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      COMMON /ACTIV/ MPA(240),MYP(max1),MYT(max1)                       
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
      DIMENSION P1(11),P2(11),P3(11),IP1(3),IP2(3),IP3(3)               
      DIMENSION V12(3),V0(3),PTEMP(11),PCT(3)                           
      DIMENSION C(11),IC(2)                                             
c      DIMENSION PMEM(max1)                                              
      N1=0                                                              
      MVOLD=0                                                           
      DO 30 IPME=1,max1                                                 
      DO 30 LPME=1,19
   30 PM(LPME,IPME)=0.
      DO 333 IPME=1,MAX1
      DO 333 LPME=1,5
333   IM(LPME,IPME)=0	                                                           
      DO 1 N=1,2                                                        
      AC(N)=A(N)                                                        
      ZC(N)=Z(N)                                                        
      E(N)=0.                                                           
      DO 1 I=1,3                                                        
      PC(N,I)=0.                                                        
   1  AMC(N,I)=0.                                                       
      KPT1=0                                                            
      KYP1=0                                                            
      KYT1=0                                                            
      KS1=0                                                             
      TINT=0.                                                           
      MV=0                                                              
      DO 24 J=1,240                                                     
      NK1P(J)=0                                                         
      TIME(J)=-.1                                                       
      MPA(J)=1                                                          
   24 CONTINUE                                                          
      IF(A(1)-1.)31,31,32                                               
   31 R(1)=0.                                                           
      SMP=A(1)
      GOTO 332
   32 SMP=.94
  332 R12=R(1)+R(2)+DEL+1./(5.07*SQRT(WT*(WT+2*SMP)))                   
      IF(R(1).GT.0..OR.R(2).GT.0.) GOTO 71
********  HADRON - HADRON INTERACTIONS  *************
   70 CALL hADhAD(Mv)
      GOTO 72
   71 B1=RNDM(-1)                                                       
      B2=6.283185*RNDM(-1)                                              
      IF(BIMP.GE.0.) R12=BIMP
      B3=R12*SQRT(B1)                                                   
      R0X=B3*COS(B2)                                                    
      R0Y=B3*SIN(B2)                                                    
      R0Z=-(R(1)*.94/(WT+.94)+R(2)+DEL)*SQRT(1.-B1)                     
      PIMPAC=SQRT(R0X**2+R0Y**2)                                        
      IF(A(1)-1.)33,33,34                                               
   33 PM(1,1)=R0X                                                       
      PM(2,1)=R0Y                                                       
      PM(3,1)=R0Z                                                       
      PM(4,1)=0.                                                        
      PM(5,1)=1.                                                        
      PM(6,1)=0.                                                        
      PM(7,1)=1.                                                        
      PM(8,1)=WT                                                        
       PM(9,1)=.94                                                      
       IM(1,1)=Z(1)                                                     
       IM(2,1)=1                                                        
      PM(10,1)=0.                                                       
      PM(19,1)=1.
      PM(16,1)=0.
      IM(3,1)=0                                                         
      IM(4,1)=0                                                         
      MYP(1)=0                                                          
      MYT(1)=1                                                          
      MK1T(1)=0                                                         
      MK2T(1)=0                                                         
      R(1)=0.                                                           
      MV=1                                                              
      AC(1)=0.                                                          
      ZC(1)=0.                                                          
   34 CONTINUE                                                          
      IF(R(1)-.1) 3,3,4                                                 
   3  OBR1=0.                                                           
      GO TO 5                                                           
   4  OBR1=0.00144*ZC(1)/R(1)                                           
    5 IF(R(2)-.1)60,60,61                                               
   60 OBR2=0.                                                           
      PM(1,1)=0.                                                        
      PM(2,1)=0.                                                        
      PM(3,1)=0.                                                        
      PM(4,1)=0.                                                        
      PM(5,1)=-1.                                                       
      PM(6,1)=0.                                                        
       PM(7,1)=1.                                                       
      PM(8,1)=0.                                                        
      PM(9,1)=.94                                                       
      PM(10,1)=0.                                                       
      PM(19,1)=1.
      PM(17,1)=0.
      IM(3,1)=0                                                         
      IM(4,1)=0                                                         
      IM(1,1)=Z(2)                                                      
      IM(2,1)=1                                                         
      MYP(1)=1                                                          
      MYT(1)=0                                                          
      MK1P(1)=0                                                         
      MK2P(1)=0                                                         
      MV=1                                                              
      AC(2)=0.                                                          
      ZC(2)=0.                                                          
      GOTO 65                                                           
   61 OBR2=0.00144*ZC(2)/R(2)                                           
   65 CONTINUE                                                          
   6  NA1=AC(1)+.1                                                      
      NA2=AC(2)+.1                                                      
      NIN=0
      nsm=0
      print*,'   *** projectile nucleus ***'
      DO 101 II=1,NA1  
  101 PRINT 103,II,(XYZ(1,J,II),J=1,3)
      print*,'   *** target nucleus ***'   
      DO 102 II=1,NA2 
  102 PRINT 103,II,(XYZ(2,J,II),J=1,3)
  103 FORMAT(5X,I5,10X,3F10.3)                   
      NIN=0
       LST=KPT1+KYP1+KYT1+KS1
*** DETERMINATION OF THE RANGE WHERE INTERACTION OCCURS 
*   AND IT'S TYPE;
* NR=1 - COLLISION OF PROJECTILE NUCLEONS WITH TARGET ONES;
* NR=2 - INTERACTION OF CASCADE PARTICLE WITH PROJECT. NUCLEON;
* NR=3 - INTERACTION OF CASCADE PARTICLE WITH TARGET NUCLEON;
* NR=4 - INTERACTION OF CASCADE PARTICLES WITH EACH OTHER.
      CALL PUNCTN(M,P1,P2,IP1,IP2,V12,U12,T12,SIG,SAB,TINT,NR,          
     *N1,N2,N3,NA1,NA2,MV,DEL,R0X,R0Y,R0Z,WT,VT,VP,DRR,TPR,PDL)        
      IF(M)15,15,7                                                      
    7 PRINT 43,TINT,NR                                            
   43 FORMAT(15X,'CASCADE EVOLUTION TIME = ',F7.3,5X,
     * 'TYPE = ',I2) 
      PRINT*,'*** CHARACTERISTICS OF INTERACTING PARTNERS ***'
      print 999
      CALL OUTPRT(N1,P1,IP1(1),IP1(2))                                 
*** IM(3,N1)=1 OR IM(4,N1)=1 - DECAY OF RESONANCE.
      IF(IM(4,N1).EQ.1.AND.NR.EQ.3) GOTO 411
      IF(IM(3,N1).EQ.1.AND.NR.EQ.2) GOTO 411
      CALL OUTPRT(N2,P2,IP2(1),IP2(2))
      IF(NR.EQ.4) GOTO 411                   
      IF(sab.GT.0) CALL OUTPRT(N3,P3,IP3(1),IP3(2))                          
  411 IF(MV) 77,77,41                                                   
*** SHIFT OF ALL CASCADE PARTICLES ACCORDING TO MINIMUM
*   INTERACTION TIME 'DRR'. 
   41 DO 44 J=1,MV                                                      
      VJ=SQRT(PM(8,J)*(PM(8,J)+2.D0*PM(9,J)))/(PM(8,J)+PM(9,J))         
      DRJ=DRR*VJ                                                        
      PM(1,J)=PM(1,J)+DRJ*PM(4,J)*PM(7,J)                               
      PM(2,J)=PM(2,J)+DRJ*PM(4,J)*PM(6,J)                               
      PM(3,J)=PM(3,J)+DRJ*PM(5,J)                                       
   44 CONTINUE                                                          
   77 IF(NR-2) 8,9,8                                                    
    8 NU=2                                                              
      IF(NR.EQ.1.OR.NR.EQ.4) GOTO 10                                         
      IF(IM(4,N1)-1) 10,20,10                                           
   20 IF(IM(2,N1)) 25,25,26                                             
   25 IF(PM(10,N1)-.1) 27,27,28                                         
   27 CALL W0DEC(0,V12,P1,IP1,N1,NP,MV)                                 
      GOTO 14                                                           
   28 CALL RODEC(0,V12,P1,IP1,N1,NP,MV)                                 
      GOTO 14                                                           
   26 CALL ISODEC(0,V12,P1,IP1,N1,NP,MV)                                
      GOTO 14                                                           
    9 NU=1                                                              
      IF(IM(3,N1)-1)10,40,10                                            
   40 IF(IM(2,N1)) 35,35,36                                             
   35 IF(PM(10,N1)-.1) 37,37,38                                         
   37 CALL W0DEC(0,V12,P1,IP1,N1,NP,MV)                                 
      GOTO 13                                                           
   38 CALL RODEC(0,V12,P1,IP1,N1,NP,MV)                                 
      GOTO 13                                                           
   36 CALL ISODEC(0,V12,P1,IP1,N1,NP,MV)                                
      GOTO 13                                                           
*** CALCULATION OF ELASTIC AND INELASTIC INTERACTIONS.
  10  CALL TIPULN(P1,IP1,P2,IP2,V12,U12,T12,SIG,SAB,MV,NP,NABS,         
     *P3,IP3,N3,NU,NIN)                                                 
       nsm=1
      PRINT 347,NP
  347 FORMAT(2X,'THE NUMBER OF PRODUCED PARTICLES IN 
     * IN HADRON-HADRON INTER.',2X,I5) 
       IF(NP) 6,6,11                                                    
  11  IF(NR-2) 12,13,14                                                 
  12  CALL PAULIC(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,TINT,         
     *R0X,R0Y,R0Z,IP,OBR1,OBR2,VT,VP,NIN)                               
      KPT1=KPT1+IP                                                      
      if(IP.EQ.0) KPTAL=KPTAL+1                                         
      IF(MV.EQ.0) KPT0=KPT0+1-IP                                        
      GOTO 50                                                            
  13  CALL PAULIB(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,TINT,         
     *R0X,R0Y,R0Z,IP,OBR1,VP,NIN)                                          
      KYP1=KYP1+IP                                            
      GOTO 50                                                            
   14 IF(NR-3) 145,145,144                                              
  145 CALL PAULIA(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,IP,OBR2,
     * VT,NIN)  
      KYT1=KYT1+IP                                               
      GOTO 50                                                            
  144 CALL PAULID(P1,P2,IP1,IP2,R0X,R0Y,R0Z,N1,N2,V12,NP,MV,IP,TINT,NIN)
      KS1=KS1+IP                                                        
      KSALL=KSALL+1                                                     
   50 IF(IP)51,51,52
   51 PRINT 53
   53 FORMAT(40X,'!!! PAULI VETO !!!') 
      GOTO 6
   52 PRINT*,'***CASCADE PARTICLES IN MEMORY ***'           
      PRINT 999
  999 FORMAT(3x,'N',5X,'X',7X,'Y',7X,'Z',7X,'PX',7X,'PY',7X,'PZ',7X,'T',
     * 8X,'M',5X,'Q',2X,'B') 
      DO 54 I=1,MV  
      DO 55 L=1,9 
   55 C(L)=PM(L,I)
      IC(1)=IM(1,I)
      IC(2)=IM(2,I)
      CALL OUTPRT(I,C,IC(1),IC(2)) 
   54 continue                                              
      GOTO 6
c   15 IF(E(1).GT..001.OR.E(2).GT..001) GOTO 17
 15   IF(R(1).EQ.0..OR.R(2).EQ.0..AND.Mv.GT.1) GOTO 17
      IF(R(1).GT.0..AND.R(2).GT.0..AND.Mv.GT.0) GOTO 17
  16  NEL=1                                                             
      GO TO 22                                                          
  17  DO 18 N=1,2                                                       
      DO 18 L=1,3                                                       
  18  AMC(N,L)=AMC(N,L)*5.06
c      GOTO 21
*** DECAY OF RESONANCES WHEN CASCADE PROCESS IS OVER.
  72  CALL RESDEC(MV)                                                   
  21  NEL=0                                                             
  22  RETURN                                                            
      END                                                                 
************************************************************************
      SUBROUTINE OUTPRT(I,A,IA,IB)
************************************************************************
* PURPOSE:  OUTPUT THE CHARACTERISTICS OF CASCADE PROCESS
************************************************************************
      DIMENSION A(11),P(3)
       T1=SQRT(A(8)*(A(8)+2.*A(9)))
      P(1)=T1*A(4)*A(7)
      P(2)=T1*A(4)*A(6)
      P(3)=T1*A(5)
      PRINT 99,I,(A(J),J=1,3),(P(K),K=1,3),A(8),A(9),IA,IB
   99 FORMAT(I4,3F8.2,5F9.3,2I3)
      RETURN 
      END                                                                    
************************************************************************
      SUBROUTINE BLNUNU(NEL,RN,DEL,TPR,MV,BIMP,PDL)                            
************************************************************************
*  PURPOSE:  THIS IS THE BASIC SUBROUTINE TO CALCULATE CASCADES;       *
*            IT IS THE SAME AS SUBROUTINE BLTREE EXCEPT PRINTING.      * 
************************************************************************    
      parameter (max1=2000)
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/LNUC/BR(2)                                                 
      common/inel/inel              
      COMMON/PMIM/pm(19,max1),IM(5,max1)                                
      COMMON /OLDN/OLDNY(max1)                                          
      COMMON/KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KS1,LST
      COMMON /KSALL/KSALL,KPTAL,KPT0                                    
     0COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     1MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      COMMON /ACTIV/ MPA(240),MYP(max1),MYT(max1)                       
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
      DIMENSION P1(11),P2(11),P3(11),IP1(3),IP2(3),IP3(3)               
      DIMENSION V12(3),V0(3),PTEMP(11),PCT(3)                           
      DIMENSION C(11),IC(2)                                             
      DIMENSION PMEM(max1)                                              
      N1=0                                                              
      MVOLD=0                                                           
      DO 30 IPME=1,max1                                                 
   30 PMEM(IPME)=0.                                                     
      DO 1 N=1,2                                                        
      AC(N)=A(N)                                                        
      ZC(N)=Z(N)                                                        
      E(N)=0.                                                           
      DO 1 I=1,3                                                        
      PC(N,I)=0.                                                        
   1  AMC(N,I)=0.                                                       
      KPT1=0                                                            
      KYP1=0                                                            
      KYT1=0                                                            
      KS1=0                                                             
      TINT=0.                                                           
      MV=0                                                              
      DO 24 J=1,240                                                     
      NK1P(J)=0                                                         
      TIME(J)=-.1                                                       
      MPA(J)=1                                                          
   24 CONTINUE                                                          
      IF(A(1)-1.)31,31,32                                               
   31 R(1)=0.                                                           
      SMP=A(1)
      GOTO 332
   32 SMP=.94
  332 R12=R(1)+R(2)+DEL+1./(5.07*SQRT(WT*(WT+2*SMP)))                   
      IF(R(1).GT.0..OR.R(2).GT.0.) GOTO 71
   70 CALL hADhAD(Mv)
      GOTO 72
   71 B1=RNDM(-1)                                                       
      B2=6.283185*RNDM(-1)                                              
      IF(BIMP.GE.0.) R12=BIMP
      B3=R12*SQRT(B1)                                                   
      R0X=B3*COS(B2)                                                    
      R0Y=B3*SIN(B2)                                                    
      R0Z=-(R(1)*.94/(WT+.94)+R(2)+DEL)*SQRT(1.-B1)                     
      PIMPAC=SQRT(R0X**2+R0Y**2)                                        
      IF(A(1)-1.)33,33,34                                               
   33 PM(1,1)=R0X                                                       
      PM(2,1)=R0Y                                                       
      PM(3,1)=R0Z                                                       
      PM(4,1)=0.                                                        
      PM(5,1)=1.                                                        
      PM(6,1)=0.                                                        
      PM(7,1)=1.                                                        
      PM(8,1)=WT                                                        
       PM(9,1)=.94                                                      
       IM(1,1)=Z(1)                                                     
       IM(2,1)=1                                                        
      PM(10,1)=0.                                                       
	PM(18,1)=1.    ! UJ 27.06.2003 to flag target and projectile protons
      PM(19,1)=1.
      PM(16,1)=0.
      IM(3,1)=0                                                         
      IM(4,1)=0                                                         
      MYP(1)=0                                                          
      MYT(1)=1                                                          
      MK1T(1)=0                                                         
      MK2T(1)=0                                                         
      R(1)=0.                                                           
      MV=1                                                              
      AC(1)=0.                                                          
      ZC(1)=0.                                                          
   34 CONTINUE                                                          
      IF(R(1)-.1) 3,3,4                                                 
   3  OBR1=0.                                                           
      GO TO 5                                                           
   4  OBR1=0.00144*ZC(1)/R(1)                                           
    5 IF(R(2)-.1)60,60,61                                               
   60 OBR2=0.                                                           
      PM(1,1)=0.                                                        
      PM(2,1)=0.                                                        
      PM(3,1)=0.                                                        
      PM(4,1)=0.                                                        
      PM(5,1)=-1.                                                       
      PM(6,1)=0.                                                        
       PM(7,1)=1.                                                       
      PM(8,1)=0.                                                        
      PM(9,1)=.94                                                       
      PM(10,1)=0.                                                       
      PM(19,1)=1.
      PM(17,1)=0.
      IM(3,1)=0                                                         
      IM(4,1)=0                                                         
      IM(1,1)=Z(2)                                                      
      IM(2,1)=1                                                         
      MYP(1)=1                                                          
      MYT(1)=0                                                          
      MK1P(1)=0                                                         
      MK2P(1)=0                                                         
      MV=1                                                              
      AC(2)=0.                                                          
      ZC(2)=0.                                                          
      GOTO 65                                                           
   61 OBR2=0.00144*ZC(2)/R(2)                                           
   65 CONTINUE                                                          
   6  NA1=AC(1)+.1                                                      
      NA2=AC(2)+.1                                                      
      NIN=0
      nsm=0
      LST=KPT1+KYP1+KYT1+KS1
      CALL PUNCTN(M,P1,P2,IP1,IP2,V12,U12,T12,SIG,SAB,TINT,NR,          
     *N1,N2,N3,NA1,NA2,MV,DEL,R0X,R0Y,R0Z,WT,VT,VP,DRR,TPR,PDL)        
      IF(M)15,15,7                                                      
    7 IF(MV) 77,77,41                                                   
   41 DO 44 J=1,MV                                                      
      VJ=SQRT(PM(8,J)*(PM(8,J)+2.D0*PM(9,J)))/(PM(8,J)+PM(9,J))         
      DRJ=DRR*VJ                                                        
      PM(1,J)=PM(1,J)+DRJ*PM(4,J)*PM(7,J)                               
      PM(2,J)=PM(2,J)+DRJ*PM(4,J)*PM(6,J)                               
      PM(3,J)=PM(3,J)+DRJ*PM(5,J)                                       
   44 CONTINUE                                                          
   77 IF(NR-2) 8,9,8                                                    
    8 NU=2                                                              
      IF(NR.EQ.1.OR.NR.EQ.4) GOTO 10                                         
      IF(IM(4,N1)-1) 10,20,10                                           
   20 IF(IM(2,N1)) 25,25,26                                             
   25 IF(PM(10,N1)-.1) 27,27,28                                         
   27 CALL W0DEC(0,V12,P1,IP1,N1,NP,MV)                                 
      GOTO 14                                                           
   28 CALL RODEC(0,V12,P1,IP1,N1,NP,MV)                                 
      GOTO 14                                                           
   26 CALL ISODEC(0,V12,P1,IP1,N1,NP,MV)                                
      GOTO 14                                                           
    9 NU=1                                                              
      IF(IM(3,N1)-1)10,40,10                                            
   40 IF(IM(2,N1)) 35,35,36                                             
   35 IF(PM(10,N1)-.1) 37,37,38                                         
   37 CALL W0DEC(0,V12,P1,IP1,N1,NP,MV)                                 
      GOTO 13                                                           
   38 CALL RODEC(0,V12,P1,IP1,N1,NP,MV)                                 
      GOTO 13                                                           
   36 CALL ISODEC(0,V12,P1,IP1,N1,NP,MV)                                
      GOTO 13                                                           
  10  CALL TIPULN(P1,IP1,P2,IP2,V12,U12,T12,SIG,SAB,MV,NP,NABS,         
     *P3,IP3,N3,NU,NIN)                                                 
       nsm=1
       IF(NP) 6,6,11                                                    
  11  IF(NR-2) 12,13,14                                                 
  12  CALL PAULIC(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,TINT,         
     *R0X,R0Y,R0Z,IP,OBR1,OBR2,VT,VP,NIN)                               
      KPT1=KPT1+IP                                                      
      if(IP.EQ.0) KPTAL=KPTAL+1                                         
      IF(MV.EQ.0) KPT0=KPT0+1-IP                                        
      GOTO 6                                                            
  13  CALL PAULIB(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,TINT,         
     *R0X,R0Y,R0Z,IP,OBR1,VP,NIN)                                          
      KYP1=KYP1+IP                                            
      GOTO 6                                                            
   14 IF(NR-3) 145,145,144                                              
  145 CALL PAULIA(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V12,NP,MV,IP,OBR2,
     * VT,NIN)  
      KYT1=KYT1+IP                                               
      GOTO 6                                                            
  144 CALL PAULID(P1,P2,IP1,IP2,R0X,R0Y,R0Z,N1,N2,V12,NP,MV,IP,TINT,NIN)
      KS1=KS1+IP                                                        
      KSALL=KSALL+1                                                     
      GO TO 6                                                           
   15 IF(E(1).GT..005.OR.E(2).GT..005) GOTO 17
      IF(R(1).EQ.0..OR.R(2).EQ.0..AND.Mv.GT.1) GOTO 17
      IF(R(1).GT.0..AND.R(2).GT.0..AND.Mv.GT.0) GOTO 17
  16  NEL=1                                                             
      GO TO 22                                                          
  17  DO 18 N=1,2                                                       
      DO 18 L=1,3                                                       
  18  AMC(N,L)=AMC(N,L)*5.06
c      GOTO 21
  72  CALL RESDEC(MV)                                                   
  21  NEL=0                                                             
  22  RETURN                                                            
      END                                                               
************************************************************************
	SUbROUTINE PARAM(DEL,TPR,PDL)
************************************************************************
*   PURPOSE:    CALCULATE PARAMETERS OF GAUSSIAN DENSITY DISTRIBUTION  *
*               FOR LIGHT NUCLEI (A<12).                               *
************************************************************************
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON /LNUC/ BR(2)                                               
      COMMON/BL1000/AM,AMF                                                
      COMMON/BL0999/RADNCL                                                    
        AM=.1
        AMF=.1
        RADNCL=1.3
        PDL=1.
	DEL=1.7
cuj        TPR=1.4	! uj 2/07/2003
        TPR=2.0		! uj 2/07/2003
	VPI=.025
	EPS(1)=.007
	EPS(2)=.007
	CR(1)=.545
	CR(2)=.545
	DR(1)=.18
	DR(2)=.18
        BR(1)=.0225*A(1)+1.37
        BR(2)=.0225*A(2)+1.37
	RETURN
	END                                                              
************************************************************************    
      SUbROUTINE hADhAD(NP)
************************************************************************
*     PURPOSE : CALCULATION OF HADRON-HADRON INTERACTIONS              *
************************************************************************    
      parameter (max1=2000)
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2)
      COMMON/PMIM/ pm(19,max1),IM(5,max1)
       DIMENSION P1(11),P2(11),IP1(3),IP2(3),v(3),P(3),PN(11)
      Mv=0
      P2(5)=1.-2.*RNDM(-1)
      P2(4)=SQRT(1.-P2(5)**2)
       P2(6)=0.
       P2(7)=1.
   21 S=3.*RNDM(-1)
      IF(RNDM(-1)-S*EXP(1.-S)) 22,22,21
   22  P2(8)=.03*S
       P2(9)=.94
        P2(11)=1.
       IP2(1)=Z(2)
        IF(A(2)) 15,15,16
   15  IP2(2)=0
       P2(9)=.14
       GOTO 17
   16  IP2(2)=1
   17  P1(4)=0.
       P1(5)=1.
       P1(6)=0.
       P1(7)=1.
       P1(8)=WT
        P1(11)=1.
       IF(A(1))1,1,2
   1   P1(9)=.14
       IP1(2)=0
       GOTO 3
   2    P1(9)=.94
       IP1(2)=1
   3   IP1(1)=Z(1)
       CALL TINvU(P1,P2,T,v,U)
       NCT=0
   10   CALL bINEL(P1,IP1,P2,IP2,V,U,T,Mv,NP,NIN)
       IF(NIN) 11,11,12
   12  NCT=NCT+1
       IF(NCT-10) 10,10,11
   11   IF(NP-2)7,7,8
   7   DO 9 L=1,10
   9   PM(L,2)=PM(L,3)
        IM(1,2)=IM(1,3)
        IM(2,2)=IM(2,3)
   8    DO 4 I=1,NP
       DO 5 J=1,3
   5   P(J)=PM(J,I)
       PN(9)=PM(9,I)
       CALL CINEMA(P,v,PN)
       DO 6 J=1,9
   6   PM(J,I)=PN(J)
   4   CONTINUE
       RETURN
        END
***********************************************************************
	SUBROUTINE MOM(k)
***********************************************************************
	PARAMETER (max1=2000)
	COMMON /PMIM/pm(19,max1),IM(5,max1)
	PMD=SQRT(PM(8,K)*(PM(8,K)+2.*PM(9,K)))
	PM(1,K)=PMD*PM(4,K)*PM(6,K)
	PM(2,K)=PMD*PM(4,K)*PM(7,K)
	PM(3,K)=PMD*PM(5,K)
	RETURN 
	END      
************************************************************************
      SUBROUTINE PUNCTN(M,P1,P2,IP1,IP2,V1,U1,TR2,SG,SAB,TINT,NR,       
     *N1,N2,N3,NA1,NA2,MV,DEL,R0X,R0Y,R0Z,WT,VT,VP,DR,TPR,PDL)         
************************************************************************
*    PURPOSE:  CHOOSES THE INTERACTING PAIR OF HADRONS OR DECAYING     *
*              RESONANCE WITH MINIMUM COLLISION OR DECAY TIME AMONG    *
*              ALL POSSIBLE ONES.                                      *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     NA1 - NUMBER OF NUCLEONS IN THE PROJECTILE NUCLEUS DURING        *
*           CASCADE EVOLUTION ( NA1=AC(1) ).                           *    
*     NA2 - NUMBER OF NUCLEONS IN THE TARGET NUCLEUS DURING            *
*           CASCADE EVOLUTION ( NA2=AC(2) ).                           *
*     MV  - CURRENT NUMBER OF CASCADE PARTICLES IN MEMORY.             *
*     DEL - PARAMETER SPECIFYING A RADIUS OF STRONG INTERACTIONS;      *
*           R0X,R0Y,R0Z - COORDINATES OF ENTRY POINT OF PROJECTILE     *
*           NUCLEUS IN LAB SYSTEM.                                     *
*     WT   - BEAM ENERGY IN LAB [GEV/NUCLEON]                          *
*     TINT - TIME OF THE EVOLUTION OF COLLISION PROCESS.               *
*     TPR  - FORMATION TIME PARAMETER.                                 *
*     PDL  - PARAMETER SPECIFYING NUCLEAR TRANSPARENCY                 *    
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES:                                                 *
*     M         M=1 -  THERE ARE INTERACTING PARTNERS OR               *
*                      OR DECAYING RESONANCE                           *
*               M=0 -  THERE ARE NOT INTERACTING PARTNERS              *
*                      OR DECAYING RESONANCES; CASCADE PROCESS IS OVER.*
*     P1,P2 -          ARRAYS OF KINEMATICAL CHARACTERISTICS OF        *
*                      INTERACTING PARTNERS OR DECAYING RESONANCE.     *
*               P1(1),P1(2),P1(3), P2(1),P2(2),P2(3) - COORDINATES;    *
*               P1(4),P1(5), P2(4),P2(5) - SIN AND COS OF ZENITH ANGLE;*
*               P1(6),P1(7), P2(6),P2(7) - SIN AND COS OF AZIMUTH ANGLE;
*               P1(8), P2(8) - KINETIC ENERGIES;                       *
*               P1(9), P2(9) - MASSES;                                 *
*               P1(10), P2(10) - WIDTHS OF RESONANCES;                 *
*               P1(11), P2(11) - FRACTIONS OF NORMAL CROSS SECTIONS.   *  
*     IP1,IP2 -        THEIR CHARGES AND BARYON NUMBERS.               *
*                      IP1(1),IP2(1) - CHARGES                         *
*                      IP1(2),IP2(2) - BARYON NUMBERS                  *
*                      IP1(3),IP2(3) -NOT USED                         *
*     V1            -  VELOCITY OF C.M.S. OF COLLIDING PARTICLES       *
*                      IN LAB.SYSTEM.                                  *
*     U1            - TOTAL ENERGY OF SYSTEM OF COLLIDING PARTICLES    *
*                     IN C.M.S.                                        *
*     TR2        - KINETIC ENERGY OF PARTICLE P1 IN THE REST FRAME     *
*                  OF PARTICLE P2.                                     *
*     SG         - TOTAL CROSS SECTION FOR THE INTERACTION P1 AND P2   *
*     SAB        - CROSS SECTION FOR THE ABSORPSION OF PION            *
*                  ON QUASIDEUTERON PAIR.                              *
*     NR         _ LABEL INDICATING THE TYPE OF THE INTERACTION        *
*                  NR=1 --> INTERACTION BETWEEN PROJECTILE AND         *
*                            TARGET NUCLEONS;                          *
*                  NR=2 --> INTERACTION BETWEEN CASCADE PARTICLE       *
*                           AND PROJECTILE NUCLEON;                    *
*                  NR=3 --> INTERACTION BETWEEN CASCADE PARTICLE       *
*                           AND TARGET NUCLEON;                        *
*                  NR=4 --> INTERACTION BETWEEN CASCADE PARTICES.      *
*     N1         - NO OF THE FIRST INTERACTING PARTICLE.               *
*     N2         - NO OF THE SECOND INTERACTING PARTICLE.              *
*     N3         - NO OF NUCLEON-PARTNER FOR PION ABSORBSION           *
*                  ON QUASIDEUTERON PAIR.                              *
*     VT         - NOT USED.                                           *
*     VP         - NOT USED.                                           *
************************************************************************     
       parameter (max1=2000)
      COMMON/PMIM/ pm(19,max1),IM(5,max1)                               
      COMMON/KYPT/ KPT,KYP,KYT,KPT1,KYP1,KYT1,KS1,LST                   
      COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     *MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      common/inel/inel              
      DIMENSION P1(11),P2(11),PT1(11),PT2(11),YT1(11),YT2(11),YP1(11)   
      DIMENSION IP1(3),IP2(3),IPT1(3),IPT2(3),IYT1(3),IYT2(3),YP2(11)   
      DIMENSION IYP1(3),IYP2(3),V1(3),VPT(3),VYP(3),VYT(3)              
      DIMENSION PIN(11),PN(11),IIN(3),IPN(3),VKS(3)                     
      TPT=-.1                                                           
      TYP=-.1                                                           
      TYT=-.1                                                           
      TYK=-.1                                                           
      TU=-.1                                                            
      TAU=-.1                                                           
      TY=-.1                                                            
      R0=SQRT(R0X**2+R0Y**2+R0Z**2)*1.01                                
      V0=SQRT(WT*(WT+1.88D0))/(WT+.94D0)                                
      RT=SQRT(R0X**2+R0Y**2+(R0Z+V0*TINT)**2)                           
      DRT=RT-R0
      IF(DRT)1,1,44
    1 IF(NA1) 3,3,52                                                    
   52 IF(NA2) 3,3,2                                                     
   2  CALL TIMEC(0,NA1,NA2,R0X,R0Y,R0Z,TINT,DEL,TPT,PT1,PT2,IPT1,       
     *IPT2,VPT,UPT,TRPT,SIGPT,SABPT,NPT1,NPT2,NPT3,VT,VP,PDL)               
    3 IF(MV-1) 7,44,4                                                   
    4 CALL TIMED(0,MV,DEL,TYK,PIN,PN,IIN,IPN,VKS,UK,TNK,SIGK,NK1,NK2,
     *TPR)
   44 IF(NA1) 6,6,5                                                     
   5  CALL TIMEB(0,MV,NA1,R0X,R0Y,R0Z,TINT,DEL,TYP,YP1,YP2,IYP1,IYP2,   
     *VYP,UYP,TRYP,SIGP,SABYP,NYP1,NYP2,NYP3,VMTP,VMPP,VLB,TPR,PDL)    
    6 IF(NA2)7,7,66                                                     
   66 CALL TIMEA(0,MV,NA2,DEL,TYT,YT1,YT2,IYT1,IYT2,VYT,UYT,TRYT,       
     *SIGT,SABYT,  NYT1,NYT2,NYT3,VMTT,VMPT,VLA,TPR,PDL)               
    7 IF(TYP) 12,12,8                                                   
    8 IF(TYT) 9,9,10                                                    
    9 TY=TYP                                                            
      NR=2                                                              
      GO TO 13                                                          
   10 IF(TYP-TYT)  9,99,11                                              
   99 IF(VLA-VLB) 11,9,9                                               
   11 TY=TYT                                                            
      NR=3                                                              
      GO TO 13                                                          
   12 IF(TYT) 13,13,11                                                  
   13 IF(TPT) 17,17,14                                                  
   14 IF(TYK) 15,15,16                                                  
   15 TU=TPT                                                            
      IF(TY)67,67,68                                                    
   68 IF(TY-TU) 18,67,67                                                
   67 NR=1                                                              
      GOTO 56                                                           
   16 IF(TPT-TYK) 15,15,54                                              
   54 TU=TYK                                                            
      IF(TY) 70,70,71                                                   
   71 IF(TY-TU) 18,18,70                                                
   70 NR=4                                                              
      GOTO 56                                                           
   17 IF(TYK) 55,55,54                                                  
   55 IF(TY) 35,35,18                                                   
   18 TAU=TY                                                            
      GOTO 41                                                           
   56 TAU=TU                                                            
      GOTO 41                                                           
  119 IF(MV) 22,22,20                                                   
   20 DO 21 K=1,MV                                                      
      PM(11,K)=PM(11,K)-TAU                                             
      PM(12,K)=PM(12,K)-TAU                                             
      PM(13,K)=PM(13,K)-TAU                                             
      PM(14,K)=PM(14,K)-TAU                                             
      PM(15,K)=PM(15,K)-TAU                                             
      PM(16,K)=PM(16,K)+TAU
      PM(17,K)=PM(17,K)+TAU
   21 CONTINUE                                                          
      DR=TAU                                                            
   22 TINT=TINT+TAU                                                     
      IF(NA1)36,36,42                                                   
   42 DO 40 K=1,NA1                                                     
   40 TIME(K)=TIME(K)-TAU                                               
      GOTO 36                                                           
   41 IF(NR-1) 23,23,27                                                 
   23 CALL TIMEC(1,NA1,NA2,R0X,R0Y,R0Z,TINT,DEL,TPT,PT1,PT2,IPT1,       
     *IPT2,VPT,UPT,TRPT,SIGPT,SABPT,NPT1,NPT2,NPT3,VT,VP,PDL)               
      DO 24 L=1,11                                                      
      P1(L)=PT1(L)                                                      
   24 P2(L)=PT2(L)                                                      
      DO 25 L=1,3                                                       
      IP1(L)=IPT1(L)                                                    
      IP2(L)=IPT2(L)                                                    
   25 V1(L)=VPT(L)                                                      
      KPT=KPT+1                                                         
      U1=UPT                                                            
      TR2=TRPT                                                          
      SG=SIGPT                                                          
      SAB=SABPT                                                         
      N1=NPT1                                                           
      N2=NPT2                                                           
      N3=NPT3                                                           
      GO TO 34                                                          
   27 IF(NR-2) 28,28,31                                                 
   28 CALL TIMEB(1,MV,NA1,R0X,R0Y,R0Z,TINT,DEL,TYP,YP1,YP2,IYP1,IYP2,   
     *VYP,UYP,TRYP,SIGP,SABYP,NYP1,NYP2,NYP3,VMTP,VMPP,VLB,TPR,PDL)    
      IF(IM(3,NYP1)-1)48,48,47                                          
   47 KYP=KYP+1                                                         
      U1=UYP                                                            
      TR2=TRYP                                                          
      VT=VMTP                                                           
      VP=VMPP                                                           
      SG=SIGP                                                           
      SAB=SABYP                                                         
   48 DO 29 L=1,11                                                      
      P1(L)=YP1(L)                                                      
   29 P2(L)=YP2(L)                                                      
      DO 30 L=1,3                                                       
      IP1(L)=IYP1(L)                                                    
      IP2(L)=IYP2(L)                                                    
   30 V1(L)=VYP(L)                                                      
      N1=NYP1                                                           
      N2=NYP2                                                           
      N3=NYP3                                                           
      GO TO 34                                                          
   31 IF(NR-3) 58,58,59                                                 
   58 CALL TIMEA(1,MV,NA2,DEL,TYT,YT1,YT2,IYT1,IYT2,VYT,UYT,TRYT,       
     *SIGT,SABYT,NYT1,NYT2,NYT3,VMTT,VMPT,VLA,TPR,PDL)                 
      IF(IM(4,NYT1)-1) 38,38,39                                         
   39 KYT=KYT+1                                                         
      U1=UYT                                                            
      TR2=TRYT                                                          
      VT=VMTT                                                           
      VP=VMPT                                                           
      SG=SIGT                                                           
      SAB=SABYT                                                         
   38 DO 32 L=1,11                                                      
      P1(L)=YT1(L)                                                      
   32 P2(L)=YT2(L)                                                      
      DO 33 L=1,3                                                       
      IP1(L)=IYT1(L)                                                    
      IP2(L)=IYT2(L)                                                    
   33 V1(L)=VYT(L)                                                      
      N1=NYT1                                                           
      N2=NYT2                                                           
      N3=NYT3                                                           
      GOTO 34                                                           
  59  CALL TIMED(1,MV,DEL,TYK,PIN,PN,IIN,IPN,VKS,UK,TNK,SIGK,NK1,NK2,   
     *TPR)
      U1=UK                                                             
      TR2=TNK                                                           
      VT=0.                                                             
      VP=0.                                                             
      SG=SIGK                                                           
      SAB=0.                                                            
      DO 60 L=1,11                                                      
      P1(L)=PIN(L)                                                      
   60 P2(L)=PN(L)                                                       
      DO 61 L=1,3                                                       
      IP1(L)=IIN(L)                                                     
      IP2(L)=IPN(L)                                                     
   61 V1(L)=VKS(L)                                                      
      N1=NK1                                                            
      N2=NK2                                                            
      N3=0                                                              
   34 M=1                                                               
      GO TO 119                                                         
   35 M=0                                                               
   36 RETURN                                                            
      END                                                               
***********************************************************************
      SUBROUTINE TIMEA(MZ,MV,NA2,DEL,TAU,P1,P2,IP1,IP2,                 
     *V1,U1,TR2,SIG,SAB,N1,N2,N3,VMT,VMP,VLI,TPR,PDL)                    
************************************************************************
*   PURPOSE : MZ=0 - CHOOSES THE INTERACTION OF CASCADE PARTICLE WITH  *
*             TARGET NUCLEON OR DECAYING RESONANCE WITH MINIMUM        *
*             TIME;                                                    *
*             MZ=1 - CALCULATES KINEMATICAL CHARACTERISTICS OF         *
*             INTERACTING PARTNERS OR DECAYING RESONANCE.              *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     MZ  - LABEL INDICATING FIRST AND SECOND RUN OF SUBROUTINE:       *
*           IN THE FIRST RUN MZ=0 - SUBROUTINE CHOOSES PARTNERS        * 
*           OR DECAYING RESONANCE WITH MINIMUM TIME;                   *
*           IN THE SECOND RUN MZ=1, - SUBROUTINE CALCULATES KINEMATICAL*
*           CHARACTERISTICS OF PARTNERS OR DECAYING RESONANCE.         *
*     MV  - CURRENT NUMBER OF CASCADE PARTICLES IN THE ARRAYS PM AND IM*
*     NA2 - CURRENT NUMBER OF NUCLEONS IN THE TARGET NUCLEUS, NA2=AC(2)*
*     DEL - PARAMETER DETERMINING A RADIUS OF STRONG INTERACTION.      *
*     TPR - PARAMETER SPECIFYING FORMATION TIME.                       *
*     PDL  - PARAMETER SPECIFYING NUCLEAR TRANSPARENCY                 *
*----------------------------------------------------------------------*
*    OUTPUT QUANTITIES:                                                *
*     TAU     - MINIMUM TIME TO INTERACTION.                           *
*     P1,P2   - ARRAYS OF KINEMATICAL CHARACTERISTICS OF INTERACTING   *
*               PARTNERS OR DECAYING RESONANCE WITH MINIMUM TIME.      *
*     IP1,IP  - ARRAYS OF THEIR CHARGES AND BARYON NUMBERS.            *
*     V1      - VELOCITY OF CMS OF INTERACTING PARTNERS OR             *
*               DECAYING RESONANCE.                                    *
*     U1      - TOTAL ENERGY OF INTERACTING PARTNERS IN CMS.           *
*     TR2     - KINETIC ENERGY OF INTERACTING CASCADE PARTICLE         *
*               IN REST FRAME OF TARGET NUCLEON.                       *
*     SIG     - TOTAL CROSS SECTION FOR THE INTERACTION.               *
*     SAB     - CROSS SECTION FOR ABSORPTION OF PION ON QUASIDEUTERON  *
*               PAIR.                                                  *  
*     N1      - NO OF INTERACTING CASCADE PARTNER IN ARRAYS PM AND IM. *
*     N2      - NO OF INTERACTING TARGET NUCLEON IN ARRAYS XYZ         *
*               (COORNINATES) AND IZ (NUCLEON CHARGES).                *
*     N3      - NO OF NUCLEON-PARTHER FOR PION ABSORBSION IN ARRAYS    *
*               XYZ AND IZ.                                                
*     VMT     - POTENTIAL FOR INTERACTING CASCADE PARTICLE.            *
*     VMP     - NOT USED.                                              *
*     VLI     - VELOCITY OF DECAYING RESONANCE.                        *
************************************************************************      
      parameter (max1=2000)
      common/inel/inel              
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
      COMMON/PMIM/ pm(19,max1),IM(5,max1)                               
     0COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     1MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      COMMON /ACTIV/ MPA(240),MYP(max1),MYT(max1)                       
      DIMENSION V0(3),PIN(11),IIN(3),P(11),IP(3),PN(11),IPN(3)          
      DIMENSION P1(11),P2(11),IP1(3),IP2(3),V1(3),V(3),B(3),IT(3)       
      DIMENSION KP(2),KT(2),SGM(MAX1),PA2(3)
      SAVE SGM,M1,M2                                                          
      TPAR=TPR 
      IF(MZ)22,23,22                                                    
   22 K=N1                                                              
      GOTO 34                                                           
   23 K=1                                                               
      TAU=-.1                                                           
    1 IF(MYT(K))70,70,2                                                 
   70 IF(PM(10,K))15,15,71                                              
   71 IF(PM(14,K))34,34,65                                              
    2 IF(PM(10,K))36,36,63
c******** correction 31 July 1999 *********                                  
   36 IF(MK1T(K))34,34,644
 644  if(pm(9,k)-.9) 645,64,64
 645  if(mk2t(k)) 34,34,64
c******************************************
   63 IF(IM(4,K)-1) 34,71,36                                            
   65 TAUK=PM(14,K)                                                     
c******** correction 31 July 1999 *********
      if(mk1t(k))34,34,66
c******************************************
   64 TAUK=PM(13,K)                                                     
   66 K1=MK1T(K)                                                        
      K2=MK2T(K)                                                        
      GOTO 33                                                           
   34 XK0=PM(1,K)                                                       
      YK0=PM(2,K)                                                       
      ZK0=PM(3,K)                                                       
      DO 3 L=1,10                                                       
    3 PIN(L)=PM(L,K)                                                    
C      VMTI=POTENN(2,A(2),PIN,IIN,AR(2),CR(2),DR(2),TR(2),EPS(2),VPI)    
C      PIN(8)=PM(8,K)+VMTI                                   
      PORG=PM(19,K)
      IIN(1)=IM(1,K)                                                    
      IIN(2)=IM(2,K)                                                    
      IIN(3)=1                                                          
      ITT=0                                                             
      IF(PIN(8).LT.TR(2)/2..AND.IIN(2).GT.0) GOTO 15
      PMO=SQRT(PIN(8)*(PIN(8)+2.D0*PIN(9)))
      EMO=PIN(8)+PIN(9)
      VLI=PMO/EMO                                                       
      GM=EMO/PIN(9)
      IF(MZ) 21,21,4                                                    
    4 K1=MK1T(K)                                                        
      K2=MK2T(K)                                                        
      IF(IM(4,K)-1)47,47,8                                              
   47 DD=PM(14,K)*VLI                                                   
      TAUK=PM(14,K)                                                     
      V(1)=VLI*PM(4,K)*PM(7,K)                                          
      V(2)=VLI*PM(4,K)*PM(6,K)                                          
      V(3)=VLI*PM(5,K)                                                  
      GOTO 24                                                           
   21 IF(PM(10,K))611,611,72                                              
c********* correction 31 July 1999 *********
 611  if(pin(8)-0.015)15,15,61
   72 IF(IM(4,K)-1)37,622,77                                             
 622  if(mk1t(k))77,77,62
c*******************************************
   62 K1=MK1T(K)                                                        
      K2=MK2T(K)                                                        
   37 CALL PLAM(PIN,DS,TRES,TPR)                                        
      PM(14,K)=TRES                                                     
      IF(MYT(K))43,43,73                                                
   73 IF(IM(4,K)-1)61,74,61                                             
   74 IF(PM(14,K)-PM(13,K))75,75,76                                     
   75 TAUK=PM(14,K)                                                     
      IM(4,K)=1                                                         
      GOTO 33                                                           
   76 TAUK=PM(13,K)                                                     
      IM(4,K)=2                                                         
      GOTO 33                                                           
   77 TRES=PM(14,K)                                                     
      DS=TRES*VLI                                                       
   61 CONTINUE                                             
   67 DLK=1./(5.07*SQRT(PIN(8)*(PIN(8)+2.*PIN(9))))+DEL                 
    6 CALL BLCENN(NA2,PIN,IIN,DLK,NC,K1,K2,0,2,ITT,RSI)                 
      IF(NC) 7,7,8                                                      
    7 IF(PM(10,K))39,39,40                                              
   40 TAUK=TRES                                                         
      IM(4,K)=1                                                         
      MYT(K)=0                                                          
      GOTO 33                                                           
   39 MYT(K)=0                                                          
      GO TO 15                                                          
C    8 IF(MK1T(K))88,89,89 
C   88 IF(K.EQ.M1.AND.K1.EQ.M2) GOTO 39                                
    8 D=((XYZ(2,1,K1)-XK0)*PIN(7)+(XYZ(2,2,K1)-YK0)*PIN(6))*            
     1PIN(4)+(XYZ(2,3,K1)-ZK0)*PIN(5)                                   
      DD=D                                                              
      TAUK=DD/VLI                                                       
   24 PIN(1)=XK0+DD*PIN(4)*PIN(7)                                       
      PIN(2)=YK0+DD*PIN(4)*PIN(6)                                       
      PIN(3)=ZK0+DD*PIN(5)                                              
      VMTI=POTENN(2,A(2),PIN,IIN,AR(2),CR(2),DR(2),TR(2),EPS(2),VPI)    
      VMPI=0.
      PIN(8)=PM(8,K)+VMTI                                    
      IF(MZ)50,50,51                                                    
   51 IF(IM(4,K)-1)50,52,50                                             
   52 DO 53 II=1,10                                                     
   53 PN(II)=PIN(II)                                                    
      DO 78 I=1,3                                                       
   78 IPN(I)=IIN(I)                                                     
      TAU=TAUK                                                          
      GOTO 54                                                           
   50 CALL PARTNN(2,K1,PN,IPN)                                          
      XR2=SQRT(XYZ(2,1,K1)**2+XYZ(2,2,K1)**2+XYZ(2,3,K1)**2)
      PR2=1.
      IF(XR2.LT.(R(2)-2.4)) PR2=PDL 
      PN(9)=.94                                                         
      CALL TINVU(PIN,PN,TIN1,V,U)                                       
      v2=SQRT(v(1)**2+v(2)**2+v(3)**2)
      IF(TIN1.EQ.0.) GOTO 6
      SAB=SIGMAN(IIN,IPN,3,PIN(8))                                      
   31 SIG=SIGMAN(IIN,IPN,0,TIN1)                                        
      IF(MZ)20,20,12                                                    
   20 SGH=SIG+SAB*Z(2)/A(2)
   92 IF(PM(17,K)) 18,18,19
   18 TAK=TAUK
      GOTO 32
   19 TAK=PM(17,K)+TAUK
   32 TKI=1.-PORG
      EPR=TAK/(TPAR*GM)
      IF(EPR.GT.50.) EPR=50.
      SGT=SGh*(1.-TKI*EXP(-EPR))*PR2
   94 TEM1=SGT/(31.4159*DLK**2) 
c      IF(RNDM(-1)-TEM1) 10,10,9           !UJ  2/07/2003
      IF(RSI-SQRT(SGT/31.4159))10,10,9           !UJ  2/07/2003
    9 MK1T(K)=0                                                         
      ITT=1                                                             
      GO TO 6                                                           
   10 MK1T(K)=K1                                                        
      MK2T(K)=K2                                                        
      PM(13,K)=TAUK                                                     
      SGM(K)=SGT/SGH
      IF(PM(10,K))44,44,42                                              
   42 IF(PM(14,K)-PM(13,K)) 43,43,44                                    
   43 IM(4,K)=1                                                         
      TAUK=TRES                                                         
      GOTO 33                                                           
   44 IM(4,K)=2                                                         
   33 IF(TAUK)15,15,112
  112 IF(TAU) 12,12,11                                                  
   11 IF(TAU-TAUK) 15,15,12                                             
   12 TAU=TAUK                                                          
      IF(MZ)27,27,28                                                    
   28 U1=U                                                              
      VMP=VMPI                                                          
      VMT=VMTI                                                          
      TR2=TIN1                                                          
   54 DO 13 L=1,10                                                      
      P1(L)=PIN(L)                                                      
   13 P2(L)=PN(L)          
      goto 139
      if(pin(9).lt..9)goto 139
C************ 11 APRIL 1999 ****
      rpt=sqrt(p1(1)**2+p1(2)**2+p1(3)**2)
      ratio=rpt/r(2)
      if(ratio.gt.1.) ratio=1.
      REX=SQRT(SIG/31.4159)
      if(rex-1.)140,140,141
 141  fex=.0001*(1.-ratio)/rex
      goto 142
 140  fex=.0001*(1.-ratio)
 142  P1(8)=PIN(8)-fex
      IF(P1(8))137,137,138
 137  P1(8)=PIN(8)
	GOTO 139
 138  PMD=SQRT(fex*(fex+1.88*AC(2)))
	PA2(1)=PMD*P1(4)*P1(6)
	PA2(2)=PMD*P1(4)*P1(7)
	PA2(3)=PMD*P1(5)
        CALL BLRECN(2,PA2(1),PA2(2),PA2(3),P1(1),P1(2),P1(3))
C*******************************
 139  P1(11)=SGM(K)
      P2(11)=pr2
      DO 14 L=1,3                                                       
      IP1(L)=IIN(L)                                                     
      V1(L)=V(L)                                                        
   14 IP2(L)=IPN(L)                                                     
   27 N1=K                                                              
      N2=K1                                                             
      N3=K2                                                             
      IF(MZ)15,15,17                                                    
   15 IF(K-MV) 16,17,17                                                 
   16 K=K+1                                                             
      GO TO 1                                                           
   17 RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE TIMEC(MZ,NA1,NA2,R0X,R0Y,R0Z,TINT,DEL,TAU,             
     *P1,P2,IP1,IP2,V1,U1,TR2,SIG,SAB,N1,N2,N3,UT,UP,PDL)                   
************************************************************************
*   PURPOSE : MZ=0 - CHOOSES THE INTERACTION OF PROJECTILE NUCLEON WITH*
*             TARGET NUCLEON WITH MINIMUM TIME IN THE CASE ON          *
*             NUCLEUS-NUCLEUS COLLISION.                               *
*             MZ=1 - CALCULATES KINEMATICAL CHARACTERISTICS OF         *
*             INTERACTING PARTNERS.                                    *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:
*     NA1 - CURRENT NUMBER OF NUCLEONS IN THE PROJ. NUCLEUS, NA1=AC(1).*    
*     NA2 - CURRENT NUMBER OF NUCLEONS IN THE TARGET NUCLEUS,NA2=AC(2).*    
*     R0X, - COORDINATES OF ENTRY POINT OF PROJECTILE NUCLEUS IN       *
*     R0Y, - THE REST FRAME OF                                         *
*     R0Z  - TARGET NUCLEUS.                                           *
*     TINT - CURRENT TIME OF COLLISION PROCESS.                        *      
*     DEL - PARAMETER DETERMINING A RADIUS OF STRONG INTERACTION.      *
*     TPR - PARAMETER SPECIFYING FORMATION TIME.                       *
*     PDL - PARAMETER SPECIFYING NUCLEAR TRANSPERANCY                  *
*----------------------------------------------------------------------*
*    OUTPUT QUANTITIES:                                                *
*     TAU     - MINIMUM TIME TO INTERACTION.                           *
*     P1,P2   - ARRAYS OF KINEMATICAL CHARACTERISTICS OF INTERACTING   *
*               PARTNERS.                                              *
*     IP1,IP  - ARRAYS OF THEIR CHARGES AND BARYON NUMBERS.            *
*     V1      - VELOCITY OF CMS OF INTERACTING PARTNERS OR             *
*               DECAYING RESONANCE.                                    *
*     U1      - TOTAL ENERGY OF INTERACTING PARTNERS IN CMS.           *
*     TR2     - KINETIC ENERGY OF INTERACTING PROJECTILE NUCLEON       *
*               IN REST FRAME OF TARGET NUCLEON.                       *
*     SIG     - TOTAL CROSS SECTION FOR THE INTERACTION.               *
*     SAB     - =0                                                     *
*     N1      - NO OF INTERACTING PROJECTILE NUCLEON IN ARRAYS XYZ     *
*               (COORNINATES) AND IZ (NUCLEON CHARGES).                *  
*     N2      - NO OF INTERACTING TARGET NUCLEON IN ARRAYS XYZ         *
*               (COORNINATES) AND IZ (NUCLEON CHARGES).                *  
*     N3      - NOT USED.                                              *
*     UT      - POTENTIAL FOR INTERACTING PROJECTILE NUCLEON INSIDE    *
*               TARGET NUCLEUS.                                        *
*     UP      - POTENTIAL FOR INTERACTING TARGET NUCLEON INSIDE        *
*               PROJECTILE NUCLEUS.                                    *  
************************************************************************      
      parameter (max1=2000)
      common/inel/inel              
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
     0COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     1MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      COMMON/ACTIV/ MPA(240),MYP(max1),MYT(max1)                        
      DIMENSION V0(3),PIN(11),IIN(3),P(11),IP(3),CS(3),PN(11),IPN(3)    
      DIMENSION P1(11),P2(11),IP1(3),IP2(3),V1(3),V(3),B(3),IT(3)       
      IF(MZ)22,23,22                                                    
   22 K=N1                                                              
      GOTO 34                                                           
   23 K=1                                                               
      TAU=-.1                                                           
    1 IF(MPA(K)) 16,16,2                                                
    2 IF(NK1P(K))34,34,35                                               
   35 TAUK=TIME(K)                                                      
      K1=NK1P(K)                                                        
      K2=0                                                              
      GOTO 33                                                           
   34 V0(1)=0.                                                          
      V0(2)=0.                                                          
      V0(3)=SQRT(WT*(WT+1.88D0))/(WT+.94D0)                             
      B(1)=0.                                                           
      B(2)=0.                                                           
      B(3)=-V0(3)                                                       
      XK0 =XYZ(1,1,K)+R0X                                               
      YK0 =XYZ(1,2,K)+R0Y                                               
      ZK0 =XYZ(1,3,K)*(.94/(WT+.94))+R0Z+V0(3)*TINT                     
      XR1=SQRT(XYZ(1,1,K)**2+XYZ(1,2,K)**2+XYZ(1,3,K)**2)
      PR1=1.
      IF(XR1.LT.(R(1)-2.4)) PR1=PDL
      CALL PARTNN(1,K,P,IP)                                             
      CMS=SQRT(P(8)*(P(8)+2.*P(9)))                                     
      CS(1)=CMS*P(4)*P(7)                                               
      CS(2)=CMS*P(4)*P(6)                                               
      CS(3)=CMS*P(5)                                                    
      PIN(9)=P(9)                                                       
      CALL CINEMA(CS,V0,PIN)                                            
      PIN(1)=XK0                                                        
      PIN(2)=YK0                                                        
      PIN(3)=ZK0                                                        
      TL=PIN(8)                                                         
      PMI=P(9)                                                          
      IIN(1)=IP(1)                                                      
      UPK=.94-P(9)                                                      
      IIN(2)=1                                                          
      IIN(3)=0                                                          
      ITT=0                                                             
      IF(MZ) 21,21,3                                                    
    3 K1=NK1P(K)                                                        
      K2=0                                                              
      GO TO 7                                                           
   21 DL1=1./(5.07*SQRT(WT*(WT+2.*1.88)))+DEL                           
    5 CALL BLCENN(NA2,PIN,IIN,DL1,NC,K1,K2,0,2,ITT,RSI)                 
      IF(NC) 6,6,7                                                      
    6 MPA(K)=0                                                          
      IPP=1                                                             
      GO TO 16                                                          
    7 PIN(3)=XYZ(2,3,K1)                                                
      DZ=PIN(3)-ZK0                                                     
      TAUK=DZ/V0(3)                                                     
      TIME(K)=TAUK                                                      
   24 VTK=POTENN(2,A(2),PIN,IIN,AR(2),CR(2),DR(2),TR(2),EPS(2),VPI)     
      PIN(8)=TL+VTK                                                     
      PIN(9)=.94                                                        
      PIN(10)=0.                                                        
      XR2=SQRT(XYZ(2,1,K1)**2+XYZ(2,2,K1)**2+XYZ(2,3,K1)**2)
      PR2=1.
      IF(XR2.LT.(R(2)-2.4)) PR2=PDL
      CALL PARTNN(2,K1,PN,IPN)                                          
      UTK=.94-PN(9)                                                     
      DO 8 I=1,10                                                       
    8 P(I)=PN(I)                                                        
      DO 9 I=1,3                                                        
    9 IP(I)=IPN(I)                                                      
      CMS=SQRT(P(8)*(P(8)+2.*P(9)))                                     
      CS(1)=CMS*P(4)*P(7)                                               
      CS(2)=CMS*P(4)*P(6)                                               
      CS(3)=CMS*P(5)                                                    
      CALL CINEMA(CS,B,PN)                                              
      TLN=PN(8)                                                         
      P(1)=P(1)-R0X                                                     
      P(2)=P(2)-R0Y                                                     
      P(3)=(P(3)-V0(3)*(TINT+TAUK)-R0Z)*(WT+.94)/.94                    
      VPK=POTENN(1,A(1),P,IP,AR(1),CR(1),DR(1),TR(1),EPS(1),VPI)        
      P(8)=TLN+VPK                                                      
      CMS=SQRT(P(8)*(P(8)+2.*P(9)))                                     
      CS(1)=CMS*PN(4)*PN(7)                                             
      CS(2)=CMS*PN(4)*PN(6)                                             
      CS(3)=CMS*PN(5)                                                   
      IF(PN(5).EQ.1.) CS(3)=-CS(3)
      PN(9)=P(9)                                                        
      CALL CINEMA(CS,V0,PN)                                             
      DO 91 I=1,3                                                       
   91 PN(I)=XYZ(2,I,K1)                                                 
      PN(9)=.94                                                         
      CALL TINVU(PIN,PN,TIN1,V,U)                                       
      IF(TIN1.LE.0.) GOTO 16
      SIG=SIGMAN(IIN,IPN,0,TIN1)                                        
      IF(SIG.EQ.0.) GOTO 16                                             
      IF(MZ)20,20,13                                                    
   20 PR12=PR1*PR2
      SG=SIG*PR12
      TEMP=.1*SG/(3.14159*DL1**2)                                       
c      IF(RNDM(-1)-TEMP) 11,11,10                                        
      IF(RSI-SQRT(SG/31.4159)) 11,11,10
   10 NK1P(K)=0                                                         
      ITT=1                                                             
      GO TO 5                                                           
   11 NK1P(K)=K1                                                        
      TIME(K)=TAUK                                                      
   33 IF(TAUK) 16,16,112
  112 IF(TAU) 13,13,12                                                  
   12 IF(TAU-TAUK) 16,16,13                                             
   13 TAU=TAUK                                                          
      IF(MZ)27,27,28                                                    
   28 DO 14 L=1,10                                                      
      P1(L)=PIN(L)                                                      
   14 P2(L)=PN(L)                                                       
      P1(11)=PR1
      P2(11)=PR2
      DO 15 L=1,3                                                       
      IP1(L)=IIN(L)                                                     
      IP2(L)=IPN(L)                                                     
   15 V1(L)=V(L)                                                        
      U1=U                                                              
      UT=UTK                                                            
      UP=UPK                                                            
      SAB=0.                                                            
      TR2=TIN1                                                          
   27 N1=K                                                              
      N2=K1                                                             
      N3=K2                                                             
      IF(MZ)16,16,18                                                    
   16 IF(K-NA1) 17,18,18                                                
   17 K=K+1                                                             
      GO TO 1                                                           
   18 RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE TIMEB(MZ,MV,NA1,R0X,R0Y,R0Z,TINT,DEL,TAU,              
     *P1,P2,IP1,IP2,V1,U1,TR2,SIG,SAB,N1,N2,N3,VMT,VMP,VLI,TPR,PDL)    
************************************************************************
*   PURPOSE : MZ=0 - CHOOSES THE INTERACTION OF CASCADE PARTICLE WITH  *
*             PROJECTILE NUCLEON OR DECAYING RESONANCE WITH MINIMUM    *
*             TIME;                                                    *
*             MZ=1 - CALCULATES KINEMATICAL CHARACTERISTICS OF         *
*             INTERACTING PARTNERS OR DECAYING RESONANCE.              *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     MZ  - LABEL INDICATING FIRST AND SECOND RUN OF SUBROUTINE:       *
*           IN THE FIRST RUN MZ=0 - SUBROUTINE CHOOSES PARTNERS        * 
*           OR DECAYING RESONANCE WITH MINIMUM TIME;                   *
*           IN THE SECOND RUN MZ=1, - SUBROUTINE CALCULATES KINEMATICAL*
*           CHARACTERISTICS OF PARTNERS OR DECAYING RESONANCE.         *
*     MV  - CURRENT NUMBER OF CASCADE PARTICLES IN THE ARRAYS PM AND IM*
*     NA1 - CURRENT NUMBER OF NUCLEONS IN POJECTILE NUCLEUS, NA1=AC(1).*
*     R0X, - COORDINATES OF ENTRY POINT OF PROJECTILE NUCLEUS IN       *
*     R0Y, - THE REST FRAME OF                                         *
*     R0Z  - TARGET NUCLEUS.                                           *
*     TINT - CURRENT TIME OF COLLISION PROCESS.                        *    
*     DEL - PARAMETER DETERMINING A RADIUS OF STRONG INTERACTION.      *
*     TPR - PARAMETER SPECIFYING DELAY TIME.                           *
*     PDL  - PARAMETER SPECIFYING NUCLEAR TRANSPERANCY                 *
*----------------------------------------------------------------------*
*    OUTPUT QUANTITIES:                                                *
*     TAU     - MINIMUM TIME TO INTERACTION.                           *
*     P1,P2   - ARRAYS OF KINEMATICAL CHARACTERISTICS OF INTERACTING   *
*               PARTNERS OR DECAYING RESONANCE WITH MINIMUM TIME.      *
*     IP1,IP  - ARRAYS OF THEIR CHARGES AND BARYON NUMBERS.            *
*     V1      - VELOCITY OF CMS OF INTERACTING PARTNERS OR             *
*               DECAYING RESONANCE.                                    *
*     U1      - TOTAL ENERGY OF INTERACTING PARTNERS IN CMS.           *
*     TR2     - KINETIC ENERGY OF INTERACTING CASCADE PARTICLE         *
*               IN REST FRAME OF PROJECTILE NUCLEON.                   *
*     SIG     - TOTAL CROSS SECTION FOR THE INTERACTION.               *
*     SAB     - CROSS SECTION FOR ABSORPTION OF PION ON QUASIDEUTERON  *
*               PAIR.                                                  *
*     N1      - NO OF INTERACTING CASCADE PARTNER IN ARRAYS PM AND IM. *
*     N2      - NO OF INTERACTING PROJECTILE NUCLEON IN ARRAYS XYZ     *
*               (COORNINATES) AND IZ (NUCLEON CHARGES).                *
*     N3      - NO OF NUCLEON-PARTHER FOR PION ABSORPSION IN ARRAYS    *
*               XYZ AND IZ.                                            *
*     VMT     - NOT USED.                                              *  
*     VMP     - POTENTIAL FOR INTERACTING CASCADE PARTICLE.            *  
*     VLI     - VELOCITY OF DECAYING RESONANCE.                        *
************************************************************************     
      parameter (max1=2000)
      common/inel/inel              
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
      COMMON/PMIM/ pm(19,max1),IM(5,max1)                               
     0COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     1MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      COMMON /ACTIV/ MPA(240),MYP(max1),MYT(max1)                       
      DIMENSION V0(3),PIN(11),IIN(3),P(11),IP(3),PN(11),IPN(3)          
      DIMENSION P1(11),P2(11),IP1(3),IP2(3),V1(3),V(3),IT(3),CS(3)      
      DIMENSION KP(2),KT(2),SGM(MAX1),PA1(3)
      SAVE SGM                                                          
      V0(1)=0.                                                          
      V0(2)=0.                                                          
      V0(3)=-SQRT(WT*(WT+1.88D0))/(WT+.94D0)
      TPAR=TPR
      IF(MZ)22,23,22                                                    
   22 K=N1                                                              
      GOTO 34                                                           
   23 K=1                                                               
      TAU=-.1                                                           
    1 IF(MYP(K))15,15,2                                                 
C   70 IF(PM(10,K)) 15,15,71
   71 IF(PM(12,K))34,34,65                                              
    2 IF(PM(10,K))36,36,63
c*******correction 31 July 1999 ***********                                   
   36 IF(MK1P(K))34,34,644
 644  if(pm(9,k)-.9) 645,64,64
 645  if(mk2p(k))34,34,64
c******************************************
   63 IF(IM(3,K)-1)34,71,36                                             
   65 TAUK=PM(12,K)                                                     
c********* correction 31 July 1999 ********
      if(mk1p(k))34,34,66
c******************************************
   64 TAUK=PM(11,K)                                                     
   66 K1=MK1P(K)                                                        
      K2=MK2P(K)                                                        
      GOTO 33                                                           
   34 XK0=PM(1,K)-R0X                                                   
      YK0=PM(2,K)-R0Y                                                   
      ZK0=(PM(3,K)+V0(3)*TINT-R0Z)*(WT+.94)/.94                         
      DO 3 I=1,10                                                       
    3 P(I)=PM(I,K)                                                      
      IP(1)=IM(1,K)                                                     
      IP(2)=IM(2,K)                                                     
      CSM=SQRT(P(8)*(P(8)+2.*P(9)))                                     
      CS(1)=CSM*PM(4,K)*PM(7,K)                                         
      CS(2)=CSM*PM(4,K)*PM(6,K)                                         
      CS(3)=CSM*PM(5,K)                                                 
      PIN(9)=PM(9,K)                                                    
      PIN(10)=PM(10,K)                                                  
      PORG=PM(19,K)
      CALL CINEMA (CS,V0,PIN)                                           
      IF(CSM.EQ.0.) PIN(5)=-1.                                          
      IIN(1)=IM(1,K)                                                    
      IIN(2)=IM(2,K)                                                    
      IIN(3)=1                                                          
      ITT=0            
      IF(PIN(8).LT.TR(1)/2..AND.IIN(2).GT.0) GOTO 15
      PIN(1)=XK0                                                        
      PIN(2)=YK0                                                        
      PIN(3)=ZK0                                                        
      VMPI=POTENN(1,A(1),PIN,IIN,AR(1),CR(1),DR(1),TR(1),EPS(1),VPI)    
      TL=PIN(8)                                                         
C      PIN(8)=TL+VMPI
      IF(TL.LE.0.) GO TO 39                          
      PMO=SQRT(PIN(8)*(PIN(8)+2.D0*PIN(9)))
      EMO=PIN(8)+PIN(9)
      VLI=PMO/EMO
      GM=EMO/PIN(9)
      TCS=(1.D0-PIN(5)*VLI*V0(3))*(WT+.94D0)/.94D0                      
      IF(MZ)21,21,4                                                     
   4  K1=MK1P(K)                                                        
      K2=MK2P(K)                                                        
      IF(IM(3,K)-1)47,47,8                                              
   47 TAUK=PM(12,K)
      TK=TAUK/TCS
      DD=TK*VLI                                               
      V(1)=VLI*PIN(4)*PIN(7)                                            
      V(2)=VLI*PIN(4)*PIN(6)                                            
      V(3)=VLI*PIN(5)                                                   
      GOTO 24                                                           
   21 IF(PM(10,K)) 611,611,72
c****** correction 31 July 1999 ***********
 611  if(pin(8)-0.015) 15,15,61
   72 IF(IM(3,K)-1)37,622,77
 622  if(mk1p(k))77,77,62
c******************************************
   62 K1=MK1P(K)                                                        
      K2=MK2P(K)                                                        
   37 CALL PLAM(PIN,DS,TRE,TPR)                                         
      TRES=TRE*TCS
      PM(12,K)=TRES
      IF(MYP(K))43,43,73                                                
   73 IF(IM(3,K)-1)61,74,61                                             
   74 TD=PM(11,K)
      IF(TRES-TD)75,75,76                                     
   75 TAUK=TRES                                                     
      IM(3,K)=1                                                         
      GOTO 33                                                           
   76 TAUK=TD                                                     
      IM(3,K)=2                                                         
      GOTO 33                                                           
   77 TRES=PM(12,K) 
      DS=TRES*VLI                                                       
   61 CONTINUE
   67 DLK=1./(5.07*SQRT(PIN(8)*(PIN(8)+2.*PIN(9))))+DEL                 
    6 CALL BLCENN(NA1,PIN,IIN,DLK,NC,K1,K2,0,1,ITT,RSI)                 
      IF(NC) 7,7,8                                                      
    7 IF(PM(10,K))39,39,40                                              
   40 TAUK=TRES                                                         
      IM(3,K)=1                                                         
      MYP(K)=0                                                          
      GOTO 33                                                           
   39 MYP(K)=0                                                          
      GO TO 15                                                          
    8 D=((XYZ(1,1,K1)-XK0)*PIN(7)+(XYZ(1,2,K1)-YK0)*PIN(6))*            
     1PIN(4)+(XYZ(1,3,K1)-ZK0)*PIN(5)                                   
      DD=D                                                              
      TAUK=(DD/VLI)*TCS
   24 PIN(1)=XK0+DD*PIN(4)*PIN(7)                                       
      PIN(2)=YK0+DD*PIN(4)*PIN(6)                                       
      PIN(3)=ZK0+DD*PIN(5)                                              
      VMPI=POTENN(1,A(1),PIN,IIN,AR(1),CR(1),DR(1),TR(1),EPS(1),VPI)    
      VMTI=0.    
      PIN(8)=TL+VMPI
      IF(MZ)50,50,51                                                    
   51 IF(IM(3,K)-1)50,52,50                                             
   52 DO 53 II=1,10                                                     
   53 PN(II)=PIN(II)                                                    
      TAU=TAUK                                                          
      DO 78 I=1,3                                                       
   78 IPN(I)=IIN(I)                                                     
      GOTO 54                                                           
   50 CALL PARTNN(1,K1,PN,IPN)                                          
      XR1=SQRT(XYZ(1,1,K1)**2+XYZ(1,2,K1)**2+XYZ(1,3,K1)**2)
      PR1=1.
      IF(XR1.LT.(R(1)-2.4)) PR1=PDL 
      PN(9)=0.94                                                        
      CALL TINVU(PIN,PN,TIN1,V,U)                                       
      v2=SQRT(v(1)**2+v(2)**2+v(3)**2)
      IF(TIN1.EQ.0.) GOTO 6
      SAB=SIGMAN(IIN,IPN,3,TIN1)                                        
   31 SIG=SIGMAN(IIN,IPN,0,TIN1)                                        
      IF(MZ)20,20,12                                                    
   20 SGH=SIG+SAB*Z(2)/A(1)
   92 IF(PM(16,K)) 18,18,19
   18 TAK=TAUK
      GOTO 32
   19 TAK=(PM(16,K)+TAUK)/TCS
   32 TKI=1.-PORG
      EPR=TAK/(TPAR*GM)
      IF(EPR.GT.50.) EPR=50.
      SGT=SGh*(1.-TKI*EXP(-EPR))*PR1
   94 TEM1=SGT/(31.4159*DLK**2)
c      IF(RNDM(-1)-TEM1)10,10,9
      IF(RSI-SQRT(SGT/31.4159))10,10,9
    9 MK1P(K)=0                                                         
      ITT=1                                                             
      GO TO 6                                                           
   10 MK1P(K)=K1                                                        
      MK2P(K)=K2                                                        
      PM(11,K)=TAUK
      SGM(K)=SGT/SGH
      IF(PM(10,K)) 44,44,42                                             
   42 IF(PM(12,K)-PM(11,K))43,43,44
   43 IM(3,K)=1                                                         
      TAUK=TRES                                                         
      GOTO 33                                                           
   44 IM(3,K)=2                                                         
   33 IF(TAUK) 15,15,112
  112 IF(TAU) 12,12,11                                                  
   11 IF(TAU-TAUK) 15,15,12                                             
   12 TAU=TAUK                                                          
      IF(MZ)27,27,28                                                    
   28 U1=U                                                              
      VMT=VMTI                                                          
      VMP=VMPI                                                          
      TR2=TIN1                                                          
   54 DO 13 L=1,10                                                      
      P1(L)=PIN(L)                                                      
   13 P2(L)=PN(L)
      goto 139
      if(pin(9).lt..9)goto 139
C************ 11 APRIL 1999 ****
      rpt=sqrt(p1(1)**2+p1(2)**2+p1(3)**2)
      ratio=rpt/r(1)
      if(ratio.gt.1.) ratio=1.
      if(rex-1.)140,140,141
 141  fex=.0001*(1.-ratio)/rex
      goto 142
 140  fex=.0001*(1.-ratio)
 142  P1(8)=PIN(8)-fex
      IF(P1(8))137,137,138
 137  P1(8)=PIN(8)
	GOTO 139
 138  PMD=SQRT(fex*(fex+1.88*AC(2)))
	PA1(1)=PMD*P1(4)*P1(6)
	PA1(2)=PMD*P1(4)*P1(7)
	PA1(3)=PMD*P1(5)
        CALL BLRECN(1,PA1(1),PA1(2),PA1(3),P1(1),P1(2),P1(3))
C*******************************
 139  P1(11)=SGM(K)
      P2(11)=pr1
      DO 14 L=1,3                                                       
      IP1(L)=IIN(L)                                                     
      IP2(L)=IPN(L)                                                     
   14 V1(L)=V(L)                                                        
   27 N1=K                                                              
      N2=K1                                                             
      N3=K2                                                             
      IF(MZ)15,15,17                                                    
   15 IF(K-MV) 16,17,17                                                 
   16 K=K+1                                                             
      GO TO 1                                                           
   17 RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE TIMED(MZ,MV,DEL,TAU,PIN,PN,IIN,IPN,V,U,TIN1,SIG,N1,N2, 
     * TPR)
************************************************************************
*   PURPOSE : MZ=0 - CHOOSES THE CASCADE-CASCADE INTERACTION           *
*             WITH MINIMUM TIME.                                       *
*             MZ=1 - CALCULATES KINEMATICAL CHARACTERISTICS OF         *
*             INTERACTING PARTNERS.                                    *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     MV  - CURRENT NUMBER OF CASCADE PARTICLES IN THE ARRAYS PM AND IM*
*     DEL - PARAMETER DETERMINING A RADIUS OF STRONG INTERACTION.      *
*     TPR - PARAMETER SPECIFYING DELAY TIME.                           *
*----------------------------------------------------------------------*
*    OUTPUT QUANTITIES:                                                *
*     TAU     - MINIMUM TIME UP TO INTERACTION.                        *
*     PIN,PN   - ARRAYS OF KINEMATICAL CHARACTERISTICS OF INTERACTING  *
*               PARTNERS WITH MINIMUM TIME.                            *
*     IIN,IPN  - ARRAYS OF THEIR CHARGES AND BARYON NUMBERS.           *
*     V       - VELOCITY OF CMS OF INTERACTING PARTNERS IN LAB.        *
*     U       - TOTAL ENERGY OF INTERACTING PARTNERS IN CMS.           *
*     TIN1     - KINETIC ENERGY OF THE FIRST CASCADE PARTICLE          *
*               IN REST FRAME OF THE SECOND ONE.                       *
*     SIG     - TOTAL CROSS SECTION FOR THE INTERACTION.               *
*     N1      - NO OF THE FIRST CASCADE PARTNER IN ARRAYS PM AND IM.   *
*     N2      - NO OF THE SECOND CASCADE PARTNER IN ARRAYS PM AND IM.  *
************************************************************************
      parameter (max1=2000)
      common/inel/inel              
      COMMON/KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KS1,LST
      COMMON/PMIM/ pm(19,max1),IM(5,max1)                               
      COMMON/CARCIL/ NK1P(240),TIME(240),MK1P(max1),MK2P(max1),         
     1 MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                       
      COMMON /ACTIV/ MPA(240),MYP(max1),MYT(max1)                       
      DIMENSION V0(3),PIN(11),IIN(3),PN(11),IPN(3),R1(3),R2(3),         
     1 PI(3),PR(3),PS(11),V(3),PL(3),V1(3),V2(3),SGM(MAX1)              
      JIM=0                                                             
      MI=0                                                              
      IF(MZ) 22,23,22                                                   
   22 K=N1                                                              
      K1=MKAS(K)                                                        
      GOTO 19                                                           
 23   KK=1                                                              
      TAU=-.1                                                           
      TU=-.1                                                            
      M=MV-1
    1 CONTINUE                                                          
      K=KK                                                              
  120 I=1                                                               
  113 K1=K+I                                                            
      MOLD=0                                                            
      MI=0                                                              
      MM=0                                                              
      IF(PM(15,K).LE.0..AND.PM(15,K1).LE.0.) GO TO 3                    
    2 IF(MKAS(K).EQ.K1.AND.MKAS(K1).EQ.K) GOTO 33                       
      GO TO 110                                                         
    3 IF(KAS(K).EQ.KAS(K1)) GOTO 110                                    
      IF(KAS(K).NE.LST.AND.KAS(K1).NE.LST) GOTO 110
  771 CONTINUE
c      IF(PM(19,K).LT..1.AND.PM(19,K1).LT..1) GOTO 110
C      IF(PM(18,K).EQ.PM(18,K1)) GOTO 110
   19 IF(PM(9,K).LE.PM(9,K1).OR.MZ.EQ.0) GO TO 8                        
      DO 43 L=1,10                                                      
      PIN(L)=PM(L,K1)                                                   
   43 PN(L)=PM(L,K)                                                     
      PR1=PM(19,K1)
      PR2=PM(19,K)
      TPIN=PM(17,K1)
      TPN=PM(17,K)
      DO 5 L=1,2                                                        
      IIN(L)=IM(L,K1)                                                   
    5 IPN(L)=IM(L,K)                                                    
      MI=1                                                              
      GOTO 99                                                           
    8 DO 90 L=1,10                                                      
      PIN(L)=PM(L,K)                                                    
   90 PN(L)=PM(L,K1)                                                    
      PR1=PM(19,K)
      PR2=PM(19,K1)
      TPIN=PM(17,K)
      TPN=PM(17,K1) 
      DO 91 L=1,2                                                       
      IIN(L)=IM(L,K)                                                    
   91 IPN(L)=IM(L,K1)                                                   
   99 IIN(3)=1                                                          
      IPN(3)=1                                                          
      CALL TINVU(PIN,PN,TIN1,V,U)                                       
      IF(TIN1.LT..250) GOTO 110
      KL=K                                                              
      K1L=K1                                                            
      IF(MZ)121,121,122                                                 
  122 VK=SQRT(PIN(8)*(PIN(8)+2.D0*PIN(9)))/(PIN(8)+PIN(9))              
      VK1=SQRT(PN(8)*(PN(8)+2.D0*PN(9)))/(PN(8)+PN(9))                  
      TAUK=PM(15,K)                                                     
      DK=TAUK*VK                                                        
      DK1=TAUK*VK1                                                      
      PIN(1)=PIN(1)+DK*PIN(4)*PIN(7)                                    
      PIN(2)=PIN(2)+DK*PIN(4)*PIN(6)                                    
      PIN(3)=PIN(3)+DK*PIN(5)                                           
      PN(1)=PN(1)+DK1*PN(4)*PN(7)                                       
      PN(2)=PN(2)+DK1*PN(4)*PN(6)                                       
      PN(3)=PN(3)+DK1*PN(5)                                             
      IF(MI) 82,82,83
   82 KF=K
      KB=K1
      GO TO 84
   83 KF=K1
      KB=K
   84 TAK=PM(17,KF)+TAUK
       GPIN=(PM(9,KF)+PM(8,KF))/PM(9,KF)
       EPR=TAK/(TPR*GPIN)
       PIN(11)=1.-(1.-PR1)*EXP(-EPR)
      TAK=PM(17,KB)+TAUK
       GPN=(PM(9,KB)+PM(8,KB))/PM(9,KB)
       EPR=TAK/(TPR*GPN)
       PN(11)=1.-(1.-PR2)*EXP(-EPR)
      GOTO 10                                                           
  121 DIST=SQRT((PIN(1)-PN(1))**2+(PIN(2)-PN(2))**2)
c      IF(DIST.GT.4.0) GOTO 110
      EPI=PIN(8)+PIN(9)                                                 
      EPN=PN(8)+PN(9)                                                   
      PPN=SQRT(PN(8)*(PN(8)+2.D0*PN(9)))                                
      PPI=SQRT(PIN(8)*(PIN(8)+2.D0*PIN(9)))                             
      GPIN=EPI/PIN(9)
      GPN=EPN/PN(9)
      V0(1)=0.D0                                                        
      V0(2)=0.D0                                                        
      V0(3)=-PPN/EPN                                                    
      PL(1)=PPN*PN(4)*PN(7)                                             
      PL(2)=PPN*PN(4)*PN(6)                                             
      PL(3)=PPN*PN(5)                                                   
      PI(1)=PPI*PIN(4)*PIN(7)                                           
      PI(2)=PPI*PIN(4)*PIN(6)                                           
      PI(3)=PPI*PIN(5)                                                  
      V2(1)=PL(1)/EPN
      V2(2)=PL(2)/EPN
      V2(3)=PL(3)/EPN
      DIREC=(PN(1)-PIN(1))*(PL(1)-PI(1))+(PN(2)-PIN(2))*(PL(2)-PI(2))   
     *+(PN(3)-PIN(3))*(PL(3)-PI(3))                                     
      IF(DIREC.GT.0.) GO TO 110                                         
      RI=PIN(1)*PN(7)+PIN(2)*PN(6)                                      
      R1(1)=RI*PN(5)-PIN(3)*PN(4)                                       
      R1(2)=-PIN(1)*PN(6)+PIN(2)*PN(7)                                  
      R1(3)=RI*PN(4)+PIN(3)*PN(5)                                       
      RP=PN(1)*PN(7)+PN(2)*PN(6)                                        
      R2(1)=RP*PN(5)-PN(3)*PN(4)                                        
      R2(2)=-PN(1)*PN(6)+PN(2)*PN(7)                                    
      R2(3)=RP*PN(4)+PN(3)*PN(5)                                        
      XK0=R1(1)-R2(1)                                                   
      YK0=R1(2)-R2(2)                                                   
      ZK0=(R1(3)-R2(3))*EPN/PN(9)                                       
      PRI=PI(1)*PN(7)+PI(2)*PN(6)                                       
      PR(1)=PRI*PN(5)-PI(3)*PN(4)                                       
      PR(2)=-PI(1)*PN(6)+PI(2)*PN(7)                                    
      PR(3)=PRI*PN(4)+PI(3)*PN(5)                                       
      PS(9)=PIN(9)                                                      
      CALL CINEMA(PR,V0,PS)                                             
      IF(PS(8).LE.0.02) GO TO 110                                        
C     COORDINATES OF INITIAL PARTICLE IN THE PARTNER FRAME              
      XY=XK0*PS(7)+YK0*PS(6)                                            
      X1=XY*PS(5)-ZK0*PS(4)                                             
      Y1=-XK0*PS(6)+YK0*PS(7)                                           
      Z1=XY*PS(4)+ZK0*PS(5)                                             
      IF(Z1)11,110,110                                                  
   11 PMS=SQRT(PS(8)*(PS(8)+2.D0*PS(9)))                                
      RSI=SQRT(X1**2+Y1**2)                                             
      IF(RSI-DEL) 12,110,110
   12 ES=PS(8)+PS(9)                                                    
      VS=PMS/ ES                                                        
      V1(1)=VS*PS(4)*PS(7)
      V1(2)=VS*PS(4)*PS(6)
      V1(3)=VS*PS(5)
      G1=1./SQRT(1.D0-VS**2)
      TUK=-Z1/VS
      TAUK=TUK*(EPN/PN(9))*(1.+V1(1)*V2(1)+V1(2)*V2(2)+V1(3)*V2(3))     
   10 IF(PM(9,K1)-PM(9,K))18,18,130                                     
   18 TPI=(U**2-(PN(9)+PIN(9))**2)/(2.*PIN(9))                          
      SIG=SIGMAN(IPN,IIN,0,TPI)                                         
      GO TO 9                                                           
  130 CONTINUE                                                          
      SIG=SIGMAN(IIN,IPN,0,TIN1)                                        
    9 CONTINUE                                                          
      IF(MZ) 20,20,52                                                   
   20 IF(TPIN) 101,101,102
  101 TAK1=TAUK
      GOTO 103
  102 TAK1=TPIN+TAUK
  103 IF(TPN) 104,104,105
  104 TAK2=TAUK
      GOTO 106
  105 TAK2=TPN+TAUK
  106 TK1=(1.-PR1)*EXP(-TAK1/(TPR*GPIN))
      TK2=(1.-PR2)*EXP(-TAK2/(TPR*GPN))           
      TK12=TK1*TK2
      SGT=SiG*(1.-TK12)
      IF(IIN(2).EQ.0.AND.IPN(2).EQ.0) SGT=.73*SGT
      RCIL=SQRT(.1*SGT/3.14159)
      IF(RSI-RCIL) 17,17,110                                            
 33   TAUK=PM(15,K)
 17   IF(TU) 13,13,14
 13   IF(TAUK) 110,110,42
 14   IF(TU-TAUK) 110,110,42
 42   TU=TAUK
      PM(15,K)=TAUK
      PM(15,K1)=TAUK
      MKAS(K)=K1
      MKAS(K1)=K
      JIM=JIM+1
  110 IF(K1.GE.Mv) GOTO 117
  111 I=I+1                                                             
      GOTO 113                                                          
   52 IF(MI)46,46,119                                                   
   46 N1=KL                                                             
      N2=K1L                                                            
      GOTO 117                                                          
  119 N1=K1L                                                            
      N2=KL                                                             
  117 IF(MZ) 53,53,48                                                   
 53   IF(KK-M)15,48,48                                                  
 15   KK=KK+1                                                           
      IF(TAU) 36,36,37
 36   IF(TU) 1,1,41
 37   IF(TU) 1,1,38
 38   IF(TAU-TU) 39,39,41
 41   TAU=TU
      N1=K
      N2=K1
 39   TU=-.1
      GOTO 1                                                            
 48   RETURN                                                            
      END                                                               
************************************************************************
      FUNCTION SIGMAN(IP,IC,I,T)                                        
************************************************************************
*    PURPOSE : CALCULATES CROSS SECTIONS OF THE INTERACTION            *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     IP(1) - CHARGE OF THE FIRST PARTNER.                             *
*     IP(2) - BARYON NUMBER OF FIRST PARTNER.                          *
*     IP(3) - NOT USED.                                                *
*     IC(1) - CHARGE OF THE SECOND PARTNER.                            *
*     IC(2) - BARYON NUMBER OF SECOND PARTNER.                         *
*     IC(3) - NOT USED.                                                *
*     T     - KINETIC ENERGY OF THE FIRST PARTNER IN REST FRAME OF     *
*             THE SECOND ONE.                                          *
*     I     - LABEL FOR TYPE OF CROSS SECTION AND REACTION.            *
************************************************************************    
      DIMENSION IS(47),IP(3),IC(3),V(3),P(10),C(10)                     
      DATA IS/210,211,220,221,120,121,122,110,111,123,214,215,224,225,  
     1114,115,126,124,125,510,511,410,411,514,515,516,112,113,116,117,  
     2127,130,131,132,133,134,135,136,137,212,213,216,217,222,223,226,  
     3227 /                                                             
      M=IP(2)+IC(2)                                                     
      IF(M.EQ.0) M=1
      N=IP(1)+IC(1)                                                     
      IF(M-1) 1,1,4                                                     
    1 IF(IP(1)) 2,8,3                                                   
    2 IF(IC(1)) 6,6,7                                                   
    3 IF(IC(1)) 7,7,6                                                   
    4 IF(M-2) 5,5,6                                                     
    5 IF(N-1)6,7,6                                                      
    6 K=1                                                               
      GO TO 9                                                           
    7 K=2                                                               
      GO TO 9                                                           
    8 K=3                                                               
    9 IK=100*M+10*K+I                                                   
      DO 11 J=1,26                                                      
      IF(IK-IS(J)) 11,10,11                                             
   10 SIGMAN=APSIG(T,J)                                                 
      RETURN                                                            
   11 CONTINUE                                                          
      DO 27 J=1,21                                                      
      IF(IK-IS(26+J)) 27,12,27                                          
   12 KN=J                                                              
     0GO TO (13,14,13,15,16,17,18,19,14,20,21,22,23,13,13,13,24,13,13,  
     125,26), KN                                                        
  13  SIGMAN=0.                                                         
      RETURN                                                            
   14 SIGMAN=APSIG(T,10)                                                
      RETURN                                                            
   15 SIGMAN=APSIG(T,15)+APSIG(T,16)                                    
      RETURN                                                            
   16 SIGMAN=APSIG(T,17)+APSIG(T,18)+APSIG(T,19)                        
      RETURN                                                            
   17 SIGMAN=(APSIG(T,5)+APSIG(T,8))/2.                                 
      RETURN                                                            
 18   SIGMAN=(APSIG(T,6)-APSIG(T,7)+APSIG(T,9))/2.                      
      RETURN                                                            
   19 SIGMAN=APSIG(T,7)                                                 
      RETURN                                                            
   20 SIGMAN=(APSIG(T,15)+APSIG(T,18))/2.                               
      RETURN                                                            
   21 SIGMAN=(APSIG(T,16)+APSIG(T,19))/2.                               
      RETURN                                                            
   22 SIGMAN=APSIG(T,17)/2.                                             
      RETURN                                                            
   23 SIGMAN=(APSIG(T,15)+APSIG(T,16)+APSIG(T,17)+APSIG(T,18)           
     1+APSIG(T,19))/2.                                                  
      RETURN                                                            
   24 SIGMAN=APSIG(T,11)+APSIG(T,12)                                    
      RETURN                                                            
   25 SIGMAN=APSIG(T,13)                                                
      RETURN                                                            
   26 SIGMAN=2.*APSIG(T,13)+APSIG(T,14)                                 
      RETURN                                                            
   27 CONTINUE                                                          
      SIGMAN=0.                                                         
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE BLELN(V,U,T,P,IP,C,IC,MV,NP)                           
************************************************************************
*     PURPOSE : CALCULATES OF CHARACTERISTICS OF SECONDARIES           *
*        IN ELASTIC AND CHARGE EXCHANGE REACTIONS                      *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     V  - VELOCITY OF CMS OF THE INTERACTING PARTNERS IN THE LAB. 
*          SYSTEM.                                                     *
*     U     - TOTAL ENERGY OF THE INERACTING PARTNERS IN THE CMS.      *
*     T     - KINETIC ENERGY OF THE FIRST PARTNER IN THE REST FRAME    *
*           - OF THE SECOND PARTNER.                                   *
*     P ,IP - KINEMATICAL CHARACTERISTICS & QUANTUM NUMBERS OF THE     *
*             FIRST PARTNER.                                           *
*     C ,IC - KINEMATICAL CHARACTERISTICS & QUANTUM NUMBERS OF THE     *
*             SECOND PARTNER.                                          *
*     MV    - CURRENT NUMBER OF CASCADE PARTICLES IN                   *
*             THE ARRAYS PM & IM.                                      *   
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES:                                                 *
*     NP  - NUMBER OF SECONDARY PARTICLES IN THE ELEMENTARY ACT.       *
************************************************************************      
      parameter (max1=2000)
      common/inel/inel              
      COMMON /PMIM/ pm(19,max1),IM(5,max1)                              
      DIMENSION V(3),P(11),IP(3),C(11),IC(3),A(3),B(3)                  
    1 I1=IP(1)                                                          
      I2=IC(1)                                                          
   14 CT=COSELN(IP,IC,T,P(9))                                           
      IF(CT.GE.1.) GOTO 14
      IF(CT.LE.-1.) GOTO 14
      F=6.283185*RNDM(-1)                                               
      CALL ABEL(P,V,U,A,B,CT,F,P(9),C(9))                               
      IF(MV-max1) 15,15,17                                              
   15 IM(1,MV+3)=I1                                                     
      IM(2,MV+3)=IP(2)                                                  
      PM(9,MV+3)=P(9)                                                   
      PM(10,MV+3)=P(10)                                                 
      PM(19,Mv+3)=P(11)
      DO 16 I=1,3                                                       
      PM(I,MV+3)=A(I)                                                   
   16 PM(I,MV+1)=B(I)                                                   
      IM(1,MV+1)=I2                                                     
      IM(2,MV+1)=IC(2)                                                  
      PM(9,MV+1)=C(9)                                                   
      PM(10,MV+1)=C(10)                                                 
      PM(19,MV+1)=C(11)
      PM(16,Mv+1)=0.
      PM(17,Mv+1)=0.
      NP=2                                                              
      RETURN                                                            
   17 NP=0                                                              
      PRINT 18                                                          
   18 FORMAT (45X,27H MEMORY IS EXCEEDED IN BLEL)                       
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE BLABSN(P,IP,C,IC,MV,NP,V,D,ID,NE,IE)                   
************************************************************************
*     PURPOSE: CALCULATES PION ABSORPTION ON QUASIDEUTERON PAIR.       *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     P, IP    - CHARACTERISTICS OF PION.                              *
*     C, IC    - CHARACTERISTICS OF THE FIRST NUCLEON OF QUASIDEUTERON.*     
*     D, ID    - CHARACTERISTICS OF THE SECOND NUCLEON OF QUASIDEUTERON*     
*     MV          - NUMBER OF PARTICLES IN THE ARRAYS PM AND IM        *
*     V        - VELOCITY OF CMS IN INTERACTING PARTICLES IN LAB.      *
*     NE          _ CHARGE OF THE FIRST SECONDARY.                     *
*     IE          - CHARGE OF THE OTHER SECONDARY.                     *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES:                                                 *
*     NP  -  NUMBER OF SECONDARY PARTICLES (HERE NP=2)                 *    
************************************************************************      
      parameter (max1=2000)
      common/inel/inel              
      COMMON /PMIM/ pm(19,max1),IM(5,max1)                              
      DIMENSION P(11),IP(3),C(11),IC(3),V(3),D(11),ID(3),X(3),Y(3),W(3) 
      PD=SQRT(D(8)*(D(8)+1.88))                                         
      PC=SQRT(C(8)*(C(8)+1.88))                                         
      X(1)=PC*C(4)*C(7)+PD*D(4)*D(7)                                    
      X(2)=PC*C(4)*C(6)+PD*D(4)*D(6)                                    
      X(3)=PC*C(5)+PD*D(5)                                              
      XM=SQRT(X(1)*X(1)+X(2)*X(2)+X(3)*X(3))                            
      D(5)=X(3)/XM                                                      
      T1=1.-D(5)*D(5)                                                   
      IF(T1) 1,1,2                                                      
    1 D(4)=0.                                                           
      D(5)=1.                                                           
      D(6)=0.                                                           
      D(7)=1.                                                           
      GO TO 3                                                           
    2 D(4)=SQRT(T1)                                                     
      T2=XM*D(4)                                                        
      D(6)=X(2)/T2                                                      
      D(7)=X(1)/T2                                                      
    3 D(8)=SQRT(XM*XM+1.88*1.88)-1.88                                   
      D(9)=1.88                                                         
      CALL TINVU (P,D,T,V,U)                                            
      IF(MV-max1) 4,4,13                                                
    4 M1=MV+1                                                           
      M3=MV+3                                                           
      IM(1,M1)=NE                                                       
      IM(1,M3)=IE                                                       
      IM(2,M1)=1                                                        
      IM(2,M3)=1                                                        
      IF(T-.455) 9,9,10                                                 
    9 CT=APCOS(19,T)                                                    
      GO TO 11                                                          
   10 CT=1.-2.*RNDM(-1)                                                 
   11 F=6.283185*RNDM(-1)                                               
      CALL ABEL(P,V,U,Y,W,CT,F,.94,.94)                                 
      PM(9,M1)=.94                                                      
      PM(9,M3)=.94                                                      
      PM(10,M1)=0.                                                      
      PM(10,M3)=0.                                                      
      PM(19,M1)=1.
      PM(16,M1)=0.
      PM(17,M1)=0.
      PM(19,M3)=1.
      PM(16,M3)=0.
      PM(17,M3)=0.
      DO 12 J=1,3                                                       
      PM(J,M1)=-Y(J)                                                    
   12 PM(J,M3)=Y(J)                                                     
      NP=2                                                              
      RETURN                                                            
   13 NP=0                                                              
      PRINT 14                                                          
   14 FORMAT(45X,28H MEMORY IS FXCEEDED IN BLABS)                       
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE TIPULN(P,IP,C,IC,V,U,T,SIG,SAB,MV,NP,NABS,C1,IC1,      
     *N3,NU,NIN)                                                        
************************************************************************
*     PURPOSE: DEFINES THE TYPE OF THE INTERACTION (ELASTIC,           *
*              INELASTIC OR PION ABSORPTION).                          *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *    
*     P ,IP - KINEMATICAL CHARACTERISTICS & QUANTUM NUMBERS OF THE     *
*             FIRST PARTNER.                                           *
*     C ,IC - KINEMATICAL CHARACTERISTICS & QUANTUM NUMBERS OF THE     *
*             SECOND PARTNER.                                          *    
*     V     - VELOCITY OF CMS OF THE INTERACTING PARTNERS IN THE LAB.  *
*             SYSTEM.                                                  *
*     U     - TOTAL ENERGY OF THE INERACTING PARTNERS IN THE CMS.      *
*     T     - KINETIC ENERGY OF THE FIRST PARTNER IN THE REST FRAME    *
*           - OF THE SECOND PARTNER.                                   *
*     SIG   - TOTAL CROSS SECTION FOR THE INTERACTION.                 *
*     SAB   - CROSS SECTION FOR ABSORPTION OF PION ON QUASIDEUTERON    *
*               PAIR.                                                  *  
*     MV    - CURRENT NUMBER OF CASCADE PARTICLES IN                   *
*             THE ARRAYS PM & IM.                                      *   
*     N3    - NO OF THE SECOND NUCLEON IN  PION ABSORPTION PROCESS.    *
*     NU    - LABEL SPECIFYING NUCLEUS WHERE INTERACTION TAKES PLACE:  *
*               NU=1 --> PROJ.NUCLEUS, NU=2 --> TARGET NUCLEUS.        *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES:                                                 *
*     NP    -  NUMBER OF SECONDARY PARTICLES IN THE ELEM. ACT.         *
*     C1,IC1 -  CHARAC. OF SECOND NUCLEON IN PION ABSORPTION PROCESS.  *
*     NABS   - LABEL (=0 -> PION ABSORPTION PROCESS, =1 -> OTHER)      *
************************************************************************      
      common/inel/inel              
      DIMENSION P(11),IP(3),C(11),IC(3),V(3),C1(11),IC1(3)              
c****************
c      sab=0.
c***************
   1  NABS=0                                                            
      IF(SIG.EQ.0..AND.SAb.EQ.0.) GOTO 20
   4  IF(N3) 15,15,5                                                    
   5  BETABS=SAB/(SIG+SAB)                                              
      IF(RNDM(-1)-BETABS) 6,6,15                                        
   6  CALL PARTNN(NU,N3,C1,IC1)                                         
      C1(9)=C(9)                                                        
      IE1=IP(1)                                                         
      IE2=IC(1)                                                         
      IE3=IC1(1)                                                        
      IF(IE1)   11,10,7                                                 
   7  IF(IE2+IE3-1) 9,8,15                                              
   8  NE1=1                                                             
      NE2=1                                                             
      GO TO 14                                                          
   9  NE1=1                                                             
      NE2=0                                                             
      GO TO 14                                                          
  10  NE1=IE2                                                           
      NE2=IE3                                                           
      GO TO 14                                                          
   11 IF(IE2+IE3-1) 15,12,13                                            
  12  NE1=0                                                             
      NE2=0                                                             
      GO TO 14                                                          
   13  NE1=0                                                            
      NE2=1                                                             
 14   CALL BLABSN(P,IP,C,IC,MV,NP,V,C1,IC1,NE1,NE2)                     
      NABS = 1                                                          
      RETURN                                                            
   15 BETAEL=(SIGMAN(IP,IC,1,T)+SIGMAN(IP,IC,2,T))/SIGMAN(IP,IC,0,T)    
      IF(RNDM(-1)-BETAEL) 16,16,17                                      
   16 IF(P(9).GT..2) GOTO 20                                            
      IF(C(9).NE..94) GO TO 20                                          
      CALL ISEX(P,IP,C,IC,T,V,U,MV,NP)                                  
      IF(NP.EQ.1) RETURN                                                
   20 CALL BLELN(V,U,T,P,IP,C,IC,MV,NP)                                 
      RETURN
   21 NCT=NCT+1
      IF(NCT-10) 22,22,18
   17 NCT=0
cuj   17 IF(P(9).GT.0.94) GOTO 22      ! UJ  27.06.2003
cuj	NCT=0			    ! UJ
   22 NIN=0
      CALL BLACT(P,IP,C,IC,V,U,T,MV,NP,NIN,L1,L2)                       
      IF(NIN.GT.0) GOTO 21                                              
      CALL CHARGE(IP,IC,MV,NP,L1,L2,NIN)                                
      IF(NIN.GT.0) GOTO 21
  18  RETURN                                                            
      END                                                               
***********************************************************************
      FUNCTION SGR(DM,GM,PT)                                            
***********************************************************************
C  CALCULATION OF RESONANT CROSS SECTION                                
*----------------------------------------------------------------------
      DM0=1.232                                                         
       HC=0.197                                                         
      HCT=(HC*PT)**2                                                    
      EPI=SQRT(0.14**2+HCT)                                             
      EN=SQRT(0.94**2+HCT)                                              
      DM=EN+EPI                                                         
      DMM0=(DM**2-DM0**2)**2                                            
       DMG=(DM0*GM)*(DM0*GM)                                            
      SGR=(251.3274*DMG)/((PT**2)*(DMG+DMM0))                           
      RETURN                                                            
      END                                                               
*************************************************************************
      SUBROUTINE ISEX(P,IP,C,IC,T,V,U,MV,NP)                            
************************************************************************
*  PURPOSE:  CALCULATES THE REACTION PI+N --> DELTA.                   *
*----------------------------------------------------------------------*    
      parameter (max1=2000)
      COMMON /PMIM/ pm(19,max1),IM(5,max1)                              
      DIMENSION P(11),C(11),IP(3),IC(3),V(3),PL(3)                      
      NT=MV+1                                                           
      NP=0                                                              
      IKS=IP(1)+IC(1)+2                                                 
      GOTO (1,2,2,1),IKS                                                
    1 AK=1.                                                             
    2 AK=.3334                                                          
      IF(IP(1).EQ.0)AK=.66667                                           
      GM=GAMD(0,NT,PT,P,V)                                              
      SIGR=SGR(DM,GM,PT)                                                
      SELX=SIGMAN(IP,IC,1,T)+SIGMAN(IP,IC,2,T)                          
      PR=SIGR*AK/SELX                                                   
      IF(RNDM(-1).GE.PR) RETURN                                         
      IF(DM.LE.1.1) RETURN                                              
      PM(1,NT)=0.                                                       
      PM(2,NT)=0.                                                       
      PM(3,NT)=0.                                                       
    3 PM(4,NT)=0.                                                       
      PM(5,NT)=1.                                                       
      PM(6,NT)=0.                                                       
      PM(7,NT)=1.                                                       
    5 PM(9,NT)=DM                                                       
      PM(8,NT)=0.                                                       
      IM(1,NT)=IKS-2                                                    
      IM(2,NT)=1                                                        
      IM(3,NT)=0                                                        
      IM(4,NT)=0                                                        
      PM(10,NT)=GM                                                      
      PM(19,NT)=1.
      PM(16,NT)=0.
      PM(17,NT)=0.
      NP=1                                                              
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE PLAM(P,DAM,TR,TPR)                                     
************************************************************************
*   PURPOSE:   CALCULATES PATHLENGTH AND LIFE-TIME OF RESONANCE.      *
*----------------------------------------------------------------------*
*   INPUT:                                                             *
*     P     - KINEMATICAL CHARAC. OF RESONANCE.                        *
*----------------------------------------------------------------------*
*   OUTPUT:                                                            *
*     DAM   - IT'S PATHLENGTH.                                        *    
*     TR    - IT'S LIFE TIME.                                          * 
************************************************************************     
      DIMENSION P(11)                                                   
      PR=SQRT(P(8)*(P(8)+2.D0*P(9)))                                    
      ER=P(8)+P(9)                                                      
      BR=PR/ER                                                          
      GAM=ER/P(9)                                                       
      IF(P(9)-1.) 1,1,2
    2 GM=ABS(P(9)-1.232)
      IF(GM.EQ.0.) GM=.0001
      GOTO 10
    1 GM=P(10)
   10 B=GAM/GM                                                          
      T=.197*B                                                         
      IF(P(9)-1.) 3,3,4
    3 TR=-T*ALOG(1.-RNDM(-1))
      GOTO 5
    4 TR=T
    5 DAM=BR*TR                                                         
      RETURN                                                            
      END                                                               
************************************************************************
      FUNCTION COSELN(IP,IC,T,W)                                        
************************************************************************
*     PURPOSE:   CALCULATES THE COS(TETA) FOR SECONDARIES              *
*                IN ELASTIC REACTIONS.                                 *
************************************************************************
*   INPUT QUANTITIES:                                                  *
*     IP - ARRAY OF QUANTUM NUMBERS OF THE FIRST PARTNER.              *
*     IC - THE SAME FOR THE SECOND PARTNER.                            *
*     T     - KINETIC ENERGY OF THE FIRST PARTNER IN THE REST FRAME    *
*             OF THE SECOND PARTNER.                                   *
*     W     - MASS OF THE FIRST PARTNER.                               *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITY:                                                   *
*     COSELN --> COS(TETA) IN THE ELASTIC INTERACTION                  *
************************************************************************     
      DIMENSION IP(3),IC(3)                                             
      M=IP(2)+IC(2)                                                     
      N=IP(1)+IC(1)                                                     
      TM=SQRT(T*(T+2.*W))
      IF(M-1) 1,1,4                                                     
    1 IF(IP(1)) 2,8,3                                                   
    2 IF(IC(1)) 6,6,7                                                   
    3 IF(IC(1)) 7,7,6                                                   
    4 IF(M-2) 5,5,6                                                     
    5 IF(N-1) 6,7,6                                                     
    6 K=1                                                               
      GO TO 9                                                           
    7 K=2                                                               
      GO TO 9                                                           
    8 K=3                                                               
    9 IF(M-2) 31,21,10                                                  
   10 IF(T-.147) 11,11,12                                               
   11 A1=25.2*T**(-.843)                                                
      GO TO 13                                                          
   12 A1=130.*T**0.0145                                                 
   13 A2=11.3*T**0.432                                                  
      IF(T-.055) 14,14,15                                               
   14 A3=.22*T**(-1.35)                                                 
      GO TO 16                                                          
   15 A3=0.000043*T**(-4.32)                                            
   16 A4=130.*T**1.33                                                   
      T1=1.-EXP(-A2*3.141592)                                           
      IF(A4-50.)101,101,100                                             
  100 T2=1.                                                             
      GO TO 102                                                         
  101 T2=1.-EXP(-A4*3.141592)                                           
  102 Z=A3*(2.-T2)/((1.+A4*A4)*(A1*(2.-T1)/(1.+A2*A2)+A3*(2.-T2)/(1.+   
     1A4*A4)))                                                          
      IF(RNDM(-1)-Z) 17,17,19                                           
   17 X=-ALOG(1.-RNDM(-1)*T2)/A4                                        
      IF(RNDM(-1)-SIN(X)) 18,18,17                                      
   18 Y=3.141592-X                                                      
      GO TO 20                                                          
   19 Y=-ALOG(1.-RNDM(-1)*T1)/A2                                        
      IF(RNDM(-1)-SIN(Y)) 20,20,19                                      
   20 COSELN=COS(Y)                                                     
      RETURN                                                            
   21 IF(K-2) 24,22,22                                                  
   22 IF(T-.97) 23,23,26                                                
   23 COSELN=APCOS(3,T)                                                 
      RETURN                                                            
   24 IF(T-.46) 25,25,26                                                
   25 COSELN=1.-2.*RNDM(-1)                                             
      RETURN                                                            
   26 IF(T-2.8) 27,27,28                                                
   27 COSELN=(1.+APCOS(1,T))/2.                                         
      RETURN                                                            
   28 IF(T-10.) 29,29,30                                                
   29 BSL=5.6+.33*TM-.0017*TM**2                                         
      GOTO 60                                        
   30 BSL=8.25+.55*(ALOG(2.*TM*W))**2 
   60 EPR=BSL*TM
      IF(EPR.GT.50.) EPR=50.
      RND=RNDM(-1)
      COSELN=1.+(2.*ALOG(1.+RND*(EXP(-EPR)-1.)))/(BSL*TM)               
      RETURN                                                            
 31   IF(K-2) 32,41,50                                                  
   32 IF(T-.08) 33,33,34                                                
   33 COSELN=APCOS(4,T)                                                 
      RETURN                                                            
   34 IF(T-.3) 35,35,36                                                 
  35  COSELN=APCOS(5,T)                                                 
      RETURN                                                            
   36 IF(T-1.) 37,37,38                                                 
   37 COSELN=APCOS(6,T)                                                 
      RETURN                                                            
   38 IF(T-2.4) 39,39,40                                                
   39 COSELN=APCOS(7,T)                                                 
      RETURN                                                            
   40 TM=SQRT(T*(T+2.*W))                                               
      EPR=7.5*TM
      IF(EPR.GT.50.) EPR=50.
      COSELN=1.+(2.*ALOG(1.+RNDM(-1)*(EXP(-EPR)-1.)))/(7.5*TM)          
      RETURN                                                            
   41 IF(T-.08) 42,42,43                                                
   42 COSELN=APCOS(8,T)                                                 
      RETURN                                                            
   43 IF(T-.3)  44,44,45                                                
   44 COSELN=APCOS(9,T)                                                 
      RETURN                                                            
   45 IF(T-1.) 46,46,47                                                 
   46 COSELN=APCOS(10,T)                                                
      RETURN                                                            
   47 IF(T-2.4) 48,48,49                                                
   48 COSELN=APCOS(11,T)                                                
      RETURN                                                            
   49 TM=SQRT(T*(T+2.*W))                                               
      EPR=7.5*TM
      IF(EPR.GT.50.) EPR=50.
      COSELN=1.+(2.*ALOG(1.+RNDM(-1)*(EXP(-EPR)-1)))/(7.5*TM)           
      RETURN                                                            
   50 IF(RNDM(-1)-.5) 32,32,41                                          
      END                                                               
************************************************************************
      SUBROUTINE PARTNN(NU,I,C,IC)                                      
************************************************************************
*     PURPOSE : CHOICE OF PARTNER FOR INTERACTION OR ABSORBSION        *
*                INSIDE THE NUCLEUS                                    *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     NU     -LABEL (NU=1--> PROJECTILE NUCLEUS; NU=2 --> TARGET NUCL.)*
*     I      - NO OF NUCLEON-PARTNER IN THE ARRAY XYZ                  *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES:                                                 *
*     C      - KINEMATICAL CHARACTERISTICS OF PARTNER IN NUCLEUS       *
*     IC     - QUANTUM NUMBERS OF THE PARTNER (IC(1) - CHARGE,         *
*              IC(2) - BARYON NUMBER)                                  *
************************************************************************     
      common/inel/inel              
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON /RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC           
      COMMON /QDEU/QD(2,3),IQU(2)                                       
      COMMON/XYZIZ/ XYZ(2,3,240),IZ(2,240)                              
      DIMENSION C(11),IC(3)                                             
      DO 1 K=1,3                                                        
 1    C(K)=XYZ(NU,K,I)                                                  
      IF(IQU(NU).NE.2) GO TO 2                                          
      ALF=1.                                                            
      IF(IZ(NU,I).EQ.0) ALF=-1.                                         
      Q2=QD(NU,1)**2+QD(NU,2)**2+QD(NU,3)**2                            
      C(8)=SQRT(Q2+0.8836)-.94                                          
      Q=SQRT(Q2)                                                        
      C(5)=ALF*QD(NU,3)/Q                                               
      C(4)=SQRT(1.-C(5)*C(5))                                           
      C(6)=ALF*QD(NU,2)/(Q*C(4))                                        
      C(7)=ALF*QD(NU,1)/(Q*C(4))                                        
      GO TO 3                                                           
 2    CONTINUE                                                          
      RR=SQRT(C(1)**2+C(2)**2+C(3)**2)                                  
      TF=TR(NU)                                            
      EPR=(RR-AR(NU))/CR(NU)
      IF(EPR.GT.10.) EPR=10.
      TFR=TF*(((1.+EXP(-AR(NU)/CR(NU)))/                                
     1(1.+EXP(EPR)))**.66667)                                           
      C(8)=TFR*(RNDM(-1)**(2./3.))                                      
      C(5)=1.-2.*RNDM(-1)                                               
      C(4)=SQRT(1.-C(5)**2)                                             
      F=6.283185*RNDM(-1)                                               
      C(6)=SIN(F)                                                       
      C(7)=COS(F)                                                       
 3    CONTINUE                                                          
      C(9)=.94                                                          
      C(10)=0.                                                          
      IC(2)=1                                                           
      IC(1)=IZ(NU,I)                                                    
      IC(3)=NU                                                          
      C(9)=C(9)-POTENN(NU,A(NU),C,IC,AR(NU),CR(NU),DR(NU),TR(NU),EPS(NU)
     *,VPI)                                                             
      C(10)=0.                                                          
      if(C(8).EQ.0.) print 100,nu,iqu(nu),tf,tfr,c(8),c(3)
  100 format(10x,' iqu=',i5,5x,'nu=',i5,5x,'tf=',f10.4,5x,'tfr=',
     *f10.5,5x,'c(8)=',f10.5,5x,'c(3)=',f10.5)
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE BLCENN(NA,P,IP,DEL,NC,K1,K2,IK,NU,ITT,RSI)             
************************************************************************
*     PURPOSE : DETERMINES POSSIBLE PARTNERS INSIDE CYLINDER OF        *
*               RADIUS 'DEL' ALONG THE PATH OF CASCADE PARTICLE        *
*               PROPAGATING IN NUCLEUS.                                *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES :                                                 *
*     NA      - CURRENT NUMBER OF NUCLEONS IN THE NUCLEUS.             *
*     P, IP   - CHARACTERISTICS OF CASCADE PARTICLE.                   *
*     DEL     - RADIUS OF CYLINDER.                                    *
*     IK      - NOT USED                                               *
*     NU      - LABEL OF THE NUCLEUS; 1 -> PROJECTILE, 2 -> TARGET.    *
*     ITT     - LABEL; =0 -> DETERMINES NUCLEAR NUCLEONS INSIDE CYLIDRE*
*                            AND CHOOSES THE NEARIST ONE;              *
*                      >0 -> CHOOSES THE NEXT NUCLEON.                 *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES:                                                 *
*     NC - THE NUMBER OF NUCLEONS INSIDE CYLINDER.                     *
*     K1 - NO OF THE NEAREST ONE.                                      *
*     K2 - NO OF THE SECOND NEAREST NUCLEON FOR PION ABSORPSION.       *
************************************************************************     
      DIMENSION KD(240),KW(240),W(240),RAD(240),RD(240)                 
      COMMON /XYZIZ/XYZ(2,3,240),IZ(2,240)                              
      DIMENSION P(11),IP(3),RI(3),RC(3)                                 
      SAVE RD,KD
    1 IF(ITT)33,33,30                                                   
   33 NC=0                                                              
      K=1                                                               
      IF(IP(3))2,2,4                                                    
   2  DO 3 I=1,3                                                        
   3  RI(I)=P(I)                                                        
      GO TO 5                                                           
   4  RI(1)=(P(1)*P(7)+P(2)*P(6))*P(5)-P(3)*P(4)                        
      RI(2)= -P(1)*P(6)+P(2)*P(7)                                       
      RI(3)=(P(1)*P(7)+P(2)*P(6))*P(4)+P(3)*P(5)                        
   5  IF(IP(3)) 6,6,8                                                   
   6  DO 7 I=1,3                                                        
   7  RC(I)=XYZ(NU,I,K)                                                 
      GO TO 9                                                           
   8  RC(1)=(XYZ(NU,1,K)*P(7)+XYZ(NU,2,K)*P(6))*P(5)-XYZ(NU,3,K)*P(4)   
      RC(2)= -XYZ(NU,1,K)*P(6)+XYZ(NU,2,K)*P(7)                         
      RC(3)=(XYZ(NU,1,K)*P(7)+XYZ(NU,2,K)*P(6))*P(4)+XYZ(NU,3,K)*P(5)   
      GO TO 11                                                          
   9  IF(NU-1) 10,10,11                                                 
  10  IF(RI(3)-RC(3)) 17,17,12                                          
  11  IF(RC(3)-RI(3)) 17,17,12                                          
  12  IF(ABS(RC(3)-RI(3))-.0001) 17,17,13                               
  13  IF((RC(1)-RI(1))**2+(RC(2)-RI(2))**2-DEL**2) 14,14,17             
  14  NC=NC+1                                                           
      W(NC)=RC(3)                                                       
      RAD(NC)=SQRT((RC(1)-RI(1))**2+(RC(2)-RI(2))**2)                   
       KW(NC)=K                                                         
   17 IF(K-NA)18,22,22                                                  
  18  K=K+1                                                             
      GO TO 5                                                           
   22 IF(NC)20,20,32                                                    
   20 K2=0                                                              
   21 RETURN                                                            
   32 KK=NC                                                             
      J=0                                                               
   39 I=1                                                               
   40 S=W(I)                                                            
      IS=I                                                              
   41 IF(I-KK)42,45,45                                                  
   42 IF(S-W(I+1))43,43,44                                              
   43 I=I+1                                                             
       GOTO 41                                                          
   44 I=I+1                                                             
       GOTO 40                                                          
   45 J=J+1                                                             
      KD(J)=KW(IS)                                                      
      RD(J)=RAD(IS)                                                     
      W(IS)=W(KK)                                                       
      RAD(IS)=RAD(KK)                                                   
       KW(IS)=KW(KK)                                                    
      KK=KK-1                                                           
      IF(KK)48,48,39                                                    
   30 IF(NC-1)15,15,34                                                  
   15 NC=0                                                              
       GOTO 21                                                          
   34 NC=NC-1                                                           
      DO 35 M=1,NC                                                      
      RD(M)=RD(M+1)                                                     
   35 KD(M)=KD(M+1)                                                     
   48 K1=KD(1)                                                          
      RSI=RD(1)                                                         
      IF(NC-1)20,20,49                                                  
   49 K2=KD(2)                                                          
       GOTO 21                                                          
      END                                                               
************************************************************************
      FUNCTION POTENN(NA,A,P,IP,AR,CR,DR,TR,EPS,VPI)                    
************************************************************************
*     PURPOSE : CALCULATION OF POTENTIAL OF PARTICLE IN THE NUCLEUS.   *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     A        - ATOMIC NUMBER OF NUCLEUS.                             *
*     P, IP    - CHARACTERISTICS OF PARTICLE.                          *
*     AR,CR    - PARAMETERS OF WOODS-SAXON DISTRIBUTION OF NUC.DENSITY.*
*     DR       - CUT PARAMETER FOR NUCLEAR DENSITY TAIL.               *
*     TR         - FERMI ENERGY OF PARTICLE IN THE CENTRE OF NUCLEUS   *
*     EPS        - MEAN BOND ENERGY OF NUCLEON INSIDE THE NUCLEUS.     *
*     VPT        - MEAN POTENTIAL OF THE PION IN THE NUCLEUS ( 25[MEB])*
************************************************************************     
      COMMON /RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC           
      DIMENSION P(11),IP(3)                                             
      TN=TR*AC(NA)/A                                                    
      IF(EPS.LT.0.002) GO TO 4                                          
      IF(IP(2)) 2,1,2                                                   
    1 POTENN=VPI*AC(NA)/A                                               
      RETURN                                                            
 2    R=AR+CR*ALOG((1.-DR)/DR)                                          
      RR=SQRT(P(1)**2+P(2)**2+P(3)**2)/R                                
      EPR=(R*RR-AR)/CR
      IF(EPR.GT.10.) EPR=10.
      IF(RR-1.5) 3,4,4                                                  
   3  POTENN=EPS+TN*(((1.+EXP(-AR/CR))/(1.+EXP(EPR)))**.66667)          
      IF(P(10).GT.0.) POTENN=POTENN+VPI                                 
      RETURN                                                            
   4  POTENN =0.                                                        
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE BLRECN(NU,PX,PY,PZ,X,Y,Z)                              
************************************************************************
*     PURPOSE : CALCULATES THE CHANGE OF MOMENTUM AND ANGULAR          *
*               MOMENTUM OF THE NUCLEUS AFTER INTERACTION OR           *
*               ABSORPTION OF PARICLE INSIDE THE NUCLEUS.              *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES:                                                  *
*     NU       -  LABEL: NU=1--> PROJECTILE NUCLEUS; NU=2 --> TARGET   *
*                        NUCLEUS.                                      *
*     PX,PY,PZ -  MOMENTUM COMPONENTS OF THE PARTICLE.                 *
*     X,Y,Z    - COORDINATES OF THE PARTICLE                           *
************************************************************************    
      COMMON /RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC           
      PC(NU,1)=PC(NU,1)+PX                                              
      PC(NU,2)=PC(NU,2)+PY                                              
      PC(NU,3)=PC(NU,3)+PZ                                              
      AMC(NU,1)=AMC(NU,1)+Z*PY-Y*PZ                                     
      AMC(NU,2)=AMC(NU,2)+X*PZ-Z*PX                                     
      AMC(NU,3)=AMC(NU,3)+Y*PX-X*PY                                     
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE TFERMN(R0,A,AR,CR,DR,TR,R,BR)                          
************************************************************************
*     PURPOSE : CALCULATION OF FERMI ENERGY OF NUCLEON IN ThE CENTRE   *
*               OF NUCLEUS AND RADIUS OF NUCLEUS.                      *
************************************************************************
*   INPUT QUANTITIES :                                                 *
*     R0 - PARAMETER FOR CALCULATION OF RADIUS OF NUCLEUS. ( 1.07[FM] )*
*     A  - NUMBER OF NUCLEONS IN THE NUCLEUS                           *
*     CR - PARAMETER OF WOODS-SAXON DISTRIBUTION OF NUCLEAR DENSITY    *
*     DR - CUT PARAMETER FOR NUCLEAR DENSITY TAIL                      *
*     BR - PARAMETER FOR BUILDING OF DISTRIBUTION OF NUCLEONS          *
*          IN ThE NUCLEUS WIDH A(<=16)                                 *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES:                                                 *
*     TR  - FERMI ENERGY OF NUCLEON IN THE CENTRE OF NUCLEUS. [GEV]    *
*     AR  - PARAMETER OF WOODS-SAXON DIST. OF NUCLEAR DENSITY          *
*     R   - RADIUS OF NUCLEUS. [FM]                                    *
************************************************************************    
      DIMENSION W(8),X(8)                                               
      DATA W /.1012285363,.2223810345,.3137066459,.362683734,           
     1        .3626837834,.3137066459,.2223810345,.1012285363/          
      DATA X / .9602898565, .7966664774, .5255324099, .1834346425,      
     2        -.1834346425,-.5255324099,-.7966664774,-.9602898565/      
      AR=R0*A**0.333333                                                 
      R=AR+CR*ALOG((1.-DR)/DR)                                          
      S=0.                                                              
      DO 1 K=1,8                                                        
      EPR=R*(X(K)+1.)/(CR+CR)
      IF(EPR.GT.10.) EPR=10.
      SK=(((X(K)+1.)**2)/(EXP(EPR)+EXP(AR/CR)))*W(K)                    
    1 S=S+SK                                                            
      S=R*R*R*S/8.                                                      
      ROZERO=A/(12.56637*(EXP(AR/CR)+1.)*S)                             
      TR=0.1985*(ROZERO/2.)**0.666667                                   
      IF(A.LE.4.) R=BR*SQRT(-ALOG(DR))                                  
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE CINEMA (A,V,P)                                         
************************************************************************
*     PURPOSE : LORENTZ TRANSFORMATION FROM C.M.S TO LAB SYSTEM.       *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES :                                                 *
*     A - COMPONENTS OF MOMENTUM OF THE PARTICLE.                      *
*     V - VELOCITY COMPONENTS OF C.M.S.                                *
*     P(9) - MASS OF THE PARTICLE.                                     *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES :                                                *
*     P(4),P(5),P(6),P(7) - SIN(TETA),COS(TETA),SIN(PHI),COS(PHI).     *
*     P(8)                - KINET.ENERGY OF PARTICLE IN THE LAB.       *
************************************************************************      
      DIMENSION A(3),V(3),P(11)                                         
      AV=A(1)*V(1)+A(2)*V(2)+A(3)*V(3)                                  
      V2=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)                                  
      IF(V2) 5,5,6                                                      
    5 P(1)=A(1)                                                         
      P(2)=A(2)                                                         
      P(3)=A(3)                                                         
      GOTO 7                                                            
    6 T1=SQRT(1.D0-V2)                                                  
      T2=AV*(1.D0/T1-1.D0)/V2                                           
      T=SQRT(A(1)*A(1)+A(2)*A(2)+A(3)*A(3)+P(9)*P(9))-P(9)              
      DO 1 I=1,3                                                        
    1 P(I)=A(I)+V(I)*(T2+(T+P(9))/T1)                                   
    7 PM=SQRT(P(1)*P(1)+P(2)*P(2)+P(3)*P(3))                            
      IF(PM.EQ.0.) GOTO 2                                               
      P(5)=P(3)/PM                                                      
      T4=1.D0-P(5)*P(5)                                                 
      IF(T4)12,12,3
   12 P(4)=0.
      P(6)=0.
      P(7)=1.
      GOTO 4                                                            
    2 P(4)=0.                                                           
      P(5)=1.                                                           
      P(6)=0.                                                           
      P(7)=1.                                                           
      GO TO 4                                                           
    3 P(4)=SQRT(T4)                                                     
      T3=PM*P(4)                                                        
      P(7)=P(1)/T3                                                      
      P(6)=P(2)/T3                                                      
    4 P(8)=SQRT(PM*PM+P(9)*P(9))-P(9)                                   
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE TINVU (P,C,T,V,U)                                      
************************************************************************
*     PURPOSE : CALCULATION OF KINEMATICS OF INTERACTION.              *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES :                                                 *
*     P -  CHARACTER. OF THE FIRST COLLIDING PARTICLE IN LAB.          *
*     C -  CHARACTER. OF PARTNER IN LAB.                               *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES :                                                *
*     T - KINETIC ENERGY OF COLLIDING PARTICLE IN THE REST             *
*         FRAME OF PARTNER.                                            *
*     V - VELOCITY COMPONENTS OF C.M.S. IN LAB.                        *
*     U - TOTAL ENERGY OF SYSTEM OF COLLIDING PARTICLES IN C.M.S.      *
************************************************************************      
      DIMENSION P(11),C(11),V(3)                                        
      A=SQRT(P(8)*(P(8)+2.D0*P(9)))                                     
      A1=A*P(4)*P(7)                                                    
      A2=A*P(4)*P(6)                                                    
      A3=A*P(5)                                                         
      B=SQRT(C(8)*(C(8)+2.D0*C(9)))                                     
      B1=B*C(4)*C(7)                                                    
      B2=B*C(4)*C(6)                                                    
      B3=B*C(5)                                                         
      E1=P(8)+P(9)                                                      
      E2=C(8)+C(9)                                                      
      V(1)=(A1+B1)/(E1+E2)                                              
      V(2)=(A2+B2)/(E1+E2)                                              
      V(3)=(A3+B3)/(E1+E2)                                              
      V2=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)                                  
      U=(E1+E2)*SQRT(1.D0-V2)                                           
      IF(P(9)-C(9)) 1,1,2                                               
    1 T=(U*U-(P(9)+C(9))*(P(9)+C(9)))/(C(9)+C(9))                       
      IF(T.LT.0.) T=0.
      RETURN                                                            
    2 T=(U*U-(P(9)+C(9))*(P(9)+C(9)))/(P(9)+P(9))                       
      IF(T.LT.0.) T=0.
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE ABEL (P,V,U,A,B,C,F,W1,W2)                             
************************************************************************
*    PURPOSE : CALCULATE MOMENTUM OF SECONDARY PARTICLES IN C.M.S.     *
*              IN REACTIONS OF ABSORBSION AND ELASTIC SCATTERING       *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES :                                                 *
*     P   -  CHARACTER. OF COLLIDING PARTICLE.                         *
*     V   -  VELOCITY OF C.M.S.                                        *
*     U     -  TOTAL ENERGY OF COLLIDING PARTICLES IN C.M.S.           *
*     C     -  COS(TETA) OF SECONDARY PARTICLE IN THE C.M.S.           *
*     F     -  AZIMUTHAL ANGLE (PHI) OF SECON.PARTICLE                 *
*     W1    -  MASS OF FIRST SECONDARY PARTICLE                        *
*     W2    -  MASS OF SECOND SECONDARY PARTICLE                       *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES :                                                *
*     A AND B -COMPONENTS OF MOMENTA OF SEC.PARTICLES IN THE C.M.S.    *
************************************************************************     
      DIMENSION P(11),V(3),A(3),B(3),D(3),R(3)                          
      E=(U*U+W1*W1-W2*W2)/(U+U)                                         
      X=SQRT(E*E-W1*W1)                                                 
      V2=V(1)*V(1)+V(2)*V(2)+V(3)*V(3)                                  
      T=SQRT(P(8)*(P(8)+2.*P(9)))                                       
      S=T*(P(4)*P(7)*V(1)+P(4)*P(6)*V(2)+P(5)*V(3))                     
      Y=S*(1./SQRT(1.-V2)-1.)/V2-(P(8)+P(9))/SQRT(1.-V2)                
      R(1)=T*P(4)*P(7)+Y*V(1)                                           
      R(2)=T*P(4)*P(6)+Y*V(2)                                           
      R(3)=T*P(5)+Y*V(3)                                                
      D(1)=X*SQRT(1.-C*C)*COS(F)                                        
      D(2)=X*SQRT(1.-C*C)*SIN(F)                                        
      D(3)=X*C                                                          
      CALL ROTOR (R,V,D,A)                                              
      DO 1 I=1,3                                                        
    1 B(I)=-A(I)                                                        
      RETURN                                                            
      END                                                               
************************************************************************
      FUNCTION APCOS(J,T)                                               
************************************************************************
*     PURPOSE : CALCULATE COS(TETA) USING APPROXIMATIONS OF            *
*               EXPERIMENTAL DATA.                                     *
************************************************************************    
      DIMENSION A(4,4,21), A1(112),A2(96),A3(128)                       
      EQUIVALENCE (A(1,1,1),A1(1)),(A(1,1,8),A2(1)),(A(1,1,14),A3(1))   
      DATA A1/                                                          
     12.7404,-9.6998,10.4,2.3882,-7.5137,44.096,-74.379,46.038,7.5479,  
     2-39.274,64.835,-41.609,-1.8369,8.6911,-13.06,7.188,-30.853,106.24,
     3-129.39,54.339,19.465,-68.102,96.358,-56.827,-3.4831,12.341,      
     4-18.592,12.024,0.18941,-0.6788,1.0665,-0.7291,0.10258,-1.0542,    
     511.389,-16.638,-0.49607,11.8,-90.857,164.76,1.5437,-33.769,251.92,
     6-450.71,-1.2021,25.336,-186.58,332.54,0.15789,2.9671,-5.5251,     
     76.8925,-7.0218,-205.34,569.51,-898.58,134.96,4872.2,-14674.,23924.
     8,-821.16,-32586.,100980.,-165530.,0.31531,-7.4981,43.295,-76.36,  
     9-6.5373,193.07,-1018.1,1742.6,46.864,-1303.,6729.1,-11075.,-95.192
     A,2637.3,-12857.,20294.,-17.953,109.72,-239.54,228.26,91.968,      
     B-519.63,1126.6,-1074.,-132.7,741.12,-1600.,1524.9,58.598,-318.74, 
     C677.51,640.11,0.42169,147.05,-653.35,915.07,-3.5198,-260.19,      
     D1225.,-1748.1,3.6373,155.92,-752.01,1079.6,-0.78041,-30.563,      
     E147.95,-212.5/                                                    
      DATA A2/                                                          
     1-0.38288,3.7587,-6.5144,6.774,103.81,-272.82,477.59,-512.22,      
     2-1788.2,4305.2,-7931.4,9347.1,7147.5,-3339.5,-4139.2,-4436.4,     
     3.24991,32.028,-118.82,150.99,-2.6994,-460.45,1895.9,-2519.,16.268,
     42138.4,-9126.2,12431.,-29.654,-3182.3,13944.,-19342.,3.9025,      
     5-91.126,323.73,-400.48,-20.619,491.7,-1715.5,2114.3,33.004,-766.84
     6,2700.3,-3352.5,-16.367,373.94,-1320.2,1642.3,19.402,-224.46,     
     7747.33,-935.7,-44.18,471.94,-1485.6,1805.5,31.567,-301.76,907.63, 
     8-1077.3,-6.8648,60.476,-175.2,203.81,0.14988,2.8753,-5.3078,6.2233
     9,-5.9558,-162.03,430.79,-625.48,128.75,3140.2,-7918.9,10983.,     
     A-851.61,-18780.,44607.,-58790.,0.53689,-13.216,81.011,-142.85,    
     B-10.55,296.29,-1695.7,2893.5,69.621,-1924.5,10620.,-17468.,       
     C-138.65,3928.1,-20293.,32058./                                    
      DATA A3/                                                          
     10.085591,5.039,-13.782,14.661,0.054284,-9.2324,36.397,-42.962,    
     2-0.051111,4.6003,-20.534,27.731,7.4514E-03,-0.62529,2.9159,-4.1101
     3,0.071622,3.096,-11.125,18.13,0.092581,-3.2186,20.273,-33.245,    
     4-0.051531,0.89886,-7.5084,13.188,5.8258E-03,-1.7288E-03,0.70224,  
     5-1.4856,0.0823,0.134,3.7716,-4.0562,0.010802,-0.33688,1.1727,     
     6-0.67476,-2.1798E-03,0.052166,-0.25816,0.32048,6.5764E-05,-1.4711E
     7-03,7.8209E-03,-0.01058,0.11138,0.60396,3.0174,-4.419,-0.017709,  
     80.23015,-1.8187,3.4518,2.0977E-03,-0.025458,0.21626,-0.40692,     
     9-5.4799E-05,5.9111E-04,-5.5552E-03,0.010647,0.17288,7.108,-17.961,
     A16.403,-0.14504,-13.032,41.781,-40.799,0.04539,8.3515,-30.26,     
     B32.882,-4.7961E-03,-1.4095,5.3505,-6.0946,0.037596,1.4331,-3.135, 
     C6.4864,0.23827,1.8253,1.7648,-16.735,-0.1541,-1.5201,-1.5692,     
     D17.185,0.025037,0.30588,0.3252,-3.5277,0.12489,1.3573,0.82338,    
     E-1.4595,-0.051577,-0.35778,-1.169,1.8078,7.4864E-03,0.032888,     
     F0.23744,-0.39802,-2.988E-04,-7.5117E-04,-0.011402,0.019505,0.1847,
     G1.9269,-3.2979,3.6843,-0.07393,0.27213,1.06,-2.3354,0.018907,     
     H-0.056473,-0.16487,0.38426,-9.2984E-04,2.5506E-03,7.3052E-03,     
     I-0.01722/                                                         
      S1=0.                                                             
      S2=0.                                                             
      B=RNDM(-1)                                                        
      DO 1 N=1,4                                                        
      DO 1 K=1,4                                                        
      S1=S1+A(N,K,J)*(T**(K-1))*(B**(N-1))                              
    1 S2=S2+A(N,K,J)*T**(K-1)                                           
      C=2.*SQRT(B)*(S1+(1-S2)*B**4)-1.                                  
      IF(ABS(C)-1.) 2,2,3                                               
    2 APCOS =C                                                          
      RETURN                                                            
    3 APCOS=SIGN(1.,C)                                                  
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE ROTOR (A,B,C,D)                                        
************************************************************************
*     PURPOSE :  ROTATION OF COORDINATE SYSTEM                         *
************************************************************************      
      DIMENSION A(3),B(3),C(3),D(3)                                     
      AM=SQRT(A(1)*A(1)+A(2)*A(2)+A(3)*A(3))                            
      IF(AM.EQ.0.) GOTO 1                                               
      B1=(A(1)*B(1)+A(2)*B(2)+A(3)*B(3))/AM                             
      BBB=     B(1)*B(1)+B(2)*B(2)+B(3)*B(3)-B1*B1                      
      IF(BBB.LE.0.) GOTO 1                                              
      B2=SQRT(BBB)                                                      
      A1=A(2)*B(3)-A(3)*B(2)                                            
      A2=A(3)*B(1)-A(1)*B(3)                                            
      A3=A(1)*B(2)-A(2)*B(1)                                            
      D(1)=C(1)*B(1)/B2+(C(3)-B1*C(1)/B2)*A(1)/AM+C(2)*A1/(B2*AM)       
      D(2)=C(1)*B(2)/B2+(C(3)-B1*C(1)/B2)*A(2)/AM+C(2)*A2/(B2*AM)       
      D(3)=C(1)*B(3)/B2+(C(3)-B1*C(1)/B2)*A(3)/AM+C(2)*A3/(B2*AM)       
      RETURN                                                            
    1 D(1)=C(1)                                                         
      D(2)=C(2)                                                         
      D(3)=C(3)                                                         
      RETURN                                                            
      END                                                               
************************************************************************
      FUNCTION APSIG(T,L)                                               
************************************************************************
*     PURPOSE : INTERPOLATION OF CROSS SEACTION USING EXPERIMENTAL DATA*
*----------------------------------------------------------------------*
*   INPUT QUANTITIES :                                                 *
*     T        - KINETIC ENERGY OF THE INTERACTING PARTICLE            *
************************************************************************     
      DIMENSION S(30,26),S1(120),S2(180),S3(270),S4(120),S5(90)         
      DIMENSION A(30,5), A1(150)                                        
     0EQUIVALENCE (S(1,1),S1(1)),(S(1,5),S2(1)),(S(1,11),S3(1)),        
     1(S(1,20),S4(1)),(S(1,24),S5(1)),(A(1,1),A1(1))                    
      DATA S1/                                                          
     117613.,330.,154.,96.,70.,51.,38.4,30.,23.6,22.4,22.2,22.6,23.4,   
     224.7,29.5,40.5,48.5,47.4,47.,46.7,46.,45.,43.,41.2,40.8,40.5,39.8,
     339.0,39.,38.5,17613.,330.,154.,96.,70.,51.,38.4,30.,23.6,22.4,22.2
     4,22.6,22.8,22.9,5*25.,22.,19.5,17.5,15.,13.7,12.,11.,9.8,8.8,8.5, 
     56.5,20357.,950.,480.,300.,200.,160.,108.,74.,50.,41.,36.5,34.,32.5
     6,32.,34.2,36.1,37.8,38.4,39.,40.,40.,40.5,41.1,42.3,42.3,42.,41.1,
     739.6,39.5,39.2,20357.,950.,480.,300.,200.,160.,108.,74.,50.,41.,  
     836.5,34.,32.5,31.2,30.8,25.1,19.1,18.,16.5,16.,15.2,14.,12.1,10.8,
     910.7,11.,9.6,8.,7.,6.2/                                           
      DATA S2/                                                          
     12*6.,6.5,7.,8.5,10.5,16.5,25.3,57.5,68.5,64.5,52.,40.5,25.7,29.,  
     232.1,45.6,38.,44.3,54.,58.,45.7,35.3,35.,34.3,34.,32.4,30.4,26.5, 
     325.,2*2.,2.25,2.3,2.4,2.5,4.5,7.7,20.,25.2,23.9,21.,16.,10.5,11.2,
     414.,20.3,16.,19.3,26.,26.5,18.7,11.7,10.5,9.6,9.5,6.8,5.9,5.,4.,  
     52*4.,4.25,4.7,6.1,8.,12.,17.6,37.5,43.3,40.6,31.,24.5,13.1,11.,9.5
     6,8.3,5.,5.5,6.8,7.,3.5,2.,1.9,1.8,1.6,.22,.15,.048,.009,1.9,2.3,  
     73.5,5.5,9.,14.,28.,60.,163.,195.,185.,145.,113.,45.,25.2,21.6,15.6
     8,15.2,19.5,22.8,24.5,27.6,36.7,41.,39.,32.3,28.9,27.7,24.9,23.5,  
     91.9,2.3,3.5,5.5,9.,14.,28.,60.,163.,195.,184.9,144.8,112.8,44.4,  
     A23.2,18.6,10.8,7.7,9.,10.2,11.3,13.5,16.9,19.,16.9,12.6,5.7,5.6,  
     B4.9,4.,20.,20.,19.,18.,17.,16.,12., 6., 2., 0.5,20*0./            
      DATA S3/                                                          
     13*0.,.05,.1,.3,.6,1.2,2.2,3.2,3.4,3.6,3.7,3.8,3*3.9,5*4.,3.9,3.8, 
     23.5,3.2,2.9,2.7,2.4,2.1,3*0.,.4,.8,1.4,2.3,4.4,8.,10.8,15.,16.1,  
     316.6,17.,17.3,17.5,17.6,17.7,17.8,17.9,17.8,17.5,16.8,16.,13.7,   
     412.6,11.6,10.8,9.2,8.,3*0.,.3,.6,1.,1.6,2.7,5.,6.3,6.9,7.2,7.4,7.5
     5,7.6,7.7,7.8,8.,3*8.1,8.,7.9,7.7,7.,6.6,6.2,5.8,4.6,3.5,4*0.,.1,.2
     6,.5,.8,1.8,2.1,2.4,2.6,2.8,3.,3.2,3.3,3.4,3.65,3.8,3.85,3.9,2*4., 
     73.9,3.5,3.2,2.7,2.1,1.,.6,4*0.,.2,.6,1.2,2.1,3.3,4.1,4.9,6.,7.7,  
     88.9,10.,10.2,10.1,2*9.8,10.1,10.3,8.5,6.9,5.5,3.9,3.3,3.,2.7,2.2, 
     92.1,4*0.,.1,.3,.5,.7,.9,1.1,1.3,1.5,1.8,2.2,2*2.6,2*2.4,3.,3.4,3.5
     A,3.1,2.7,2.4,2.1,1.9,1.8,1.7,1.5,1.3,3*0.,.2,.7,1.4,3.,4.9,5.1,5.,
     B4.4,3.6,3.1,2.7,2.4,2.2,2.,1.7,1.6,1.4,1.3,1.1,.9,.8,3*.6,.5,.4,.3
     C,3*0.,.2,.5,.8,1.3,2.1,3.4,4.7,6.4,6.8,6.2,5.2,7.1,8.,8.9,9.8,10.3
     D,10.5,10.6,9.9,7.7,5.5,3.4,2.8,2.7,2.6,2.2,1.9,0.,.1,.5,1.2,3.,4.2
     E,5.,5.2,6.1,7.9,8.9,9.6,9.9,10.,2*10.1,10.,9.8,9.5,8.8,8.,6.8,5.9,
     F5.2,4.2,2.9,2.7,2.5,2.2,1.9/                                      
      DATA S4/                                                          
     1750.,310.,270.,213.,172.,156.,132.,116.,108.,102.,98.,96.,96.,98.,
     2102.,109.,116.,122.,128.,132.,136.,138.,3*140.,141.,140.,138.,134.
     3,132.,740.,204.,165.,117.,87.,75.,55.,40.,30.,24.,22.,3*20.,21.,  
     423.,26.,28.,31.,33.,34.,7*36.,34.,32.,532.,266.,228.,178.,148.,   
     5136.,115.,102.,94.,86.,79.,2*76.,78.,81.,85.,91.,100.,107.,112.,  
     6115.,5*118.,116.,112.,110.,106.,520.,133.,105.,68.,46.,38.,28.,   
     722.,19.,15.,11.,10.,9.5,10.,11.,12.,15.,17.,19.,21.,23.,3*24.,23.,
     82*22.,21.,2*20./                                                  
      DATA S5/                                                          
     13*0.,.2,.8,1.6,3.2,6.2,29.,32.,34.,3*35.,35.,34.,32.8,31.3,30.4   
     2,29.0,27.3,26.2,24.8,24.0,23.0,22.4,22.,21.5,21.,18.,3*0.,.3,1.2  
     3,2.2,3.4,4.9,9.,11.2,13.4,15.6,18.1,20.3,22.3,25.4,27.1,27.9,28.3,
     428.6,28.8,29.,29.1,29.2,6*29.3,0.,1.,30.,41.,43.2,43.6,43.4,43.,  
     541.,39.,32.,23.5,16.,12.,9.4,5.6,3.8,2.9,2.2,1.8,1.6,1.3,1.,.8,   
     7.7,.5,.4,.4,.3,0./                                                
      DATA A1/                                                          
     10.,.01,.02,.03,.04,.05,.07,.1,.15,.2,.25,.3,.35,.4,.5,.65,.85,.95,
     21.1,1.3,1.5,2.,3.,4.,5.,7.,10.,16.,22.,30.,0.,0.01,.02,.03,.04,   
     3.05,.075,.1,.15,.175,.2,.225,.25,.35,.45,.5,.6,.7,.8,.85,.9,1.,1.2
     4,1.3,1.4,1.6,3.,4.,10.,20.,.2,.25,.3,.35,.4,.45,.5,.55,.6,.65,.7, 
     5.75,.8,.85,.9,.95,1.,1.1,1.2,1.3,1.4,1.6,1.8,2.,2.4,2.6,2.8,3.,   
     63.5,4.,.02,.05,.06,.08,.1,.11,.13,.15,.17,.2,.25,.3,.35,.4,.45,.5,
     7.55,.6,.65,.7,.75,.8,.85,.9,1.,2.,3.,5.,10.,30.,.0185,.02,.0224,  
     8.025,.0282,.0316,.0355,.04,.05,.056,.063,.071,.08,.09,.1,.126,    
     9.157,.2,.25,.316,.4,.5,.63,.8,1.,1.26,1.57,2.,2.5,30./            
      IF(L-4) 1,1,2                                                     
    1 I=1                                                               
      GO TO 9                                                           
    2 IF(L-10) 3,3,4                                                    
    3 I=2                                                               
      GO TO 9                                                           
    4 IF(L-19) 5,5,6                                                    
    5 I=3                                                               
      GO TO 9                                                           
    6 IF(L-23) 7,7,8                                                    
    7 I=4                                                               
      GO TO 9                                                           
    8 I=5                                                               
    9 J=1                                                               
   10 IF(T-A(J,I)) 11,16,17                                             
   11 IF(J-1) 12,12,13                                                  
   12 APSIG=0.                                                          
      RETURN                                                            
   13 IF(J-29) 15,14,14                                                 
   14 APSIG=S(30,L)
      RETURN
   15 H1=S(J-1,L)                                                       
      H2=S(J,L)                                                         
      H3=S(J+1,L)                                                       
      P1=A(J-1,I)                                                       
      P2=A(J,I)                                                         
      P3=A(J+1,I)                                                       
      GO TO 18                                                          
   16 APSIG=S(J,L)                                                      
      RETURN
   17 IF(J-29) 21,21,14
   21 J=J+1                                                             
      GO TO 10                                                          
   18 D=(P2-P3)*P1**2+(P3-P1)*P2**2+(P1-P2)*P3**2                       
      X=H1*(P2-P3)+H2*(P3-P1)+H3*(P1-P2)                                
      Y=(H2-H3)*P1**2+(H3-H1)*P2**2+(H1-H2)*P3**2                       
      Z=(P2*H3-P3*H2)*P1**2+(P3*H1-P1*H3)*P2**2+(P1*H2-P2*H1)*P3**2     
      APSIG =(X/D)*T**2+(Y/D)*T+Z/D                                     
      RETURN                                                            
      END  
************************************************************************
      SUBROUTINE BLPINN(JTY,JK,RN)                                             
************************************************************************
*  PURPOSE : NUCLEONS COORDINATES GENERATION                           *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES : 
*     JTY       - PARAMETER SPECIFYING NUCLEAR MODEL:
*                 0 - DILUTE FERMI GAS MODEL;                          *
*                 1 - LATTICE MODEL                                    *
*     JK        - SERVICE INDEX                                        *
*     RN        - PARAMETER SPECIFYING MINIMUM DISTANCE BETWEEN        
*                 CENTERS OF NUCLEONS IN THE EACH NUCLEUS              *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES :                                                *
*     XXX(2,3,240)  - COORDINATES OF NUCLEONS IN THE NUCLEI.           *
*     JZ(2,240)     - ELECTRIC CHARGE OF NUCLEONS IN THE NUCLEI.       *
*********************************************************************    
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
       COMMON/LNUC/ BR(2)                                               
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON /QDEU/QD(2,3),IQU(2)                                       
      COMMON/XYZIZ/XXX(2,3,240),JZ(2,240)                               
      T=(2.*RN)**2                                                      
      DO 17 NU=1,2                                                      
      IF(A(NU)-1.) 17,17,5                                              
   5  II=A(NU)                                                          
      NZ=Z(NU)                                                          
      IQU(NU)=II                                                        
 24   CONTINUE                                                          
      DO 16 I=1,II                                                      
      IF(II.NE.2.OR.I.NE.2) GO TO 19                                    
      XXX(NU,1,2)=-XXX(NU,1,1)                                          
      XXX(NU,2,2)=-XXX(NU,2,1)                                          
      XXX(NU,3,2)=-XXX(NU,3,1)                                          
      JZ(NU,2)=0                                                        
      GO TO 17                                                          
   19 CONTINUE                                                          
   6  B1=RNDM(-1)                                                       
      RI=R(NU)*(B1**(1./3.))                                            
      IF(II.NE.3) GO TO 3                                               
      FB=EXP(-RI*RI/BR(NU)**2)                                          
      GO TO 9                                                           
 3    CONTINUE                                                          
      IF(A(NU)-2.5)1,1,2                                                
 1    IF(RNDM(-1)-0.0645) 32,32,31                                      
   31 DSTATE=-1.                                                        
      RI=8.*RNDM(-1)                                                    
      FB=.000823*EXP(-RI*RI/198.950)+.011710*EXP(-RI*RI/57.725)+        
     *.048600*EXP(-RI*RI/13.805)+.108672*EXP(-RI*RI/2.945)              
      IF(RI.LE.3.0) FB=FB-.194806*EXP(-RI*RI/0.388)                     
      FB=FB**2*(RI**2)                                                  
      RI=0.5*RI                                                         
      FB=100.*FB                                                        
      IF(2.45*RNDM(-1)-FB)10,10,31                                      
   32 DSTATE=1.                                                         
      RI=5.1*RNDM(-1)                                                   
      FB=.07383*EXP(-RI*RI/2.457)+.0017656*EXP(-RI*RI/7.798)+           
     *.0042485*EXP(-RI*RI/8.192)+.0001859*EXP(-RI*RI/32.040)            
      IF(RI.LT.4.4) FB=FB+.4422*EXP(-RI*RI/0.8045)                      
      FB=100.*FB**2*(RI**6)                                             
      RI=0.5*RI                                                         
      IF(4.5*RNDM(-1)-FB) 10,10,32
c**************FCC**************                                      
 2    IF(JTY-1) 220,221,221
221	CALL FCC(JK,NU,A(NU),Z(NU))
	GOTO 17                                                        
c*******************************
220   IF(A(NU)-12.) 7,8,8                                               
  7   BUB=BR(NU)**2                                                     
      BSS=BUB*(1.-1./A(NU))+.64                                               
      GM1=(1.-BUB/BSS)/2.                                            
      GM2=BUB*RI**2/BSS**2                                              
      GM3=EXP(-RI**2/BSS)                                               
      FSS=(1.+(Z(NU)-2.)/3.*(3.*GM1+GM2))*GM3                           
      FB=FSS/(1.+(Z(NU)-2.)*GM1)                                         
      GO TO 9                                                           
   8  FB=(1.+EXP(-AR(NU)/CR(NU)))/(1.+EXP((RI-AR(NU))/CR(NU)))          
   9  IF(RNDM(-1)-FB) 10,10,6                                           
  10  CT=1.-2.*RNDM(-1)                                                 
      IF(II.NE.3) GO TO 20                                              
      GO TO (21,22,23),I                                                
 21   R1=RI                                                             
      GO TO 20                                                          
 22   R2=RI                                                             
      XXX(NU,3,I)=RI*CT                                                 
      GO TO 16                                                          
 23   R3=RI                                                             
      R11=XXX(NU,1,1)**2+XXX(NU,2,1)**2                                 
      AD=0.5*(R3**2-R1**2-R2**2)-XXX(NU,3,1)*XXX(NU,3,2)                
      DIS=R11*(R2**2-XXX(NU,3,2)**2)-AD**2                              
      IF(DIS.LT.0.) GO TO 24                                            
      XXX(NU,1,2)=(AD*XXX(NU,1,1)+XXX(NU,2,1)*SQRT(DIS))/R11            
      XXX(NU,2,2)=(AD*XXX(NU,2,1)-XXX(NU,1,1)*SQRT(DIS))/R11            
      DO 25 J=1,3                                                       
 25   XXX(NU,J,3)=-XXX(NU,J,1)-XXX(NU,J,2)                              
      JZ(NU,1)=1                                                        
      JZ(NU,2)=NZ-1                                                     
      JZ(NU,3)=0                                                        
      DO 26 J1=1,II                                                     
      DO 27 J2=1,II                                                     
      IF(J2.EQ.J1) GO TO 27                                             
      IF((XXX(NU,1,J1)-XXX(NU,1,J2))**2+(XXX(NU,2,J1)-XXX(NU,2,J2))**2+ 
     *(XXX(NU,3,J1)-XXX(NU,3,J2))**2-T)24,24,27                         
 27   CONTINUE                                                          
 26   CONTINUE                                                          
      GO TO 16                                                          
 20   CONTINUE                                                          
      FI=6.283185*RNDM(-1)                                              
      ST=SQRT(1.-CT*CT)                                                 
      RS=RI*ST                                                          
      CF0=COS(FI)                                                       
      SF0=SIN(FI)                                                       
      XXX(NU,1,I)=RS*CF0                                                
      XXX(NU,2,I)=RS*SF0                                                
      XXX(NU,3,I)=RI*CT                                                 
      IF(II.NE.2) GO TO 18                                              
      IF(XXX(NU,1,1)**2+XXX(NU,2,1)**2+XXX(NU,3,1)**2-RN**2)1,18,18     
 18   CONTINUE                                                          
      IF(I-NZ)11,11,12                                                  
   11 JZ(NU,I)=1                                                        
      GOTO 13                                                           
   12 JZ(NU,I)=0                                                        
  13  IF(I-1) 115,115,14                                                
  14  KM=I-1                                                            
      DO 15 K=1,KM                                                      
      IF((XXX(NU,1,I)-XXX(NU,1,K))**2+(XXX(NU,2,I)-XXX(NU,2,K))**2      
     *+(XXX(NU,3,I)-XXX(NU,3,K))**2-T)6,15,15                           
  15  CONTINUE                                                          
 115  CONTINUE                                                          
      IF(II.NE.2) GO TO 16                                              
 28   B2=RNDM(-1)                                                       
      IF(DSTATE.GT.0.) GO TO 33                                         
      QI=.23*B2                                                         
      FQ=20.714*EXP(-370.595*QI*QI)+10.054*EXP(-88.625*QI*QI)+          
     *2.215*EXP(-18.904*QI*QI)-.194*EXP(-2.494*QI*QI)                   
      IF(QI.LE.0.13) FQ=FQ+9.312*EXP(-1277.260*QI*QI)                   
      FQ=FQ**2*(QI**2)                                                  
      IF(.88*RNDM(-1)-FQ)29,29,28                                       
  33  QI=.75*B2                                                         
      FQ=5.343*EXP(-5.165*QI*QI)+44.416*EXP(-15.774*QI*QI)              
      IF(QI.LE..69) FQ=FQ+60.496*EXP(-50.065*QI*QI)+                    
     *172.971*EXP(-52.592*QI*QI)                                        
      IF(QI.LE.0.36) FQ=FQ+895.549*EXP(-205.697*QI*QI)                  
      FQ=FQ**2*(QI**6)                                                  
      IF(0.2181*RNDM(-1)-FQ)29,29,28                                    
 29   IF(DSTATE.GT.0.) GO TO 34                                         
      CTNEW=CT                                                          
      STNEW=ST                                                          
      CFNEW=CF0                                                         
      SFNEW=SF0                                                         
      SNAN=RNDM(-1)                                                     
      IF(SNAN.LT..5) GO TO 35                                           
      CTNEW=-CTNEW                                                      
      SFNEW=-SFNEW                                                      
      CFNEW=-CFNEW                                                      
      GO TO 35                                                          
 34   RIQI=RI*QI                                                        
      QIRI=2.45/(2.*5.06*RIQI)                                          
      IF(QIRI.GT.1.) GO TO 24                                           
      STOLD=ASIN(QIRI)                                                  
      SNAN=RNDM(-1)                                                     
      CTOLD=COS(STOLD)                                                  
      IF(SNAN.LT..5) CTOLD=-CTOLD                                       
      ANGLE=6.28318*RNDM(-1)                                            
      SFOLD=SIN(ANGLE)                                                  
      CFOLD=COS(ANGLE)                                                  
      CALL SHPROT(CT,SF0,CF0,CTOLD,SFOLD,CFOLD,                         
     *CTNEW,SFNEW,CFNEW,STNEW)                                          
 35   CONTINUE                                                          
      QS=QI*STNEW                                                       
      QD(NU,1)=QS*CFNEW                                                 
      QD(NU,2)=QS*SFNEW                                                 
      QD(NU,3)=QI*CTNEW                                                 
 16   CONTINUE                                                          
  17  CONTINUE                                                          
      RETURN
      END
      SUBROUTINE FCC(JK,NA,A,Z)
      COMMON/XYZIZ/ XXX(2,3,240),IZ(2,240)
      DIMENSION XYZ(2,3,240),LZ(2,240)
      INTEGER NU(200,3,2),NX(9),NY(9),NZ(9),NW(9),PP(3),NN(3),PN(3)
      DATA NZ,NY/0,0,0,4,4,2,2,2,2,2,4,0,0,2,2,0,2,4/
      DATA NX,NW/2,0,4,0,2,0,2,4,2,1,2,2,2,3,1,1,3,3/
      DATA E,R,D,B1,B2,B3/.7638,.415,.71637,9.0403,-6.7589,2.91/
      SAVE XYZ,LZ
C	JK=0
C	NA=1
C	A=197
C	Z=82
	PR=Z
	NE=A-Z
	NBAR=A
	IF(JK-1)1,1,390
C
C   FCC COORDINATION ROUTINE (I=1,2 IS PROTONS,NEUTRONS)
C
 1     DO 240 I=1,2
          KOUNT=0
C
C  N=1 TO 7 ARE THE PRINCIPAL EIGENVALUE SHELLS (0 THROUGH 6 ).
C 
       DO 230 N=1,7
          JJ=-1
C
C    THE J-LOOP IS THE J-SUBSHELLS WITHIN THE N-SHELLS.
C
         DO 220 J=N*2-1,1,-2
            JJ=JJ+2
            JZ=JJ
            CALL ZVALUE(N,I,J,JZ)
            AJ=FLOAT(J)
            IF((AJ-1.0)/4.0.NE.FLOAT(INT((AJ-1.0)/4.0))) GO TO 190
            JX=(J+1)/2
            JY=JX
            CALL XVALUE(J,I,JX)
            CALL STORE(NU,KOUNT,JX,JY,JZ,I)
            GO TO 200
  190           JX=(J+4)/2
                JY=JX-2
                CALL XVALUE(J,I,JX)
                CALL STORE(NU,KOUNT,JX,JY,JZ,I)
                CALL STORE(NU,KOUNT,JX,JY,JZ,I)
  200       DO 210 M=1,INT((FLOAT(J)-1.0)/4.0)
               JX=ABS(JX)+2
               JY=ABS(JY)-2
               CALL XVALUE(J,I,JX)
               CALL STORE(NU,KOUNT,JX,JY,JZ,I)
               CALL STORE(NU,KOUNT,JX,JY,JZ,I)
  210       CONTINUE
  220     CONTINUE
  230   CONTINUE
  240 CONTINUE
C MAIN ROUTINE: first calculate 1st-3rd neighbor bonds (p-p,n-n,p-n)
      DO 260 K=1,3
         PP(K)=0
         NN(K)=0
         PN(K)=0
  260 CONTINUE
C
      DO 330 KK=1,9
        DO 280 N=1,PR-1
          DO 270 M=N+1,PR
            NZDIF=ABS(NU(N,3,1)-NU(M,3,1))
            NYDIF=ABS(NU(N,2,1)-NU(M,2,1))
            NXDIF=ABS(NU(N,1,1)-NU(M,1,1))
            IF(NZDIF.EQ.NZ(KK).AND.NYDIF.EQ.NY(KK).AND.NXDIF.EQ.
     1      NX(KK)) PP(NW(KK))=PP(NW(KK))+1
  270     CONTINUE
  280   CONTINUE
        DO 300 N=1,NE-1
          DO 290 M=N+1,NE
            NZDIF=ABS(NU(N,3,2)-NU(M,3,2))
            NYDIF=ABS(NU(N,2,2)-NU(M,2,2))
            NXDIF=ABS(NU(N,1,2)-NU(M,1,2))
            IF(NZDIF.EQ.NZ(KK).AND.NYDIF.EQ.NY(KK).AND.NXDIF.EQ.
     1      NX(KK)) NN(NW(KK))=NN(NW(KK))+1
  290     CONTINUE
  300   CONTINUE
        DO 320 N=1,NE
          DO 310 M=1,PR
          NZDIF=ABS(NU(N,3,2)-NU(M,3,1))
          NYDIF=ABS(NU(N,2,2)-NU(M,2,1))
          NXDIF=ABS(NU(N,1,2)-NU(M,1,1))
          IF(NZDIF.EQ.NZ(KK).AND.NYDIF.EQ.NY(KK).AND.NXDIF.EQ.
     1    NX(KK)) PN(NW(KK))=PN(NW(KK))+1
  310     CONTINUE
  320   CONTINUE
  330 CONTINUE
C Calculate the total Coulomb repulsion when neighbors = 0.7638 Mev
c
      Q=0.0
      DO 350 N=1,PR-1
        DO 340 M=N+1,PR
          Q=Q+E/(SQRT(FLOAT((NU(N,1,1)-NU(M,1,1))**2+(NU(N,2,1)-
     1    NU(M,2,1))**2+(NU(N,3,1)-NU(M,3,1))**2))/SQRT(8.0))**2
  340   CONTINUE
  350 CONTINUE
C

C Find the total binding energy, using only 1st and 2nd neighbor bonds

c
      BE=B1*(PP(1)+NN(1))+B2*(PP(2)+NN(2))+B3*PN(1)-Q
C
C. Calculate the RMS radius (D is used to change from lattice units to fm)
c
      S=0.0
      DO 360 N=1,PR
	XP=D*FLOAT(NU(N,1,1))
	YP=D*FLOAT(NU(N,2,1))
	ZP=D*FLOAT(NU(N,3,1))
        S=S+SQRT(XP**2+YP**2+ZP**2)+R
	XYZ(NA,1,N)=XP
	XYZ(NA,2,N)=YP
	XYZ(NA,3,N)=ZP
	LZ(NA,N)=1
  360 CONTINUE
      DO 370 N=1,NE
	XN=D*FLOAT(NU(N,1,2))
	YN=D*FLOAT(NU(N,2,2))
	ZN=D*FLOAT(NU(N,3,2))
        S=S+SQRT(XN**2+YN**2+ZN**2)+R
	NEU=PR+N
	XYZ(NA,1,NEU)=XN
	XYZ(NA,2,NEU)=YN
	XYZ(NA,3,NEU)=ZN
	LZ(NA,NEU)=0
  370 CONTINUE
	NBAR=PR+NE
      S=S/FLOAT(NBAR)
c      goto 999
  390 CT=1.-2.*RNDM(-1)
      TT=ACOS(CT)
	PI=3.14159
      FY=2.*PI*RNDM(-1)                                                    
      DO 385 L=1,NBAR
	XY=XYZ(NA,1,L)*COS(TT)+XYZ(NA,3,L)*SIN(TT)                                           
      XXX(NA,1,L)=XY*COS(FY)-XYZ(NA,2,L)*SIN(FY)                     
      XXX(NA,2,L)=XY*SIN(FY)+XYZ(NA,2,L)*COS(FY)                     
      XXX(NA,3,L)=-XYZ(NA,1,L)*SIN(TT)+XYZ(NA,3,L)*COS(TT)                                   
385	IZ(NA,L)=LZ(NA,L)
 999    continue
	GOTO 2
C
C    Print coordinates of nucleons
c
	RM=0.
	PRINT*,'  X     Y     Z'
	DO 771 K=1,NBAR
	PRINT 777,XYZ(NA,1,K),XYZ(NA,2,K),XYZ(NA,3,K)
	RM=RM+SQRT(XYZ(NA,1,K)**2+XYZ(NA,2,K)**2+XYZ(NA,3,K)**2)+R
777	FORMAT(3F6.3)
771	CONTINUE
	RMS=RM/FLOAT(NBAR)	
C
C.   Print out results.
c
      PRINT*,'       NEIBOURS   P-P BONDS   N-N BONDS   P-N BONDS'
      DO 380 K=1,3
         PRINT*,K,PP(K),NN(K),PN(K)
  380 CONTINUE
C 
      PRINT*,' The RMS radus = ',S
      PRINT*,' Total repulsion = ',Q,'  Total BE = ',BE
c      go to 250
 2	RETURN
      END
      SUBROUTINE XVALUE(J,I,JX)
      JX=JX*(-1)**((J-1)/2+I-1)
      RETURN
      END
      SUBROUTINE ZVALUE(N,I,J,JZ)
      JZ=JZ*((-1)**(N+I-2+(J-1)/2))
      RETURN
      END
      SUBROUTINE STORE(NU,KOUNT,JX,JY,JZ,I)
      INTEGER NU(200,3,2)
      KOUNT=KOUNT+2
      NU(KOUNT-1,1,I)=JX
      NU(KOUNT-1,2,I)=JY
      NU(KOUNT-1,3,I)=JZ
      NU(KOUNT,1,I)=-JX
      NU(KOUNT,2,I)=-JY
      NU(KOUNT,3,I)=JZ
      ITEMP=JX
      JX=JY
      JY=ITEMP
      RETURN
      END
      SUBROUTINE SHPROT(CT,SF,CF,CTP,SFP,CFP,CTR,SFR,CFR,STR)           
      ST=SQRT(1.-CT**2)                                                 
      STP=SQRT(1.-CTP**2)                                               
      CTR=CTP*CT-STP*CFP*ST                                             
      TEMP=1.-CTR**2                                                    
      IF(TEMP)1,1,2                                                     
 1    SFR=0.                                                            
      CFR=1.                                                            
      TEMP1=ABS(CTR)                                                    
       CTR=CTR/TEMP1                                                    
      STR=0.                                                            
      GO TO 3                                                           
 2    CONTINUE                                                          
      STR=SQRT(TEMP)                                                    
      SFR=(STP*CFP*CT*SF+STP*SFP*CF+CTP*ST*SF)/STR                      
      CFR=(STP*CFP*CT*CF-STP*SFP*SF+CTP*ST*CF)/STR                      
 3    CONTINUE                                                          
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE PAULID(P1,P2,IP1,IP2,R0X,R0Y,R0Z,N1,N2,V,NP,MV,IP,TINT,
     * NIN)
************************************************************************
*     PURPOSE: CHECKS ON PAULI PRINCIPLE FOR THE CASCADE-CASCADE       *
*              INTERACTION, TRANSFORMS KINEMATICAL CHARACTERISTICS     *
*              OF PRODUCED SECONDARIES INTO LAB SYSTEM AND REARRANGES  *
*              ARRAYS PM AND IM.                                       *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES :                                                 *
*     P1,IP1  - ARRAYS OF KINEMATICAL CHARACTERISTICS AND QUANTUM      *
*               NUMBERS OF THE FIRST INTERACTING CASCADE PARTICLE.     *
*     P2,IP2  - ARRAYS OF KINEM.CHARACTERISTICS & QUANTUM NUMBERS      *
*               OF THE SECOND INTERACTING CASCADE PARTICLE.            *
*               ABSORPSION OF MESON).                                  *
*     R0X,    - COORDINATES OF ENTRY POINT OF PROJECTILE NUCLEUS IN    *
*     R0Y,    - THE REST FRAME OF                                      *
*     R0Z     - TARGET NUCLEUS.                                        *    
*     N1      - NO OF THE FIRST CASCADE PARTNER IN ARRAYS PM AND IM.   *
*     N2      - NO OF THE SECOND CASCADE PARTNER IN ARRAYS PM AND IM.  *
*     V       - VELOCITY OF CMS OF INTERACTING PAIR OR DECAYING        *
*               RESONANCE IN THE REST FRAME OF TARGET NUCLEUS.         *
*     NP      - NUMBER OF PRODUCED SECONDARY PARTICLES                 *
*               IN ELEMENTARY  ACT.                                    *
*     TINT    - CURRENT TIME OF CASCADE PROCESS.                       *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES :                                                *
*     MV      - CURRENT NUMBER OF CASCADE PARTICLES IN THE ARRAYS      *
*               PM & IM.                                               *
*     IP      - LABEL;                                                 *
*             IP=1 - PAULI PRINCIPLE LETS THE INTERACTION;             *
*             IP=0 - PAULI PRINCIPLE  FORBIDS THE COLLISION.           *    
************************************************************************       
      parameter (max1=2000)
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
      COMMON /KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KS1,LST                   
      COMMON/PMIM/ pm(19,max1),IM(5,max1)                               
     0COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     1MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      COMMON /ACTIV/ MPA(240),MYP(max1),MYT(max1)                       
      DIMENSION P1(11),P2(11),IP1(3),IP2(3)                             
      DIMENSION B(3),V0(3),PS(3),V(3),P(3),PT(3),PN1(11),PN2(11),PL(11) 
      IR=0                                                              
      IP=0                                                              
      IF(NIN.GT.0) GOTO 22
      V0(1)=0.                                                          
      V0(2)=0.                                                          
      V0(3)=-SQRT(WT*(WT+1.88D0))/(WT+.94D0)                            
      B(1)=0.                                                           
      B(2)=0.                                                           
      B(3)=-V0(3)                                                       
C      PI + N = DELTA                                                   
      IF(NP.EQ.1) GO TO 23                                              
      PM(18,Mv+3)=PM(18,N1)
      PM(18,Mv+1)=PM(18,N2)
      IF(P1(11)-1.) 171,172,172
  171 PM(16,MV+3)=PM(16,N1)
      PM(17,MV+3)=PM(17,N1)
  172 IF(P2(11)-1.) 173,174,174
  173 PM(16,MV+1)=PM(16,N2)
      PM(17,MV+1)=PM(17,N2) 
  174 IF(NP-2) 1,1,3                                                    
   1  DO 2 L=1,19                                                       
   2  PM(L,MV+2)=PM(L,MV+3)                                             
      IM(1,MV+2)=IM(1,MV+3)                                             
      IM(2,MV+2)=IM(2,MV+3)                                             
      GO TO 233                                                         
    3 DO 4 L=1,19                                                       
      TEMP=PM(L,MV+2)                                                   
      PM(L,MV+2)=PM(L,MV+3)                                             
    4 PM(L,MV+3)=TEMP                                                   
      DO 5 L=1,2                                                        
      ITEMP=IM(L,MV+2)                                                  
      IM(L,MV+2)=IM(L,MV+3)                                             
   5  IM(L,MV+3)=ITEMP                                                  
 233  XK0=P1(1)-R0X                                                     
      YK0=P1(2)-R0Y                                                     
      ZK0=(P1(3)+V0(3)*TINT-R0Z)*(WT+.94)/.94                           
      R1=SQRT(XK0**2+YK0**2+ZK0**2)                                     
      EPR1=(R1-AR(1))/CR(1)
      IF(EPR1.GT.10.) EPR1=10.
      IF(R1-R(1)) 455,455,456                                           
 456  TFR1=0.                                                           
      GO TO 457                                                         
 455  TFR1=TR(1)*(((1.+EXP(-AR(1)/CR(1)))/(1.+EXP(EPR1)))*              
     **(2./3.))                                                         
 457  CONTINUE                                                          
      R2=SQRT(P1(1)**2+P1(2)**2+P1(3)**2)                               
      EPR2=(R2-AR(2))/CR(2)
      IF(EPR2.GT.10.) EPR2=10.
      IF(R2-R(2)) 555,555,556                                           
 556  TFR2=0.                                                           
      GO TO 557                                                         
 555  TFR2=TR(2)*(((1.+EXP(-AR(2)/CR(2)))/(1.+EXP(EPR2)))*              
     **(2./3.))                                                         
 557  CONTINUE                                                          
      IF(PM(10,MV+1))101,101,11                                         
 101  IF(PM(9,MV+1).NE..94) GO TO 11                                    
      DO 7 I=1,3                                                        
 7    PS(I)=PM(I,MV+1)                                                  
      PN1(9)=.94                                                        
      CALL CINEMA(PS,V,PN1)                                             
      T1=PN1(8)                                                         
      IF(T1-TFR2) 22,22,11                                              
 11   IF(PM(10,MV+2)) 102,102,19                                        
 102  IF(PM(9,MV+2).NE..94) GO TO 19                                    
      DO 12 I=1,3                                                       
 12   PS(I)=PM(I,MV+2)                                                  
      PN2(9)=.94                                                        
      CALL CINEMA(PS,V,PN2)                                             
      T2=PN2(8)                                                         
      IF(T2-TFR2) 22,22,16                                              
 16   T=T2                                                              
 17   PSM=SQRT(T*(T+1.88))                                              
      PS(1)=PSM*PN2(4)*PN2(7)                                           
      PS(2)=PSM*PN2(4)*PN2(6)                                           
      PS(3)=PSM*PN2(5)                                                  
      PN2(9)=.94                                                        
      CALL CINEMA(PS,V0,PN2)                                            
      T2=PN2(8)                                                         
      IF(T2-TFR1) 22,22,19                                              
 19   IF(PM(10,MV+1)) 103,103,23                                        
 103  IF(PM(9,MV+1).NE..94) GO TO 23                                    
      T=T1                                                              
 21   PSM=SQRT(T*(T+1.88))                                              
      PS(1)=PSM*PN1(4)*PN1(7)                                           
      PS(2)=PSM*PN1(4)*PN1(6)                                           
      PS(3)=PSM*PN1(5)                                                  
      PN1(9)=.94                                                        
      CALL CINEMA(PS,V0,PN1)                                            
      T1=PN1(8)                                                         
      IF(T1-TFR1) 22,22,23                                              
 22   IP=0                                                              
      MKAS(N1)=-1                                                       
      MKAS(N2)=-1                                                       
      PM(15,N1)=0.                                                      
      PM(15,N2)=0.                                                      
      RETURN                                                            
 23   CONTINUE                                                          
      DO 210 M=1,MV                                                     
      IF(MKAS(M)-N1) 212,211,212
  211 MKAS(M)=0
      PM(15,M)=0.                                                       
  212 IF(MKAS(M)-N2) 210,213,210
  213 MKAS(M)=0                                                         
      PM(15,M)=0.
  210 CONTINUE
      L=1                                                               
      NP1=NP
  226 M=MV+L                                                            
      DO 227 I=1,3                                                      
  227 P(I)=PM(I,M)                                                      
      IF(P(3)) 31,31,32
   31 PM(18,M)=PM(18,N2)
      GOTO 33
   32 PM(18,M)=PM(18,N1)
   33 PL(9)=PM(9,M)                                                     
      CALL CINEMA(P,V,PL)                                               
  228 MK1P(M)=0                                                         
      MK2P(M)=0                                                         
      MK1T(M)=0                                                         
      MK2T(M)=0                                                         
      MYP(M)=1                                                          
      MYT(M)=1                                                          
      MKAS(M)=0                                                         
      PM(15,M)=0.                                                       
      IM(3,M)=0                                                         
      IM(4,M)=0                                                         
      IM(5,M)=2
      IF(NP.EQ.1) IM(5,M)=1
      IF(NP.EQ.2.AND.L.EQ.1) IM(5,M)=IM(5,N2)
      IF(NP.EQ.2.AND.L.EQ.2) IM(5,M)=IM(5,N1)
      KAS(M)=KPT1+KYP1+KYT1+KS1+1                                       
      IF(P(3)) 111,111,112                                               
  111 DO 113 I=1,3                                                      
  113 PM(I,M)=P2(I)                                                     
      GOTO 114                                                          
  112 DO 43 I=1,3                                                       
   43 PM(I,M)=P1(I)                                                     
  114 DO 44 I=4,8                                                       
   44 PM(I,M)=PL(I)                                                     
      IF(L-NP1) 45,46,46                                                
   45 L=L+1                                                             
      GO TO 226                                                         
   46 IP=1                                                              
      MV=MV+NP1-1                                                       
      DO 51 K=1,19                                                      
  51  PM(K,N1)=PM(K,MV+1)                                               
      IM(1,N1)=IM(1,MV+1)                                               
      IM(2,N1)=IM(2,MV+1)                                               
      IM(3,N1)=IM(3,MV+1)                                               
      IM(4,N1)=IM(4,MV+1)                                               
      IM(5,N1)=IM(5,MV+1)
      MK1P(N1)=MK1P(MV+1)                                               
      MK2P(N1)=MK2P(MV+1)                                               
      MK1T(N1)=MK1T(MV+1)                                               
      MK2T(N1)=MK2T(MV+1)                                               
      NOM=MKAS(MV+1)                                                    
      MKAS(N1)=NOM                                                      
      IF(NOM.LE.0) GO TO 125                                            
      MKAS(NOM)=N1                                                      
      IF(N1.GT.NOM.OR.PM(9,N1).NE.PM(9,NOM)) GO TO 125                  
 125  CONTINUE                                                          
      KAS(N1)=KAS(MV+1)                                                 
      MYP(N1)=MYP(MV+1)                                                 
      MYT(N1)=MYT(MV+1)                                                 
      MV=MV-1                                                           
      DO 30 K=1,19                                                      
   30 PM(K,N2)=PM(K,MV+1)                                               
      IM(1,N2)=IM(1,MV+1)                                               
      IM(2,N2)=IM(2,MV+1)                                               
      IM(3,N2)=IM(3,MV+1)                                               
      IM(4,N2)=IM(4,MV+1)                                               
      IM(5,N2)=IM(5,MV+1)
      MK1P(N2)=MK1P(MV+1)                                               
      MK2P(N2)=MK2P(MV+1)                                               
      MK1T(N2)=MK1T(MV+1)                                               
      MK2T(N2)=MK2T(MV+1)                                               
      NOM=MKAS(MV+1)                                                    
      MKAS(N2)=NOM                                                      
      IF(NOM.LE.0) GO TO 126                                            
      MKAS(NOM)=N2                                                      
 126  CONTINUE                                                          
      KAS(N2)=KAS(MV+1)                                                 
      MYP(N2)=MYP(MV+1)                                                 
      MYT(N2)=MYT(MV+1)                                                 
  108 RETURN                                                            
      END                                                               
***********************************************************************
      SUBROUTINE PAULIA(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V,NP,MV,IP,OBR2,  
     * VT,NIN)
************************************************************************
*     PURPOSE: CHECKS ON PAULI PRINCIPLE FOR THE INTERACTION AND       *
*              RESONANCE DECAY IN TARGET NUCLEUS, TRANSFORMS           *
*              KINEMATICAL CHARACTERISTICS OF PRODUCED SECONDARIES INTO*
*              LAB SYSTEM AND REARRANGES ARRAYS PM,IM AND XYZ,IZ.      *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES :                                                 *
*     P1,IP1  - ARRAYS OF KINEMATICAL CHARACTERISTICS AND QUANTUM      *
*               NUMBERS OF THE INTERACTING CASCADE PARTICLE OR         *
*               DECAYING RESONANCE.                                    *
*     P2,IP2  - ARRAYS OF KINEM.CHARACTERISTICS & QUANTUM NUMBERS      *
*               OF THE TARGET NUCLEON.                                 *
*     P3,IP3  - ARRAYS OF KINEM.CHARAC. & QUANTUM NUMBERS.             *
*               OF THE SECOND NUCLEON-PARTNER (IN THE CASE             *
*               ABSORPSION OF MESON).                                  *
*     N1      - NO OF INTERACTING CASCADE PARTNER IN ARRAYS PM AND IM. *
*     N2      - NO OF INTERACTING TARGET NUCLEON IN ARRAYS XYZ         *
*               (COORNINATES) AND IZ (NUCLEON CHARGES).                *
*     N3      - NO OF NUCLEON-PARTHER FOR PION ABSORBSION IN ARRAYS    *
*               XYZ AND IZ.                                            *    
*     V       - VELOCITY OF CMS OF INTERACTING PAIR OR DECAYING        *
*               RESONANCE IN THE REST FRAME OF TARGET NUCLEUS.         *
*     NP      - NUMBER OF PRODUCED SECONDARY PARTICLES                 *
*               IN ELEMENTARY  ACT OR IN RESONANCE DECAY.              *
*     OBR2    - CUT OFF ENERGY IN THE TARGET NUCLEUS.                  *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES :                                                *
*     MV      - CURRENT NUMBER OF CASCADE PARTICLES IN THE ARRAYS 
*               PM & IM.                                               *
*     IP      - LABEL;                                                 *
*             IP=1 - PAULI PRINCIPLE LETS THE INTERACTION;             *
*             IP=0 - PAULI PRINCIPLE  FORBIDS THE COLLISION.           *    
************************************************************************       
      parameter (max1=2000)
      common/inel/inel              
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
      COMMON/PMIM/ pm(19,max1),IM(5,max1)                               
      COMMON /KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KS1,LST                   
     0COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     1MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      COMMON /ACTIV/ MPA(240),MYP(max1),MYT(max1)                       
      DIMENSION P1(11),P2(11),P3(11),IP1(3),IP2(3),IP3(3)               
      DIMENSION B(3),V0(3),PS(3),V(3),P(3),PT(3),PN1(11),PN2(11),PL(11) 
      IR=0                                                              
      IPP=0                                                             
      IF(NIN.GT.0) GOTO 23
      PM(18,Mv+1)=2.
      IbAR=1
C      PI + N = DELTA                                                   
      IF(NP.EQ.1) GOTO 106                                              
      PM(18,Mv+3)=PM(18,N1)
      IF(P1(11)-1.) 171,172,172
  171 PM(16,MV+3)=PM(16,N1)
      PM(17,MV+3)=PM(17,N1)
  172 PM(18,Mv+3)=PM(18,N1)
      V0(1)=0.                                                          
      V0(2)=0.                                                          
      V0(3)=-SQRT(WT*(WT+1.88D0))/(WT+.94D0)                            
      B(1)=0.                                                           
      B(2)=0.                                                           
      B(3)=-V0(3)                                                       
      IF(IM(4,N1).EQ.1) IR=1                                            
      R2=SQRT(P1(1)**2+P1(2)**2+P1(3)**2)                               
      EPR2=(R2-AR(2))/CR(2)
      IF(EPR2.GT.10.) EPR2=10.
      TR2=TR(2)*AC(2)/A(2)                                              
      IF(R2-2.*R(2))55,55,56                                               
   56 TFR2=0.                                                           
      GOTO 57                                                           
   55 TFR2=TR2*(((1.+EXP(-AR(2)/CR(2)))/(1.+EXP(EPR2)))*                
     **(2./3.))                                                         
   57 IF(NP-2) 1,1,3                                                    
   1  DO 2 L=1,19                                                       
   2  PM(L,MV+2)=PM(L,MV+3)                                             
      IM(1,MV+2)=IM(1,MV+3)                                             
      IM(2,MV+2)=IM(2,MV+3)                                             
c      if(r2.gt.r(2)) goto 53      
C     DECAY OF MESON RESONANCES                                         
      IbAR=IM(2,MV+1)+IM(2,MV+2)                                        
      IF(IbAR)53,53,6
    3 IF(NP-3) 15,15,16                                                 
   15 IBAR=IM(2,MV+1)+IM(2,MV+2)+IM(2,MV+3)                             
      IF(IBAR) 53,53,16                                                 
   16 DO 4 L=1,19                                                       
      TEMP=PM(L,MV+2)                                                   
      PM(L,MV+2)=PM(L,MV+3)                                             
    4 PM(L,MV+3)=TEMP                                                   
      DO 5 L=1,2                                                        
      ITEMP=IM(L,MV+2)                                                  
      IM(L,MV+2)=IM(L,MV+3)                                             
   5  IM(L,MV+3)=ITEMP                                                  
    6 IF(PM(10,MV+1)) 13,13,8                                           
  13  DO 7 L=1,3                                                        
   7  PS(L)=PM(L,MV+1)                                                  
      PN1(9)= PM(9,MV+1)                                                
      CALL CINEMA(PS,V,PN1)                                             
      IF(PN1(8)-TFR2) 23,23,8                                           
    8 IF(PM(10,MV+2)) 14,14,11                                          
   14 IF(IM(2,MV+2))  11,11,9                                           
   9  DO 10 L=1,3                                                       
  10  PS(L)=PM(L,MV+2)                                                  
      PN2(9)=PM(9,MV+2)                                                 
      CALL CINEMA(PS,V,PN2)                                             
      IF(PN2(8)-TFR2) 23,23,11                                          
   11 IF(IM(4,N1)-1) 106,53,106                                         
  106 IF(MV-1)224,224,12                                                
12    J2=AC(2)+.1                                                       
      J22=J2-1                                                          
      J12=J2                                                            
      IF(IP1(2).EQ.0.AND.IM(2,MV+2).EQ.1)IPP=1                          
      DO 93 L=1,MV                                                      
      IF(MK1T(L)-N2) 61,60,61                                           
60    MK1T(L)=-1                                                        
      MK2T(L)=0                                                         
   61 IF(IPP.EQ.1.AND.J12.EQ.N3) J2=J22                                 
      IF(MK1T(L)-J2) 63,62,63                                           
62    MK1T(L)=N2                                                        
63    IF(MK2T(L)-N2) 65,64,65                                           
64    MK2T(L)=-1                                                        
   65 IF(MK2T(L)-J2) 75,66,75                                           
66    MK2T(L)=N2                                                        
   75 IF(NP.EQ.1) GOTO 93                                               
      IF(IP1(2)) 76,76,93                                               
76    IF(IM(2,MV+2)) 93,93,77                                           
77    IF(MK1T(L)-N3) 79,78,79                                           
78    MK1T(L)=-1                                                        
      MK2T(L)=0                                                         
   79 IF(N3.EQ.J12) GOTO 93                                             
      IF(MK1T(L)-J22) 81,80,81                                          
80    MK1T(L)=N3                                                        
81    IF(MK2T(L)-N3) 83,82,83                                           
82    MK2T(L)=-1                                                        
   83 IF(MK2T(L)-J22) 93,84,93                                          
84    MK2T(L)=N3                                                        
93    CONTINUE                                                          
  224 NA1=AC(1)+.1                                                      
      J2=AC(2)+.1                                                       
      J22=J2-1                                                          
      DO 101 L=1,NA1                                                    
      IF(NK1P(L)-N2)69,68,69                                            
   68 NK1P(L)=0                                                         
      TIME(L)=-.1                                                       
   69 IF(IPP.EQ.1.AND.J12.EQ.N3) J2=J22                                 
      IF(NK1P(L)-J2) 85,70,85                                           
   70 NK1P(L)=N2                                                        
   85 IF(NP.EQ.1) GOTO 101                                              
      IF(IP1(2))94,94,101                                               
   94 IF(IM(2,MV+2))101,101,95                                          
   95 IF(NK1P(L)-N3)87,86,87                                            
86    NK1P(L)=0                                                         
      TIME(L)=-.1                                                       
   87 IF(N3.EQ.J12) GOTO 101                                            
      IF(NK1P(L)-J22)101,88,101                                         
88    NK1P(L)=N3                                                        
  101 CONTINUE                                                          
      GO TO 24                                                          
  23  IP=0                                                              
      IF(IM(4,N1)-1) 103,104,103                                        
  104 PM(14,N1)=-.1                                                     
      PM(12,N1)=-.1                                                     
      GOTO 105                                                          
  103 MK1T(N1)=-1                                                        
      MK2T(N1)=-1                                                        
  105 IM(4,N1)=0
      RETURN                                                            
  24  E2=IP2(1)                                                         
      J2=AC(2)+.1                                                       
      IF(IPP.EQ.1.AND.N3.EQ.J12) J2=J22                                 
      DO 25 L=1,3                                                       
  25  XYZ(2,L,N2)=XYZ(2,L,J2)                                           
      IZ(2,N2)=IZ(2,J2)                                                 
      ZC(2)=ZC(2)-E2                                                    
      AC(2)=AC(2)-1.                                                    
      P2M=SQRT(P2(8)*(P2(8)+1.88))                                      
      PT(1)=-P2M*P2(4)*P2(7)                                            
      PT(2)=-P2M*P2(4)*P2(6)                                            
      PT(3)=-P2M*P2(5)                                                  
      CALL BLRECN(2,PT(1),PT(2),PT(3),P2(1),P2(2),P2(3))                
      E(2)=E(2)+TFR2-P2(8)                                              
      IF(IP1(2))47,47,53                                                
   47 IF(NP.EQ.1) GOTO 53                                               
      IF(IM(2,MV+2))53,53,48                                            
   48 E(2)=E(2)+POTENN(2,A(2),P3,IP3,AR(2),CR(2),DR(2),TR2,EPS(2),VPI)
     *-EPS(2)-P3(8)                                                     
      E3=IP3(1)                                                         
      P3M=SQRT(P3(8)*(P3(8)+1.88))                                      
      P(1)=-P3M*P3(4)*P3(7)                                             
      P(2)=-P3M*P3(4)*P3(6)                                             
      P(3)=-P3M*P3(5)                                                   
      CALL BLRECN(2,P(1),P(2),P(3),P3(1),P3(2),P3(3))                   
      J2=AC(2)+.1                                                       
      IF(N3.EQ.J12) GOTO 225                                            
      DO 49 I=1,3                                                       
   49 XYZ(2,I,N3)=XYZ(2,I,J2)                                           
      IZ(2,N3)=IZ(2,J2)                                                 
  225 AC(2)=AC(2)-1.                                                    
      ZC(2)=ZC(2)-E3                                                    
   53 L=1                                                               
      NP1=NP                                                            
      VPI1=VPI                                                          
      NABSN=0                                                           
      NABS=0                                                            
  26  M=MV+L                                                            
      EL=IM(1,M)                                                        
      QL=IM(2,M)                                                        
      IF(IM(2,M).EQ.2.AND.PM(10,M).GT.0.) QL=1                          
      DO 27 I=1,3                                                       
  27  P(I)=PM(I,M)                                                      
      IF(IM(4,N1).EQ.1) GOTO 72
      IF(P(3)) 71,71,72
   71 PM(18,M)=2.
      GOTO 73
   72 PM(18,M)=PM(18,N1)
   73 PL(9)=PM(9,M)                                                     
      CALL CINEMA(P,V,PL)                                               
      TL=PL(8)                                                          
      IF(NP-2) 120,120,121                                              
  121 IF(L-2)120,120,122                                                
  122 TL0=TL                                                            
      GOTO 124                                                          
  120 TL0=TL-QL*(TFR2+EPS(2))-(1.-QL)*VPI1                              
      IF(IR.EQ.1) TL0=TL-QL*VT-(1.-QL)*VPI1
  124 IF(TL0.LE.0.) TL0=0.001                                              
      CUT=EPS(2)*QL+OBR2*EL                                             
      MK1T(M)=0                                                         
      MK2T(M)=0                                                         
      PM(15,M)=0.                                                       
      MKAS(M)=0                                                         
      KAS(M)=KPT1+KYP1+KYT1+KS1+1                                       
      IF(EL) 28,29,29                                                   
  28  CUT=0                                                             
  29  MK1P(M)=0                                                         
      MK2P(M)=0                                                         
      IF(IR) 58,58,59                                                   
   59 MYP(M)=1                                                          
      MYT(M)=0                                                          
      GOTO 50                                                           
   58 MYP(M)=1                                                          
      MYT(M)=1                                                          
   50 IM(3,M)=0                                                         
      IM(5,M)=2
      IF(NP.EQ.1) IM(5,M)=1
      IF(NP.EQ.2.AND.L.EQ.2) IM(5,M)=IM(5,N1)   
      IF(NP.EQ.1) GOTO 42                                               
      IM(4,M)=0                                                         
      NABS=NABS+IM(2,M)                                                 
      IF(PM(10,M).GT.0.) GOTO 42                                        
      IF(TL0) 31,31,30                                                  
  30  IF(TL0-CUT) 31,31,42                                              
 31   IF(AC(2).LT.2..OR.ZC(2).LT.1.) GO TO 42                           
      IF(QL) 136,136,135                                                
 136  IF(ZC(2).EQ.AC(2).AND.EL.GT.0.) GO TO 42                          
      IF(ZC(2).EQ.0..AND.EL.LT.0.) GO TO 42                             
      NA2=AC(2)                                                         
      DO 130 I=1,NA2                                                    
      IF(EL) 131,135,132                                                
 131  IF(IZ(2,I)) 130,130,133                                           
 133  IZ(2,I)=0                                                         
      GO TO 135                                                         
 132  IF(IZ(2,I)) 134,134,130                                           
 134  IZ(2,I)=1                                                         
      GO TO 135                                                         
 130  CONTINUE                                                          
 135  CONTINUE                                                          
      E(2)=E(2)+TL0+EPS(2)*QL+(.14+VPI)*(1.-QL)                         
      CALL BLRECN(2,PL(1),PL(2),PL(3),P2(1),P2(2),P2(3))                
      NABSN=NABSN+IM(2,M)                                               
      AC(2)=AC(2)+QL                                                    
      ZC(2)=ZC(2)+EL                                                    
      IF(IM(2,M)-1) 38,32,38                                            
   32 IF(NABSN-1) 33,33,35                                              
  33  J2=AC(2)+.1                                                       
      DO 34 I=1,3                                                       
  34  XYZ(2,I,J2)=P2(I)                                                 
      GO TO 37                                                          
  35  J2=AC(2)+.1                                                       
      DO 36 I=1,3                                                       
  36  XYZ(2,I,J2)=P1(I)                                                 
  37  IZ(2,J2)=IM(1,M)                                                  
  38  IF(L-NP1) 39,41,41                                                
  39  DO 40 K=1,19                                                      
  40  PM(K,M)=PM(K,MV+NP1)                                              
      IM(1,M)=IM(1,MV+NP1)                                              
      IM(2,M)=IM(2,MV+NP1)                                              
      IM(3,M)=IM(3,MV+NP1)                                              
      IM(4,M)=IM(4,MV+NP1)                                              
      IM(5,M)=IM(5,MV+NP1)
      NP1=NP1-1                                                         
      GO TO 26                                                          
  41  NP1=NP1-1                                                         
      GO TO 46                                                          
  42  IF(NP.EQ.1) GOTO 111
      IF(IbAR.EQ.0) GOTO 112
      IF(P(3))111,111,112
  111 DO 113 I=1,3                                                      
  113 PM(I,M)=P2(I)                                                     
      GOTO 114                                                          
  112 DO 43 I=1,3                                                       
  43  PM(I,M)=P1(I)                                                     
  114 DO 44 I=4,7                                                       
   44 PM(I,M)=PL(I)                                                     
      PM(8,M)=TL0                                                       
      IF(L-NP1) 45,46,46                                                
  45  L=L+1                                                             
      GO TO 26                                                          
   46 IP=1                                                              
      MV=MV+NP1-1                                                       
      DO 51 K=1,19                                                      
  51  PM(K,N1)=PM(K,MV+1)                                               
      IM(1,N1)=IM(1,MV+1)                                               
      IM(2,N1)=IM(2,MV+1)                                               
      IM(3,N1)=IM(3,MV+1)                                               
      IM(4,N1)=IM(4,MV+1)                                               
      IM(5,N1)=IM(5,MV+1)
      MK1P(N1)=MK1P(MV+1)                                               
      MK2P(N1)=MK2P(MV+1)                                               
      MK1T(N1)=MK1T(MV+1)                                               
      MK2T(N1)=MK2T(MV+1)                                               
      KAS(N1)=KAS(MV+1)                                                 
      NOM=MKAS(MV+1)                                                    
      MKAS(N1)=NOM                                                      
      IF(NOM.LE.0) GO TO 125                                            
      MKAS(NOM)=N1                                                      
 125  CONTINUE                                                          
      MYP(N1)=MYP(MV+1)                                                 
      MYT(N1)=MYT(MV+1)                                                 
107   RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE PAULIC(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V,NP,MV,TINT,     
     *R0X,R0Y,R0Z,IP,OBR1,OBR2,VT,VP,NIN)                               
************************************************************************
*     PURPOSE: CHECKS ON PAULI PRINCIPLE FOR THE INTERACTION           *
*              OF PROJECTILE NUCLEON WITH TARGET NUCLEON, TRANSFORMS   *
*              KINEMATICAL CHARACTERISTICS OF PRODUCED SECONDARIES INTO*
*              LAB SYSTEM AND REARRANGES ARRAYS PM,IM AND XYZ,IZ.      *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES :                                                 *
*     P1,IP1  - ARRAYS OF KINEMATICAL CHARACTERISTICS AND QUANTUM      *
*               NUMBERS OF THE INTERACTING CASCADE PARTICLE OR         *
*               DECAYING RESONANCE.                                    *
*     P2,IP2  - ARRAYS OF KINEM.CHARACTERISTICS & QUANTUM NUMBERS      *
*               OF THE TARGET NUCLEON.                                 *
*     P3,IP3  - NOT USED                                               *
*     N1      - NO OF INTERACTING CASCADE PARTNER IN ARRAYS PM AND IM. *
*     N2      - NO OF INTERACTING TARGET NUCLEON IN ARRAYS XYZ         *
*               (COORNINATES) AND IZ (NUCLEON CHARGES).                *
*     N3      - NOT USED                                               *
*     V       - VELOCITY OF CMS OF INTERACTING PAIR OF NUCLEONS        *
*               IN THE REST FRAME OF TARGET NUCLEUS.                   *
*     NP      - NUMBER OF PRODUCED SECONDARY PARTICLES                 *
*               IN ELEMENTARY  ACT.                                    *
*     R0X,    - COORDINATES OF ENTRY POINT OF PROJECTILE NUCLEUS IN    *
*     R0Y,    - THE REST FRAME OF                                      *
*     R0Z     - TARGET NUCLEUS.                                        *    
*     TINT    - CURRENT TIME OF COLLISION PROCESS.                     *   
*     OBR1    - CUT OFF ENERGY IN THE PROJECTILE NUCLEUS.              * 
*     OBR2    - CUT OFF ENERGY IN THE TARGET NUCLEUS.                  * 
*     VT      - POTENTIAL FOR INTERACTING PROJECTILE NUCLEON INSIDE    *
*               TARGET NUCLEUS.                                        *
*     VP      - POTENTIAL FOR INTERACTING TARGET NUCLEON INSIDE        *
*               PROJECTILE NUCLEUS.                                    *  
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES :                                                *
*     MV      - CURRENT NUMBER OF CASCADE PARTICLES IN THE ARRAYS      *
*               PM & IM.                                               *
*     IP      - LABEL;                                                 *
*             IP=1 - PAULI PRINCIPLE LETS THE INTERACTION;             *
*             IP=0 - PAULI PRINCIPLE  FORBIDS THE COLLISION.           *
************************************************************************     
      parameter (max1=2000)
      common/inel/inel              
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
      COMMON/PMIM/ pm(19,max1),IM(5,max1)                               
      COMMON /KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KS1,LST                   
     0COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     1MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      COMMON /ACTIV/ MPA(240),MYP(max1),MYT(max1)                       
      DIMENSION P1(11),P2(11),P3(11),IP1(3),IP2(3),IP3(3)               
      DIMENSION B(3),V0(3),PS(3),V(3),PT(3),PN1(11),PN2(11),PL(11),P(11)
      IF(NIN.GT.0) GOTO 22
      V0(1)=0.                                                          
      V0(2)=0.                                                          
      V0(3)= -SQRT(WT*(WT+1.88D0))/(WT+.94D0)                           
      B(1)=0.                                                           
      B(2)=0.                                                           
      B(3)=-V0(3)                                                       
      TFR1=VP-EPS(1)                                                    
      TFR2=VT-EPS(2)                                                    
      PM(18,Mv+1)=2.
      PM(18,Mv+3)=1.
      IF(NP-2) 1,1,3                                                    
    1 DO 2 L=1,19                                                       
    2 PM(L,MV+2)=PM(L,MV+3)                                             
      IM(1,MV+2)=IM(1,MV+3)                                             
      IM(2,MV+2)=IM(2,MV+3)                                             
      GO TO 6                                                           
    3 DO 4 L=1,19                                                       
      TEMP=PM(L,MV+2)                                                   
      PM(L,MV+2)=PM(L,MV+3)                                             
    4 PM(L,MV+3)=TEMP                                                   
      DO 5 L=1,2                                                        
      ITEMP=IM(L,MV+2)                                                  
      IM(L,MV+2)=IM(L,MV+3)                                             
    5 IM(L,MV+3)=ITEMP                                                  
    6 IF(PM(10,MV+1))101,101,11                                         
  101 DO 7 I=1,3                                                        
    7 PS(I)=PM(I,MV+1)                                                  
      PN1(9)=PM(9,MV+1)                                            
      CALL CINEMA(PS,V,PN1)                                             
      T1=PN1(8)                                                         
      DO 8 I=1,3                                                        
    8 PS(I)=PN1(I)                                                      
      P(9)=PN1(9)                                                          
      CALL CINEMA(PS,V0,P)                                              
      T=P(8)                                                            
      T=T-VP                                                            
      IF(T) 10,9,9                                                      
    9 PSM=SQRT(T*(T+1.88))                                              
      PS(1)=PSM*P(4)*P(7)                                               
      PS(2)=PSM*P(4)*P(6)                                               
      PS(3)=PSM*P(5)                                                    
      P(9)=PN1(9)                                                          
      CALL CINEMA(PS,B,P)                                               
      T=P(8)                                                            
      IF(T-TFR2) 22,22,11                                               
   10 IF(T+EPS(1)) 22,22,11                                             
   11 IF(PM(10,MV+2)) 102,102,19                                        
  102 DO 12 I=1,3                                                       
   12 PS(I)=PM(I,MV+2)                                                  
      PN2(9)=PM(9,MV+2)                                               
      CALL CINEMA(PS,V,PN2)                                             
      T2=PN2(8)                                                         
      DO 13 I=1,3                                                       
   13 PS(I)=PN2(I)                                                      
      P(9)=PN2(9)                                                          
      CALL CINEMA(PS,V0,P)                                              
      T=P(8)                                                            
      T=T-VP                                                            
      IF(T) 15,14,14                                                    
   14 PSM=SQRT(T*(T+1.88))                                              
      PS(1)=PSM*P(4)*P(7)                                               
      PS(2)=PSM*P(4)*P(6)                                               
      PS(3)=PSM*P(5)                                                    
      P(9)=PN2(9)                                                          
      CALL CINEMA(PS,B,P)                                               
      T=P(8)                                                            
      IF(T-TFR2) 22,22,16                                               
   15 IF(T+EPS(1)) 22,22,16                                             
   16 T=T2-VT                                                           
      IF(T) 18,18,17                                                    
   17 PSM=SQRT(T*(T+1.88))                                              
      PS(1)=PSM*PN2(4)*PN2(7)                                           
      PS(2)=PSM*PN2(4)*PN2(6)                                           
      PS(3)=PSM*PN2(5)                                                  
      CALL CINEMA(PS,V0,PN2)                                            
      T2=PN2(8)                                                         
      IF(T2-TFR1) 22,22,19                                              
   18 IF(T+EPS(2)) 22,22,19                                             
   19 IF(PM(10,MV+1)) 103,103,23                                        
  103 T=T1-VT                                                           
      IF(T) 20,20,21                                                    
   20 IF(T+EPS(2)) 22,22,23                                             
   21 PSM=SQRT(T*(T+1.88))                                              
      PS(1)=PSM*PN1(4)*PN1(7)                                           
      PS(2)=PSM*PN1(4)*PN1(6)                                           
      PS(3)=PSM*PN1(5)                                                  
      CALL CINEMA(PS,V0,PN1)                                            
      T1=PN1(8)                                                         
      IF(T1-TFR1) 22,22,23                                              
   22 IP=0                                                              
      NK1P(N1)=-1                                                        
      TIME(N1)=-.1                                                      
      RETURN                                                            
   23 IF(MV) 37,37,24                                                   
24    J1=AC(1)+.1                                                       
      J2=AC(2)+.1                                                       
      DO 85 L=1,MV                                                      
      IF(MK1P(L)-N1) 71,70,71                                           
70    MK1P(L)=-1                                                        
      MK2P(L)=0                                                         
71    IF(MK1P(L)-J1) 73,72,73                                           
72    MK1P(L)=N1                                                        
73    IF(MK1T(L)-N2) 75,74,75                                           
74    MK1T(L)=-1                                                        
      MK2T(L)=0                                                         
75    IF(MK1T(L)-J2) 77,76,77                                           
   76 MK1T(L)=N2                                                        
77    IF(MK2P(L)-N1) 79,78,79                                           
78    MK2P(L)=-1                                                        
79    IF(MK2P(L)-J1) 81,80,81                                           
80    MK2P(L)=N1                                                        
81    IF(MK2T(L)-N2) 83,82,83                                           
82    MK2T(L)=-1                                                        
83    IF(MK2T(L)-J2) 85,84,85                                           
84    MK2T(L)=N2                                                        
85    CONTINUE                                                          
   37 NA1=AC(1)+.1                                                      
      J2=AC(2)+.1                                                       
      DO 32 L=1,NA1                                                     
      IF(NK1P(L)-N2) 26,25,26                                           
25    NK1P(L)=0                                                         
      TIME(L)=-.1                                                       
26    IF(NK1P(L)-J2) 32,27,32                                           
27    NK1P(L)=N2                                                        
32    CONTINUE                                                          
      E1=IP1(1)                                                         
      E2=IP2(1)                                                         
      J2=AC(2)+.1                                                       
      DO 38 I=1,3                                                       
   38 XYZ(2,I,N2)=XYZ(2,I,J2)                                           
      IZ(2,N2)=IZ(2,J2)                                                 
      ZC(2)=ZC(2)-E2                                                    
      AC(2)=AC(2)-1.                                                    
      AM1=P1(9)-VP                                                      
      VT1=POTENN(2,A(2),P1,IP1,AR(2),CR(2),DR(2),TR(2),EPS(2),VPI)      
      PSM=SQRT((P1(8)-VT1)*(P1(8)-VT1+2.*AM1))                          
      PS(1)=PSM*P1(4)*P1(7)                                             
      PS(2)=PSM*P1(4)*P1(6)                                             
      PS(3)=PSM*P1(5)                                                   
      P(9)=AM1                                                          
      CALL CINEMA(PS,V0,P)                                              
      T=P(8)                                                            
      E(1)=E(1)+.94-AM1-EPS(1)-T                                        
      DO 39 I=1,3                                                       
   39 P(I)=-P(I)                                                        
      CALL BLRECN(1,P(1),P(2),P(3),XYZ(1,1,N1),XYZ(1,2,N1),XYZ(1,3,N1)) 
      AM2=P2(9)-VT                                                      
      PSM=SQRT(P2(8)*(P2(8)+2.*AM2))                                    
      PS(1)=PSM*P2(4)*P2(7)                                             
      PS(2)=PSM*P2(4)*P2(6)                                             
      PS(3)=PSM*P2(5)                                                   
      P(9)=AM2                                                          
      CALL CINEMA(PS,V0,P)                                              
      T=P(8)                                                            
      DO 40 I=1,10                                                      
   40 P3(I)=P2(I)                                                       
      P3(1)=P3(1)-R0X                                                   
      P3(2)=P3(2)-R0Y                                                   
      P3(3)=(P3(3)-B(3)*TINT-R0Z)*(WT+.94)/.94                          
      VP2=POTENN(1,A(1),P3,IP2,AR(1),CR(1),DR(1),TR(1),EPS(1),VPI)      
      PSM=SQRT((T-VP2)*(T-VP2+2.*AM2))                                  
      PS(1)=PSM*P(4)*P(7)                                               
      PS(2)=PSM*P(4)*P(6)                                               
      PS(3)=PSM*P(5)                                                    
      P(9)=AM2                                                          
      CALL CINEMA(PS,B,P)                                               
      T=P(8)                                                            
      E(2)=E(2)+.94-AM2-EPS(2)-T                                        
      DO 41 I=1,3                                                       
   41 P(I)=-P(I)                                                        
      CALL BLRECN(2,P(1),P(2),P(3),P2(1),P2(2),P2(3))                   
      J1=AC(1)+.1                                                       
      DO 42 I=1,3                                                       
   42 XYZ(1,I,N1)=XYZ(1,I,J1)                                           
      IZ(1,N1)=IZ(1,J1)                                                 
      ZC(1)=ZC(1)-E1                                                    
      AC(1)=AC(1)-1.                                                    
      NK1P(N1)=NK1P(J1)                                                 
      MPA(N1)=MPA(J1)                                                   
      TIME(N1)=TIME(J1)                                                 
      L=1                                                               
      NP1=NP                                                            
   43 M=MV+L                                                            
      EL=IM(1,M)                                                        
      QL=IM(2,M)                                                        
      IF(QL.GT.1.) QL=1.                                                
      DO 44 I=1,3                                                       
   44 PS(I)=PM(I,M)                                                     
      PSQ=PS(3)
      IF(PS(3)) 97,97,95
   97 PM(18,M)=2.
      GOTO 96
   95 PM(18,M)=1.
   96 PL(9)=PM(9,M)                                                     
      CALL CINEMA(PS,V,PL)                                              
      TL=PL(8)                                                          
      DO 45 I=1,3                                                       
   45 PS(I)=PL(I)                                                       
      P(9)=PM(9,M)                                                      
      CALL CINEMA(PS,V0,P)                                              
      T=P(8)                                                            
      T=T-VP*QL-(1.-QL)*VPI                                             
      IF(T.LE.0.) T=.001
   46 PSM=SQRT(T*(T+2.*PM(9,M)))                                        
      PS(1)=PSM*P(4)*P(7)                                               
      PS(2)=PSM*P(4)*P(6)                                               
      PS(3)=PSM*P(5)                                                    
      PN1(9)=PM(9,M)                                                    
      CALL CINEMA(PS,B,PN1)                                             
      TL1=PN1(8)                                                        
      TL0=TL1-VT*QL-(1.-QL)*VPI                                         
      IF(TL0.LE.0.) TL0=.002                                            
      IF(EL) 47,48,48                                                   
   47 CUT=0.                                                            
      GO TO 49                                                          
   48 CUT=EPS(2)*QL+OBR2*EL                                             
   49 IF(TL0-CUT) 50,50,53                                              
   50 IF(AC(2).LT.2..OR.ZC(2).LT.1.) GOTO 64                            
      IF(QL) 136,136,135                                                
 136  IF(ZC(2).EQ.AC(2).AND.EL.GT.0.) GO TO 64                          
      IF(ZC(2).EQ.0..AND.EL.LT.0.) GO TO 64                             
      NA2=AC(2)                                                         
      DO 130 I=1,NA2                                                    
      IF(EL) 131,135,132                                                
 131  IF(IZ(2,I)) 130,130,133                                           
 133  IZ(2,I)=0                                                         
      GO TO 135                                                         
 132  IF(IZ(2,I)) 134,134,130                                           
 134  IZ(2,I)=1                                                         
      GO TO 135                                                         
 130  CONTINUE                                                          
 135  CONTINUE                                                          
      E(2)=E(2)+TL0+EPS(2)*QL+(.14+VPI)*(1.-QL)                         
      CALL BLRECN(2,PN1(1),PN1(2),PN1(3),P1(1),P1(2),P1(3))             
      AC(2)=AC(2)+QL                                                    
      ZC(2)=ZC(2)+EL                                                    
      IF(IM(2,M)-1) 60,51,60                                            
   51 J2=AC(2)+.1                                                       
      DO 52 I=1,3                                                       
   52 XYZ(2,I,J2)=P2(I)                                                 
      IZ(2,J2)=IM(1,M)                                                  
      GO TO 60                                                          
   53 CONTINUE                                                          
      T=TL-VT*QL-(1.-QL)*VPI                                            
      IF(T) 50,50,54                                                    
  54  PSM=SQRT(T*(T+2.*PM(9,M)))                                        
      PS(1)=PSM*PL(4)*PL(7)                                             
      PS(2)=PSM*PL(4)*PL(6)                                             
      PS(3)=PSM*PL(5)                                                   
      P(9)=PM(9,M)                                                      
      CALL CINEMA(PS,V0,P)                                              
      T=P(8)                                                            
      T=T-VP*QL-(1.-QL)*VPI                                             
      IF(EL) 55,56,56                                                   
  55  CUT=0.                                                            
      GO TO 57                                                          
  56  CUT=EPS(1)*QL+OBR1*EL                                             
  57  IF(T-CUT) 58,58,64                                                
   58 IF(AC(1).LT.2..OR.ZC(1).LT.1.) GOTO 64                            
      IF(PM(10,M).GT.0.) GO TO 64                                       
      IF(QL) 146,146,145                                                
  146 IF(ZC(1).EQ.AC(1).AND.EL.GT.0.) GO TO 64                          
      IF(ZC(1).EQ.0..AND.EL.LT.0.) GO TO 64                             
      NA1=AC(1)                                                         
      DO 140 I=1,NA1                                                    
      IF(EL) 141,145,142                                                
 141  IF(IZ(1,I)) 140,140,143                                           
 143  IZ(1,I)=0                                                         
      GO TO 145                                                         
 142  IF(IZ(1,I)) 144,144,140                                           
 144  IZ(1,I)=1                                                         
      GO TO 145                                                         
 140  CONTINUE                                                          
 145  CONTINUE                                                          
      E(1)=E(1)+T+EPS(1)*QL+(.14+VPI)*(1.-QL)                           
      X=P1(1)-R0X                                                       
      Y=P1(2)-R0Y                                                       
      ZT=(P1(3)-B(3)*TINT-R0Z)*(WT+.94)/.94                             
      CALL BLRECN(1,P(1),P(2),P(3),X,Y,ZT)                              
      AC(1)=AC(1)+QL                                                    
      ZC(1)=ZC(1)+EL                                                    
      IF(IM(2,M)-1) 60,59,60                                            
  59  J1=AC(1)+.1                                                       
      XYZ(1,1,J1)=X                                                     
      XYZ(1,2,J1)=Y                                                     
      XYZ(1,3,J1)=ZT                                                    
      IZ(1,J1)=IM(1,M)                                                  
      NK1P(J1)=0                                                        
      TIME(J1)=-.1                                                      
      MPA(J1)=1                                                         
   60 IF(L-NP1) 61,63,63                                                
   61 DO 62 K=1,19                                                      
   62 PM(K,M)=PM(K,MV+NP1)                                              
      IM(1,M)=IM(1,MV+NP1)                                              
      IM(2,M)=IM(2,MV+NP1)                                              
      IM(3,M)=IM(3,MV+NP1)                                              
      IM(4,M)=IM(4,MV+NP1)                                              
      IM(5,M)=IM(5,MV+NP1)
      NP1=NP1-1                                                         
      GO TO 43                                                          
   63 NP1=NP1-1                                                         
      GO TO 68                                                          
   64 IF(PSQ)69,69,91
   69 DO 65 I=1,3                                                       
   65 PM(I,M)=P2(I)                                                     
      GOTO 92                                                           
   91 DO 93 I=1,3                                                       
   93 PM(I,M)=P1(I)                                                     
   92 DO 66 I=4,7                                                       
   66 PM(I,M)=PN1(I)                                                    
      PM(8,M)=TL0                                                       
      MK1P(M)=0                                                         
      MK2P(M)=0                                                         
      MK1T(M)=0                                                         
      MK2T(M)=0                                                         
      MYP(M)=1                                                          
      MYT(M)=1                                                          
      MKAS(M)=0                                                         
      KAS(M)=KPT1+KYP1+KYT1+KS1+1                                       
      PM(15,M)=0.                                                       
      IM(3,M)=0                                                         
      IM(4,M)=0                                                         
      IM(5,M)=2
      IF(L-NP1) 67,68,68                                                
   67 L=L+1                                                             
      GO TO 43                                                          
   68 IP=1                                                              
      MV=MV+NP1                                                         
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE PAULIB(P1,P2,P3,IP1,IP2,IP3,N1,N2,N3,V,NP,MV,          
     *TINT,R0X,R0Y,R0Z,IP,OBR1,VP,NIN)                                     
************************************************************************
*     PURPOSE: CHECKS ON PAULI PRINCIPLE FOR THE INTERACTION AND       *
*              RESONANCE DECAY INSIDE PROJECTILE NUCLEUS, TRANSFORMS   *
*              KINEMATICAL CHARACTERISTICS OF PRODUCED SECONDARIES INTO*
*              LAB SYSTEM AND REARRANGES ARRAYS PM,IM AND XYZ,IZ.      *
*----------------------------------------------------------------------*
*   INPUT QUANTITIES :                                                 *
*     P1,IP1  - ARRAYS OF KINEMATICAL CHARACTERISTICS AND QUANTUM      *
*               NUMBERS OF THE INTERACTING CASCADE PARTICLE OR         *
*               DECAYING RESONANCE.                                    *
*     P2,IP2  - ARRAYS OF KINEM.CHARACTERISTICS & QUANTUM NUMBERS      *
*               OF THE PROJECTILE NUCLEON.                             *
*     P3,IP3  - ARRAYS OF KINEM.CHARAC. & QUANTUM NUMBERS.             *
*               OF THE SECOND NUCLEON-PARTNER (IN THE CASE             *
*               ABSORPSION OF MESON).                                  *
*     N1      - NO OF INTERACTING CASCADE PARTNER IN ARRAYS PM AND IM. *
*     N2      - NO OF INTERACTING PROJECTILE NUCLEON IN ARRAYS XYZ     *
*               (COORNINATES) AND IZ (NUCLEON CHARGES).                *
*     N3      - NO OF NUCLEON-PARTHER FOR PION ABSORBSION IN ARRAYS    *
*               XYZ AND IZ.                                            *    
*     V       - VELOCITY OF CMS OF INTERACTING PAIR OR DECAYING        *
*               RESONANCE IN THE REST FRAME OF PROJECTILE NUCLEUS.     *
*     NP      - NUMBER OF PRODUCED SECONDARY PARTICLES                 *
*               IN ELEMENTARY  ACT OR IN RESONANCE DECAY.              *
*     R0X,    - COORDINATES OF ENTRY POINT OF PROJECTILE NUCLEUS IN    *
*     R0Y,    - THE REST FRAME OF                                      *
*     R0Z     - TARGET NUCLEUS.                                        *
*     TINT    - CURRENT TIME OF COLLISION PROCESS.                     *   
*     OBR1    - CUT OFF ENERGY IN THE PROJECTILE NUCLEUS.              *
*----------------------------------------------------------------------*
*   OUTPUT QUANTITIES :                                                *
*     MV      - CURRENT NUMBER OF CASCADE PARTICLES IN THE ARRAYS      *
*               PM & IM.                                               *
*     IP      - LABEL;                                                 *
*             IP=1 - PAULI PRINCIPLE LETS THE INTERACTION;             *
*             IP=0 - PAULI PRINCIPLE  FORBIDS THE COLLISION.           *    
************************************************************************     
      parameter (max1=2000)
      common/inel/inel              
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
      COMMON/PMIM/ pm(19,max1),IM(5,max1)                               
      COMMON /KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KS1,LST                   
     0COMMON/CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),          
     1MK1T(max1),MK2T(max1),MKAS(max1),KAS(max1)                        
      COMMON/ACTIV/MPA(240),MYP(max1),MYT(max1)                         
      DIMENSION P1(11),P2(11),P3(11),IP1(3),IP2(3),IP3(3)               
      DIMENSION B(3),V0(3),PS(3),V(3),P(3),PT(3),PN1(11),PN2(11),PL(11) 
      IR=0                                                              
      IPP=0                                                             
      IF(NIN.GT.0) GOTO 23
      IbAR=1
      PM(18,Mv+1)=1.
      V0(1)=0.                                                          
      V0(2)=0.                                                          
      V0(3)=-SQRT(WT*(WT+1.88D0))/(WT+.94D0)                            
      B(1)=0.                                                           
      B(2)=0.                                                           
      B(3)=-V0(3)                                                       
C     PI + N = DELTA                                                    
      IF(NP.EQ.1) GOTO 106                                              
      IF(P1(11)-1.) 171,172,172
  171 PM(16,MV+3)=PM(16,N1)
      PM(17,MV+3)=PM(17,N1)
  172 PM(18,Mv+3)=PM(18,N1)
      R1=SQRT(P1(1)**2+P1(2)**2+P1(3)**2)                               
      EPR1=(R1-AR(1))/CR(1)
      IF(EPR1.GT.10.) EPR1=10.
      TR1=TR(1)*AC(1)/A(1)                                              
      IF(IM(3,N1).EQ.1) IR=1                                            
      IF(R1-2.*R(1))55,55,56                                               
   56 TFR1=0.                                                           
      GOTO 57                                                           
   55 TFR1=TR1*(((1.+EXP(-AR(1)/CR(1)))/(1.+EXP(EPR1)))*                
     **(2./3.))                                                         
   57 IF(NP-2) 1,1,3                                                    
   1  DO 2 L=1,19                                                       
   2  PM(L,MV+2)=PM(L,MV+3)                                             
      IM(1,MV+2)=IM(1,MV+3)                                             
      IM(2,MV+2)=IM(2,MV+3)                                             
c      if(r1.gt.r(1)) goto 53
C     DECAY OF MESON RESONANCES                                         
      IbAR=IM(2,Mv+1)+IM(2,Mv+2)
      IF(IbAR)53,53,6
    3 IF(NP-3)15,15,16                                                  
   15 IBAR=IM(2,MV+1)+IM(2,MV+2)+IM(2,MV+3)                             
      IF(IBAR) 53,53,16                                                 
   16 DO 4 L=1,19                                                       
      TEMP=PM(L,MV+2)                                                   
      PM(L,MV+2)=PM(L,MV+3)                                             
   4  PM(L,MV+3)=TEMP                                                   
      DO 5 L=1,2                                                        
      ITEMP=IM(L,MV+2)                                                  
      IM(L,MV+2)=IM(L,MV+3)                                             
   5  IM(L,MV+3)=ITEMP                                                  
    6 IF(PM(10,MV+1)) 13,13,8                                           
   13 DO 7 L=1,3                                                        
    7 PS(L)=PM(L,MV+1)                                                  
      PN1(9)=PM(9,MV+1)                                                 
      CALL CINEMA(PS,V,PN1)                                             
      IF(PN1(8)-TFR1) 23,23,8                                           
    8 IF(PM(10,MV+2)) 14,14,11                                          
   14 IF(IM(2,MV+2)) 11,11,9                                            
    9 DO 10 L=1,3                                                       
  10  PS(L)=PM(L,MV+2)                                                  
      PN2(9)=PM(9,MV+2)                                                 
      CALL CINEMA(PS,V,PN2)                                             
      IF(PN2(8)-TFR1) 23,23,11                                          
   11 IF(IM(3,N1)-1)106,53,106                                          
  106 IF(MV-1) 24,24,12                                                 
12    J1=AC(1)+.1                                                       
      J11=J1-1                                                          
      J12=J1                                                            
      IF(IP1(2).EQ.0.AND.IM(2,MV+2).EQ.1) IPP=1                         
      DO 77 L=1,MV                                                      
      IF(MK1P(L)-N2) 61,60,61                                           
60    MK1P(L)=-1                                                        
      MK2P(L)=0                                                         
   61 IF(IPP.EQ.1.AND.J12.EQ.N3) J1=J11                                 
      IF(MK1P(L)-J1) 63,62,63                                           
62    MK1P(L)=N2                                                        
63    IF(MK2P(L)-N2) 65,64,65                                           
64    MK2P(L)=-1                                                        
65    IF(MK2P(L)-J1) 67,66,67                                           
66    MK2P(L)=N2                                                        
   67 IF(NP.EQ.1) GOTO 77                                               
      IF(IP1(2)) 68,68,77                                               
68    IF(IM(2,MV+2)) 77,77,69                                           
69    IF(MK1P(L)-N3) 71,70,71                                           
70    MK1P(L)=-1                                                        
      MK2P(L)=0                                                         
   71 IF(N3.EQ.J12) GOTO 77                                             
      IF(MK1P(L)-J11) 73,72,73                                          
72    MK1P(L)=N3                                                        
73    IF(MK2P(L)-N3) 75,74,75                                           
74    MK2P(L)=-1                                                        
75    IF(MK2P(L)-J11) 77,76,77                                          
76    MK2P(L)=N3                                                        
77    CONTINUE                                                          
      GO TO 24                                                          
  23  IP=0                                                              
      IF(IM(3,N1)-1)103,104,103                                         
  104 PM(12,N1)=-.1                                                     
      PM(14,N1)=-.1                                                     
      IM(3,N1)=0
      GOTO 105                                                          
 103  MK1P(N1)=-1                                                        
      MK2P(N1)=-1
      IM(3,N1)=0                                                        
 105  RETURN                                                            
   24 E2=IP2(1)                                                         
      J1=AC(1)+.1                                                       
      IF(IPP.EQ.1.AND.N3.EQ.J12) J1=J11                                 
      DO 25 I=1,3                                                       
  25  XYZ(1,I,N2)=XYZ(1,I,J1)                                           
      IZ(1,N2)=IZ(1,J1)                                                 
      ZC(1)=ZC(1)-E2                                                    
      AC(1)=AC(1)-1.                                                    
      NK1P(N2)=NK1P(J1)                                                 
      TIME(N2)=TIME(J1)                                                 
      MPA(N2)=MPA(J1)                                                   
      P2M=SQRT(P2(8)*(P2(8)+1.88))                                      
      PT(1)=-P2M*P2(4)*P2(7)                                            
      PT(2)=-P2M*P2(4)*P2(6)                                            
      PT(3)=-P2M*P2(5)                                                  
      CALL BLRECN(1,PT(1),PT(2),PT(3),P2(1),P2(2),P2(3))                
      E(1)=E(1)+TFR1-P2(8)                                              
      IF(NP.EQ.1) GOTO 53                                               
      IF(IP1(2))47,47,53                                                
   47 IF(IM(2,MV+2))53,53,48                                            
   48 E(1)=E(1)+POTENN(1,A(1),P3,IP3,AR(1),CR(1),DR(1),TR1,EPS(1),VPI)
     *-EPS(1)-P3(8)                                                     
      E3=IP3(1)                                                         
      P3M=SQRT(P3(8)*(P3(8)+1.88))                                      
      P(1)=-P3M*P3(4)*P3(7)                                             
      P(2)=-P3M*P3(4)*P3(6)                                             
      P(3)=-P3M*P3(5)                                                   
      CALL BLRECN(1,P(1),P(2),P(3),P3(1),P3(2),P3(3))                   
      J2=AC(1)+.1                                                       
      IF(N3.EQ.J12) GOTO 225                                            
      DO 49 I=1,3                                                       
   49 XYZ(1,I,N3)=XYZ(1,I,J2)                                           
      IZ(1,N3)=IZ(1,J2)                                                 
      NK1P(N3)=NK1P(J2)                                                 
      TIME(N3)=TIME(J2)                                                 
      MPA(N3)=MPA(J2)                                                   
  225 AC(1)=AC(1)-1.                                                    
      ZC(1)=ZC(1)-E3                                                    
   53 L=1                                                               
      NP1=NP                                                            
      NABSN=0                                                           
      NABS=0                                                            
  26  M=MV+L                                                            
      EL=IM(1,M)                                                        
      QL=IM(2,M)                                                        
      IF(QL.EQ.2..AND.PM(10,M).GT.0.)QL=1.                              
      DO 27 I=1,3                                                       
  27  P(I)=PM(I,M)                                                      
      IF(IM(4,N1).EQ.1) GOTO 92
      IF(P(3)) 91,91,92
   91 PM(18,M)=1.
      GOTO 93
   92 PM(18,M)=PM(18,N1)
   93 PL(9)=PM(9,M)                                                     
      CALL CINEMA(P,V,PL)                                               
      TL=PL(8)                                                          
      IF(NP-2) 120,120,121                                              
  121 IF(L-2)120,120,122                                                
  122 TL0=TL                                                            
      GOTO 124                                                          
  120 TL0=TL-QL*(TFR1+EPS(1))-(1.-QL)*VPI                               
      IF(IR.EQ.1) TL0=TL-QL*VP-(1-QL)*VPI
  124 IF(TL0.LE.0.) TL0=0.001                                              
      CUT=EPS(1)*QL+OBR1*EL                                             
      MK1P(M)=0                                                         
      MK2P(M)=0                                                         
      MKAS(M)=0                                                         
      PM(15,M)=0.                                                       
      KAS(M)=KPT1+KYP1+KYT1+KS1+1                                       
      IF(EL) 28,29,29                                                   
  28  CUT=0.                                                            
  29  MK1T(M)=0                                                         
      MK2T(M)=0                                                         
      IF(IR) 51,51,52                                                   
   52 MYP(M)=0                                                          
      MYT(M)=1                                                          
      GOTO 58                                                           
   51 MYP(M)=1                                                          
      MYT(M)=1                                                          
   58 IM(3,M)=0                                                         
      IM(4,M)=0                                                         
      IM(5,M)=2
      IF(NP.EQ.1) IM(5,M)=1
      IF(NP.EQ.2.AND.L.EQ.2) IM(5,M)=IM(5,N1)
      IF(NP.EQ.1) GOTO 42                                               
      NABS=NABS+IM(2,M)                                                 
      IF(PM(10,M).GT.0.) GOTO 42                                        
      IF(TL0) 31,31,30                                                  
  30  IF(TL0-CUT) 31,31,42                                              
 31   IF(AC(1).LT.2..OR.ZC(1).LT.1.) GO TO 42                           
      IF(QL) 146,146,145                                                
 146  IF(ZC(1).EQ.AC(1).AND.EL.GT.0.) GO TO 42                          
      IF(ZC(1).EQ.0..AND.EL.LT.0.) GO TO 42                             
      NA1=AC(1)                                                         
      DO 140 I=1,NA1                                                    
      IF(EL) 141,145,142                                                
 141  IF(IZ(1,I)) 140,140,143                                           
 143  IZ(1,I)=0                                                         
      GO TO 145                                                         
 142  IF(IZ(1,I)) 144,144,140                                           
 144  IZ(1,I)=1                                                         
      GO TO 145                                                         
 140  CONTINUE                                                          
 145  CONTINUE                                                          
      E(1)=E(1)+TL0+EPS(1)*QL+(.14+VPI)*(1.-QL)                         
      CALL BLRECN(1,PL(1),PL(2),PL(3),P2(1),P2(2),P2(3))                
      NABSN=NABSN+IM(2,M)                                               
      AC(1)=AC(1)+QL                                                    
      ZC(1)=ZC(1)+EL                                                    
      IF(IM(2,M)-1) 38,32,38                                            
  32  IF(NABSN-1) 33,33,35                                              
  33  J1=AC(1)+.1                                                       
      DO 34 I=1,3                                                       
  34  XYZ(1,I,J1)=P2(I)                                                 
      GO TO 37                                                          
  35  J1=AC(1)+.1                                                       
      DO 36 I=1,3                                                       
  36  XYZ(1,I,J1)=P3(I)                                                 
  37  IZ(1,J1)=IM(1,M)                                                  
      NK1P(J1)=0                                                        
      TIME(J1)=-.1                                                      
      MPA(J1)=1                                                         
  38  IF(L-NP1) 39,41,41                                                
  39  DO 40 K=1,19                                                      
  40  PM(K,M)=PM(K,MV+NP1)                                              
      IM(1,M)=IM(1,MV+NP1)                                              
      IM(2,M)=IM(2,MV+NP1)                                              
      IM(5,M)=IM(5,MV+NP1)
      NP1=NP1-1                                                         
      GO TO 26                                                          
  41  NP1=NP1-1                                                         
      GO TO 45                                                          
  42  PSM=SQRT(TL0*(TL0+2.*PM(9,M)))                                    
      PS(1)=PSM*PL(4)*PL(7)                                             
      PS(2)=PSM*PL(4)*PL(6)                                             
      PS(3)=PSM*PL(5)                                                   
      PN2(9)=PM(9,M)                                                    
      CALL CINEMA(PS,B,PN2)                                             
      TTL=PN2(8)
      IF(NP.EQ.1) GOTO 80
      IF(IbAR.EQ.0) GOTO 81
      IF(P(3))80,80,81                                                   
   80 PM(1,M)=P2(1)+R0X                                                 
      PM(2,M)=P2(2)+R0Y                                                 
      PM(3,M)=P2(3)*(.94/(WT+.94))+R0Z+B(3)*TINT                        
      GOTO 82                                                           
   81 PM(1,M)=P1(1)+R0X                                                 
      PM(2,M)=P1(2)+R0Y                                                 
      PM(3,M)=P1(3)*(.94/(WT+.94))+R0Z+B(3)*TINT                        
   82 DO 43 I=4,7                                                       
  43  PM(I,M)=PN2(I)                                                    
      PM(8,M)=TTL                                                       
      IF(L-NP1) 44,45,45                                                
  44  L=L+1                                                             
      GO TO 26                                                          
   45 IP=1                                                              
      MV=MV+NP1-1                                                       
      DO 50 K=1,19                                                      
  50  PM(K,N1)=PM(K,MV+1)                                               
      IM(1,N1)=IM(1,MV+1)                                               
      IM(2,N1)=IM(2,MV+1)                                               
      IM(3,N1)=IM(3,MV+1)                                               
      IM(4,N1)=IM(4,MV+1)                                               
      IM(5,N1)=IM(5,MV+1)
      MK1P(N1)=MK1P(MV+1)                                               
      MK2P(N1)=MK2P(MV+1)                                               
      MK1T(N1)=MK1T(MV+1)                                               
      MK2T(N1)=MK2T(MV+1)                                               
      KAS(N1)=KAS(MV+1)                                                 
      NOM=MKAS(MV+1)                                                    
      MKAS(N1)=NOM                                                      
      IF(NOM.LE.0) GO TO 125                                            
      MKAS(NOM)=N1                                                      
 125  CONTINUE                                                          
      MYP(N1)=MYP(MV+1)                                                 
      MYT(N1)=MYT(MV+1)                                                 
      RETURN                                                            
      END                                                               
