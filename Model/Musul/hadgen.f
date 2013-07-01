************************************************************************
      SUBROUTINE BINEL(P1,IP1,P2,IP2,V,U,TIN1,MV,NP,NIN)  
************************************************************************
*   PURPOSE:    BLOCK OF INELASTIC SCATTERING IN THE CASE OF           *
*               HADRON-HADRON COLLISION.                               *     
************************************************************************
      DIMENSION P1(11),IP1(3),IP2(3),V(3),P2(11)          
    1 NIN = 0                                                           
      IK = 0                                                            
   21 CALL BLACT(P1,IP1,P2,IP2,V,U,TIN1,MV,NP,NIN,        
     1 LAB1,LAB2)                                                       
      IF(NIN.GT.0) RETURN
   18 CALL CHARGE(IP1,IP2,MV,NP,LAB1,LAB2,NIN)                    
      RETURN                                                            
         END                                                            
************************************************************************
      SUBROUTINE CHARGE(IP1,IP2,MV,NP,LAB1,LAB2,NIN)              
************************************************************************
*   PURPOSE:  DETERMINES SECONDARY PARTICLES CHARGES IN INELASTIC
*     SCATTERING.                                                       
*----------------------------------------------------------------------*
*  INPUT QUANTATIES:                                                   *
*      IP1,IP2  - QUANTUM NUMBERS OF INTERACTING PARTICLES             *
*      MV       - CURRENT NUMBER OF CASCADE PARTICLES IN MEMORY.       *
*      NP       - THE NUMBER OF SECONDARIES DIRECTLY PRODUCED          *
*                  IN THE LAST INTERACTION.                            *
*      LAB1     - LABEL INDICATING THE TYPE OF THE FIRST PARTNER.      *
*      LAB2       THE SAME OF THESECOND PARTNER.                       *
************************************************************************
      parameter (max1=2000)
      DIMENSION IP1(3),IP2(3),PM(19,max1),IM(5,MAX1)        
      COMMON/PMIM/PM,IM                                           
      ME=IP1(1)+IP2(1)                                            
      MER=IM(1,MV+1)+IM(1,MV+3)                                   
      MR=ME-MER                                                         
      NS=NP-2                                                           
      IF(NP.EQ.2) GOTO 40                                               
   21 LAMBDA = 2                                                        
      NW=0                                                              
   30 IF (LAMBDA-3) 31,36,31                                            
   31 TEMP2=RNDM(-1)                                                    
       MTEMP=MV+LAMBDA                                                  
      IF(PM(9,MTEMP)-.783)11,12,12                                   
   12 IM(1,MTEMP)=0                                                  
      NW=NW+1                                                           
      GOTO 36                                                           
   11 IF (TEMP2-1./3.) 33,32,32                                         
   32 IF (TEMP2-2./3.) 34,35,35                                         
   33 IM(1,MTEMP) = 1                                                
       GO TO 36                                                         
   34 IM(1,MTEMP) = 0                                                
       GO TO 36                                                         
   35 IM(1,MTEMP) = -1                                               
       GO TO 36                                                         
   37 LAMBDA = LAMBDA+1                                                 
       GO TO 30                                                         
   36 IF (LAMBDA-NP) 37,38,38                                           
   38 SIGQ = 0.                                                         
      IF(NW-NS)15,22,22                                                 
   22 IF(MR)23,40,23                                                    
   23 NIN=1                                                             
      RETURN                                                            
   15 DO 39 I=1,NP                                                      
      ITEMP = MV+I                                                      
      SIGQ = SIGQ+IM(1,ITEMP)                                        
   39 CONTINUE                                                          
      IF(ME-SIGQ)21,40,21                                               
   40 CONTINUE                                                          
      RETURN                                                            
         END                                                            
************************************************************************
      SUBROUTINE BLACT(P1,IP1,P2,IP2,V,UU,TIN1,MV,NP,NIN,LAB1,LAB2) 
************************************************************************
*   PURPOSE: THE MAIN SUBROUTINE GENERATING INELASTIC INTERACTION;     *
*            CHARACTERISTICS OF GENERATED SECONDARIES ARE WRITTEN IN   *
*            ARRAYS 'PM' AND 'IM'.                                     *
*----------------------------------------------------------------------*
*    INPUT QUANTATIES:                                                 *
*     P1,IP1    - CHARACTERISTICS OF THE FIRST PARTNER                 *
*     P2,IP2    - CHARACTERISTICS OF THE SECOND PARTNER                *
*     V         - VELOCITY OF C.M.S IN LAB. SYSTEM                     *
*     UU        - TOTAL ENERGY OF PARTNERS IN C.M.S.                   *
*     TIN1      - KINETIC ENERGY OF THE FIRST PARTNER IN THE REST      *
*                 FRAME OF THE SECOND ONE.                             *
*     MV        - CURRENT NUMBER OF CASCADE PARTICLES IN MEMORY.       *
*----------------------------------------------------------------------*
*    OUTPUT QUANTATIES:                                                *
*     NP        - THE NUMBER OF DIRECTLY PRODUCED SECONDARIES          *
*     LAB1,LAB2 - LABELS INDICATING TYPES OF PARTNERS.                 *
************************************************************************
      parameter (max1=2000)
      COMMON/PMIM/PM,IM                                           
      COMMON/KNEY/ TMNEY,MET                                            
      DIMENSION PM(19,max1),IM(5,MAX1),P1(11),IP1(3),V(3), 
     *PL(3),PLST(3),SIGMA(3),P2(11),IP2(3),XX(max1)              
      MET=0                                                             
      NTRY=0                                                            
      MIR=0                                                             
      IF(IP1(2))19,19,18                                             
   18 PK=.45
      GOTO 25                                                           
   19 IF(IP2(2)) 118,118,119
  118 PK=.67
      GOTO 25
  119 PK=.58
   25 A=PK*6.
      B=6.-A                                                         
      SK=(A-1.)/(A+B-2.)                                                
      SNR=SK**(A-1.)*(1.-SK)**(B-1.) 
      PIPRM=.20+.8*EXP(-.25*TIN1)                                       
      PIKF=.25+.15*ALOG(UU/2.)                                          
      MQ=IP1(1)+IP2(1)                                            
      MB=IP1(2)+IP2(2)
      FMQ=FLOAT(MQ)                                                     
      IF(IP1(2))3,3,4                                                
    3 IF(P1(9).GT..14.OR.P2(9).GT..94) MIR=1                    
      GOTO 20                                                           
    4 IF(P1(9).GT..94.OR.P2(9).GT..94)MIR=1                     
   20 IF(UU-P1(9)-P2(9)-.16) 6,6,7                              
    6 IF(MIR)1,1,7                                                      
    7 CALL BLLABL(LAB1,LAB2,P1,IP1,P2,IP2)                
  335 CALL CLEAD(0,LAB1,LAB2,P1,P2,IP1,IP2,MV,NIN,UU)     
      U=UU-PM(9,MV+1)-PM(9,MV+3)                                  
      IF(NIN.GT.0) GOTO 1                                               
      IF(U)335,335,336                                                  
  336 IF(U-.16)339,338,338                                              
  339 CMB=PM(9,MV+1)                                                 
      CMF=PM(9,MV+3)                                                 
      GOTO 79                                                           
  338 THR=U/UU                                                          
   35 SS=RNDM(-1)                                                       
      IF(RNDM(-1)-SS**(A-1.)*(1.-SS)**(B-1.)/SNR) 36,36,35              
   36 TEMP2=SS                                                          
      IF(TEMP2-THR) 337,35,35                                           
  337 NTRY=1                                                            
      U2=(1.-TEMP2)*UU                                                  
       EKL=TEMP2*UU                                                     
      UMX=EKL
      IF(TEMP2.GT..6)UMX=(.6-TEMP2/5.)*UU
      IF(IP2(2).EQ.0.AND.TEMP2.GT..5) UMX=(.5-TEMP2/5.)*UU
      NP=0                                                              
   24 L1=1                                                              
       L2=3                                                             
      CMB=PM(9,MV+1)                                                 
      CMF=PM(9,MV+3)                                                 
      IF(U2-CMB-CMF)35,35,700                                           
  700 IF(EKL-.16) 77,77,70                                              
   70 EF=(U2**2+CMF**2-CMB**2)/(2.*U2)                                  
      EB=(U2**2+CMB**2-CMF**2)/(2.*U2)                                  
      if(ef.le.cmf) goto 1
      PPT=SQRT(EF**2-CMF**2)                                            
      LL=0
  680 IF(LL-10)68,335,335
   68 SS=RNDM(-1)                                                       
      SSN=5.*SS                                                         
      IF(RNDM(-1)-SSN*EXP(1.-SSN))69,69,68                        
   69 POP=PIKF*SSN                                                   
	LL=LL+1
      IF(POP.GT.PPT) GOTO 680                                           
      STP=POP/PPT                                                       
       CTP=SQRT(1.-STP**2)                                              
      CTT=-CTP                                                          
       STT=SQRT(1.-CTT**2)                                              
      PI=3.14159                                                        
       FIN=2.*PI*RNDM(-1)                                               
       FTA=FIN+PI                                                       
   62 PM(1,MV+1)=PPT*STT*COS(FTA)                                    
      PM(2,MV+1)=PPT*STT*SIN(FTA)                                    
      PM(3,MV+1)=PPT*CTT                                             
      PM(4,MV+1)=PPT                                                 
      PM(1,MV+3)=-PM(1,MV+1)                                      
      PM(2,MV+3)=-PM(2,MV+1)                                      
      PM(3,MV+3)=-PM(3,MV+1)                                      
      PM(4,MV+3)=PPT                                                 
      X=2.*PM(4,MV+3)/UU                                             
      IF(NP)23,23,26                                                    
   26 TEMP2=1.-X                                                        
      GOTO 40                                                           
   79 PMS=PM(9,MV+1)+PM(9,MV+3)                                   
      IF(IP1(2)-1) 50,51,51                                          
   51 IF(MQ) 52,53,54                                                   
   52 IF(PMS-1.88) 35,35,53                                             
   54 IF(MQ-3) 53,52,23                                                 
   50 IF(MQ) 57,53,58                                                   
   57 IF(ABS(FMQ)-1.) 59,59,61                                          
   59 IF(PMS-1.723) 53,335,53                                           
   61 IF(PM(9,MV+1)-.94) 335,335,63                                  
   63 IF(PM(9,MV+3)-.783) 53,335,53                                  
   58 IF(MQ-2) 53,59,61                                                 
   53 CALL CLEAD(1,LAB1,LAB2,P1,P2,IP1,IP2,MV,NIN,UU)     
      IF(NIN.GT.0) GOTO 1
      NP=2                                                              
      U2=UU                                                             
      CMB=PM(9,MV+1)
      CMF=PM(9,MV+3)
      IF(U2-CMB-CMF) 335,335,70           
   98 CONTINUE                                                          
   23 U1=EKL                                                            
      NTRY=NTRY+1                                                       
      IF(NTRY.GT.20) GOTO 1                                             
      NP=0                                                              
      KF=0                                                              
       KB=0                                                             
      LAMBDA=1                                                          
    2 LAMBDA=LAMBDA+1                                                   
       IF(LAMBDA.EQ.3) LAMBDA=4                                         
      LTEMP=MV+LAMBDA                                                   
      IF(max1.LE.LTEMP) GO TO 15                                        
   80 TEM=RNDM(-1)                                                      
      IF(TEM-PIPRM) 71,71,72                                            
   71 PM(9,LTEMP)=.14                                                
      PM(10,LTEMP)=0.                                                
      GOTO 75                                                           
   72 IF(TEM-.86)73,73,74                                               
   73 PM(9,LTEMP)=.77                                                
      PM(10,LTEMP)=.112                                              
      GOTO 75                                                           
   74 PM(9,LTEMP)=.783                                               
      PM(10,LTEMP)=.013                                              
   75 IM(2,LTEMP)=0                                                  
      IF(U1-PM(9,LTEMP)-.005)76,76,21                                
   76 IF(PM(9,LTEMP).GT..14.AND.TIN1.LT.2.) GOTO 80                  
      GOTO 5                                                            
   77 IF(CMB.NE.P2(9).OR.CMF.NE.P1(9)) GOTO 79                  
      GOTO 335                                                          
   21 S=RNDM(-1)                                                        
      SI=5.*S                                                           
      IF(RNDM(-1)-SI*EXP(1.-SI))22,22,21                          
   22 P=PIKF*SI                                                         
      SMT=SQRT(P**2+PM(9,LTEMP)**2)                                  
      FI=6.283185*RNDM(-1)                                              
      PT1=P*SIN(FI)                                                     
      PT2=P*COS(FI)                                                     
      RU=UMX/SMT
      IF(RU.LT.1.) goto 5
      IF(U1-SMT)21,21,81                                                
   81 YC=ALOG(RU)                                                       
   82 YMS=(1.-2.*RNDM(-1))*YC                                           
      PZ=SMT*SINH(YMS)                                                  
      EM=SMT*COSH(YMS)                                                  
      PM(1,LTEMP)=PT1                                                
      PM(2,LTEMP)=PT2                                                
      PM(3,LTEMP)=PZ                                                 
      PM(4,LTEMP)=SQRT(PT1**2+PT2**2+PZ**2)                          
      XX(LAMBDA)=2.*PZ/UU
********
      XM=2.*EM/UU
      IF(XM.GE.1.) GOTO 335
********      
      U1=U1-EM                                                          
      IF(U1)83,83,17                                                    
   17 IF(YMS.GT.0.) KF=1                                                
      IF(YMS.LE.0.) KB=1                                                
      GOTO 2                                                            
   83 U1=U1+EM                                                          
      XPI=2.*EM/UU                                                      
    5 IF(LAMBDA-2)77,77,55                                              
   55 NP=LAMBDA-1                                                       
      NPI=NP-2                                                          
      NP1=NPI+1                                                         
      IF(NP-3)77,8,10                                                   
    8 ND=2                                                              
      GOTO 11                                                           
   10 ND=NP                                                             
   11 UPI=SQRT(PM(4,MV+ND)**2+PM(9,MV+ND)**2)                     
      UPI=UPI+U1                                                        
      PM(4,MV+ND)=SQRT(UPI**2-PM(9,MV+ND)**2)                     
   91 PZD=SQRT(PM(4,MV+ND)**2-PM(1,MV+ND)**2-PM(2,MV+ND)**2)   
      IF(PM(3,MV+ND))92,92,93                                        
   92 PM(3,MV+ND)=-PZD                                               
      GOTO 94                                                           
   93 PM(3,MV+ND)=PZD                                                
   94 X=2.*PM(4,MV+L2)/UU                                            
      XX(ND)=2.*PM(3,MV+ND)/UU
      IF(TIN1.LT.1.) GO TO 95                                           
      IF((U2/UU)-.2) 95,95,96                                                 
   96 CALL ACIMET(P1,IP1,P2,IP2,L1,L2,U2,IND,MV,NP,PIKF   
     *, TIN1,TEMP2)                                                     
      IF(IND)97,40,35                                                   
   97 MET=1                                                             
      GOTO 40                                                           
   95 SIGMA(1)=0.                                                       
       SIGMA(2)=0.                                                      
       SIGMA(3)=0.                                                      
      IF(NP-3) 60,60,64                                                 
   60 NS=1                                                              
       NF=3                                                             
       N1=2                                                             
       N2=2                                                             
      GOTO 87                                                           
   64 NS=NP-1                                                           
       NF=NP                                                            
       N1=1                                                             
       N2=NP-2                                                          
   87 DO 12 I3=N1,N2                                                    
      SIGMA(1)=SIGMA(1)+PM(1,MV+I3)                                  
      SIGMA(2)=SIGMA(2)+PM(2,MV+I3)                                  
      SIGMA(3)=SIGMA(3)+PM(3,MV+I3)                                  
   12 CONTINUE                                                          
    9 PM1=PM(4,MV+NS)                                                 
      PM2=PM(4,MV+NF)                                                 
      P=SQRT(SIGMA(1)**2+SIGMA(2)**2+SIGMA(3)**2)                       
   29 IF(P.GE.PM1+PM2) GO TO 98                                           
      IF(P.LE.ABS(PM1-PM2)) GO TO 98                                      
      CT1=(PM2**2-PM1**2-P**2)/(2.*P*PM1)                                  
      CT2=(PM1**2-PM2**2-P**2)/(2.*P*PM2)                                  
      ST1=SQRT(1.-CT1**2)                                               
       ST2=SQRT(1.-CT2**2)                                              
      PI=3.14159                                                        
       FI1=2.*PI*RNDM(-1)                                               
       FI2=FI1+PI                                                       
      PL(1)=PM1*ST1*COS(FI1)                                             
       PL(2)=PM1*ST1*SIN(FI1)                                            
       PL(3)=PM1*CT1                                                     
      CALL ROTOR(SIGMA,V,PL,PLST)                                       
      DO 13 I4=1,3                                                      
   13 PM(I4,MV+NS)=PLST(I4)                                          
      PL(1)=PM2*ST2*COS(FI2)                                             
       PL(2)=PM2*ST2*SIN(FI2)                                            
       PL(3)=PM2*CT2                                                     
      CALL ROTOR(SIGMA,V,PL,PLST)                                       
   46 DO 14 I5=1,3                                                      
   14 PM(I5,MV+NF)=PLST(I5)                                          
   40 CONTINUE                                                          
      MQR=IM(1,MV+1)+IM(1,MV+3)
      DQ=ABS(MQ-MQR)
      IF(DQ.GE.2.) GOTO 335
      XX(1)=2.*PM(3,MV+1)/UU
      XX(3)=2.*PM(3,MV+3)/UU
      CALL VELOC(V,CTV,STV,CFV,SFV)                                     
      L=0                                                               
      DO 41 J=1,NP                                                      
      L=L+1                                                             
      IF(NP.EQ.2.AND.J.EQ.2) L=3                                        
      LM=MV+L                                                           
      EPM=SQRT(PM(4,LM)**2+PM(9,LM)**2) 
      IF(XX(L)) 27,27,28
   27 PM(19,LM)=P2(11)*ABS(XX(L))
c      IF(IP2(2).GT.0.AND.L.NE.1) PM(19,LM)=1.5*PM(19,LM)
      GOTO 38
   28 PM(19,LM)=P1(11)*XX(L)
c      IF(IP1(2).GT.0.AND.L.NE.3) PM(19,LM)=1.5*PM(19,LM)
   38 CONTINUE
      IF(PM(19,LM).GT.1.) PM(19,LM)=1.
      PM(16,LM)=0.
      PM(17,LM)=0.
      DO 47 JJ=1,3                                                      
   47 PL(JJ)=PM(JJ,LM)                                               
      CALL ROTV(LM,PL,CTV,STV,CFV,SFV)                                  
   41 CONTINUE                                                           
      TMNEY=TEMP2                                                       
      NIN=0                                                             
      RETURN                                                            
    1 NIN=2                                                             
      NP=0
      RETURN                                                            
   15 PRINT 16                                                          
   16 FORMAT(10X,29H MEMORY IS EXCEEDED IN  bLACT)                      
      NP=0                                                              
      RETURN                                                            
      END                                                                   
************************************************************************
      SUBROUTINE ACIMET(F1,I1,F2,I2,L1,L2,U2,IND,MV,NP,PIKF,TIN1,TM2)   
************************************************************************
*   PURPOSE:   RECALCULATES OF KINEMATICS OF LEADING PARTICLES         *
*              (REMNANTS OF PARTNERS) AT HIGH ENERGIES IN ORDER        *
*              TO TAKE INTO ACCOUNT INELASTIC DIFFRACTION.             *
*----------------------------------------------------------------------*
*   INPUT QUANTATIES:                                                  *
*     F1,I1    - CHARACTERISTICS OF THE FIRST LEADING PARICLE.         *
*     F2,I2    - CHARACTERISTICS OF THE SECOND LEADING PARTICLE.       *
*     L1,L2    - LABELS INDICATING TYPES OF LEADING PARTICES.          *
*     U2       - TOTAL ENERGY OF PARTNERS IN C.M.S.                    *
*     MV       - CURRENT NUMBER OF CASCADE PARTICES IN MEMORY          *
*     NP       - THE NUMBER OF DIRECT SECONDARIES GENERATED IN BLACT.  *
*     PIKF     - PARAMETER FOR TRANSVERSE MOMENTUM DISTR.              *
*     TIN1     - KINETIC ENERGY OF THE FIRST PARTNER IN THE REST       *
*                FRAME OF THE SECOND ONE.                              *
*     TM2      - INELASTICITY.                                         *
************************************************************************
C      DOUbLE PRECISION V
      parameter (max1=2000)
      COMMON /PMIM/ PM(19,max1),IM(5,MAX1)                              
      DIMENSION V(3),P1(3),P2(11)                                       
      DIMENSION F1(11),F2(11),I1(3),I2(3),T1(3),T2(3)                   
      PAR=.003*TIN1+.14                                                 
      N1=MV+L1                                                          
      N2=MV+L2                                                          
      NPI=NP/2                                                          
      NP2=NP-NPI                                                        
      MB=0                                                              
       MF=0                                                             
      PX=0.                                                             
       PY=0.                                                            
       PZ=0.                                                            
      KB=0                                                              
       KF=0                                                             
      DO 10 J=2,NP                                                      
      IF(J.EQ.3) GOTO 10                                                
      PX=PX+PM(1,MV+J)                                                  
      PY=PY+PM(2,MV+J)                                                  
      PZ=PZ+PM(3,MV+J)                                                  
      IF(PM(3,MV+J))4,4,5                                               
    4 KB=1                                                              
      GOTO 10                                                           
    5 KF=1                                                              
   10 CONTINUE                                                          
c      IF(RNDM(-1)-EXP(-(TM2**2/PAR))) 17,17,12                          
c   17 IF(KF)6,6,7                                                       
c    6 PDF=F1(10)                                                        
c      PDB=PM(10,N1)                                                     
c      CMF=F1(9)                                                         
c      CMB=PM(9,N1)                                                      
c      IMF=I1(1)                                                         
c      IMB=IM(1,N1)                                                      
c      IND=-1                                                            
c       GOTO 8                                                           
c    7 IF(KB)9,9,12                                                      
c    9 IF(I1(2).GT.0) GOTO 13                                            
c      IF(NPI-NP2) 12,13,12                                              
c   13 PDB=F2(10)                                                        
c      CMB=F2(9)                                                         
c      IMB=I2(1)                                                         
c      IF(I1(2)) 14,14,15                                                
c   14 PDF=F1(10)                                                        
c      CMF=F1(9)                                                         
c      IMF=I1(1)                                                         
c      GOTO 16                                                           
c   15 CMF=PM(9,N2)                                                      
c      PDF=PM(10,N2)                                                     
c      IMF=IM(1,N2)                                                      
c   16 IND=-1                                                            
c      GOTO 8                                                           
   12 IND=0                                                             
      CMB=PM(9,N1)                                                      
      CMF=PM(9,N2)                                                      
      PDB=PM(10,N1)                                                     
      PDF=PM(10,N2)                                                     
      IMB=IM(1,N1)                                                      
      IMF=IM(1,N2)                                                      
    8 PX=-PX                                                            
      PY=-PY                                                            
      PZ=-PZ                                                            
      PMOD=SQRT(PX**2+PY**2+PZ**2)                                      
      S2=U2**2-PMOD**2                                                  
      SI=(CMB+CMF)**2                                                   
      IF(S2-SI)35,35,1                                                  
    1 SM=SQRT(S2)                                                       
      V(1)=PX/U2                                                        
      V(2)=PY/U2                                                        
      V(3)=PZ/U2                                                        
      EF=(SM**2+CMF**2-CMB**2)/(2.*SM)                                  
      EB=SM-EF                                                          
      PPT=SQRT(EF**2-CMF**2)                                            
      IF(PPT) 22,22,2
   22 P1(1)=0.
      P1(2)=0.
      P1(3)=0.
      GOTO 23
    2 SS=RNDM(-1)                                                       
      SSN=3.*SS                                                         
      IF(RNDM(-1)-EXP(-SSN**2))3,3,2                                    
    3 POP=PIKF*SSN                                                  
      IF(POP.GT.PPT) GOTO 2                                             
      STP=POP/PPT                                                       
       CTP=SQRT(1.-STP**2)                                              
      CTT=-CTP                                                          
       STT=SQRT(1.-CTT**2)                                              
      PI=3.14159                                                        
       FIN=2.*PI*RNDM(-1)                                               
       FTA=FIN+PI                                                       
      P1(1)=PPT*STT*COS(FTA)                                            
      P1(2)=PPT*STT*SIN(FTA)                                            
      P1(3)=PPT*CTT                                                     
   23 P2(9)=CMB                                                         
      CALL CINEMA(P1,V,P2)                                              
      PB=SQRT(P2(8)*(P2(8)+2.*P2(9)))                                   
      T1(1)=PB*P2(4)*P2(7)                                              
      T1(2)=PB*P2(4)*P2(6)                                              
      T1(3)=PB*P2(5)                                                    
      P1(1)=PPT*STP*COS(FIN)                                            
      P1(2)=PPT*STP*SIN(FIN)                                            
      P1(3)=PPT*CTP                                                     
      P2(9)=CMF                                                         
      CALL CINEMA(P1,V,P2)                                              
      PF=SQRT(P2(8)*(P2(8)+2.*P2(9)))                                   
      T2(1)=PF*P2(4)*P2(7)                                              
      T2(2)=PF*P2(4)*P2(6)                                              
      T2(3)=PF*P2(5)                                                    
      IF(T1(3).GT.0..OR.T2(3).LT.0.) GOTO 35                            
      DO 11 I=1,3                                                       
      PM(I,N1)=T1(I)                                                    
      PM(I,N2)=T2(I)                                                    
   11 CONTINUE                                                          
      PM(4,N1)=PB                                                       
      PM(4,N2)=PF                                                       
      PM(9,N1)=CMB                                                      
      PM(9,N2)=CMF                                                      
      PM(10,N1)=PDB                                                     
      PM(10,N2)=PDF                                                     
      IM(1,N1)=IMB                                                      
      IM(1,N2)=IMF                                                      
      RETURN                                                            
   35 IND=1                                                             
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE BLLABL(LAB1,LAB2,P1,IP1,P2,IP2)                        
************************************************************************
*   PURPOSE: SPECIFY LABELS OF INERACTING PARTNERS.                    *
************************************************************************
      DIMENSION P1(11),IP1(3),P2(11),IP2(3)                             
      IF(P1(9)-.14)41,41,42                                             
   41 IF(IP1(1))43,44,45                                                
   43 LAB1=5                                                            
      GOTO 50                                                           
   44 LAB1=6                                                            
      GOTO 50                                                           
   45 LAB1=7                                                            
      GOTO 50                                                           
   42 IF(P1(9)-.783)47,48,49                                            
   48 LAB1=1                                                            
      GOTO 50                                                           
   47 IF(IP1(1))51,52,53                                                
  51  LAB1=2                                                            
      GOTO 50                                                           
   52 LAB1=3                                                            
      GOTO 50                                                           
   53 LAB1=4                                                            
      GOTO 50                                                           
   49 IF(P1(9)-.94)54,54,55                                             
   54 IF(IP1(1))56,56,57                                                
   56 LAB1=2                                                            
      GOTO 50                                                           
   57 LAB1=1                                                            
      GOTO 50                                                           
   55 IF(IP1(1)) 58,59,30                                               
   58 LAB1=6                                                            
      GOTO 50                                                           
   59 LAB1=5                                                            
      GOTO 50                                                           
   30 IF(IP1(1)-1) 31,31,32                                             
   31 LAB1=4                                                            
      GOTO 50                                                           
   32 LAB1=3                                                            
   50 CONTINUE                                                          
      IF(P2(9)-.14) 61,61,62
   61 IF(IP2(1)) 63,64,65
   63 LAB2=5
      GOTO 20
   64 LAB2=6
      GOTO 20
   65 LAB2=7
      GOTO 20
   62 IF(P2(9)-.783) 67,68,24
   68 LAB2=1
      GOTO 20
   67 IF(IP2(1)) 71,72,73
   71 LAB2=2
      GOTO 20
   72 LAB2=3
      GOTO 20
   73 LAB2=4
      GOTO 20
   24 IF(P2(9)-.94) 4,4,5                                               
    4 IF(IP2(1))6,6,7                                                   
    6 LAB2=2                                                            
      GOTO 20                                                           
    7 LAB2=1                                                            
      GOTO 20                                                           
    5 IF(IP2(1))8,9,10                                                  
    8 LAB2=6                                                            
      GOTO 20                                                           
    9 LAB2=5                                                            
      GOTO 20                                                           
   10 IF(IP2(1)-1)11,11,12                                              
   11 LAB2=4                                                            
      GOTO 20                                                           
   12 LAB2=3                                                            
   20 CONTINUE                                                          
      RETURN                                                            
      END                                                               
************************************************************************
      FUNCTION RESMAS(E0,GM)                                            
************************************************************************
*   PURPOSE:  CALCULATES DELTA-ISOBAR MASS.                            *
************************************************************************
    2 X=1.-2.*RNDM(-1)                                                  
      DE=3.*X                                                           
      IF(RNDM(-1)-(1./(DE**2+1.))) 1,1,2                                
    1 RESMAS=1.5*GM*X+E0                                                
      IF(RESMAS-1.1) 2,2,3                                              
    3 CONTINUE                                                          
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE ISODEC(IZ,V,P,IP,N1,NP,MV)                             
************************************************************************
*   PURPOSE:  CALCULATION OF ISOBAR DECAY                              *
************************************************************************
C      DOUbLE PRECISION V
      parameter (max1=2000)
      COMMON /PMIM/ PM(19,max1),IM(5,MAX1)                              
      COMMON /KSI/ KSI(3)                                               
      DIMENSION PI(3),V(3),P(11),IP(3),IS(3),IX(3)                      
      IF(MV.GT.max1) GOTO 7                                             
      VM=SQRT(V(1)**2+V(2)**2+V(3)**2)                                  
      IF(VM.EQ.0.) GOTO 111                                             
      CALL VELOC(V,CTV,STV,CFV,SFV)                                     
  111 IF(IZ) 11,11,12                                                   
   12 NM1=MV+1                                                          
      NM2=N1                                                            
      GOTO 13                                                           
   11 NM1=MV+1                                                          
      NM2=MV+3                                                          
   13 NP=2                                                              
      EPI=(P(9)**2+.0196-.883)/(2.*P(9))                                
      EN=P(9)-EPI                                                       
      PPI=SQRT(EPI**2-.0196)                                            
      FI=6.283185*RNDM(-1)                                              
      T=P(9)-.14-.94                                                    
      IF(T.LE.0) T=0.                                                   
      IKS=IP(1)+2                                                       
      IM(2,NM1)=1                                                       
      IM(2,NM2)=0                                                       
      GOTO(1,2,4,6),IKS                                                 
 1    IM(1,NM1)=0                                                       
       IM(1,NM2)=-1                                                     
       GO TO 9                                                          
 2    IF(RNDM(-1).LT.0.66667) GO TO 3                                   
      IM(1,NM1)=1                                                       
       IM(1,NM2)=-1                                                     
       GO TO 9                                                          
 3    IM(1,NM1)=0                                                       
       IM(1,NM2)=0                                                      
       GO TO 9                                                          
 4    IF(RNDM(-1).LT.0.66667) GO TO 5                                   
      IM(1,NM1)=0                                                       
       IM(1,NM2)=1                                                      
       GO TO 9                                                          
 5    IM(1,NM1)=1                                                       
       IM(1,NM2)=0                                                      
       GO TO 9                                                          
 6    IM(1,NM1)=1                                                       
       IM(1,NM2)=1                                                      
    9 DO 10 K=1,3                                                       
      IS(K)=IM(K,NM1)                                                   
   10 IX(K)=IM(K,NM2)                                                   
 14   S=1.-2.*RNDM(-1)                                                  
      SN=.25*(1.+3.*S**2)                                               
      IF(RNDM(-1)-SN) 15,15,14                                          
 15   CT=S                                                              
      ST=SQRT(1.-CT**2)                                                 
      R=PPI*ST                                                          
      PI(1)=R*COS(FI)                                                   
      PI(2)=R*SIN(FI)                                                   
      PI(3)=PPI*CT                                                      
c      IF(CT) 21,22,22
c   21 IF(IM(5,N1)-1)23,23,25
c   23 PI(1)=-PI(1)
c      PI(2)=-PI(2)
c      PI(3)=-PI(3)
c      GOTO 25
c   22 IF(IM(5,N1)-1) 25,25,23
c   25 IF(VM-.0) 112,112,113                                   
c  112 DO 114 L=1,3                                         
c      PM(L,NM2)=PI(L)                                     
c  114 PM(L,NM1)=-PI(L)                                    
c      GOTO 115                                         
  113 CALL ROTV(NM2,PI,CTV,STV,CFV,SFV)                                 
      PI(1)=-PI(1)                                                      
      PI(2)=-PI(2)                                                      
      PI(3)=-PI(3)                                                      
      CALL ROTV(NM1,PI,CTV,STV,CFV,SFV)                                 
  115 PM(4,NM1)=PPI                                                     
       PM(4,NM2)=PPI                                                    
      PM(9,NM1)=0.94                                                    
       PM(9,NM2)=0.14                                                   
      PM(10,NM1)=0.                                                     
       PM(10,NM2)=0.                                                    
      PM(19,NM1)=PM(19,N1)*.66667
      PM(19,NM2)=PM(19,N1)*.33333
      PM(16,NM1)=PM(16,N1)
      PM(17,NM1)=PM(17,N1)
      PM(16,NM2)=PM(16,N1)
      PM(17,NM2)=PM(17,N1)
      PM(19,NM1)=PM(19,N1)
      RETURN                                                            
 7    NP=-1                                                             
       PRINT 8,MV                                                       
 8    FORMAT(20X,'MEMORY IS EXCEDED IN ISODEC N=',I5)                   
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE RESDEC(MV)                                             
************************************************************************
*  PURPOSE:   MAKES RESONANCES DECAY WHEN CASCADING PROCESS IS OVER.   *
************************************************************************
C      DOUbLE PRECISION V,VR
      parameter (max1=2000)
      COMMON /PMIM/ PM(19,max1),IM(5,MAX1)                              
      COMMON /CARCIL/NK1P(240),TIME(240),MK1P(max1),MK2P(max1),         
     *MK1T(max1),MK2T(max1)/ACTIV/MPA(240),MYP(max1),MYT(max1)          
      DIMENSION P(11),IP(3),V(3),PS(3),PN(11)                           
      M=MV                                                              
      DO 1 I=1,M                                                        
      IF(PM(10,I))1,1,2                                                 
    2 DO 3 L=1,10                                                       
    3 P(L)=PM(L,I)                                                      
      IP(1)=IM(1,I)                                                     
      IP(2)=IM(2,I)                                                     
      VR=SQRT(P(8)*(P(8)+2.*P(9)))/(P(8)+P(9))                          
      V(1)=VR*P(4)*P(7)                                                 
      V(2)=VR*P(4)*P(6)                                                 
      V(3)=VR*P(5)                                                      
      IF(IP(2))4,4,5                                                    
    5 CALL ISODEC(1,V,P,IP,I,NP,MV)                                     
      GOTO 6                                                            
    4 IF(P(9)-.77) 11,11,12                                             
   11 CALL RODEC(1,V,P,IP,I,NP,MV)                                      
      GOTO 6                                                            
   12 CALL W0DEC(1,V,P,IP,I,NP,MV)                                      
    6 DO 7 L=1,3                                                        
    7 PS(L)=PM(L,I)                                                     
      PN(9)=PM(9,I)                                                     
      CALL CINEMA(PS,V,PN)                                              
      DO 8 L=4,8                                                        
    8 PM(L,I)=PN(L)                                                     
      II=NP-1                                                           
      DO 13 K=1,II                                                      
      MJ=MV+K                                                           
      DO 9 L=1,3                                                        
    9 PS(L)=PM(L,MJ)                                                    
      PN(9)=PM(9,MJ)                                                    
      CALL CINEMA(PS,V,PN)                                              
      DO 10 L=4,8                                                       
   10 PM(L,MJ)=PN(L)                                                    
      MYP(MJ)=0                                                         
      MYT(MJ)=0                                                         
      MK1P(MJ)=0                                                        
      MK2P(MJ)=0                                                        
      MK1T(MJ)=0                                                        
      MK2T(MJ)=0                                                        
   13 CONTINUE                                                          
      MV=MV+II                                                          
    1 CONTINUE                                                          
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE RODEC(IZ,V,P,IP,N1,NP,MV)                              
************************************************************************
*   PURPOSE:      CALCULATION OF RO-MESON DECAY                        *
************************************************************************
C      DOUbLE PRECISION V
      parameter (max1=2000)
      COMMON /PMIM/ PM(19,max1),IM(5,MAX1)                              
      DIMENSION PI(3),V(3),P(11),IP(3)                                  
      IF(MV.GT.max1) GOTO 7                                             
      CALL VELOC(V,CTV,STV,CFV,SFV)                                     
      IF(IZ)11,11,12                                                    
   12 NM1=N1                                                            
      NM2=MV+1                                                          
      GOTO 13                                                           
   11 NM1=MV+1                                                          
      NM2=MV+3                                                          
   13 NP=2                                                              
      EPI=P(9)/2.                                                       
      PPI=SQRT(EPI**2-.0196)                                            
      FI=6.283185*RNDM(-1)                                              
    9 SS=RNDM(-1)                                                       
      SSN=3.*SS                                                         
      IF(RNDM(-1)-EXP(-SSN**2)) 10,10,9                                 
   10 POP=0.5*SSN                                                       
      IF(POP.GE.PPI) GOTO 9                                             
      ST=POP/PPI                                                        
      CT=SQRT(1.-ST**2)                                                 
      IF(RNDM(-1).LT..5) CT=-CT                                         
      R=PPI*ST                                                          
      PI(1)=R*COS(FI)                                                   
      PI(2)=R*SIN(FI)                                                   
      PI(3)=PPI*CT                                                      
      CALL ROTV(NM1,PI,CTV,STV,CFV,SFV)                                 
      DO 1 I=1,3                                                        
    1 PI(I)=-PI(I)                                                      
      CALL ROTV(NM2,PI,CTV,STV,CFV,SFV)                                 
      PM(4,NM1)=PPI                                                     
      PM(4,NM2)=PPI                                                     
      PM(9,NM1)=.14                                                     
       PM(9,NM2)=.14                                                    
      PM(10,NM1)=0.                                                     
       PM(10,NM2)=0.                                                    
      PM(19,NM1)=PM(19,N1)*.5
      PM(19,NM2)=PM(19,N1)*.5
      PM(16,NM1)=PM(16,N1)
      PM(17,NM1)=PM(17,N1)
       PM(16,NM2)=PM(16,N1)
      PM(17,NM2)=PM(17,N1)
      IF(IP(1)) 2,3,4                                                   
    2 IM(1,NM1)=-1                                                      
      IM(1,NM2)=0                                                       
      GOTO 5                                                            
    3 IM(1,NM1)=-1                                                      
      IM(1,NM2)=1                                                       
      GOTO 5                                                            
    4 IM(1,NM1)=1                                                       
       IM(1,NM2)=0                                                      
    5 IM(2,NM1)=0                                                       
       IM(2,NM2)=0                                                      
      RETURN                                                            
    7 NP=0                                                              
      PRINT 8,MV                                                        
    8 FORMAT(20X,'MEMORY IS EXEEDED IN RODEC N =',I5)                   
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE W0DEC(IZ,V,P,IP,N1,NP,MV)                              
************************************************************************
*   PURPOSE:      CALCULATION OF W-MESON DECAY                         *
************************************************************************
C      DOUbLE PRECISION V
      parameter (max1=2000)
      COMMON /PMIM/ PM(19,max1),IM(5,MAX1)                              
      DIMENSION P(11),IP(3),RI(3),PR(3),V(3),SIG(3)                     
      IF(MV.GT.max1) GOTO 7                                             
      CALL VELOC(V,CTV,STV,CFV,SFV)                                     
      IF(IZ) 11,11,12                                                   
   12 NM1=N1                                                            
       NM2=MV+1                                                         
      NM3=MV+2                                                          
      GOTO 13                                                           
   11 NM1=MV+1                                                          
      NM2=MV+2                                                          
      NM3=MV+3                                                          
   13 NP=3                                                              
    5 U1=P(9)                                                           
      LM=1                                                              
    1 S=RNDM(-1)                                                        
      SI=4.*S                                                           
      IF(RNDM(-1)-SI**2*EXP(1.-SI**2))2,2,1                             
    2 PPI=.35*SI                                                        
      EPI=SQRT(PPI**2+.0196)                                            
      U1=U1-EPI                                                         
      IF(U1)8,8,9                                                       
    8 IF(LM-2)5,5,10                                                    
   10 U1=U1+EPI                                                         
      IF(U1-.16)5,5,15                                                  
    9 IF(LM-2)6,14,10                                                   
   14 P1=PPI                                                            
       E1=EPI                                                           
      GOTO 16                                                           
    6 FI=6.283285*RNDM(-1)                                              
      CT=1.-2.*RNDM(-1)                                                 
      ST=SQRT(1.-CT**2)                                                 
      R=PPI*ST                                                          
      RI(1)=R*COS(FI)                                                   
      RI(2)=R*SIN(FI)                                                   
      RI(3)=PPI*CT                                                      
      SIG(1)=RI(1)                                                      
      SIG(2)=RI(2)                                                      
      SIG(3)=RI(3)                                                      
      PM(4,NM1)=PPI                                                     
      PM(9,NM1)=.14                                                     
      PM(10,NM1)=0.                                                     
   16 LM=LM+1                                                           
      GOTO 1                                                            
   15 E2=U1                                                             
       P2=SQRT(E2**2-.0196)                                             
      IF(PM(4,NM1).GE.(P1+P2)) GOTO 5                                   
      IF(PM(4,NM1).LE.ABS(P1-P2)) GOTO 5                                
      CT1=(P2**2-P1**2-PM(4,NM1)**2)/(2.*PM(4,NM1)*P1)                  
      CT2=(P1**2-P2**2-PM(4,NM1)**2)/(2.*PM(4,NM1)*P2)                  
      ST1=SQRT(1.-CT1**2)                                               
       ST2=SQRT(1.-CT2**2)                                              
      PI=3.14159                                                        
       FI1=2.*PI*RNDM(-1)                                               
       FI2=FI1+PI                                                       
      CALL ROTV(NM1,RI,CTV,STV,CFV,SFV)                                 
      RI(1)=P1*ST1*COS(FI1)                                             
      RI(2)=P1*ST1*SIN(FI1)                                             
      RI(3)=P1*CT1                                                      
      CALL ROTOR(SIG,V,RI,PR)                                           
      CALL ROTV(NM2,PR,CTV,STV,CFV,SFV)                                 
      PM(4,NM2)=P1                                                      
      RI(1)=P2*ST2*COS(FI2)                                             
      RI(2)=P2*ST2*SIN(FI2)                                             
      RI(3)=P2*CT2                                                      
      CALL ROTOR(SIG,V,RI,PR)                                           
      CALL ROTV(NM3,PR,CTV,STV,CFV,SFV)                                 
      PM(4,NM3)=P2                                                      
      PM(9,NM2)=.14                                                     
       PM(9,NM3)=.14                                                    
      PM(10,NM2)=0.                                                     
       PM(10,NM3)=0.                                                    
      IM(1,NM1)=-1                                                      
      IM(1,NM2)=0                                                       
      IM(1,NM3)=1                                                       
      IM(2,NM1)=0                                                       
       IM(2,NM2)=0                                                      
       IM(2,NM3)=0                                                      
      PM(19,NM1)=PM(19,N1)*.3333
      PM(19,NM2)=PM(19,N1)*.3333
      PM(19,NM3)=PM(19,N1)*.3333
      PM(16,NM1)=PM(16,N1)
      PM(17,NM1)=PM(17,N1)
      PM(16,NM2)=PM(16,N1)
      PM(17,NM2)=PM(17,N1)
      PM(16,NM3)=PM(16,N1)
      PM(17,NM3)=PM(17,N1)
      RETURN                                                            
    7 NP=0                                                              
      PRINT 20,MV                                                       
   20 FORMAT(20X,'MEMORY IS EXEEDED IN W0DEC N = ',I5)                  
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE ROTV(M,P,CT,ST,CF,SF)                                  
************************************************************************
*   PURPOSE:   ROTATION OF COORDINATE SYSTEM.                          *
************************************************************************
      parameter (max1=2000)
      COMMON /PMIM/ PM(19,max1),IM(5,MAX1)                              
      DIMENSION P(3)                                                    
      PXY=P(1)*CT+P(3)*ST                                               
      PM(1,M)=PXY*CF-P(2)*SF                                            
      PM(2,M)=PXY*SF+P(2)*CF                                            
      PM(3,M)=-P(1)*ST+P(3)*CT                                          
      RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE VELOC(V,CT,ST,CF,SF)                                   
************************************************************************
*   PURPOSE: CALCULATION OF ANGLES FOR CMS VELOCITY VECTOR.            *
************************************************************************
C      DOUbLE PRECISION V
      DIMENSION V(3)                                                    
      VM=SQRT(V(1)**2+V(2)**2+V(3)**2)                                  
      IF(VM-.0001) 4,4,5                                                
    4 ST=0.                                                             
      CT=1.                                                             
      SF=0.                                                             
      CF=1.                                                             
      RETURN                                                            
    5 CT=V(3)/VM                                                        
      ST=SQRT(1.-CT**2)                                                 
      IF(ST) 1,1,2                                                      
    1 ST=0.                                                             
      CF=1.                                                             
      SF=0.                                                             
      GOTO 3                                                            
    2 CF=V(1)/(VM*ST)                                                   
      SF=V(2)/(VM*ST)                                                   
    3 RETURN                                                            
      END                                                               
************************************************************************
      SUBROUTINE CLEAD(MZZ,LAB1,LAB2,P1,P2,IP1,IP2,MV,NIN,UU)           
************************************************************************
*   PURPOSE:   DETERMINATION MASSES AND CHARGES OF LEADERS             *
************************************************************************
      parameter (max1=2000)
      COMMON /PMIM/ PM(19,max1),IM(5,MAX1)                        
      DIMENSION IP1(3),IP2(3),P1(11),P2(11),V(3),PP(11)                        
      DIMENSION CMASS(7,2),QCH(7,2),PBB(6,6),PMB(7,7),ROPIOM(3),        
     1 GAMMA(7,2)                                                       
       DATA  GAMMA/                                                     
     1            .013,.11,.11,.11,0.,0.,0.,                            
     2           0.,0.,.11,.11,.11,.11,  0./                            
       DATA  CMASS/                                                     
     1            .783,  .77,.77,.77,.14,.14,.14,                       
     2           .94,.94,1.232,1.232,1.232,1.232,  0./                  
       DATA  QCH/                                                       
     1          0.,-1., 0.,+1.,-1., 0.,+1.,                             
     2         +1., 0.,+2.,+1., 0.,-1., 0./                             
       DATA  ROPIOM/0.4,0.8,1./                                         
       DATA  PBB/                                                       
     1    0.167, 0.311, 0.261, 0.174, 0.087, 0.   ,                     
     2    0.311, 0.167, 0.   , 0.087, 0.174, 0.261,                     
     3    0.083, 0.   , 0.550, 0.367, 0.   , 0.   ,                     
     4    0.055, 0.028, 0.367, 0.061, 0.489, 0.   ,                     
     5    0.028, 0.055, 0.   , 0.489, 0.061, 0.367,                     
     6    0.   , 0.083, 0.   , 0.   , 0.367, 0.550/                     
       DATA  PMB/                                                       
     1   0.250,0.250,0.250,0.250,0.0  ,0.0  ,0.0  ,                     
     2   0.250,0.500,0.250,0.000,0.0  ,0.0  ,0.000,                     
     3   0.250,0.250,0.250,0.250,0.0  ,0.0  ,0.0  ,                     
     4   0.250,0.000,0.250,0.500,0.000,0.0  ,0.0  ,                     
     5   0.200,0.400,0.200,0.000,0.133,0.067,0.000,                     
     6   0.200,0.200,0.200,0.200,0.067,0.067,0.067,                     
     7   0.200,0.000,0.200,0.400,0.000,0.067,0.133/                     
      MTRY=0                                                            
       NIN=0                                                            
      MZ=MZZ                                                            
      MB=IP1(2)+IP2(2)
      MQ=IP1(1)+IP2(1)                                                  
      FMQ=FLOAT(MQ)                                                     
   20 SUM=0.                                                            
      MTRY=MTRY+1                                                       
      IF(MTRY-10)40,41,41                                               
   40 TEMP1=RNDM(-1)                                                    
      IF(IP1(2)-1) 1,2,2                                                
    1 IF(MZ)10,10,12                                                    
   12 DO 13 I=1,5                                                       
      IF(CMASS(I,1)-PM(9,MV+3))13,14,13                              
   13 CONTINUE                                                          
   14 IF(I-1)6,6,15                                                     
   15 II=I+2                                                            
      DO 16 L=I,II                                                      
   16 SUM=SUM+PMB(L,LAB1)                                               
      SS=0.                                                             
      DO 17 L=I,II                                                      
      PMES=PMB(L,LAB1)/SUM                                                
      SS=SS+PMES                                                          
      IF(TEMP1-SS) 18,18,17                                             
   17 CONTINUE                                                          
   18 IM(1,MV+3)=QCH(L,1)                                            
      GOTO 6                                                            
   10 DO 3 I=1,7                                                        
      SUM=SUM+PMB(I,LAB1)                                               
      IF(TEMP1-SUM) 4,4,3                                               
    4 PM(9,MV+3)=CMASS(I,1)                                          
      E01=CMASS(I,1)                                                    
      PM(10,MV+3)=GAMMA(I,1)                                         
      IM(2,MV+3)=0.                                                  
      IM(1,MV+3)=QCH(I,1)                                            
      GOTO 6                                                            
    3 CONTINUE                                                          
    2 IF(MZ) 21,21,22                                                   
   22 IF(PM(9,MV+3)-.94)23,23,24                                     
   23 N1=1                                                              
       N2=2                                                             
      GOTO 25                                                           
   24 N1=3                                                              
       N2=6                                                             
   25 DO 26 L=N1,N2                                                     
   26 SUM=SUM+PBB(L,LAB1)                                               
      SS=0.                                                             
      DO 27 L=N1,N2                                                     
      PBAR=PBB(L,LAB1)/SUM                                                
      SS=SS+PBAR                                                          
      IF(TEMP1-SS)28,28,27                                              
   27 CONTINUE                                                          
   28 IM(1,MV+3)=QCH(L,2)                                            
      GOTO 6                                                            
   21 DO 5 I=1,6                                                        
      SUM=SUM+PBB(I,LAB1)                                               
      IF(TEMP1-SUM) 7,7,5                                               
    7 E01=CMASS(I,2)                                                    
      PM(9,MV+3)=E01                                                 
      IF(E01-.94) 71,71,70                                              
   71 PM(10,MV+3)=0.                                                 
      GOTO 73                                                           
   70 PM(9,MV+3)=RESMAS(E01,GAMMA(I,2))                              
   72 NM1=MV+3                                                          
      PM(10,MV+3)=GAMD(1,NM1,PI,PP,V)                                
   73 IM(2,MV+3)=1                                                   
      IM(1,MV+3)=QCH(I,2)                                            
      GOTO 6                                                            
    5 CONTINUE                                                          
    6 CONTINUE                                                          
      SUM=0.                                                            
      TEMP2=RNDM(-1)                                                    
      IF(MZ)29,29,30                                                    
   30 IF(PM(9,MV+1)-.9) 233,33,33                                    
C    TARGET PARTICLE IS MESON
  233 DO 213 I=1,5
	IF(CMASS(I,1)-PM(9,MV+1)) 213,214,213
  213 CONTINUE
  214 IF(I-1) 112,112,215
  215 II=I+2
      DO 216 L=I,II
  216 SUM=SUM+PMB(L,LAB2)
      SS=0.
      DO 217 L=I,II
      PMES=PMB(L,LAB2)/SUM
      SS=SS+PMES
      IF(TEMP2-SS) 218,218,217
  217 CONTINUE
  218 IM(1,MV+1)=QCH(L,1)
      GOTO 112
   33 IF(PM(9,MV+1)-.94) 333,333,34
  333 N1=1                                                              
       N2=2                                                             
      GOTO 35                                                           
   34 N1=3                                                              
       N2=6                                                             
   35 DO 36 L=N1,N2                                                     
   36 SUM=SUM+PBB(L,LAB2)                                               
      SS=0.                                                             
      DO 37 L=N1,N2                                                     
      PBAR=PBB(L,LAB2)/SUM                                                
      SS=SS+PBAR                                                          
      IF(TEMP2-SS)38,38,37                                              
   37 CONTINUE                                                          
   38 IM(1,MV+1)=QCH(L,2)                                            
      GOTO 11                                                           
   29 IF(MB) 299,299,300
  299 DO 203 I=1,7
      SUM=SUM+PMB(I,LAB2)
      IF(TEMP2-SUM) 204,204,203
  204 PM(9,MV+1)=CMASS(I,1)
      E02=CMASS(I,1)
      PM(10,MV+1)=GAMMA(I,1)
      IM(2,MV+1)=0
      IM(1,MV+1)=QCH(I,1)
      GOTO 11
  203 CONTINUE
  300 DO 8 I=1,6                                                        
      SUM=SUM+PBB(I,LAB2)                                               
      IF(TEMP2-SUM) 9,9,8                                               
    9 E02=CMASS(I,2)                                                    
      PM(9,MV+1)=E02                                                 
      IF(E02-.94) 81,81,80                                              
   81 PM(10,MV+1)=0.                                                 
      GOTO 83                                                           
   80 PM(9,MV+1)=RESMAS(E02,GAMMA(I,2))                              
   82 NM1=MV+1                                                          
      PM(10,MV+1)=GAMD(1,NM1,PI,PP,V)                                
   83 IM(2,MV+1)=1                                                   
      IM(1,MV+1)=QCH(I,2)                                            
      GOTO 11                                                           
    8 CONTINUE                                                          
   11 IF(MZ)112,112,113                                                 
  112 DU=UU-PM(9,MV+1)-PM(9,MV+3)                                 
      IF(DU) 42,42,113                                                  
   42 IF(MB) 20,20,242
  242 IF(E01-.94) 51,43,44                                              
   51 IF(E01-.77) 53,52,59                                              
   53 IF(UU-1.24) 55,55,60                                              
   60 MZ=1                                                              
      GOTO 50                                                           
   55 PM(9,MV+1)=.94                                                 
      PM(10,MV+1)=0.                                                 
      MZ=1                                                              
      GOTO 33                                                           
   52 IF(UU-1.24) 54,54,56                                              
   54 PM(9,MV+3)=.14                                                 
      PM(10,MV+3)=0.                                                 
      PM(9,MV+1)=.94                                                 
      PM(10,MV+1)=0.                                                 
      MZ=1                                                              
      GOTO 12                                                           
   56 IF(UU-1.71) 57,57,58                                              
   57 PM(9,MV+3)=.14                                                 
      PM(10,MV+3)=0.                                                 
      IF(MQ) 65,62,66                                                   
   66 IF(MQ-3) 62,63,63                                                 
   65 IF(ABS(FMQ)-2.) 62,63,63                                          
   63 IF(PM(9,MV+1)-1.0) 64,64,62                                    
   64 PM(9,MV+1)=P2(9)                                               
      PM(10,MV+1)=.11                                                
   62 MZ=1                                                              
      GOTO 12                                                           
   58 IF(UU-1.88) 61,61,60                                              
   61 IF(MQ) 67,55,68                                                   
   68 IF(MQ-3) 55,75,75                                                 
   67 IF(ABS(FMQ)-2.) 55,75,75                                          
   75 PM(9,MV+3)=.14                                                 
      PM(10,MV+3)=0.                                                 
      MZ=1                                                              
      GOTO12                                                            
   59 IF(UU-1.24) 54,54,76                                              
   76 IF(UU-1.723) 57,57,77                                             
   77 IF(UU-1.883) 78,78,90                                             
   78 IF(MQ) 79,55,85                                                   
   79 IF(ABS(FMQ)-2.) 86,87,87                                          
   87 PM(9,MV+3)=.14                                                 
      PM(10,MV+3)=0.                                                 
      MZ=1                                                              
      GOTO12                                                            
   86 PM(9,MV+3)=.77                                                 
      PM(10,MV+3)=.11                                                
      PM(9,MV+1)=.94                                                 
      PM(10,MV+1)=0.                                                 
      MZ=1                                                              
      GOTO 12                                                           
   85 IF(MQ-2) 55,86,87                                                 
   90 IF(MQ) 91,60,92                                                   
   91 IF(ABS(FMQ)-2.) 60,93,93                                          
   92 IF(MQ-2) 50,50,93                                                 
   93 PM(9,MV+3)=.77                                                 
      PM(10,MV+3)=.11                                                
      MZ=1                                                              
      GOTO 50                                                           
   44 IF(E02-.94)45,45,20                                               
   45 IF(UU-2.041) 97,97,47                                             
   47 PM(9,MV+3)=RESMAS(E01,.11)                                     
      IF(UU-PM(9,MV+1)-PM(9,MV+3)) 47,47,114                      
  114 NM1=MV+3                                                          
      PM(10,MV+3)=GAMD(1,NM1,PI,PP,V)                                
      GOTO 111                                                          
   43 IF(UU-2.041) 96,96,50                                             
   96 PM(9,MV+1)=.94                                                 
      PM(10,MV+1)=0.                                                 
      MZ=1                                                              
      GOTO 33                                                           
   97 PM(9,MV+3)=.94                                                 
      PM(10,MV+3)=0.                                                 
      MZ=1                                                              
      GOTO 23                                                           
   50 PM(9,MV+1)=RESMAS(E02,.11)                                     
      IF(UU-PM(9,MV+1)-PM(9,MV+3)) 50,50,115                      
  115 NM1=MV+1                                                          
      PM(10,MV+1)=GAMD(1,NM1,PI,PP,V)                                
      GOTO 111                                                          
  111 MZ=1                                                              
  113 MQR=IM(1,MV+1)+IM(1,MV+3)                                   
      DMR=MQ-MQR                                                        
      MODMR=ABS(DMR)                                                    
      IF(MODMR.GE.2) GOTO 20                                            
      IF(MZ.GT.0.AND.MODMR.GT.0) GOTO 20                                
      RETURN                                                            
   41 NIN=1                                                             
       RETURN                                                           
      END                                                               
************************************************************************
      FUNCTION GAMD(I3,NM1,PT,P,V)                                      
************************************************************************
*     CALCULATION OF ISOBAR WIDTH                                      *
************************************************************************
C      DOUbLE PRECISION V
      parameter (max1=2000)
      COMMON /PMIM/PM(19,max1),IM(5,MAX1)                               
      DIMENSION PX(3),P(11),T(11),V(3)                                  
      ROT(A,B,C)=C/(1.+A**2+(B**2)*(B**2))                              
      IF(I3.EQ.0) GOTO 1                                                
C     CALCULATION OF ISOBAR WIDTH IN INELASTIC EVENT                    
      PXM=SQRT((PM(9,NM1)**2-.9032)**2-.06927)/(2.*PM(9,NM1))           
      GOTO 2                                                            
C     CALCULATION OF ISOBAR WIDTH IN ELASTIC EVENT                      
    1 CALL CMS(P,1,V,PX)                                                
      PXM=SQRT(SP(PX,PX))                                               
    2 PT0=1.15                                                          
       GO=120.                                                          
      PT=5.067*PXM                                                      
       A1=.83*PT                                                        
       B1=.62*PT                                                        
       C1=(PT*PT)*PT                                                    
      RT=ROT(A1,B1,C1)                                                  
      A2=.83*PT0                                                        
       B2=.62*PT0                                                       
       C2=(PT0*PT0)*PT0                                                 
       R0=ROT(A2,B2,C2)                                                 
      GT=GO*RT/R0                                                       
    3 GAMD=.001*GT                                                      
      RETURN                                                            
      END                                                               
      SUBROUTINE CMS(P,L,V,PX)                                          
C     TRANSFORMATION OF THE PARTICLE MOMENTUM FROM THE LAB. SYSTEM      
C     TO THE CMS (L=1) AND BACKWARD (L=-1)                              
      DIMENSION PX(3),V(3),P(11)                                        
      IF(L .GT. 0) GO TO 1                                              
C     DETERMINATION OF THE INITIAL MOMENTUM COMPONENTS IN CMS           
    1 CALL MOMLAB(P,PX)                                                 
      E=P(9)+P(8)                                                       
      GOTO 2                                                            
 2    VV=SP(V,V)                                                        
      IF(VV.EQ.0.) GO TO 4                                              
      A=SQRT(ABS(1.-VV))                                                
      IF(A .EQ. 0.) A=9.E-18                                            
      B=SP(PX,V)*(1./A-1.)/VV-L*E/A                                     
      DO 3 I=1,3                                                        
 3    PX(I)=PX(I)+V(I)*B                                                
      RETURN                                                            
 4    PX(1)=PX(1)                                                       
       PX(2)=PX(2)                                                      
       PX(3)=PX(3)                                                      
      RETURN                                                            
       END                                                              
      SUBROUTINE  MOM LAB(P,PM)                                         
C     CALCULATION OF THE MOMENTUM COMPONENTS                            
      DIMENSION PM(3),P(11)                                             
      PMOD=PP(P(8),P(9))                                                
       PN=PMOD*P(4)                                                     
      PM(1)=PN*P(7)                                                     
       PM(2)=PN*P(6)                                                    
       PM(3)=PMOD*P(5)                                                  
      RETURN                                                            
       END                                                              
      SUBROUTINE MOM CMS(PM,CT,FI,PC)                                   
C     CALCULATION OF THE PARTICLE MOMENTUM IN CMS                       
      DIMENSION PC(3)                                                   
      ST=SQRT(ABS(1.-CT*CT))                                            
       R=PM*ST                                                          
      PC(1)=R*COS(FI)                                                   
       PC(2)=R*SIN(FI)                                                  
       PC(3)=PM*CT                                                      
      RETURN                                                            
       END                                                              
      FUNCTION PP(T,W)                                                  
C     CALCULATION OF THE MOMENTUM                                       
      PP=SQRT(T*(T+2.*W))                                               
      RETURN                                                            
       END                                                              
      FUNCTION SP(A,B)                                                  
C     SCALAR PRODUCT OF THREE DIMENSIONAL VECTORS                       
      DIMENSION A(3),B(3)                                               
      SP=A(1)*B(1)+A(2)*B(2)+A(3)*B(3)                                  
      RETURN                                                            
       END                                                              
