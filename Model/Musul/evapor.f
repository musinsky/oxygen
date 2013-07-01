************************************************************************
      BLOCK DATA                                                          
************************************************************************
*   PURPOSE:     DATA FOR EVAPORATION.                                 * 
************************************************************************
      COMMON /BL1016/CC(6) /BL1017/VK(6)                                  
     */BL1005/ AJ(6)  /BL1006/ ZJ(6)  /BL1008/ DLM(6)                     
     */BL1014/ GAM(6)  /BL1001/ T1Y(130)  /BL1002/ T2XY(200)              
      DATA AJ/1.,1.,2.,3.,3.,4./                                          
      DATA ZJ/0.,1.,1.,1.,2.,2./,                                         
     *  DLM/          8.368,7.569,13.835,15.835,15.817,3.607/,            
     * CC/0.,0.2,0.1,0.07,0.13,0.1/,                                      
     * VK/0.,.7,.77,.8,.8,.83/,                                           
     *  GAM          /1.,1.,3.,3.,3.,2./                                  
      DATA  T1Y/                                                          
     *20.8,15.8,21.,16.8,19.8,16.5,18.8,16.5,18.5,17.2,18.26,15.05,       
     *16.01,12.04,13.27,11.09,12.17,10.26,11.04,8.41,9.79,7.36,8.15,      
     *5.63,5.88,3.17,3.32,.82,1.83,.97,2.33,1.27,2.92,1.61,2.91,1.35,     
     *2.4,.89,1.74,.36,.95,-.65,-.04,-1.73,-.96,-2.87,-2.05,-4.05,-3.4,   
     *-5.72,-3.75,-4.13,-2.42,-2.85,-1.01,-1.33,.54,-.02,1.74,.75,2.24,   
     *1.,1.98,.79,1.54,.39,1.08,0.,.78,-.35,.58,-.55,.59,-.61,.59,-.35,   
     *.32,-.96,-.52,-2.08,-2.46,-3.64,-1.55,-.96,.97,.88,2.37,1.75,2.72,  
     *1.9,2.55,1.46,1.93,.86,1.17,.08,.39,-.76,-.39,-1.51,-1.17,-2.36,    
     *-1.95,-3.06,-2.62,-3.55,-2.95,-3.75,-3.07,-3.79,-3.06,-3.77,-3.05,  
     *-3.78,-3.12,-3.9,-3.35,-4.24,-3.86,-4.92,-5.06,-6.77,-7.41,-9.18,   
     *-10.16,-11.12,-9.76,-9.23,-7.96,-7.65/                              
      DATA   T2XY/                                                        
     *-8.4,-12.9,-8.,-11.9,-9.2,-12.5,-10.8,-13.6,-11.2,-12.2,-12.81,     
     *-15.4,-13.07,-15.8,-13.81,-14.98,-12.63,-13.76,-11.37,-12.38,       
     *-9.23,-9.65,-7.64,-9.17,-8.05,-9.72,-8.87,-10.76,-8.64,-8.89,-6.6,  
     *-7.13,-4.77,-5.33,-3.06,-3.79,-1.72,-2.79,-.93,-2.19,-.52,-1.9,     
     *-.45,-2.2,-1.22,-3.07,-2.42,-4.37,-3.94,-6.08,-4.49,-4.50,-3.14,    
     *-2.93,-1.04,-1.36,.69,.21,2.11,1.33,3.29,2.46,4.3,3.32,4.79,3.62,   
     *4.97,3.64,4.63,3.07,4.06,2.49,3.3,1.46,2.06,.51,.74,-1.18,-1.26,    
     *-3.54,-3.97,-5.26,-4.18,-3.71,-2.1,-1.7,-.08,-.18,.94,.27,1.13,     
     *.08,.91,-.31,.49,-.78,.08,-1.15,-.23,-1.41,-.42,-1.55,-.55,-1.66,   
     *-.66,-1.73,-.75,-1.74,-.78,-1.69,-.78,-1.6,-.75,-1.46,-.67,-1.26,   
     *-.51,-1.04,-.53,-1.84,-2.42,-4.52,-4.76,-6.33,-6.76,-7.81,-5.8,     
     *-5.37,-3.63,-3.35,-1.75,-1.88,-.61,-.9,.09,-.32,.55,-.13,.7,-.06,   
     *.49,-.2,.4,-.22,.36,-.09,.58,.12,.75,.15,.7,.17,1.11,.89,1.85,      
     *1.62,2.54,2.29,3.2,2.91,3.84,3.53,4.48,4.15,5.12,4.78,5.75,5.39,    
     *6.31,5.91,6.87,6.33,7.13,6.61,7.3,6.31,6.27,4.83,4.49,2.85,2.32,    
     *.58,-.11,-.98,.81,1.77,3.37,4.13,5.6,6.15,7.29,7.35,7.95,7.67,      
     *8.16,7.83,8.31,8.01,8.53,8.27/                                      
       END                                                                
************************************************************************
      SUBROUTINE EVAPOR(ENEXT,ATWGHT,CHARGE,PNX,PNY,PNZ,        
     *AM,AMF,RADNCL,FU,KSTART)                                        
************************************************************************
*  PURPOSE:   
      parameter (max1=2000,max2=max1+500)
      COMMON /BL1001/T1Y(130) /BL1002/T2XY(200)                           
      COMMON/BL1011/VJ(6) /BL1015/RJ(7)                                   
     */BL1005/AJ(6)  /BL1006/ZJ(6) /BL1008/DLM(6)                         
     */BL1014/GAM(6) /BL1016/CC(6) /BL1017/VK(6)                          
     */BL1003/U,AI,Z/BL1009/AFJ(7) /BL1010/ZFJ(6)                         
     */WRPR/ WP(9)                                                        
      COMMON/BLANGL/ANGL(4)
      COMMON/OUTRES/W(6,max2)                               
      DIMENSION GJ(7),BJ(7)                                               
      FU=-1.                                                          
      U=ENEXT*1000.                                                       
       AI=ATWGHT                                                          
       Z=CHARGE                                                           
       REMN=940.*AI                                                       
      VNX=(PNX/REMN)*1000.                                                
       VNY=(PNY/REMN)*1000.                                               
       VNZ=(PNZ/REMN)*1000.                                               
      KST1 = KSTART                                                       
      max=KST1+AI
      DO 20 K=KST1,MAX                                                    
       if(ai.le.z) goto 102
      IF(AI-4.0)101,101,100                                               
 100  IF(Z-2.0)101,101,5                                                  
 101  IF(AI-1.0)102,102,103                                               
 102  RETURN                                                              
 103  EP1=U/AI                                                            
      U=U-EP1                                                             
      IF(Z-1.0)104,105,105                                                
 104  EP2=0.0                                                             
      GO TO 106                                                           
 105  EP2=1.0                                                             
 106  Z=Z-EP2                                                             
      EP3=940.0                                                           
      AI=AI-1.0                                                           
       GO TO 107                                                          
    5 DL=DELTAM(AI,Z)                                                     
      DO 6 I=1,6                                                          
      AFJ(I)=AI-AJ(I)                                                     
       ZFJ(I)=Z-ZJ(I)                                                     
      VJ(I) = COLOMB(I,RADNCL,AM,AMF)                                     
      BJ(I)=DELTAM(AFJ(I),ZFJ(I) )-(DL-DLM(I) )                           
 6    RJ(I)=U-(BJ(I)+VJ(I))                                               
      A=AI                                                                
      N1=(A-1.)/2                                                         
        N2=Z/2                                                            
      IF(N1-(A-1.)/2.) 90,91,91                                           
  90  RJ(1)=RJ(1)-12./SQRT (A)                                            
       GO TO 93                                                           
  91  IF(N2-Z/2.) 93,92,92                                                
  92  RJ(1)=RJ(1)-2.*(12./SQRT (A))                                       
  93  CONTINUE                                                            
      X=Z**2/A                                                            
      XT=((X-33.5)**2)**0.333333333                                       
      IF(X-33.5) 71,71,72                                                 
  71  BJ(7)=12.5+4.7*XT**(9./8.)                                          
       GO TO 73                                                           
 72   BJ(7)=12.5-2.7*XT                                                   
  73  CONTINUE                                                            
      L1=(A-Z)/2.                                                         
       L2=Z/2.                                                            
      IF((A-Z)/2.-L1) 75,75,74                                            
  74  BJ(7)=BJ(7)+1.0                                                     
  75  CONTINUE                                                            
      IF(Z/2.-L2) 76,76,77                                                
  76  BJ(7)=BJ(7)-0.5                                                     
  77  CONTINUE                                                            
      IZ=Z                                                                
       IA=A                                                               
      TZ=T1Y(IZ)                                                          
       TN=T2XY(IA-IZ)                                                     
      BJ(7)=BJ(7)-(TZ+TN)                                                 
      BJ(7)=BJ(7)/(1.+SQRT (U/(2.*A)))                                    
      AFJ(7)=A                                                            
       RJ(7)=U-BJ(7)                                                      
      M1=A/2.                                                             
       M2=Z/2.                                                            
      IF(M1-A/2.) 80,81,81                                                
 80   RJ(7)=RJ(7)-1./(1.+SQRT(U/(2.*A)))                                  
       GO TO 83                                                           
  81  IF(M2-Z/2.) 83,82,82                                                
 82   RJ(7)=RJ(7)-2./(1.+SQRT(U/(2.*A)))                                  
  83  CONTINUE                                                            
      CALL ARFA11(PER,AM,AMF)                                             
      DO 7 I=1,7                                                          
      IF(RJ(I))8,8,9                                                      
   8  GJ(I)=0.                                                            
       GO TO 7                                                            
 9    GJ(I)=PAMMA(I,PER,AM,AMF,RADNCL)                                    
   7  CONTINUE                                                            
      G=0.                                                                
       DO 10 I=1,7                                                        
  10  G=G+GJ(I)                                                           
      IF (G)11,11,12                                                      
  11  RETURN                                                              
  12  DO 13 J=2,7                                                         
  13  GJ(J)=GJ(J-1)+GJ(J)                                                 
      B = RNDM(-1)*G                                                      
      DO 14 J=1,7                                                         
      IF(B-GJ(J))15,14,14                                                 
  15  LM=J                                                                
       GO TO 16                                                           
  14  CONTINUE                                                            
  16  IF(LM-7)18,17,18                                                    
  17  FU=1.                                                           
      RETURN                                                              
  18  EP1=TKIN(LM,AM,AMF)                                                 
       EP2=ZJ(LM)                                                         
       EP3=940.*AJ(LM)                                                    
      U=(RJ(LM)-EP1)+VJ(LM)                                               
  19  AI=AFJ(LM)                                  
       Z=ZFJ(LM)                                                          
 107  CONTINUE                                                            
      VPM=SQRT ((2.*EP1)/EP3)                                             
         CALL ISANGL                                                      
      VPX=VPM*ANGL(4)*ANGL(3)                                             
      VPY=VPM*ANGL(4)*ANGL(2)                                             
      VPZ=VPM*ANGL(1)                                                     
      VX=VNX+VPX                                                          
       VY=VNY+VPY                                                         
       VZ=VNZ+VPZ                                                         
      VM=SQRT (VX**2+VY**2+VZ**2)                                         
      COT = VZ/VM                                                         
      SIT=SQRT(1.-COT**2)
      VXY=SQRT(VX**2+VY**2)                                               
      CPHI=VX/VXY                                                    
      SPHI=VY/VXY                                                       
      W3=(EP3*VM**2)/2000.                                            
      W(5,K)=EP2                                                         
      W(4,K)=EP3/1000.                                                   
      PMOD=SQRT(W3*(W3+2.*W(4,k)))
      W(1,K)=PMOD*SIT*CPHI
      W(2,K)=PMOD*SIT*SPHI
      W(3,K)=PMOD*COT
      KSTART=KSTART+1                                                     
      IF(KSTART-max2)20,20,110                                            
  20  CONTINUE                                                            
      RETURN                                                              
  110 PRINT 111                                                           
  111 FORMAT(5X,'MEMORY OVERFLOW IN EVAPOR')                              
      RETURN                                                              
         END                                                              
      SUBROUTINE FUSION(A,Z,EVGEV,PX,PY,PZ,KSTART)                         
      COMMON/BL1000/AM,AMF                                                
      COMMON/BL0999/RADNCL                                                    
      COMMON/BLANGL/ANGL(4)                                               
      COMMON /BL1003/URES,ARES,ZRES                                       
      COMMON /RESFIS/ARFIS1,ZRFIS1,ARFIS2,ZRFIS2                          
      DIMENSION V(3),PA(3),PS(3),P(3)                                     
      DIMENSION A1(10),Z1(10),D1(10),D2(10),CK(10),DM(10),VER(10),UF(10)  
      AMX=0.                                                              
      OST=0.383/(2.0*(A/2.)**0.3333333)        
       EV=EVGEV*1000.0                                                   
      DMA=DELTAM(A,Z)                                                     
      DO 1 K=1,10                                                         
      AK = K-1                                                            
      A1(K)=A/2.+AK*5.                                                    
       A13=A1(K)**0.3333333                                               
      A23=(A-A1(K))**.33333333                                           
      A124=1.0-1.24/(A13**2)                                              
      A224=1.0-1.24/(A23**2)                                              
      Y1=(88.4712/A13)*A124                                               
      Y2=(88.4712/A23)*A224                                               
      E1=125.8024/A1(K)-(156.9424/(A13**4))*A124+                         
     *(0.779/A13)*(1.-1.58/(A13**2))-OST                                  
      E2=125.8024/(A-A1(K))-(156.9424/(A23**4))*A224+                     
     *(0.779/A23)*(1.0-1.58/(A23**2))-OST                                 
      Z1(K)=((Y2-Y1)+Z*E2)/(E1+E2)                                        
       Z2=Z-Z1(K)                                                         
      AL1=6.8*A13*A13 - 0.142*Z1(K)*Z1(K)/A13                             
       AL2=6.8*A23*A23 - 0.142*Z2*Z2/A23                                  
       AM1=12.138*A13*A13 - 0.14484*Z1(K)*Z1(K)/A13                       
      AM2=12.138*A23*A23-0.14484*Z2*Z2/A23                                
          ETA=A23/A13                                                     
       ET2=(1.+ETA)**2                                                    
       ET4=ET2**2                                                         
       EF=1.+1.464*AL1/AM1+ETA*ETA*(AL1/AL2+1.464*AL1/AM2)                
      R1=5.76*Z1(K)*Z2*0.6147*(1.+ETA)*EF                                 
      R0=1.35                                                             
      R2=AL1*A13*R0                                                       
      PD=SQRT(ET4+R1/R2)                                                  
      ALF21=(PD-ET2)/(.4*31.4506*(1.+ETA)*0.784)                          
      ALF22=ETA*AL1*ALF21/AL2                                             
      ALF31=1.21*AL1*ALF21/AM1                                            
      ALF32=ETA*1.21*AL1*ALF21/AM2                                        
      D1(K)=AL1*ALF21*ALF21+AM1*ALF31*ALF31                               
      D2(K)=AL2*ALF22*ALF22+AM2*ALF32*ALF32                               
      CK(K)=(1.07*Z1(K)*Z2)/(A13*(1.+ALF21+ALF31)+A23*(1.+ALF22+ALF32))   
      DM(K)=DMA-DELTAM(A1(K),Z1(K))-DELTAM(A-A1(K),Z-Z1(K))               
      UF(K)=EV+DM(K)-D1(K)-D2(K)-CK(K)                                    
      IF(UF(K))2,2,4                                                      
    2 KMAX=K-1                                                            
       IF(KMAX)190,190,3                                                  
  4   VER(K)=((4.*Z1(K)*Z2)/(Z*Z))*((UF(K)/UF(1))**(11./4.))*             
     1 EXP(2.0*(SQRT(AM*A))*(SQRT(UF(K))-SQRT(UF(1))))                    
      IF(VER(K).GE.AMX) AMX=VER(K)                                        
  1   CONTINUE                                                            
      KMAX=10                                                             
  3    DEA=(A1(KMAX)-A1(1))*RNDM(-1)                                      
  5   AP=A/2.+DEA*RNDM(-1)                                                
      W=SUB(AP,A1,VER,KMAX)                                               
      BR=AMX*RNDM(-1)                                                     
      IF(BR-W) 6,6,3                                                      
  6   AP=AINT(AP+RNDM(-1))                                                
      AS=A-AP                                                             
      ZP=SUB(AP,A1,Z1,KMAX)                                               
      ZP=AINT(ZP+RNDM(-1))                                                
      ZS=Z-ZP                                                             
      TKUL=SUB(AP,A1,CK,KMAX)                                             
      DF1=SUB(AP,A1,D1,KMAX)                                              
      DF2=SUB(AP,A1,D2,KMAX)                                              
      EVT=SUB(AP,A1,UF,KMAX)                                              
      EVT=ABS(EVT)                                                        
      TK1=AP*TKUL/A                                                       
       TK2=TKUL-TK1                                                       
      EV1=DF1+AP*EVT/A                                                    
       EV2=DF2+AS*EVT/A                                     
      T1G=TK1/1000.0                                                     
      P1=SQRT(T1G*(T1G+2.*AP))                                            
      CALL ISANGL                                                         
      PA(1)=P1*ANGL(4)*ANGL(3)                                            
       PA(2)=P1*ANGL(4)*ANGL(2)                                           
      PA(3)=P1*ANGL(1)                                                    
      PS(1)=-PA(1)                                                        
       PS(2)=-PA(2)                                                       
       PS(3)=-PA(3)                                                       
      EVA=EV1/1000.                                                       
      P(1)=PA(1)+PX                                                       
       P(2)=PA(2)+PY                                                      
       P(3)=PA(3)+PZ                                                      
      CALL EVAPOR(EVA,AP,ZP,P(1),P(2),P(3),AM,AMF,RADNCL,FU,KSTART)   
      ARFIS1=ARES                                                         
      ZRFIS1=ZRES                                                         
      EVA=EV2/1000.0                                                      
      P(1)=PS(1)+PX                                                       
       P(2)=PS(2)+PY                                                      
       P(3)=PS(3)+PZ                                                      
      CALL EVAPOR(EVA,AS,ZS,P(1),P(2),P(3),AM,AMF,RADNCL,FU,KSTART)   
      ARFIS2=ARES                                                         
      ZRFIS2=ZRES                                                         
190   RETURN                                                              
      END                                                                 
      FUNCTION DELTAM(X,Y)                                                
C     CALCULATION OF MASS DEFECT                                          
      COMMON  /BL1001/T1Y(130)  /BL1002/T2XY(200)                         
	if(x-y)1,1,2
1	deltam=0.
	return
2      ES=(25.8357-44.2355*((X-2.*Y)/X)**2)*                               
     *(((1.-.62025/X**.66667)**2)*(X**.66667))                            
      EC=.779*((Y*(Y-1.))/X**.33333)*                                     
     *((1.-1.5849/X**.66667)+(1.2273/X+(1.5772/X**1.33333)))              
      EALFA=(-.4323*((Y**1.33333)/(X**.33333)))*                          
     *((1.+.49597/X)-(.57811/X**.33333+(.14518/X**.66667)))               
      EDOB=((8.367*X+31.4506*((X-2.*Y)**2/X))-(.783*Y+17.0354*X))         
      I=X                                                                 
       J=Y                                                                
       L=I-J                                                              
       T1=T1Y(J)                                                          
       T2=T2XY(L)                                                         
      DELTAM=((ES+EC)+(EALFA+EDOB))+(T1+T2)                               
      RETURN                                                              
       END                                                                
      FUNCTION TKIN(L,AM,AMF)                                             
C     KINETIC ENERGY FOR PARTICLES IN EQUILIBRIUM DECAY                   
      COMMON /BL1009/ AFJ(7) /BL1015/RJ(7) /BL1011/ VJ(6)                 
      RB=4.*AM*AFJ(L)*RJ(L)                                               
    5 B1=RNDM(-1)                                                         
      RK=1.+(1./SQRT (RB))*ALOG(B1+(1.-B1)*EXP (-SQRT (RB)))              
      IF(L-1)1,2,1                                                        
    2 BETA=(2.12/AFJ(1)**0.66667-0.05)/(0.76+2.2/AFJ(1)**0.33333)         
      Q1=1.+BETA/RJ(1)                                                    
       Q2=Q1*SQRT (Q1)                                                    
      FRK=(((3.*SQRT (3.))/2.)/Q2)*(Q1*RK-RK**3)                          
       GO TO 3                                                            
    1 FRK=((3.*SQRT (3.))/2.)*(RK-RK**3)                                  
       GO TO 3                                                            
    3 B2=RNDM(-1)                                                         
       IF(B2-FRK)4,4,5                                                    
    4 TKIN=RJ(L)*(1.-RK**2)+VJ(L)                                         
      RETURN                                                              
      END                                                                 
      FUNCTION PAMMA(J,PER,AM,AMF,RADNCL)                                 
C     PROBABILITIES FOR EQUILIBRIUM PARTICLE EMISSION                     
      COMMON/BL1009/AFJ(7) /BL1015/RJ(7)                                  
     */BL1014/GAM(6)/BL1016/CC(6)                                         
      IF(J-7)4,10,4                                                       
    4 IF(J-1)1,2,1                                                        
    2 ALFA=.76+2.2/AFJ(1)**.33333                                         
      BETA=(2.12/AFJ(1)**.666 7-.05)/ALFA                                 
       GO TO 3                                                            
    1 ALFA=1.+CC(J)                                                       
       BETA=0.                                                            
       GO TO 3                                                            
    3 Q1=AM*AFJ(J)                                                        
       Q2=Q1*RJ(J)                                                        
      Q3=(GAM(J)*AFJ(J)**.66667)*(ALFA/Q1**2)                             
      Q4=(2.*BETA*Q1-3.)/2.+Q2                                            
      Q5=(2.*BETA*Q1-3.)*(SQRT (Q2)-.5)+2.*Q2                             
      PAMMA=Q3*(Q4*EXP (-PER)+Q5*EXP (2.*SQRT (Q2)-PER))                  
      RETURN                                                              
   10 Q1=2.*SQRT (AMF*AFJ(7)*RJ(7))                                       
      Q2=33.4/(3.1415926536*RADNCL**2)                                    
      PAMMA=(Q2/(AMF*AFJ(7)))*((Q1-1.)*EXP (Q1-PER)+EXP (-PER))           
      RETURN                                                              
      END                                                                 
      SUBROUTINE ARFA11(PER,AM,AMF)                                       
C     RENORMALISATION OF EMISSION PROBABILITIES                           
      COMMON /BL1009/AFJ(7) /BL1015/RJ(7)                                 
      SFIX=30.                                                            
      SMX = 0.                                                            
      DO 101 K=1,7                                                        
      IF(RJ(K))2,2,3                                                      
    2 Q8=0.                                                               
       GO TO 1                                                            
    3 IF(K-7)4,5,4                                                        
    4 Q8=2.*SQRT (AM*AFJ(K)*RJ(K))                                        
       GO TO 1                                                            
    5 Q8=2.*SQRT (AMF*AFJ(7)*RJ(7))                                       
       GO TO 1                                                            
    1 IF(SMX-Q8) 100,100,101                                              
  100 SMX=Q8                                                              
  101 CONTINUE                                                            
      IF(SMX-SFIX)6,6,7                                                   
    6 PER=0.                                                              
       RETURN                                                             
    7 PER=SMX-SFIX                                                        
       RETURN                                                             
      END                                                                 
      FUNCTION COLOMB (L,RADNCL,AM,AMF)                                   
C     CALCULATION OF COLOMB ENERGY                                        
      COMMON /BL1006/ZJ(6) /BL1010/ZFJ(6)                                 
     */BL1009/AFJ(7) /BL1005/AJ(6) /BL1017/VK(6) /BL1003/U,AI,Z           
      IF(L-1)1,1,2                                                        
    1 COLOMB=0.                                                           
       RETURN                                                             
    2 TEMP1=VK(L)*(1.44/RADNCL)                                           
      COLOMB=TEMP1*((ZJ(L)*ZFJ(L))/(AJ(L)**.33333+AFJ(L)**.33333))        
      COLOMB=COLOMB*(1.-U/(81.*AI*AM))                                    
      IF(COLOMB)3,3,4                                                     
    3  COLOMB=0.                                                          
       RETURN                                                             
    4 RETURN                                                              
      END                                                                 
      SUBROUTINE ISANGL                                                   
C     CHOOSE ISOTROPIC DISTRIBUTED ANGLE                                  
      COMMON /BLANGL/ ANGL(4)                                             
      ANGL(1)=1.-2.*RNDM(-1)                                              
       ANGL(4)=SQRT(1.-ANGL(1)**2)                                        
      F=2.*3.1415926536*RNDM(-1)                                          
       ANGL(2)=SIN(F)                                                     
       ANGL(3)=COS(F)                                                     
      RETURN                                                              
       END                                                                
      FUNCTION SUB(U,E,F,N)                                               
      DIMENSION E(N),F(N)                                                 
      IF(U-E(1))1,1,2                                                     
    1 X1=E(1)                                                             
       X2=E(2)                                                            
      X3=E(3)                                                             
      Y1=F(1)                                                             
      Y2=F(2)                                                             
      Y3=F(3)                                                             
      GOTO 7                                                              
    2 IF(U-E(N-1))3,4,4                                                   
    4 X1=E(N-2)                                                           
      X2=E(N-1)                                                           
      X3=E(N)                                                             
      Y1=F(N-2)                                                           
      Y2=F(N-1)                                                           
      Y3=F(N)                                                             
      GOTO 7                                                              
    3 DO 5 J=2,N                                                          
      IF(U-E(J))6,5,5                                                     
    6 X1=E(J-1)                                                           
      X2=E(J)                                                             
      X3=E(J+1)                                                           
      Y1=F(J-1)                                                           
      Y2=F(J)                                                             
      Y3=F(J+1)                                                           
      GOTO 7                                                              
    5 CONTINUE                                                            
    7 SUB=Y1*(((U-X2)*(U-X3))/((X1-X2)*(X1-X3)))+                         
     * Y2*(((U-X1)*(U-X3))/((X2-X1)*(X2-X3)))+                            
     * Y3*(((U-X1)*(U-X2))/((X3-X1)*(X3-X2)))                             
      RETURN                                                              
      END                                                                 
