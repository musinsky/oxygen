************************************************************************
      program CASCADE                                                  
************************************************************************
*  AUTHOR:  GENIS MUSULMANBEKOV                                        *
*   UPDATED  1 JANUARY 2001                                            *
*  E-MAIL:  genis@jinr.ru                                              *
*                                                                      *
*  ORGANIZATION: JOINT INSTITUTE FOR NUCLEAR RESEARCH, LIT             *
*                141980 DUBNA, RUSSIA                                  *   
************************************************************************
*    INPUT QUANTATIES:                                                 *
*       NSTAR    - CONTROL PARAMETER: =-1 - PRINT THE TREE OF CASCADES;*
*                  =0 - WRITE SIMULATED EVENTS IN FILE;                *
*                  >0 - APPEND OF CONTINUATED SIMULATION TO EXISTING   *
*                   FILE; THIS TIME NSTAR IS THE NUMBER OF SKIPPED     *
*                   EVENTS.                                            *   
*       JTY      - CONTROL PARAMETER FOR NUCLEUS MODEL:                *    
*                  =0 - FERMI-GAS MODEL;                               *
*                  =1 - FCC LATTICE MODEL                              *
*       A        - ARRAY OF ATOMIC NUMBERS OF COLLIDING NUCLEI         *
*       Z        _ ARRAY OF THEIR CHARGES                              *
*       WT       - KINETIC ENERGY OF PROJECTILE NUCLEUS IN LAB. [GeV]  *   
*       BIMP     - >=0 - INPUT VALUE OF IMPACT PARAMETER [FM];         *
*                  <0 - RANDOM VALUE OF IMPACT PARAMETER THAT IS       *
*                  GIVEN BY PROGRAM.                                   *
*----------------------------------------------------------------------*
*    OUTPUT QUANTITIES:                                                *
*       NCAS     - =5*KEV, WHERE KEV IS THE MULTIPLICITY IN EVENT.     *
*       PIMPAC   - IMPACT PARAMETER IN COLLISION [FM].                 *
*       E(1)     - EXCITATION ENERGY OF PROJECTILE NUCLEUS.            *
*       E(2)     - EXCITATION ENERGY OF TARGET NUCLEUS.                *
*       MVS1     - J=MVS1 IS THE INDEX NUMBER IN ARRAY W(6,J)          *
*                  FROM WHICH PARTICES EMITTED BY EXCITED              *
*                  PROJECTILE NUCLEUS ARE WRITTEN.                     *
*       MVS2     - J=MVS2 IS THE SAME ONE FOR PARTICLES EMITTED        *
*                  BY TARGET NUCLEUS.                                  *
*       W(1,J)   -                                                     *
*       W(2,J)   -   MOMENTUM COMPONENTS OF PRODUCED PARTICLE.         *
*       W(3,J)   -                           *
*       W(4,J)   - IT'S MASS [GEV].                                    *
*       W(5,J)   - IT'S CHARGE.                                        *
************************************************************************ 
      parameter (max1=2000,max2=max1+500)
      common/inel/inel
      COMMON/DATIN/A(2),Z(2),WT,VPI,EPS(2),AR(2),CR(2),DR(2),R(2),TR(2) 
      COMMON/BL1003/UI,AI,ZI                                            
      COMMON/BL1000/AM,AMF                                                
      COMMON/BL0999/RADNCL                                                    
      COMMON /LNUC/ BR(2)                                               
      COMMON/XYZIZ/XYZ(2,3,240),IZ(2,240)                               
      COMMON/RESCAS/AC(2),ZC(2),E(2),PC(2,3),AMC(2,3),PIMPAC            
      COMMON/KYPT/KPT,KYP,KYT,KPT1,KYP1,KYT1,KS1,LST                    
      COMMON/KSALL/KSALL,KPTALL,KPT0
      COMMON /PMIM/PM(19,max1),IM(5,max1)                               
      COMMON/OUTRES/W(6,max2)                                            
      DIMENSION WW(max2)                                                
      DIMENSION PCE(3),P(9),V0(3)                                       
      EQUIVALENCE (W(1,1),WW(1))
*  INPUT PARAMETERS FOR NUCLEI
      CALL PARAM(DEL,tpr,PDL)
    1 CONTINUE                                                          
      KPT=0                                                             
      KYP=0                                                             
      KYT=0                                                             
      KPT1=0                                                            
      KYP1=0                                                            
      KYT1=0
      KS1=0
      KPTALL=0
      open(unit=2,file='oxygen.in')
      open(unit=1,file='musul_model.asc',status='REPLACE')
c free format 9/July/2003
cuj      READ(2,11)NSTAR,JTY                                                
cuj  11  FORMAT(2I5)                                                     
cuj      READ (2,2) A,Z,WT,LIMC,BIMP                                        
cuj   2  FORMAT(4F4.0,F8.3,I6,F7.3)                                       
      READ(2,*)NSTAR,JTY                                                
      READ (2,*) A,Z,WT,LIMC,BIMP      
      PRINT*,'              ****** FULL COMPONENT CASCADE MODEL ******'
   88 PRINT 3,A(1),Z(1)                                                 
    3 FORMAT(/10X,'PROJECTILE   A = ',F4.0,'  Z =',F4.0/)                 
      PRINT 4, WT                                                       
    4 FORMAT(25X,'ENERGY OF PROJECTILE WT = ',F10.3,' GEV/NUCLEON')    
      PRINT 5,A(2),Z(2)                                                 
    5 FORMAT(/14X,'TARGET   A = ',F4.0,'  Z =',F4.0/)                     
      IF(NSTAR) 56,52,55
cuj   52 WRITE (1) A,Z,WT
   52 WRITE (1, *) NSTAR,JTY
      WRITE (1, *) A,Z,WT,LIMC,BIMP
      GOTO 56
cuj   55 READ(1)A,Z,WT                                                     
   55 READ (1) NSTAR,JTY
      READ (1) A,Z,WT,LIMC,BIMP
      DO 57 II=1,NSTAR                                                  
      RSTART=RNDM(-1)
cuj   57 READ(1)NCAS,PIMPAC,(WW(K),K=1,NCAS)                     
   57 READ(1)NCAS,PIMPAC,MVS1,MVS2,(WW(K),K=1,NCAS)                     
   56 CONTINUE                                                          
      IF(A(1).LE.1..AND.A(2).LE.1.) GOTO 70
      IF(A(1)-1.1) 37,37,6                                              
   37 R(1)=0.
      TR(1)=0.
	GOTO 77
*** CALCULATION OF FERMI ENERGY 'TR(1)' IN PROJECTILE NUCLEUS ***
   6  CALL TFERMN(1.12,A(1),AR(1),CR(1),DR(1),TR(1),R(1),BR(1))         
   77 IF(A(2)-1.)78,78,7                                                
   78 TR(2)=0.                                                          
      R(2)=0.
      GOTO 79                                                           
   70 R(1)=0.
      R(2)=0.
	R12=0.
      GOTO 79
***    CALCULATION OF FERMI ENERGY 'TR(2)' IN TARGET NUCLEUS   ***
   7  CALL TFERMN(1.12,A(2),AR(2),CR(2),DR(2),TR(2),R(2),BR(2))         
   79 RN=.4                                                             
	R12=R(1)+R(2)
      AIN=A(1)+A(2)
      V0(1)=0.D0                                                        
      V0(2)=0.D0                                                        
      V0(3)=SQRT(WT*(WT+1.88D0))/(WT+.94D0)                             
          INEL=0                                                        
      IEL=0                                                             
      MM=0        
      ncentr=0                                                      
*
***    LOOP OVER 'LIMC' INELASTIC NUCLEUS-NUCLEUS COLLISIONS  ***
* 
      DO 1000 M=1,10000000                                               
      IF(R12.LE.0.) GOTO 31
***     CALCULATION OF NUCLEON COORDINATES IN NUCLEI   ***              
   32 CALL BLPINN(JTY,M,RN)                                                   
   31 DO 20 I=1,5                                                       
      DO 20 J=1,max1                                                    
   20 W(I,J)=0.                                                         
      EX1=0.
      EX2=0.
      caspr=0.
      IF(R12.LE.0.) GOTO 62
      DO 34 I=1,2                                                       
      IF(A(I)-1.) 51,51,34                                              
   51 IF(I-2)34,33,33                                                 
   33 XYZ(2,1,1)=0.                                                     
      XYZ(2,2,1)=0.                                                     
      XYZ(2,3,1)=0.                                                     
      IZ(2,1)=Z(2)                                                         
   34 CONTINUE                                                          
*    
***   main subroutine  for cascade simulation ***
*
c      if(m.gt.162) goto 620
      IF(NSTAR) 620,62,62
  620 CALL BLTREE(NEL,RN,DEL,TPR,MV,BIMP,PDL)
      GOTO 69
   62 CALL BLNUNU(NEL,RN,DEL,TPR,Mv,BIMP,PDL)                               
   69 IF(NEL) 9,9,8                                                     
   8  IEL=IEL+1                                                         
      GO TO 1000                                                        
   9  MM=MM+1                                                           
      ztot=0.
      atot=0.
      DO 22 K=1,MV
      CALL MOM(K)
      W(1,K)=PM(1,K)                                                    
      W(2,K)=PM(2,K)                                                    
      W(3,K)=PM(3,K)                                                    
      W(4,K)=PM(9,K)                                                    
      W(5,K)=IM(1,K)
      ztot=ztot+w(5,k)
      if(w(4,k).gt..9) atot=atot+1.
   22 CONTINUE
      KEv=Mv
      MVS1=0
      MVS2=0
      IF(R12.LE.0.) GOTO 40
      KEV=KEV+1                                                         
      MVS2=KEV
      IF(E(2).LT..001.AND.AC(2).GT.0.) GOTO 93
   42 IF(AC(2)-1.) 98,93,41                                              
   98 MVS2=0
      GOTO 90
   93 AI=AC(2)
      ZI=ZC(2)
      GOTO 95
*
***  EVAPORATION OF EXCITED RESIDUAL TARGET NUCLEUS  ***
*
  41  CALL EVAPOR(E(2),AC(2),ZC(2),PC(2,1),PC(2,2),PC(2,3),AM,AMF,   
     *RADNCL,FU,KEV)                                                           
      IF(AC(2).LT.100.) GOTO 99
      IF(FU.GT.0.) THEN
       CALL FUSION(AC(2),ZC(2),E(2),PC(2,1),PC(2,2),PC(2,3),KEV)
      ENDIF                                                              
 99   IF(MVS2-KEV) 955,93,93
 955  IF(AI) 90,90,95
   95 DO 96 I=1,3
   96 PCE(I)=PC(2,I)*(AI/AC(2))
      W(1,KEV)=PCE(1)
      W(2,KEV)=PCE(2)
      W(3,KEV)=PCE(3)
      W(4,KEV)=AI*.94
      W(5,KEV)=ZI
      KEV=KEV+1
   90 MVS1=KEV                                                           
      IF(E(1).LT..001.AND.AC(1).GT.0.) GOTO 81
      IF(AC(1)-1.)43,81,83                                              
   43 KEV=KEV-1
      MVS1=0
      GOTO 40
   81 AI=AC(1)                                                          
      ZI=ZC(1)                                                          
      GOTO 85                                                           
*
***   EVAPORATION OF EXCITED RESIDUAL PROJECTILE NUCLEUS   ***
*
 83   CALL EVAPOR(E(1),AC(1),ZC(1),PC(1,1),PC(1,2),PC(1,3),AM,AMF,
     *RADNCL,FU,KEV)                                       
      IF(AC(1).LT.100.) GOTO 455
      IF(FU.GT.0.) THEN
       CALL FUSION(AC(1),ZC(1),E(1),PC(1,1),PC(1,2),PC(1,3),KEV)
      ENDIF                                                              
  455 IF(MVS1-KEV)456,81,81
 456  IF(AI) 45,45,457
   45 KEV=KEV-1                                                         
 457  DO 91 L=MVS1,KEV                                                   
      PCE(1)=W(1,L)                                     
      PCE(2)=W(2,L)                                     
      PCE(3)=W(3,L)                                                 
      P(9)=W(4,L)                                                       
      CALL CINEMA(PCE,V0,P)                                             
      PMOD=SQRT(P(8)*(P(8)+2.*P(9)))
      W(1,L)=PMOD*P(4)*P(7)                   
      W(2,L)=PMOD*P(4)*P(6)                
      W(3,L)=PMOD*P(5)                                                        
   91 CONTINUE                                                          
   92 IF(AI.LT.1.) GOTO 40                                              
c      KEV=KEV+1
C   85 CONTINUE
   85 DO 84 I=1,3                                                       
   84 PCE(I)=PC(1,I)*(AI/AC(1))                                         
      P(9)=AI*.94                                                       
      CALL CINEMA(PCE,V0,P)                                             
c       IF(P(7).LT.-1.)P(7)=-.999                                       
c       IF(P(7).GT.1.) P(7)=.999                                         
      PMOD=SQRT(P(8)*(P(8)+2.*P(9)))
       W(1,KEV)=PMOD*P(4)*P(7)                   
       W(2,KEV)=PMOD*P(4)*P(6)                
       W(3,KEV)=PMOD*P(5)                           
       W(4,KEV)=P(9)                                                      
       W(5,KEV)=ZI                                                      
*       W(5,KEV)=ZI*10.0 + PM(18,KEV)
   40 CONTINUE
      PROT=0.
      nprot=0
      ztot=0.
      DO 400 K=1,KEV
      isel=0     
      iprot=0                                
      ztot=ztot+w(5,k)
      LPR=W(4,K)/.94
      cpr=lpr
      PROT=PROT+CPR
      PMOD=SQRT(W(1,K)**2+W(2,K)**2+W(3,K)**2)
      if(pmod.gt..15.and.pmod.lt..75) isel=1
      if(w(4,k).eq..94.and.w(5,k).eq.1) iprot=1
      if(iprot.gt.0.and.isel.gt.0) nprot=nprot+1
  400 CONTINUE
c      IF(AIN.NE.PROT) GOTO 1000
c***********************
      zdif=z(1)+z(2)-ztot
c      if(zdif.ne.0.) print*,'AFTER EVAPOR zdif=',zdif
C**********************
      IF(KEV.EQ.0) GOTO 1000
      IF(A(1).LE.1..AND.KEV.LE.2) GOTO 1000
      IF(A(2).LE.1..AND.KEV.LE.2) GOTO 1000
      IF(A(1).LE.1..AND.A(2).LE.1.) PIMPAC=0.
      INEL=INEL+1
*      NCAS=5*KEV                                                        
      NCAS=6*KEV
	do i = 1, KEV
            W(6,KEV)= PM(18,KEV)
	enddo
      write(1, *)ncas,pimpac,mvs1,mvs2,(ww(k),k=1,ncas)
      if(nprot-5) 63,631,631
  631 ncentr=ncentr+1
c       WRITE(1)NCAS,(WW(K),K=1,NCAS)
c      print 9999,mrcntr,kev
   63 CONTINUE                                                          
      MR=INEL+NSTAR                                                     
      mrcntr=ncentr+nstar
***    END OF LOOP   ****
      IF(mr.GE.LIMC) GOTO 2000
c      PRINT 9999,m,KEV
 9999 FORMAT(10X,'EVENT NO = ',I10,10X,'MULTIPLICITY = ',I5)
 1000 CONTINUE                                                          
 2000 CONTINUE                                                          
      DLAM=1./(5.06*SQRT(WT*(WT+1.88)))+del                             
       SIG=31.459*((R(1)+R(2)+DLAM)**2)                                 
      SIGIN=M        
      tinel=inel
      centrfr=ncentr/tinel                                                   
*   INELASTIC CROSS SECTION 
      SIGIN=(SIG*INEL)/SIGIN                                            
      PRINT 10, INEL,centrfr,SIGIN                                              
   10 FORMAT(/10X,'TOTAL EvENTS = ',I10/10X,
     * 'FRACTION OF CENTRAL EVENTS = ',F8.5/10X,
     * 'INEL.CROSS SECTION = ', F9.4/)
c      call histdo
      ENDFILE 1
      STOP
      END
