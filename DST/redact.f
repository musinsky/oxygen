c     last changes: 26 oct 2003 (01 jul 2013)
      program TransformDST
      common/bc/b(500),c(1000)/nw/nw,nwc,nk,nkon,nkan,p0(5)
      open(11,file='oxygen_DST.asc',status='REPLACE')
      nev=0
      nk=5
      nkan=1
      nkon=0
 1    continue
      call redact
      if(nkon.eq.1) go to 77
      nev=nev+1
      write(11, *) nwc,(c(i),i=1,nwc)
      go to 1
 77   continue
      write (*,*) 'Number of events = ', nev
      stop
      end
cccccccccccccccccccccccccccccccccccccccc
      subroutine redact
      common/bc/b(500),c(1000)/nw/nw,nwc,nk,nkon,nkan,p0(5)
      data p0/3.25,3.34,3.19,3.3,3.31/
      open(1,file='kis57.fbin',form='UNFORMATTED',status='OLD')
      open(2,file='dub57.fbin',form='UNFORMATTED',status='OLD')
      open(3,file='kis61.fbin',form='UNFORMATTED',status='OLD')
      open(4,file='kis65.fbin',form='UNFORMATTED',status='OLD')
      open(5,file='dub61.fbin',form='UNFORMATTED',status='OLD')
 1    ik=nkan
      read(nkan,end=77) nw,(b(i),i=1,nw)
      pn=p0(nkan)
      do 2 i=1,1000
         c(i)=0.
 2    continue
      b(7)=-b(7)
      b(8)=b(8)-3.1415926
      sina0=sin(b(7))
      cosa0=cos(b(7))
      nch=b(2)/100
      nm=(b(2)-nch*100)/10
      npl=b(1)/10000.
      nkad=(b(1)-npl*10000.)/10.
      c(18)=b(6)                !!!
      c(19)=b(9)                !!!
      c(1)=b(1)-float(int(b(1)/10.)*10)
      c(2)=int(b(2)/100.+0.01)
      c(3)=b(3)
      c(4)=b(4)
      c(5)=b(5)
      c(15)=nm
      c(20)=nch
      nz=21
      do 3 i=11,nw,5
         k=b(i+4)/1000
         if(k.ge.2.and.k.le.8) c(6)=c(6)+1.
         if(k.eq.1) c(7)=c(7)+1.
         if(k.eq.2) c(8)=c(8)+1.
         if(k.eq.3) c(9)=c(9)+1.
         if(k.eq.4) c(10)=c(10)+1.
         if(k.eq.5) c(11)=c(11)+1.
         if(k.eq.6) c(12)=c(12)+1.
         if(k.eq.7) c(13)=c(13)+1.
         if(k.eq.8) c(14)=c(14)+1.
         if(k.eq.10) c(16)=c(16)+1.
         if(k.eq.11) c(17)=c(17)+1.
         dl=b(i+4)-k*1000
cccccccccccccccccccccccccccccccccccccccc
c     opravy na hlbyne uhly fragmentov v zavislosti na ich naboj,
c     dlzky zmeranej drahy a cisla zalivky
         if(ik.eq.1.and.k.eq.11) b(i+1)=b(i+1)+(8.8-0.16*dl)*0.001
         if(ik.eq.2.and.k.eq.11) b(i+1)=b(i+1)+(13.2-0.24*dl)*0.001
         if(ik.eq.3.and.k.eq.11) b(i+1)=b(i+1)+0.008
         if(ik.eq.4.and.k.eq.11) b(i+1)=b(i+1)+0.010
         if(ik.eq.5.and.k.eq.11) b(i+1)=b(i+1)+0.009
         if(ik.eq.1.and.(k.ge.2.and.k.le.8)) b(i+1)=b(i+1)+
     +        (12.605-0.53*(k-2)+(0.0076*(k-2)-0.177)*dl)*0.001
         if(ik.eq.1.and.(k.ge.2.and.k.le.8))b(i+2)=b(i+2)-0.001
         if(ik.eq.2.and.(k.ge.2.and.k.le.8)) b(i+1)=b(i+1)+
     +        (14.0406+0.244*(k-2)-(0.008*(k-2)+0.192)*dl)*0.001
         if(ik.eq.2.and.(k.ge.2.and.k.le.8))b(i+2)=b(i+2)-0.0007
         if(ik.eq.3.and.(k.ge.2.and.k.le.8))
     *        b(i+1)=b(i+1)+(6.95+0.3*(k-2))*0.001
         if(ik.eq.3.and.(k.ge.2.and.k.le.8))b(i+2)=b(i+2)+.00015
         if(ik.eq.4.and.(k.ge.2.and.k.le.8))
     *        b(i+1)=b(i+1)+(11.275+0.2*(k-2))*0.001
         if(ik.eq.4.and.(k.ge.2.and.k.le.8))b(i+2)=b(i+2)+0.00022
         if(ik.eq.5.and.(k.ge.2.and.k.le.8)) b(i+1)=b(i+1)+0.0074
         if(ik.eq.5.and.(k.ge.2.and.k.le.8))b(i+2)=b(i+2)+0.0002
cccccccccccccccccccccccccccccccccccccccc
         sina=sin(b(i+1))
         cosa=cos(b(i+1))
         db=b(i+2)-b(8)
         bl=cosa*sin(db)
         bm=cosa0*cosa*cos(db)+sina0*sina
         bn=-sina0*cosa*cos(db)+cosa0*sina
c     bln=sqrt(bl*bl+bn*bn)
c     pci=acos(bl/bln)
c     if(bn.lt.0.) pci=6.2831852-pci
         if(k.ge.2.and.k.le.8) b(i)=3.25*b(i)/pn
         if(k.eq.11) b(i)=3.25*b(i)/pn
         c(nz)=k
         c(nz+1)=dl
         c(nz+2)=bl*b(i)
         c(nz+3)=bm*b(i)
         c(nz+4)=bn*b(i)
         c(nz+5)=b(i+3)
         nz=nz+6
 3    continue
      nwc=nz-1
      return
 77   continue
      nkan=nkan+1
      if(nkan.le.nk)go to 1
      nkon=1
      return
      end
