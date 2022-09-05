      parameter(ns=100000,ni=30)
      integer ID(ns,ni),nimage,nstar,ncoef,filnum(ni),c1,c2,cif(ni,2)
      integer fn(5),m
      real*8 xc(ns),yc(ns),mag(ns,ni),merr(ns,ni)
      real*8 trancf,dzero,t(ni,8),coef(100,6)
      real*8 scol,new(ns,ni),NC
      character fil(ni)*7,p1*14,filter(ni)*7,ch(7)*1
     

      open(21,file='crossid.dat')
      open(22,file='out1.fs')
      open(23,file='ext.cmd')
      open(24,file='trancf.cmd')
      read(21,1) p1,(fil(i),i=1,30)
1     format(a14,30a7)
      read(21,2) p1
2     format(/////a14)
      nimage=0

      do i=1,30
       read(fil(i),'(7a1)') (ch(j),j=1,7)
       k=j-1
       do j=1,k
        if(ch(j).eq.'I') then
         write(filter(i),'(a1)') 'I'
         nimage=nimage+1
         filnum(i)=5
         fn(5)=nimage
         cif(5,2)=nimage
        elseif(ch(j).eq.'V') then
         nimage=nimage+1
         write(filter(i),'(a1)') 'V'
         filnum(i)=3
         fn(3)=nimage
         cif(5,1)=nimage
         cif(2,2)=nimage
         cif(3,2)=nimage
        elseif(ch(j).eq.'B') then
         nimage=nimage+1
         write(filter(i),'(a1)') 'B'
         filnum(i)=2
         fn(2)=nimage
         cif(2,1)=nimage
         cif(3,1)=nimage
         cif(1,2)=nimage
        elseif(ch(j).eq.'U') then
         nimage=nimage+1
         write(filter(i),'(a1)') 'U'
         fn(1)=nimage
         filnum(i)=1
         cif(1,1)=nimage
        elseif(ch(j).eq.'R') then
         nimage=nimage+1
         write(filter(i),'(a1)') 'R'
         fn(4)=nimage
         filnum(i)=4
        elseif(ch(j).eq.'H') then
         nimage=nimage+1
         write(filter(i),'(a1)') 'H'
        elseif(ch(j).eq.'G') then
         nimage=nimage+1
         write(filter(i),'(a1)') 'G'
        elseif(ch(j).eq.'J') then
         nimage=nimage+1
         write(filter(i),'(a1)') 'J'
        elseif(ch(j).eq.'K') then
         nimage=nimage+1
         write(filter(i),'(a1)') 'K'
        elseif(ch(j).eq.'L') then
         nimage=nimage+1
         write(filter(i),'(a1)') 'L'
        endif
       enddo
      enddo

      print *,nimage,' images found.'
3     format(t9,f2.0,t16,f6.4,t24,f6.4,t32,f8.6,t46,f7.3,t56,f11.9,
     *t68,f7.4,t78,f2.0)
      read(23,'(//a14)') p1
      do i=1,100
       read(23,3,end=10) (t(i,j),j=1,8)
      enddo

10    nstar=0
11    format(t8,f2.0,t21,f2.0,t26,f8.2,t34,f8.2,t46,f9.4,t56,f9.4) 
      read(24,'(//a14)') p1
      do i=1,100
       read(24,11,end=20) (coef(i,j),j=1,6) 
      enddo

20    ncoef=i-1

 





      do i=1,ns
       read(21,*,end=100) xc(i),yc(i),(ID(i,j),j=1,nimage)
       read(21,*,end=100) (mag(i,j),j=1,nimage)
       read(21,*,end=100) (merr(i,j),j=1,nimage)
       nstar=nstar+1
      enddo




100   do i=1,nstar
       do j=1,nimage
        call selectcolor(filnum(j),j,icol,c1,c2)
        scol=mag(i,fn(c1))-mag(i,fn(c2))
        if(mag(i,fn(c1)).eq.0.d0.or.mag(i,fn(c2)).eq.0.d0) then
         goto 101 
        endif
        call tr(filter,scol,coef,ncoef,trancf,dzero)
        NC=trancf*scol+dzero
        do k=1,nimage
         if(t(k,1).eq.filnum(j)) then
          new(i,j) = mag(i,j) - (t(k,2)-t(k,3)*scol)*t(k,4)+NC+t(k,5)*
     *    t(k,6)+t(k,7)  
c          print *,i,j,new(i,j)
         endif
        enddo
       enddo
101    nstar=nstar
      enddo


      write(22,111) (fil(i),i=1,nimage)
111   format(//30a10)
      write(22,*) ""

201   format()
202   format()
203   format()
      do i=1,nstar
       write(22,'(30I10)') (ID(i,j),j=1,nimage)
       write(22,'(30f10.3)') (new(i,j),j=1,nimage)
       write(22,'(30f10.3)') (merr(i,j),j=1,nimage)
       write(22,*) ""
      enddo
      close(21)
      close(22)
      close(23)
      close(24)
      stop
      end
 
      subroutine selectcolor(ifilter,nimage,icol,c1,c2)
      parameter(ns=100000,ni=30)
      integer ifilter,icol,nimage,c1,c2
      do i=1,nimage
       if(ifilter.eq.1) then
        icol=1
        c1=1
        c2=2
        elseif(ifilter.eq.2) then 
        icol=2
        c1=2 
        c2=3
        elseif(ifilter.eq.3) then
        icol=2
        c1=2
        c2=3
        elseif(ifilter.eq.5) then
        icol=5   
        c1=3
        c2=5
c        elseif(ifilter.eq.1) then
c        icol=1
c       elseif(ifilter.eq.1) then
c         icol=1
       endif
      enddo
      return
      end         
  
 
      subroutine tr(filter,scol,coef,ncoef,trancf,dzero)
      parameter(ns=100000,ni=30)
      integer ifilter,icol,ncoef
      real*8 scol,trancf,dzero
      real*8 coef(100,6)
      character filter*1
      if(filter.eq.'U') then
       ifilter=1
       elseif(filter.eq.'B') then
       ifilter=2
       elseif(filter.eq.'V') then
       ifilter=3
       elseif(filter.eq.'R') then
       ifilter=4
       elseif(filter.eq.'I') then
       ifilter=5
      endif
      do i=1,ncoef
       if(ifilter.eq.coef(i,1)) then
        if(scol.ge.coef(i,3).and.scol.lt.coef(i,4)) then
         trancf=coef(i,5)
         dzero=coef(i,6)
        endif
       endif
      enddo
 
      return   
      end
 
      
      
