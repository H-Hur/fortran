!     For fortran95.
!     To extract coordinates from als file.
!     18,Jan,2008, V1.00, apco.f, apco.exe, by Hyeon-oh Hur.
!     f95 -o apco ~/fortran/subroutine/daotools/apco.f
!     gfortran -o apco ~/work/fortran/subroutine/daotools/apco.f
      parameter (ns=100000)
      integer ID(ns),nap,nstar,napstar
      real*8 xc(ns),yc(ns),fwhm,aprad,d,mag(ns),dat(8)
      real*8 axc(ns),ayc(ns),dx,dy
      character input*30,test*10       
      print *,'Input psfs file=?'
      read(*,'(a30)') input    
c      open(101,file='apc.coo')
      dx=0.d0
      dy=0.d0
c      print *,'dx,dy=?'
c      read(*,*) dx,dy
     
      write(test,'(a10)') 'apc.coo        '
      open(101,file=test)

       
      call reading(xc,yc,mag,ID,dat,nstar) 
      call readapc(input,axc,ayc,nap)

      napstar=0
      do i=1,nstar
       do j=1,nap
        d=sqrt((xc(i)-dx-axc(j))**2.d0+(yc(i)-dy-ayc(j))**2.d0)
        if(d.le.dat(8)) then
         write(101,'(2f9.3)') xc(i),yc(i)
         napstar=napstar+1
        endif
       enddo
      enddo
      print *,napstar,'stars found.'
      close(101)
      stop
      end



      subroutine readapc(input,axc,ayc,nap)
      parameter (ns=100000,ni=30)
      character input*30,temp(20)*1,title*7
      integer nhead,nap
      real*8 ann,dann,skip,axc(ns),ayc(ns)
      open(31,file=input)
      do i=1,100
       read(31,'(t1,12a1)') (temp(j),j=1,12)
       if(temp(1).eq.'#') then
        if(temp(2).eq.'K') then
         write(title,'(7a1)') (temp(j),j=4,10)
         if(title.eq.'ANNULUS') then
          backspace(31)
          read(31,'(t17,f7.3)') ann
          elseif(title.eq.'DANNULU') then
          backspace(31)
          read(31,'(t17,f7.3)') dann
         endif
        endif
       else
        backspace(31)
        goto 100
       endif
      enddo
100   do i=1,ns
       read(31,*,end=500) skip,axc(i),ayc(i)
      enddo
500   close(31)
      nap=i-1
      return
      end







      subroutine reading(xc,yc,mag,ID,dat,nstar)
      parameter(ns=100000)
      real*8 xc(ns),yc(ns),mag(ns),merr(ns),a1
      real*8 sharp(ns),chi(ns),sec
      character filename*30,filter*16,input*30,H*1
      character rhead*10,INimage*16,Pimage*16,cut*16
      integer ID(ns),nstar,ninput,i1,hour,mini
      real*8 dat(8),ut
      print *,'Input als file name=?'
c      read(*,'(a30)') input
      write(filename,'(a30)') 'als.1                                   '
      open(21,file=filename,ERR=22)
      goto 30
22    print *,'Failed to open ',filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
24    format(T17,A16)                               !character format 
25    format (T17,F7.1)                         ! data min,max format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
30    print *,'Running subroutine reading...'
      do i=1,100
       read(21,'(t1,a1)',err=31) H
c       print *,H
       if(H.eq.'#') then
        backspace(21)
        read(21,'(t4,a10)',err=31) rhead
c        print *,rhead
c        if(rhead.eq.'IMAGE     ') then
c         backspace(21)
c         read(21,24,err=31) INimage
c         elseif(rhead.eq.'DATAMIN   ') then
c         backspace(21)
c         read(21,25,err=31) dat(1)

c         elseif(rhead.eq.'DATAMAX   ') then
c         backspace(21)
c         read(21,25,err=31) dat(2)

c         elseif(rhead.eq.'OTIME     ') then
c         backspace(21)
c         read(21,24,err=31) cut
c         read(cut,'(I2,a1,I2,a1,f6.3)',err=31) hour,H,mini,H,sec
c         ut=real(hour)+real(mini)/60.d0+real(sec)/3600.d0
c         elseif(rhead.eq.'GAIN      ') then
c         backspace(21)
c         read(21,25,err=31) dat(3)

c         elseif(rhead.eq.'READNOISE ') then
c         backspace(21)
c         read(21,25,err=31) dat(4)

c         elseif(rhead.eq.'XAIRMASS  ') then
c         backspace(21)
c         read(21,25,err=31) dat(5)

c         elseif(rhead.eq.'IFILTER   ') then
c         backspace(21)
c         read(21,24,err=31) filter

c         elseif(rhead.eq.'PSFMAG    ') then
c         backspace(21)
c         read(21,25,err=31) dat(6)

c         elseif(rhead.eq.'PSFRAD    ') then
c         backspace(21)
c         read(21,25,err=31) dat(7)

         if(rhead.eq.'FITRAD    ') then
         backspace(21)
         read(21,25,err=31) dat(8)
        endif

       else
        backspace(21)
c        print *,H
        goto 32
       endif
      enddo
     
31    print *,'Failed to start reading from ',filename
32    do i=1,ns
       read(21,*,end=34,err=35) ID(i),xc(i),yc(i),mag(i),merr(i)
       read(21,*,end=34,err=36) sharp(i),chi(i)
      enddo
      goto 34
35    print *,i,id(i),xc(i),yc(i),mag(i),merr(i)
      goto 34
36    print *,i,id(i),sharp(i),chi(i)   
34    nstar=i-1
      print *,dat(8)
      close(21)
      return
      end

       
