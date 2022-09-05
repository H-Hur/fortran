!     For fortran95.
!     To extract coordinates from als file.
!     18,Jan,2008, V1.00, apco.f, apco.exe, by Hyeon-oh Hur.
      parameter (ns=100000)
      integer ID(ns),nap,nstar,napstar
      real*8 xc(ns),yc(ns),fwhm,aprad,d,mag(ns),dat(8)
      real*8 axc(ns),ayc(ns),dx,dy
      character input*16       
      print *,'Input apc file=?'
      read(*,'(a16)') input    
      open(101,file='apc.coo')
      dx=0.d0
      dy=0.d0
      print *,'dx,dy=?'
      read(*,*) dx,dy
     
       
      call reading(xc,yc,mag,ID,dat,nstar) 
      call readapc(input,aprad,axc,ayc,nap)

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



      subroutine readapc(input,aprad,axc,ayc,nap)
      parameter (ns=100000,ni=30)
      character input*16,temp(20)*3,title*7,cmer*5
      integer nimage,nhead,naperture,numap,nbad,nindef,nap
      real*8 fitrad,ann,dann,ap(16),corap,readap,skip,aprad,axc(ns)
      real*8 ayc(ns)
      common/ap/amerr(ns),amag(ns),pmag(ns)
      corap=aprad
      open(31,file=input)
      nhead=0
      do i=1,100
       read(31,'(t1,12a1)') (temp(j),j=1,12)
       if(temp(1).eq.'#') then
        nhead=nhead+1
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
100   fitrad=ann+dann

      read(31,'(t1,a3)') temp(1)
      do i=2,20
       read(31,'(t1,a3)') temp(i)
       if(temp(i).eq.temp(1)) then
        goto 200
       endif
      enddo
200   naperture=i-5

      do i=1,naperture+1
       backspace(31)
      enddo
      do i=1,naperture
       read(31,*) ap(i)
       if(ap(i).eq.corap) then
        numap=i
        goto 250
       endif
      enddo
250   rewind(31)
      nbad=0
      nindef=0
      do i=1,nhead
       read(31,'(a7)') title
      enddo

      do 300 i=1,ns
       read(31,'(a7)',end=500) title
       read(31,*) axc(i),ayc(i)
       read(31,'(a7)') title
       read(31,'(a7)') title
       read(31,*,end=500) readap,skip,skip,skip,pmag(i)
       do j=2,naperture
        if(j.eq.numap) then
         read(31,'(t52,a5)',end=500) cmer
         if (cmer.eq.'INDEF') then
          nindef=nindef+1
          goto 300
          else
          backspace(31)
         endif
         read(31,*,end=500) readap,skip,skip,skip,amag(i),amerr(i)
         if(amerr(i).gt.0.01d0) then
          nbad=nbad+1
         endif
         else
         read(31,'(a7)',end=500) title
        endif
       enddo
300   continue
500   close(31)
      nap=i-1
      if(nbad.gt.0) then
       print *,'Warring. ',input,'has ',nbad,' of high merr stars higher
     * than 0.01.'
      endif
      return
      end








      subroutine reading(xc,yc,mag,ID,dat,nstar)
      parameter(ns=100000)
      real*8 xc(ns),yc(ns),mag(ns),merr(ns),a1
      real*8 sharp(ns),chi(ns),sec
      character filename*16,filter*16,input*16,H*1
      character rhead*10,INimage*16,Pimage*16,cut*16
      integer ID(ns),nstar,ninput,i1,hour,mini
      real*8 dat(8),ut
      print *,'Input als file name=?'
      read(*,'(a16)') input
      write(filename,'(a16)') input
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
       if(H.eq.'#') then
        backspace(21)
        read(21,'(t4,a10)',err=31) rhead
        if(rhead.eq.'IMAGE     ') then
         backspace(21)
         read(21,24,err=31) INimage
         elseif(rhead.eq.'DATAMIN   ') then
         backspace(21)
         read(21,25,err=31) dat(1)

         elseif(rhead.eq.'DATAMAX   ') then
         backspace(21)
         read(21,25,err=31) dat(2)

         elseif(rhead.eq.'OTIME     ') then
         backspace(21)
         read(21,24,err=31) cut
         read(cut,'(I2,a1,I2,a1,f6.3)',err=31) hour,H,mini,H,sec
         ut=real(hour)+real(mini)/60.d0+real(sec)/3600.d0
         elseif(rhead.eq.'GAIN      ') then
         backspace(21)
         read(21,25,err=31) dat(3)

         elseif(rhead.eq.'READNOISE ') then
         backspace(21)
         read(21,25,err=31) dat(4)

         elseif(rhead.eq.'XAIRMASS  ') then
         backspace(21)
         read(21,25,err=31) dat(5)

         elseif(rhead.eq.'IFILTER   ') then
         backspace(21)
         read(21,24,err=31) filter

         elseif(rhead.eq.'PSFMAG    ') then
         backspace(21)
         read(21,25,err=31) dat(6)

         elseif(rhead.eq.'PSFRAD    ') then
         backspace(21)
         read(21,25,err=31) dat(7)

         elseif(rhead.eq.'FITRAD    ') then
         backspace(21)
         read(21,25,err=31) dat(8)
        endif

        else
        backspace(21)
        goto 32
       endif
      enddo

     
31    print *,'Failed to start reading from ',filename
32    do i=1,ns
       read(21,*,end=34) ID(i),xc(i),yc(i),mag(i),merr(i)
       read(21,*,end=34) sharp(i),chi(i)
      enddo
34    nstar=i-1
      close(21)
      return
      end

       
