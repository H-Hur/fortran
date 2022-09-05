c     gfortran -o make_psfs ~/work/fortran/subroutine/daotools/make_psfs.f
c     To build a psf list file from coordinates and aperture photometry file
c     31,32
      character input*32
      real aprad

      call psfs_head()


      print *,'Input aperture photometry file=?'
      read(*,'(t1,a32)') input
c      write(input,'(a32)') 'i3s_3.ap                                '
      print *,'FWHM of the image=?'
c      aprad=3.2  
      read(*,*) aprad
      call readapc(input,aprad)

      print *,'Input psf coordinate list=?       '
c      write(input,'(a32)') 'pst.5                                  ' 
      read(*,'(t1,a32)') input 
      call read_coo(input)

      call build_psfs(aprad)


      stop
      end


      subroutine build_psfs(aprad)
      parameter (ns=100000)
      common/coo/xc(ns),yc(ns),ncoo
      common/ap/id(ns),amerr(ns),amag(ns),axc(ns),ayc(ns),sky(ns),nap
      common/head/l1,l20,l21,l24
      character l1(19)*61,l20*1,l21(3)*80,l24*1      
      integer nn
      real dist,bmag,aprad
      open(33,file='psfs.s')
      do i=1,19
       write(33,'(t1,a61)') l1(i)
      enddo
      write(33,'(t1,a1)') l20
      do i=1,3
       write(33,'(t1,a80)') l21(i)
      enddo
      write(33,'(t1,a1)') l24
      do i=1,ncoo
       bmag=100.
       nn=0 
       do j=1,nap
        dist=sqrt((xc(i)-axc(j))**2+(yc(i)-ayc(j))**2)
        if(dist.le.aprad.and.amag(j).le.bmag) then
         bmag=amag(i)
         nn=j
        endif
       enddo
       if(nn.ne.0) then
        write(33,'(t1,i5,t10,f8.3,t20,f8.3,t30,f6.3,t42,f8.3)') 
     &  id(nn),axc(nn),ayc(nn),amag(nn),sky(nn)
       endif
      enddo 
      close(33)
      print *,'PSF star list has written to psfs.s'
      return
      end
      subroutine read_coo(input)
      parameter (ns=100000)
      character input*32
      common/coo/xc(ns),yc(ns),ncoo
      open(32,file=input)
      do i=1,ns
       read(32,*,end=99) xc(i),yc(i)
      enddo
99    ncoo=i-1
      close(32)
      return
      end
      subroutine readapc(input,aprad)
      parameter (ns=100000)
      character input*32,temp(20)*3,title*7,cmer*5
      integer nhead,naperture,numap,nindef,nap
      real ap(16),corap,readap,skip,aprad
      common/ap/id(ns),amerr(ns),amag(ns),axc(ns),ayc(ns),sky(ns),nap
      corap=aprad
      open(31,file=input)
      nhead=0
      do i=1,100
       read(31,'(t1,12a1)') (temp(j),j=1,12)
       if(temp(1).eq.'#') then
        nhead=nhead+1
       else
        backspace(31)
        goto 100
       endif
      enddo
100   read(31,'(t1,a3)') temp(1)
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
       read(31,'(t44,i6)',end=500) id(i)
       read(31,*) axc(i),ayc(i)
       read(31,*) sky(i)
       read(31,'(a7)') title
       do j=1,naperture
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
      return
      end
      subroutine psfs_head()
      common/head/l1,l20,l21,l24
      character l1(19)*61,l20*1,l21(3)*80,l24*1
      write(l1(1),'(a61)')'#K IRAF       = NOAO/IRAFV2.16.1        versi
     &on    %-23s                                                      '
      write(l1(2),'(a61)')'#K USER       = make_psfs.f             name
     &      %-23s            '
      write(l1(3),'(a61)')'#K HOST       = make_psfs.f             compu
     &ter   %-23s     '
      write(l1(4),'(a61)')'#K DATE       = 2014-04-04              yyyy-
     &mm-dd %-23s     '
      write(l1(5),'(a61)')'#K TIME       = 18:53:20                hh:mm
     &:ss   %-23s     '
      write(l1(6),'(a61)')'#K PACKAGE    = daophot                 name 
     &      %-23s                '
      write(l1(7),'(a61)')'#K TASK       = make_psfs.f             name
     &      %-23s                '
      write(l1(8),'(a61)')'#K IMAGE      = make_psfs.f             image
     &name  %-23s     '
      write(l1(9),'(a61)')'#K PHOTFILE   = make_psfs.f             filen
     &ame   %-23s     '
      write(l1(10),'(a61)')'#K PSTFILE    = make_psfs.f             file
     &name   %-23s      '
      write(l1(11),'(a61)')'#K PSFIMAGE   = make_psfs.f             imag
     &ename  %-23s     '
      write(l1(12),'(a61)')'#K GRPSFILE   = make_psfs.f             file
     &name   %-23s     '
      write(l1(13),'(a61)')'#K OPSTFILE   = make_psfs.f             file
     &name   %-23s     '
      write(l1(14),'(a61)')'#K SCALE      = 1.                      unit
     &s/pix  %-23.7g   '
      write(l1(15),'(a61)')'#K OTIME      = 01:01:01.000            time
     &unit   %-23s     '
      write(l1(16),'(a61)')'#K IFILTER    = make_psfs.f             filt
     &er     %-23s     '
      write(l1(17),'(a61)')'#K XAIRMASS   = 0.000000                numb
     &er     %-23.7g   '
      write(l1(18),'(a61)')'#K PSFRAD     = 1.0                     scal
     &eunit  %-23.7g   '
      write(l1(19),'(a61)')'#K FITRAD     = 1.0                     scal
     &eunit  %-23.7g   '
      write(l20,'(a1)')'#'
      write(l21(1),'(a80)')'#N ID    XCENTER   YCENTER   MAG         MSK
     &Y                                   \'
      write(l21(2),'(a80)')'#U ##    pixels    pixels    magnitudes  cou
     &nts                                 \'
      write(l21(3),'(a80)')'#F %-9d  %-10.3f   %-10.3f   %-12.3f     %-1
     &5.7g                                 '
      write(l24,'(a1)')'#'
      return
      end
