c     gfortran -o del_als_badpix ~/work/fortran/subroutine/daotools/del_als_badpix.f
      character input*32, als*32, badpix*32,H*1
      open(32,file='del_als_badpix.input')
      read(32,'(t1,a1)') H
      do i=1,1000
       read(32,'(t1,2a32)',end=99) als,badpix
       call read_als(als)
       call read_badpix(badpix)
       call write_output()
      enddo
99    close(32)
      stop
      end

      subroutine write_output()
      parameter (ns=100000)
      common/badpix/x1(ns),x2(ns),y1(ns),y2(ns),nbad
      common/als/id(ns),xc(ns),yc(ns),line1,line2,ndel(ns),header,
     &input,nstar,nheader,fwhm
      character line1(ns)*80,line2(ns)*80,output*32,input*32,ch(32)*1
      character header(ns)*80,H*1
      real dist,mag,merr,msky,sharp,chi,r
      integer niter,n,n_del
      read(input,'(32a1)') (ch(i), i=1,32)
      do i=1,32
       if(ch(i).eq.'.') goto 10
      enddo
10    write(ch(i+1),'(a1)') 'd'
      write(ch(i+2),'(a1)') ' '
      write(ch(i+3),'(a1)') ' '
      write(output,'(a32)') ' '
      write(output,'(32a1)') (ch(j),j=1,i),ch(i+1),ch(i+2),ch(i+3)
      open(21,file=output)
      open(31,file='badpix.rjt')
      do i=1,ns
       read(31,'(t1,a1)',end=99) H
      enddo
99    backspace(31)
      do i=1,nheader
       write(21,'(t1,a80)') header(i)
      enddo
      n_del=0
      do i=1,nstar
       do j=1,nbad
        if(xc(i).gt.x1(j)-0.5*fwhm.and.xc(i).lt.x2(j)+0.5*fwhm.and.
     &  yc(i).gt.y1(j)-0.5*fwhm.and.yc(i).lt.y2(j)+0.5*fwhm) then
         write(*,'(2f10.3)') xc(i),yc(i)
         read(line1(i),*) n,r,r,mag,merr,msky,niter
         read(line2(i),*) sharp,chi
         write(31,'(t1,i6,2f9.3,2f7.3,f10.3,2f9.3)')
     &   id(i),xc(i),yc(i),mag,merr,msky,sharp,chi
         n_del=n_del+1 
         goto 100
        endif
       enddo
       write(21,'(t1,a80)') line1(i)
       write(21,'(t1,a80)') line2(i)
100   i100=i100
      enddo
      close(21)
      write(31,'(t1,a1)') ' '
      close(31)
      print *,n_del,' stars were deleted from ',input
      print *,'the result was wrote on ',output
      return
      end
      subroutine read_badpix(input)
      parameter (ns=100000)
      common/badpix/x1(ns),x2(ns),y1(ns),y2(ns),nbad
      character input*32,H*1
      open(22,file=input)
      read(22,'(t1,a1)') H
      do i=1,ns
       read(22,*,end=99) x1(i),x2(i),y1(i),y2(i)
      enddo
99    nbad=i-1
      close(22)
      return
      end

      subroutine read_als(als_in)
      parameter (ns=100000)
      common/als/id(ns),xc(ns),yc(ns),line1,line2,ndel(ns),header,
     &input,nstar,nheader,fwhm
      character line1(ns)*80,line2(ns)*80,input*32,H*1,header(ns)*80
      character als_in*32
      write(input,'(a32)') als_in
      open(23,file=input)
      do i=1,ns
       read(23,'(t1,a1,t1,a80)') H,header(i)
       if(i.eq.26) then
        read(header(i),'(t17,f3.1)') fwhm
        write(*,'(a6,f3.1)') 'FWHM= ',fwhm
       endif
       if(H.ne.'#') then
        nheader=i-1
        backspace(23)
        goto 99
       endif
      enddo
99    do i=1,ns
       read(23,'(t1,a80)',end=199) line1(i)
       read(23,'(t1,a80)',end=199) line2(i)
       read(line1(i),*) id(i),xc(i),yc(i)
      enddo
199   nstar=i-1
      close(23)
      return
      end
