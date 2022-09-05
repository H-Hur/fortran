c     gfortran -o del_als ~/work/fortran/subroutine/daotools/del_als.f
      integer i
      call read_als()
      call read_coo() 
      call write_output()


      stop
      end

      subroutine write_output()
      parameter (ns=100000)
      common/coo/xdel(ns),ydel(ns),ncoo
      common/als/id(ns),xc(ns),yc(ns),line1,line2,ndel(ns),header,
     &input,nstar,nheader
      character line1(ns)*80,line2(ns)*80,output*16,input*16,ch(16)*1
      character header(ns)*80,H*1
      real dist,mag,merr,msky,sharp,chi,r
      integer n_del,num_del(ns),niter,n
      read(input,'(16a1)') (ch(i), i=1,16)
      do i=1,16
       if(ch(i).eq.'.') goto 10
      enddo
10    write(ch(i+1),'(a1)') 'd'
      write(ch(i+2),'(a1)') ' '
      write(ch(i+3),'(a1)') ' '
      write(output,'(a16)') ' '
      write(output,'(16a1)') (ch(j),j=1,i),ch(i+1),ch(i+2),ch(i+3)
      open(21,file=output)
      open(31,file='del_als.rjt')
      do i=1,ns
       read(31,'(t1,a1)',end=99) H
      enddo
99    backspace(31)
      do i=1,nheader
       write(21,'(t1,a80)') header(i)
      enddo
      do i=1,ncoo
       num_del(i)=0
      enddo
      n_del=0
      do i=1,nstar
       do j=1,ncoo
        dist=sqrt((xc(i)-xdel(j))**2+(yc(i)-ydel(j))**2)
        if(dist.lt.1.5) then 
         n_del=n_del+1
         num_del(j)=1
         print *,xc(i),yc(i)
         read(line1(i),*) n,r,r,mag,merr,msky,niter
         read(line2(i),*) sharp,chi
         write(31,'(t1,i6,2f9.3,2f7.3,f10.3,2f9.3)') 
     &   id(i),xc(i),yc(i),mag,merr,msky,sharp,chi
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
c      do i=1,ncoo
c       if(num_del(i).eq.0) print *,i,'st star are not deleted!'
c      enddo
      return
      end

      subroutine read_coo()
      parameter (ns=100000)
      common/coo/xdel(ns),ydel(ns),ncoo
      character input*50
      print *,'input coo file to delete?'
      read(*,'(t1,a50)') input
      open(22,file=input)
      do i=1,ns
       read(22,*,end=99) xdel(i),ydel(i)
      enddo
99    ncoo=i-1
      close(22)
      return
      end

      subroutine read_als()
      parameter (ns=100000)
      common/als/id(ns),xc(ns),yc(ns),line1,line2,ndel(ns),header,
     &input,nstar,nheader
      character line1(ns)*80,line2(ns)*80,input*16,H*1,header(ns)*80
      print *,'input als file?'
      read(*,'(t1,a16)') input
      open(23,file=input)
      do i=1,ns
       read(23,'(t1,a1,t1,a80)') H,header(i)
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
