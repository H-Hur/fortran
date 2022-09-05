c     f95 -o dstar_coo /home/gjgusdh/fortran/subroutine/fs/dstar_coo.f
      parameter (ns=50000)
      integer id(ns),nstar,nheader,ndel,inf_del(ns),num_del
      real*8 dxc(ns),dyc(ns),xc(ns),yc(ns),dist,dfit,rh
      character input*30,H*1,ch(30)*1,header(100)*80,l1(ns)*80,l2(ns)*80

      dfit=1.5d0

      print *,'input als file to delete=?'
      read(*,'(a30)') input
      open(21,file=input)
      read(input,'(30a1)') (ch(j),j=1,30)
      do i=1,30
       if(ch(i).eq.'.') goto 10
      enddo
10    write(ch(i+3),'(a1)') 'n'
      write(input,'(30a1)') (ch(j),j=1,i+3)
      open(22,file=input)
      do i=1,100
       read(21,'(t1,a80)') header(i)
       read(header(i),'(t1,a1)') H
       if(H.ne.'#') goto 20
      enddo
20    nheader=i-1
      backspace(21)
      do i=1,nheader
       write(22,'(t1,a80)') header(i)
      enddo

      do i=1,ns*2
       read(21,'(t1,a80)',end=99) l1(i)
       read(l1(i),*) id(i),xc(i),yc(i)
       read(21,'(t1,a80)') l2(i) 
       inf_del(i)=0
      enddo
99    nstar=i-1
      close(21)

      print *,'input coo file to delete=?' 
      read(*,'(a30)') input
      open(23,file=input)
      do i=1,ns
       read(23,*,end=199) dxc(i),dyc(i)
      enddo
199   close(23)
      ndel=i-1
      print *,nstar,' stars read, ',ndel,' stars to delete'
      num_del=0
      do i=1,nstar
       do j=1,ndel
        dist=sqrt((dxc(j)-xc(i))**2+(dyc(j)-yc(i))**2)
        if(dist.le.dfit) inf_del(i)=1
        if(dist.le.dfit) num_del=num_del+1
       enddo
       if(inf_del(i).eq.0) then
        write(22,'(t1,a80)') l1(i)
        write(22,'(t1,a80)') l2(i)
       endif
      enddo
      close(22)
      print *, num_del,' stars deleted'

      stop
      end

      
