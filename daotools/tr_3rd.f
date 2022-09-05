c     gfortran -o tr_3rd ~/work/fortran/subroutine/daotools/tr_3rd.f ~/work/fortran/subroutine/caltools/mat.f
      parameter (ns=100000,nx=1000)
      character image1*32,image2*32,image1_all*32,output*32
      real*8 x1(ns),x2(ns),y1(ns),y2(ns),xq(nx),yq(nx)
      real*8 x(ns),y(ns),xx,yy
      integer nc1,nc2,nstar
      print *,'coo file 1 for cid ? '
      read(*,*) image1
      print *,'coo file 2 for cid ? '
      read(*,*) image2
      print *,'coo file of stars ? '
      read(*,*) image1_all
      print *,'output coo file name? '
      read(*,*) output
   
      call read_coo1(image1,x1,y1,nc1)
      call read_coo1(image2,x2,y2,nc2)
      call read_coo2(image1_all,x,y,nstar)
      if(nc1.ne.nc2) then
       print *,image1,' ',nc1
       print *,image2,' ',nc2
       print *,'# of stars from  ',image1,' and ',image2,' is different
     &!!! check stars!!'
       goto 900
      endif
      print *,1
      call make_3xy(10,x1,y1,x2,nc1,xq,1.0d0)
      print *,2
      call make_3xy(10,x1,y1,y2,nc1,yq,1.0d0)
      print *,3
      open(22,file='tr_3rd.out')
      write(22,'(a84)') '#   xc1       yc1        xc2       yc2       xc
     &al      ycal         xdif     ydif                               '
      do i=1,nc1
       xx=xq(1) +xq(2)*x1(i) +xq(3)*y1(i) +xq(4)*x1(i)*x1(i) 
     & +xq(5)*x1(i)*y1(i) +xq(6)*y1(i)*y1(i) + xq(7)*x1(i)*x1(i)*x1(i) 
     & +xq(8)*x1(i)*x1(i)*y1(i) +xq(9)*x1(i)*y1(i)*y1(i) 
     & +xq(10)*y1(i)*y1(i)*y1(i)
       yy=yq(1) +yq(2)*x1(i) +yq(3)*y1(i) +yq(4)*x1(i)*x1(i) 
     & +yq(5)*x1(i)*y1(i) +yq(6)*y1(i)*y1(i) + yq(7)*x1(i)*x1(i)*x1(i) 
     & +yq(8)*x1(i)*x1(i)*y1(i) +yq(9)*x1(i)*y1(i)*y1(i) 
     & +yq(10)*y1(i)*y1(i)*y1(i)
       write(22,'(2f10.3,1x,2f10.3,1x,2f10.3,1x,2f10.3)') 
     & x1(i),y1(i),x2(i),y2(i),xx,yy,xx-x2(i),yy-y2(i) 
      enddo
      close(22)
      open(23,file=output)
      do i=1,nstar
       xx=xq(1) +xq(2)*x(i) +xq(3)*y(i) +xq(4)*x(i)*x(i) 
     & +xq(5)*x(i)*y(i) +xq(6)*y(i)*y(i) + xq(7)*x(i)*x(i)*x(i)
     & +xq(8)*x(i)*x(i)*y(i) +xq(9)*x(i)*y(i)*y(i)
     & +xq(10)*y(i)*y(i)*y(i)
       yy=yq(1) +yq(2)*x(i) +yq(3)*y(i) +yq(4)*x(i)*x(i)
     & +yq(5)*x(i)*y(i) +yq(6)*y(i)*y(i) + yq(7)*x(i)*x(i)*x(i)
     & +yq(8)*x(i)*x(i)*y(i) +yq(9)*x(i)*y(i)*y(i)
     & +yq(10)*y(i)*y(i)*y(i)
       if(xx.gt.2.5.and.xx.lt.2045.5.and.yy.gt.2.5.and.yy.lt.4093.5)
     & write(23,'(2f10.3)') xx,yy
      enddo
      close(23)

900   stop
      end

      subroutine read_coo1(filename,xc,yc,nc)
      parameter (ns=1000)
      integer nc
      real*8 xc(nc),yc(nc)
      character filename*32,H*1
      open(21,file=filename)
      do i=1,ns
       read(21,'(t1,a1)') H
       if(H.ne.'#') then
        backspace(21)
        goto 99
       endif
      enddo
99    do i=1,ns
       read(21,*,end=199) xc(i),yc(i)
      enddo
199   nc=i-1
      write(*,'(t1,i5,a16,a32)') nc,' stars read from ',filename
      close(21)
      return
      end
      subroutine read_coo2(filename,xc,yc,nc)
      parameter (ns=100000)
      integer nc
      real*8 xc(nc),yc(nc)
      character filename*32,H*1
      open(21,file=filename)
      do i=1,ns
       read(21,'(t1,a1)') H
       if(H.ne.'#') then
        backspace(21)
        goto 99
       endif
      enddo
99    do i=1,ns
       read(21,*,end=199) xc(i),yc(i)
      enddo
199   nc=i-1
      write(*,'(t1,i5,a16,a32)') nc,' stars read from ',filename
      close(21)
      return
      end
