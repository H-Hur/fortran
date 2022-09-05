      real*8 c_bv(21),c_u(21),bv0,ub0,V(34),VI(34),BV(34),UB(34)
      real*8 cn_bv(236),cn_u(236),f
      integer nzams
      character H*1
      data c_bv/-.35d0,-.33d0,-0.3d0,-0.2d0,-0.1d0, 0.0d0, 0.1d0, 0.2d0,
     &0.3d0, 0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,1.2d0,1.4d0,1.6d0
     &,1.8d0,2.0d0/
      data c_u/-.132d0,-.132d0,-.126d0,-.074d0,-.030d0,0.002d0,0.024d0,
     &0.026d0,0.010d0,-.017d0,-.06d0,-.104d0,-.132d0,-.132d0,-.132d0,
     &-.132d0,-.132d0,-.132d0,-.132d0,-.132d0,-.132d0/

      open(21,file='Sung2001.ZAMS.cc')
      open(22,file='U-f.B-V.ZAMS.cc')
      read(21,'(a1)') H
      do i=1,100
       read(21,*,end=100) V(i),VI(i),BV(i),UB(i)
      enddo
100   close(21)
      nzams=i-1
      print *,nzams
      write(22,'(t1,a16)') '#   B-V    U-B-f'
      do i=1,236
       cn_bv(i)=real(i)/100.d0-0.36d0
       call lint(cn_bv(i),c_bv,c_u,f,21,21)
       call lint(cn_bv(i),BV,UB,ub0,nzams,nzams)
       cn_u(i)=ub0-f
       write(22,'(t1,2f8.3)') cn_bv(i),cn_u(i)
      enddo
      close(22)
 
      stop
      end


      subroutine lint(x,xt,yt,y,n,m)
      parameter(nmax=10000)
      real*8 xt(n),yt(n),xt1(nmax),yt1(nmax)
      real*8 x,y
      integer i,k,m
      if(xt(m).lt.xt(1)) then
       do 100 i=1,m
        xt1(i)=xt(m+1-i)
        yt1(i)=yt(m+1-i)
100    continue
      else
       do 200 i=1,m
        xt1(i)=xt(i)
        yt1(i)=yt(i)
200    continue
      endif
      do 10 i=2,m
       k=i-1
       if(x-xt1(i).lt.0.d0) then
        goto 11
        elseif(x-xt1(i).eq.0.d0) then
        goto 12
        elseif (x-xt1(i).gt.0.d0) then
        goto 10
       endif
12     y=yt1(i)
       goto 14
10    continue
11    y=yt1(k)+(yt1(k+1)-yt1(k))*((x-xt1(k))/(xt1(k+1)-xt1(k)))
14    return
      end

