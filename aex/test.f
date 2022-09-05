      real  ub_s0,bv_s,bv_int,ub_int,u_final
      character input*16


      write(input,'(a16)') 'Sung2001.ZAMS.cc'
      call read_ZAMS(input)
      ub_s0=-1.081 + 0.1267
      bv_s=-0.271
      call f_iter(ub_s0,bv_s,bv_int,ub_int)
 
      stop
      end

      subroutine f_iter(ub_s0,bv_s,x,y)
      parameter (ns=100000)
      real bv_s,f0,f,c_bv(21),c_u(21)
      real x,x0,y,ub_s0
      integer niter
      common /zams/zMv(34),zVI(34),zBV(34),zUB(34)
      data c_bv/-.35d0,-.33d0,-0.3d0,-0.2d0,-0.1d0, 0.0d0, 0.1d0, 0.2d0,
     &0.3d0, 0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,1.2d0,1.4d0,1.6d0
     &,1.8d0,2.0d0/
      data c_u/-.132d0,-.132d0,-.126d0,-.074d0,-.030d0,0.002d0,0.024d0,
     &0.026d0,0.010d0,-.017d0,-.06d0,-.104d0,-.132d0,-.132d0,-.132d0,
     &-.132d0,-.132d0,-.132d0,-.132d0,-.132d0,-.132d0/
      call lint(bv_s,c_bv,c_u,f,21,21)
      x=bv_s
      do i=1,100000
       f0=f
       x0=x
       ub_s=ub_s0+f
       call cor_red(bv_s,ub_s,x,y)
       call lint(x,c_bv,c_u,f,21,21)
       if(sqrt((f-f0)**2.d0).le.0.00001.and.sqrt((x-x0)*2.d0).le.
     & 0.00001) then
        goto 500
       endif
       print *,i,x,x0,f,f0
      enddo
500   ub_int=y
      bv_int=x
c      print *,x,y,f
 
      return
      end


      subroutine cor_red(sBV,sUB,unredx,unredy)
      parameter (nf=10,ns=1000,maxiter=1000)
      common /zams/zMv(34),zVI(34),zBV(34),zUB(34)
      real fx,fy,dy,unredx,unredy,sBV,sUB,slope,eBV,eUB,Rv
      slope=0.72
      Rv=3.1
      call lint(sUB,zUB,zBV,fx,34,34)
      dx=fx-sBV
      call lint(sBV,zBV,zUB,fy,34,34)
      dy=sUB-fy
      unredx=sBV
      unredy=sUB
      do j=1,1000
       call lint(unredy,zUB,zBV,fx,34,34)
       dx=fx-unredx
       dy=dx*slope
       unredx=unredx+dx
       unredy=unredy+dy
        if(abs(dx).lt.0.00001.and.abs(dy).lt.0.00001) then
        goto 100
       endif
100    eBV=sBV-unredx
      enddo
      
c       write(*,'(a7,f7.4)') 'e(B-V)=  ',eBV
c       write(*,'(a7,f7.4)') 'slope= ',slope

      return
      end


      subroutine read_ZAMS(input)
      parameter (ns=100000)
      character H*1,input*16
      common /zams/zMv(34),zVI(34),zBV(34),zUB(34)
      open(31,file=input)
      read(31,'(a1)') H
      do i=1,34
       read(31,*,end=100) zMv(i),zVI(i),zBV(i),zUB(i)
      enddo
100   nzams=i-1
      close(31)
      return
      end


      subroutine lint(x,xt,yt,y,n,m)
      parameter(nmax=10000)
      real xt(n),yt(n),xt1(nmax),yt1(nmax)
      real x,y
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

