      subroutine sort2(ns,E1,i2,E2,E3)
      real*8 C1,C2,C3,E1(ns),E2(ns),E3(ns)
      integer i,j,i2(ns),ic2
      do i=1,ns
       do j=1,ns-1
        if (E1(j+1).lt.E1(j)) then
         C1=E1(j)
         C2=E2(j)
         C3=E3(j)
         ic2=i2(j)
         E1(j)=E1(j+1)
         E2(j)=E2(j+1)
         E3(j)=E3(j+1)
         i2(j)=i2(j+1)
         E1(j+1)=C1
         E2(j+1)=C2
         E3(j+1)=C3
         i2(j+1)=ic2
        endif
       enddo
      enddo
      return
      end
c     역행렬 구하기                                        
      subroutine make_inverse(n,ref,inv)
      parameter (nx=1000)
      integer n
      real*8 mat(nx,nx),inv(nx,nx),pre(nx),top(nx),ref(nx,nx),col(nx)
      do i=1,n
       do j=1,n
        inv(i,j)=0.d0
        mat(i,j)=ref(i,j)
       enddo
       inv(i,i)=1.d0
      enddo
      do k=1,n
       do i=k,n
        pre(i)=mat(i,k)
        do j=1,n
         mat(i,j)=mat(i,j)/pre(i)
         if(i.eq.k) then
          top(j)=mat(k,j)
         endif
         inv(i,j)=inv(i,j)/pre(i)
        enddo
        if(i.ne.k) then
         do j=1,n
          mat(i,j)=mat(i,j)-top(j)
          inv(i,j)=inv(i,j)-inv(k,j)
         enddo
        endif
       enddo
      enddo
      do i=1,n
       do j=1,n
        if(j.lt.i.and.mat(i,j).lt.1.E-30) then
         mat(i,j) =0.d0
        endif
       enddo
      enddo
      do i=n,1,-1
       do j=n,i+1,-1
        if(j.gt.i) then
         do j2=1,n
          top(j2)=inv(j,j2)
          col(j2)=mat(i,j2)
         enddo
         do j2=1,n
          inv(i,j2)=inv(i,j2)-top(j2)*mat(i,j)
         enddo
         mat(i,j)=0.d0
        endif
       enddo
      enddo
      return
      end
c     행렬식 결정
      subroutine cal_det(n,mat,det)
      parameter (nx=1000)
      integer n,n1,n2
      real*8 mat(nx,nx),cd1(nx),cd2(nx),det,d1,d2
      d1=0.d0
      d2=0.d0
      if(n.eq.2) then
       det=mat(1,1)*mat(2,2)-mat(1,2)*mat(2,1) 
      else
       do i=1,n
        cd1(i)=1.d0
        cd2(i)=1.d0
        do j=1,n
         n1=j
         n2=i+j
         if(n2.ge.n+1) then
          n2=n2-n
         endif
         cd1(i)=cd1(i)*mat(n1,n2)
         n2=i-j
         if(n2.le.0) then
          n2=n2+n
         endif
         cd2(i)=cd2(i)*mat(n1,n2)
        enddo
        d1=d1+cd1(i)
        d2=d2+cd2(i)
       enddo
       det=d1-d2
      endif
      return
      end
C     행렬을 복사 
      subroutine matcopy(n,mat1,mat2)
      parameter (nx=1000)
      integer n
      real*8 mat1(nx,nx),mat2(nx,nx)
      do i=1,n
       do j=1,n
        mat2(i,j)=mat1(i,j)
       enddo
      enddo
      return
      end

      subroutine matcopy2(n,m1,m2)
      integer n
      real*8 m1(n),m2(n)
      do i=1,n
       m2(i)=m1(i)
      enddo
      return
      end
c 행렬과 행렬식을 print.
      subroutine printing(n,mat,s)
      parameter (nx=1000)
      integer n
      real*8 mat(nx,nx),s(n),meq
      call cal_det(n,mat,meq)
      do i=1,n
       write(*,*)(mat(i,j),j=1,n),s(i)
      enddo
      print *,'meq=',meq
      print *,''
      return
      end
c 크래머 적용(사용 후 반드시 recall문으로 원래 행렬로 돌릴 것.)
      subroutine cramer(n,ncol,mat,rgt)
      parameter (nx=1000)
      integer n,ncol
      real*8 mat(nx,nx),rgt(nx),slot(nx)
      do i=1,n
       slot(i)=mat(i,ncol)
       mat(i,ncol)=rgt(i)
       rgt(i)=slot(i)
      enddo
      return
      end

      subroutine make_zero(a,b,c,d,e,f,g,h)
      real*8 a,b,c,d,e,f,g,h
      a=0.d0
      b=0.d0
      c=0.d0
      d=0.d0
      e=0.d0
      f=0.d0
      g=0.d0
      h=0.d0
      return
      end
  
      subroutine cal_eq(n,x,y,eq,aaa,sc)
      parameter (nx=1000)
      real*8 x,y,aaa,eq(nx),sc
      integer n
      aaa=0.d0
!!!!! aaa = eq(1) + eq(2)*x + eq(3)*y + eq(4)*xx + eq(5)*xy + eq(6)*yy + !eq(7)*xxx + eq(8)*xxy + eq(9)*xyy + eq(10)*yyy
      if(n.ge.11) n=10
      do i=1,n
       if(i.eq.1) aaa=aaa+eq(1)
       if(i.eq.2) aaa=aaa+eq(2)*x/sc
       if(i.eq.3) aaa=aaa+eq(3)*y/sc
       if(i.eq.4) aaa=aaa+eq(4)*x*x/sc/sc
       if(i.eq.5) aaa=aaa+eq(5)*x*y/sc/sc
       if(i.eq.6) aaa=aaa+eq(6)*y*y/sc/sc
       if(i.eq.7) aaa=aaa+eq(7)*x*x*x/sc/sc/sc
       if(i.eq.8) aaa=aaa+eq(8)*x*x*y/sc/sc/sc
       if(i.eq.9) aaa=aaa+eq(9)*x*y*y/sc/sc/sc
       if(i.eq.10) aaa=aaa+eq(10)*y*y*y/sc/sc/sc
      enddo
      return
      end

      subroutine cal_eq_s(n,x,y,eq,aaa,sc)
      parameter (nx=1000)
      real*8 x,y,aaaa,eq(nx),sc
      real aaa
      integer n
      aaaa=0.d0
!!!!! aaa = eq(1) + eq(2)*x + eq(3)*y + eq(4)*xx + eq(5)*xy + eq(6)*yy + !eq(7)*xxx + eq(8)*xxy + eq(9)*xyy + eq(10)*yyy
      if(n.ge.11) n=10
      do i=1,n
       if(i.eq.1) aaaa=aaaa+eq(1)
       if(i.eq.2) aaaa=aaaa+eq(2)*x/sc
       if(i.eq.3) aaaa=aaaa+eq(3)*y/sc
       if(i.eq.4) aaaa=aaaa+eq(4)*x*x/sc/sc
       if(i.eq.5) aaaa=aaaa+eq(5)*x*y/sc/sc
       if(i.eq.6) aaaa=aaaa+eq(6)*y*y/sc/sc
       if(i.eq.7) aaaa=aaaa+eq(7)*x*x*x/sc/sc/sc
       if(i.eq.8) aaaa=aaaa+eq(8)*x*x*y/sc/sc/sc
       if(i.eq.9) aaaa=aaaa+eq(9)*x*y*y/sc/sc/sc
       if(i.eq.10) aaaa=aaaa+eq(10)*y*y*y/sc/sc/sc
      enddo
      aaa=real(aaaa)
      return
      end

      subroutine make_3xy(n,xc,yc,aaa,nstar,eq,sc)
      parameter (nx=1000)
      integer nstar
      real*8 x,y,xx,xy,yy,xxx,xxy,xyy,yyy,xxxx,xxxy,xxyy,xyyy,yyyy
      real*8 x5,x4y1,x3y2,x2y3,x1y4,y5,x6,x5y1,x4y2,x3y3,x2y4,x1y5,y6
      real*8 ax,ay,axx,axy,ayy,a,sc,ax3,ax2y,axy2,ay3,aaa(nstar)
      real*8 c(nx,nx),t(nx),eq(nx),inv(nx,nx),xc(nstar),yc(nstar)
      call make_zero(ax,ay,axx,axy,ayy,a,x,y)
      call make_zero(xx,xy,yy,xxx,xxy,xyy,yyy,xxxx)
      call make_zero(xxxy,xxyy,xyyy,yyyy,x5,x4y1,x3y2,x2y3)
      call make_zero(ax3,ax2y,axy2,ay3,x1y4,y5,x6,x5y1)
      call make_zero(x4y2,x3y3,x2y4,x1y5,y6,ax,ay,axx)
!!!!! aaa = eq(1) + eq(2)*x + eq(3)*y + eq(4)*xx + eq(5)*xy + eq(6)*yy + eq(7)*xxx + eq(8)*xxy + eq(9)*xyy + eq(10)*yyy
      do i=1,nstar
       x    = x    + xc(i)/ sc
       y    = y    + yc(i)/ sc
       xx   = xx   + (xc(i)/sc)**2
       yy   = yy   + (yc(i)/sc)**2
       xy   = xy   + xc(i)*yc(i)/sc/sc
       xxx  = xxx  + (xc(i)/sc)**3
       xxy  = xxy  + (yc(i)/sc)*(xc(i)/sc)**2
       xyy  = xyy  + (xc(i)/sc)*(yc(i)/sc)**2
       yyy  = yyy  + (yc(i)/sc)**3
       xxxx = xxxx + (xc(i)/sc)**4
       xxxy = xxxy + (yc(i)/sc)*(xc(i)/sc)**3
       xxyy = xxyy + (xc(i)*yc(i)/sc/sc)**2
       xyyy = xyyy + (xc(i)/sc)*(yc(i)/sc)**3
       yyyy = yyyy + (yc(i)/sc)**4
       x5   = x5   + (xc(i)/sc)**5
       x4y1 = x4y1 + ((xc(i)/sc)**4) * (yc(i)/ sc)
       x3y2 = x3y2 + ((xc(i)/sc)**3) * ((yc(i)/sc)**2)
       x2y3 = x2y3 + ((xc(i)/sc)**2) * ((yc(i)/sc)**3)
       x1y4 = x1y4 + ( xc(i)/ sc   ) * ((yc(i)/sc)**4)
       y5   = y5   + (yc(i)/sc)**5
       x6   = x6   + (xc(i)/sc)**6
       x5y1 = x5y1 + ((xc(i)/sc)**5) * (yc(i)/ sc)
       x4y2 = x4y2 + ((xc(i)/sc)**4) * ((yc(i)/sc)**2)
       x3y3 = x3y3 + ((xc(i)/sc)**3) * ((yc(i)/sc)**3)
       x2y4 = x2y4 + ((xc(i)/sc)**2) * ((yc(i)/sc)**4)
       x1y5 = x1y5 + ((xc(i)/sc)   ) * ((yc(i)/sc)**5)
       y6   = y6   + (yc(i)/sc)**6

       a    = a    + aaa(i)
       ax   = ax   + aaa(i)*xc(i)/sc
       ay   = ay   + aaa(i)*yc(i)/sc
       axx  = axx  + aaa(i)*(xc(i)/sc)**2
       ayy  = ayy  + aaa(i)*(yc(i)/sc)**2
       axy  = axy  + aaa(i)*yc(i)*xc(i)/sc/sc
       ax3  = ax3  + aaa(i)*(xc(i)/sc)**3
       ax2y = ax2y + aaa(i)*((xc(i)/sc)**2)*(yc(i)/sc)
       axy2 = axy2 + aaa(i)*(xc(i)/sc)*((yc(i)/sc)**2)
       ay3  = ay3  + aaa(i)*(yc(i)/sc)**3
      enddo
      c(1,1) = real(nstar)
      c(1,2) = x
      c(1,3) = y
      c(1,4) = xx
      c(1,5) = xy
      c(1,6) = yy
      c(1,7) = xxx
      c(1,8) = xxy
      c(1,9) = xyy
      c(1,10)= yyy

      c(2,1) = x
      c(2,2) = xx
      c(2,3) = xy
      c(2,4) = xxx
      c(2,5) = xxy
      c(2,6) = xyy
      c(2,7) = xxxx
      c(2,8) = xxxy
      c(2,9) = xxyy
      c(2,10)= xyyy

      c(3,1) = y
      c(3,2) = xy
      c(3,3) = yy
      c(3,4) = xxy
      c(3,5) = xyy
      c(3,6) = yyy
      c(3,7) = xxxy
      c(3,8) = xxyy
      c(3,9) = xyyy
      c(3,10)= yyyy

      c(4,1) = xx
      c(4,2) = xxx
      c(4,3) = xxy
      c(4,4) = xxxx
      c(4,5) = xxxy
      c(4,6) = xxyy
      c(4,7) = x5
      c(4,8) = x4y1
      c(4,9) = x3y2
      c(4,10)= x2y3

      c(5,1) = xy
      c(5,2) = xxy
      c(5,3) = xyy
      c(5,4) = xxxy
      c(5,5) = xxyy
      c(5,6) = xyyy
      c(5,7) = x4y1
      c(5,8) = x3y2
      c(5,9) = x2y3
      c(5,10)= x1y4

      c(6,1) = yy
      c(6,2) = xyy
      c(6,3) = yyy
      c(6,4) = xxyy
      c(6,5) = xyyy
      c(6,6) = yyyy
      c(6,7) = x3y2
      c(6,8) = x2y3
      c(6,9) = x1y4
      c(6,10)= y5

      c(7,1) = xxx
      c(7,2) = xxxx
      c(7,3) = xxxy
      c(7,4) = x5
      c(7,5) = x4y1
      c(7,6) = x3y2
      c(7,7) = x6
      c(7,8) = x5y1
      c(7,9) = x4y2
      c(7,10)= x3y3

      c(8,1) = xxy
      c(8,2) = xxxy
      c(8,3) = xxyy
      c(8,4) = x4y1
      c(8,5) = x3y2
      c(8,6) = x2y3
      c(8,7) = x5y1
      c(8,8) = x4y2
      c(8,9) = x3y3
      c(8,10)= x2y4

      c(9,1) = xyy
      c(9,2) = xxyy
      c(9,3) = xyyy
      c(9,4) = x3y2
      c(9,5) = x2y3
      c(9,6) = x1y4
      c(9,7) = x4y2
      c(9,8) = x3y3
      c(9,9) = x2y4
      c(9,10)= x1y5

      c(10,1) = yyy
      c(10,2) = xyyy
      c(10,3) = yyyy
      c(10,4) = x2y3
      c(10,5) = x1y4
      c(10,6) = y5
      c(10,7) = x3y3
      c(10,8) = x2y4
      c(10,9) = x1y5
      c(10,10)= y6
      t(1)   = a
      t(2)   = ax
      t(3)   = ay
      t(4)   = axx
      t(5)   = axy
      t(6)   = ayy
      t(7)   = ax3
      t(8)   = ax2y
      t(9)   = axy2
      t(10)  = ay3

      call make_inverse(n,c,inv)
      do i=1,n
       eq(i)=0.d0
       do j=1,n
        eq(i)=eq(i)+t(j)*inv(i,j)
       enddo
      enddo
c      write(*,*)(eq(i),i=1,n)
c      print *,' '
      return
      end

      subroutine make_2xy(n,xc,yc,aaa,nstar,eq,sc)
      parameter (nx=1000)
      integer nstar
      real*8 x,y,xx,xy,yy,xxx,xxy,xyy,yyy,xxxx,xxxy,xxyy,xyyy,yyyy
      real*8 x5,x4y1,x3y2,x2y3,x1y4,y5,x6,x5y1,x4y2,x3y3,x2y4,x1y5,y6
      real*8 ax,ay,axx,axy,ayy,a,sc,ax3,ax2y,axy2,ay3,aaa(nstar)
      real*8 c(nx,nx),t(nx),eq(nx),inv(nx,nx),xc(nstar),yc(nstar)
      sc=1.d0
      call make_zero(ax,ay,axx,axy,ayy,a,x,y)
      call make_zero(xx,xy,yy,xxx,xxy,xyy,yyy,xxxx)
      call make_zero(xxxy,xxyy,xyyy,yyyy,x5,x4y1,x3y2,x2y3)
      call make_zero(ax3,ax2y,axy2,ay3,x1y4,y5,x6,x5y1)
      call make_zero(x4y2,x3y3,x2y4,x1y5,y6,ax,ay,axx)
!!!!! aaa = ep(1)*x + eq(2)*y + eq(3)*xy + eq(4)*xx + eq(5)*yy + eq(6)
      do i=1,nstar
       x    = x    + xc(i)/ sc
       y    = y    + yc(i)/ sc
       xx   = xx   + (xc(i)/sc)**2
       yy   = yy   + (yc(i)/sc)**2
       xy   = xy   + xc(i)*yc(i)/sc/sc
       xxx  = xxx  + (xc(i)/sc)**3
       xxy  = xxy  + (yc(i)/sc)*(xc(i)/sc)**2
       xyy  = xyy  + (xc(i)/sc)*(yc(i)/sc)**2
       yyy  = yyy  + (yc(i)/sc)**3
       xxxx = xxxx + (xc(i)/sc)**4
       xxxy = xxxy + (yc(i)/sc)*(xc(i)/sc)**3
       xxyy = xxyy + (xc(i)*yc(i)/sc/sc)**2
       xyyy = xyyy + (xc(i)/sc)*(yc(i)/sc)**3
       yyyy = yyyy + (yc(i)/sc)**4

       a    = a    + aaa(i)
       ax   = ax   + aaa(i)*xc(i)/sc
       ay   = ay   + aaa(i)*yc(i)/sc
       axx  = axx  + aaa(i)*(xc(i)/sc)**2
       ayy  = ayy  + aaa(i)*(yc(i)/sc)**2
       axy  = axy  + aaa(i)*yc(i)*xc(i)/sc/sc
      enddo
      c(1,1) = xx
      c(1,2) = xy
      c(1,3) = xxy
      c(1,4) = xxx
      c(1,5) = xyy
      c(1,6) = x

      c(2,1) = xy
      c(2,2) = yy
      c(2,3) = xyy
      c(2,4) = xxy
      c(2,5) = yyy
      c(2,6) = y
      c(3,1) = xxy
      c(3,2) = xyy
      c(3,3) = xxyy
      c(3,4) = xxxy
      c(3,5) = xyyy
      c(3,6) = xy

      c(4,1) = xxx
      c(4,2) = xxy
      c(4,3) = xxxy
      c(4,4) = xxxx
      c(4,5) = xxyy
      c(4,6) = xx
      c(5,1) = xyy
      c(5,2) = yyy
      c(5,3) = xyyy
      c(5,4) = xxyy
      c(5,5) = yyyy
      c(5,6) = yy

      c(6,1) = x
      c(6,2) = y
      c(6,3) = xy
      c(6,4) = xx
      c(6,5) = yy
      c(6,6) = nstar

      t(1)   = ax
      t(2)   = ay
      t(3)   = axy
      t(4)   = axx
      t(5)   = ayy
      t(6)   = a

      call make_inverse(n,c,inv)
      do i=1,n
       eq(i)=0.d0
       do j=1,n
        eq(i)=eq(i)+t(j)*inv(i,j)
       enddo
      enddo
c      write(*,*)(eq(i),i=1,n)
c      print *,' '
      return
      end
      subroutine cal_eq_select(x,y,eq,n,sc,aaa,
     &n0,                   ! C
     &n11,n12,              ! X Y
     &n21,n22,n23,          ! X2 XY Y2
     &n31,n32,n33,n34,      ! X3 X2Y XY2 Y3
     &n41,n42,n43,n44,n45)  ! X4 X3Y X2Y2 XY3 Y4
      integer nstar,n,n0,n11,n12,n21,n22,n23,n31,n32,n33,n34
      integer n41,n42,n43,n44,n45,nc(15),num_c(15)
      real*8 x,y,eq(15),sc,aaa,t(15)
      nc(1)=n0
      nc(2)=n11
      nc(3)=n12
      nc(4)=n21
      nc(5)=n22
      nc(6)=n23
      nc(7)=n31
      nc(8)=n32
      nc(9)=n33
      nc(10)=n34
      nc(11)=n41
      nc(12)=n42
      nc(13)=n43
      nc(14)=n44
      nc(15)=n45
      t(1)  = 1.d0
      t(2)  = x
      t(3)  = y 
      t(4)  = x * x
      t(5)  = x * y
      t(6)  = y * y
      t(7)  = x * x * x
      t(8)  = x * x * y
      t(9)  = x * y * y
      t(10) = y * y * y
      t(11) = x * x * x * x
      t(12) = x * x * x * y
      t(13) = x * x * y * y
      t(14) = x * y * y * y
      t(15) = y * y * y * y !* y! * y * y* y * y * y
      n=0
      do i=1,15
       if(nc(i).eq.1) then
        n=n+1
        num_c(n)=i
       endif
      enddo
      aaa=0.d0 
      do i=1,n
       aaa = aaa + t(num_c(i))*eq(i)
      enddo
      return
      end
      subroutine make_3xy_select(xc,yc,aaa,nstar,n,eq,sc,
     &n0,                   ! C
     &n11,n12,              ! X Y
     &n21,n22,n23,          ! X2 XY Y2
     &n31,n32,n33,n34,      ! X3 X2Y XY2 Y3
     &n41,n42,n43,n44,n45)  ! X4 X3Y X2Y2 XY3 Y4
      parameter (nx=1000,ns=100000)
      integer nstar,n,n0,n11,n12,n21,n22,n23,n31,n32,n33,n34
      integer n41,n42,n43,n44,n45,nc(15),num_c(15)
      real*8 ax4,ax3y,ax2y2,axy3,ay4,ax,ay,axx,axy,ayy,a,sc,ax3
      real*8 ax2y,axy2,ay3,aaa(nstar),c(nx,nx),t(nx),eq(nx),inv(nx,nx)
      real*8 xc(nstar),yc(nstar),temp(15,15,ns),t1,t2
      nc(1)=n0
      nc(2)=n11
      nc(3)=n12
      nc(4)=n21
      nc(5)=n22
      nc(6)=n23
      nc(7)=n31
      nc(8)=n32
      nc(9)=n33
      nc(10)=n34
      nc(11)=n41
      nc(12)=n42
      nc(13)=n43
      nc(14)=n44
      nc(15)=n45
      call make_zero(ax,ay,axx,axy,ayy,a,ay4,ax4)
      call make_zero(ax4,ax3y,ax2y2,axy3,ax3,ax2y,axy2,ay3)
!!!!! aaa = eq(1) + eq(2)*x + eq(3)*y + eq(4)*xx + eq(5)*xy + eq(6)*yy + eq(7)*xxx + eq(8)*xxy + eq(9)*xyy + eq(10)*yyy + eq(11)*xxxx + eq(12)*xxxy + eq(13)*xxyy + eq(14)*xyyy + eq(15)*yyyy
      do i=1,15
       if(nc(i).eq.1) then
        n=n+1
        num_c(n)=i
       endif
      enddo
      do i=1,nstar
       do j=1,n
        if(num_c(j).eq.1)  t1=1.d0
        if(num_c(j).eq.2)  t1=xc(i)/ sc
        if(num_c(j).eq.3)  t1=yc(i)/ sc
        if(num_c(j).eq.4)  t1=(xc(i)/sc)**2
        if(num_c(j).eq.5)  t1=xc(i)*yc(i)/sc/sc
        if(num_c(j).eq.6)  t1=(yc(i)/sc)**2
        if(num_c(j).eq.7)  t1=(xc(i)/sc)**3
        if(num_c(j).eq.8)  t1=(yc(i)/sc)*(xc(i)/sc)**2
        if(num_c(j).eq.9)  t1=(xc(i)/sc)*(yc(i)/sc)**2
        if(num_c(j).eq.10) t1=(yc(i)/sc)**3
        if(num_c(j).eq.11) t1=(xc(i)/sc)**4
        if(num_c(j).eq.12) t1=(yc(i)/sc)*(xc(i)/sc)**3
        if(num_c(j).eq.13) t1=(xc(i)*yc(i)/sc/sc)**2
        if(num_c(j).eq.14) t1=(xc(i)/sc)*(yc(i)/sc)**3
        if(num_c(j).eq.15) t1=(yc(i)/sc)**4 
        do k=1,n
         if(num_c(k).eq.1)  t2=1.d0
         if(num_c(k).eq.2)  t2=xc(i)/ sc
         if(num_c(k).eq.3)  t2=yc(i)/ sc
         if(num_c(k).eq.4)  t2=(xc(i)/sc)**2
         if(num_c(k).eq.5)  t2=xc(i)*yc(i)/sc/sc
         if(num_c(k).eq.6)  t2=(yc(i)/sc)**2
         if(num_c(k).eq.7)  t2=(xc(i)/sc)**3
         if(num_c(k).eq.8)  t2=(yc(i)/sc)*(xc(i)/sc)**2
         if(num_c(k).eq.9)  t2=(xc(i)/sc)*(yc(i)/sc)**2
         if(num_c(k).eq.10) t2=(yc(i)/sc)**3
         if(num_c(k).eq.11) t2=(xc(i)/sc)**4
         if(num_c(k).eq.12) t2=(yc(i)/sc)*(xc(i)/sc)**3
         if(num_c(k).eq.13) t2=(xc(i)*yc(i)/sc/sc)**2
         if(num_c(k).eq.14) t2=(xc(i)/sc)*(yc(i)/sc)**3
         if(num_c(k).eq.15) t2=(yc(i)/sc)**4 
         temp(j,k,i)=t1*t2
        enddo
       enddo
       a    = a    + aaa(i)
       ax   = ax   + aaa(i)*xc(i)/sc
       ay   = ay   + aaa(i)*yc(i)/sc
       axx  = axx  + aaa(i)*(xc(i)/sc)**2
       ayy  = ayy  + aaa(i)*(yc(i)/sc)**2
       axy  = axy  + aaa(i)*yc(i)*xc(i)/sc/sc
       ax3  = ax3  + aaa(i)*(xc(i)/sc)**3
       ax2y = ax2y + aaa(i)*((xc(i)/sc)**2)*(yc(i)/sc)
       axy2 = axy2 + aaa(i)*(xc(i)/sc)*((yc(i)/sc)**2)
       ay3  = ay3  + aaa(i)*(yc(i)/sc)**3
       ax4  = ax4  + aaa(i)*(xc(i)/sc)**4
       ax3y = ax3y + aaa(i)*((xc(i)/sc)**3)*(yc(i)/sc)
       ax2y2= ax2y2+ aaa(i)*((xc(i)/sc)**2)*((yc(i)/sc)**2)
       axy3 = axy3 + aaa(i)*((xc(i)/sc))*((yc(i)/sc)**3)
       ay4  = ay4  + aaa(i)*(yc(i)/sc)**4 
      enddo
      do i=1,15
       do j=1,15
        c(i,j)=0.d0
       enddo
      enddo

      do i=1,n  
       if(num_c(i).eq.1)  t(i)  = a
       if(num_c(i).eq.2)  t(i)  = ax 
       if(num_c(i).eq.3)  t(i)  = ay
       if(num_c(i).eq.4)  t(i)  = axx
       if(num_c(i).eq.5)  t(i)  = axy
       if(num_c(i).eq.6)  t(i)  = ayy
       if(num_c(i).eq.7)  t(i)  = ax3
       if(num_c(i).eq.8)  t(i)  = ax2y
       if(num_c(i).eq.9)  t(i)  = axy2
       if(num_c(i).eq.10) t(i)  = ay3
       if(num_c(i).eq.11) t(i)  = ax4
       if(num_c(i).eq.12) t(i)  = ax3y
       if(num_c(i).eq.13) t(i)  = ax2y2
       if(num_c(i).eq.14) t(i)  = axy3
       if(num_c(i).eq.15) t(i)  = ay4
      enddo   
      do j=1,n
       do k=1,n
        do i=1,nstar
         c(j,k)=c(j,k)+temp(j,k,i)
        enddo
c        print *,j,k,c(j,k),t(j)
       enddo
      enddo

      call make_inverse(n,c,inv)
      do i=1,n
       eq(i)=0.d0
       do j=1,n
        eq(i)=eq(i)+t(j)*inv(i,j)
c        print *,i,j,num_c(i),num_c(j),t(j)
       enddo
      enddo
      return
      end

