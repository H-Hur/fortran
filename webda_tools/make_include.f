c for only fortran77 to draw color-magnitude daigram.
c f77 -o 'make_include' 'make_include.f' -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11
c g77 -o make_include make_include.f -L/Applications/scisoft/i386/Packages/sm/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 -L/Users/apple/utils/sm_mac/scisoft/lib -laquaterm

      parameter (ns=100000)
      integer nstar,no(ns),nobs(ns,4)
      real v(ns),bv(ns),ub(ns),pm(ns),vi(ns)
      character H*1
      call read_ZAMS()
      open(21,file='combine.dat')
      read(21,'(t1,a1)') H
      do i=1,ns
       read(21,*,end=99) no(i),v(i),bv(i),ub(i),vi(i),(nobs(i,j),j=1,4)
     & ,pm(i)
      enddo
99    nstar=i-1
      close(21)
      call make_include(v,bv,ub,pm,no,nobs,nstar)
      stop
      end

      subroutine make_include(v,bv,ub,pm,no,nobs,nstar)
      parameter (ns=100000)
      integer nstar,no(ns),nobs(ns,4)
      real v(ns),bv(ns),ub(ns),pm(ns),ebv0,ebv1
      common /zams/zam(41,4)
      real pp(1),ebv,eub,ubz0,ubz1,zbv1(41),zub1(41),zbv0(41),zub0(41),q
      call sm_device ('postlandfile include.ps')
      call sm_graphics
      call sm_erase()
      call sm_location(3000,30000,3000,30000)
      call sm_limits(-0.5,0.4,1.,-1.2)
      call sm_box(1,2,0,0)
      call sm_xlabel('B-V')
      call sm_ylabel('U-B')
      print *,'E(B-V)0, E(B-V)1 =?'
      read(*,*) ebv0,ebv1
      open(21,file='include.dat')
      write(21,'(t1,a4)') '# No'
      do i=1,41
       zbv0(i)=zam(i,2)+ebv0
       zub0(i)=zam(i,3)+ebv0*0.72
       zbv1(i)=zam(i,2)+ebv1
       zub1(i)=zam(i,3)+ebv1*0.72
      enddo
      call sm_ctype('black')
      call sm_expand(1.)
      do i=1,nstar
       call lint(bv(i),zbv0,zub0,ubz0,41,41)
       call lint(bv(i),zbv1,zub1,ubz1,41,41)
       q=ub(i)-0.72*bv(i)
       if(nobs(i,2).ne.0.and.nobs(i,3).ne.0.and.
     & ubz1.le.ub(i).and.ubz0.ge.ub(i).and.
     & ((bv(i).gt.0.5.and.q.lt.-0.5).or.ub(i).lt.-0.1)) then
        pp(1)=403.+0.7
        write(21,*) no(i)
       else
        pp(1)=401.+0.1
       endif
       call sm_ptype(pp,1)
       if(nobs(i,2).ne.0.and.nobs(i,3).ne.0) then
        call sm_points(bv(i),ub(i),1)
       endif
      enddo
      call draw_zams(2,3,0.,0.,3.,0)
      call draw_zams(2,3,ebv0,ebv0*0.72,1.,0)
      call draw_zams(2,3,ebv1,ebv1*0.72,1.,0)
      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      close(21)
      return
      end


      subroutine draw_zams(nx,ny,dx,dy,weight,ltype)
      common /zams/zam(41,4)
      integer nx,ny,ltype
      real dx,dy,x(41),y(41),weight
      call sm_lweight(weight)
      call sm_ltype(ltype)
      call sm_ctype('black')
      do i=1,41
       x(i)=zam(i,nx)+dx
       y(i)=zam(i,ny)+dy
      enddo
      call sm_conn(x,y,41)
      call sm_lweight(0.5)
      call sm_ltype(0)
      return
      end


      subroutine read_ZAMS()
      character H*1
      common /zams/zam(41,4)
      open(31,file='Sung2013.ZAMS.cc')
      read(31,'(a1)') H
      do i=1,41
       read(31,*,end=100) zam(i,1),zam(i,4),zam(i,2),zam(i,3)
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
12    y=yt1(i)
      goto 14
10    continue
11    y=yt1(k)+(yt1(k+1)-yt1(k))*((x-xt1(k))/(xt1(k+1)-xt1(k)))
14    return
      end


