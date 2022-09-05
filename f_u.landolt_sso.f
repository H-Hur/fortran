! To calculate the f(B-V)0, which is a non-linear factor caused by
! balmer discontinuty on the Landolt standard stars observed using SSO
! U filter system.

! This program is made for Landolt standard stars, not for target stars.
! First made date : 24th, Apr, 2008 by Hyeon-Oh Hur
      parameter(ns=1000)
      integer nobs(ns),date(ns),ncollum,nc_star,nc_mag,nc_merr,nc_V
      integer nc_bv,nc_ub,nc_vr,nc_vi,nc_ut,nc_air,nc_wei,nc_date
      integer nc_nobs,nstar
      real*8 c_bv(21),c_u(21),mag(ns),bv(ns),ub(ns),vr(ns),vi(ns)
      real*8 ut(ns),air(ns),weight(ns),fbv(ns),dat(20,ns),v(ns)
      real*8 merr(ns),wei(ns),f_col,bv_col
      character input*16,header*200,coll*3,star(ns)*16,line*200
      character H(200)*1
      data c_bv/-.35d0,-.33d0,-0.3d0,-0.2d0,-0.1d0, 0.0d0, 0.1d0, 0.2d0,
     &0.3d0, 0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,1.2d0,1.4d0,1.6d0
     &,1.8d0,2.0d0/
      data c_u/-.132d0,-.132d0,-.126d0,-.074d0,-.030d0,0.002d0,0.024d0,
     &0.026d0,0.010d0,-.017d0,-.06d0,-.104d0,-.132d0,-.132d0,-.132d0,
     &-.132d0,-.132d0,-.132d0,-.132d0,-.132d0,-.132d0/

      print *,'Input file name=?'
      read(*,'(a16)') input
      open(31,file=input)
      
      read(31,'(a200)') header
      ncollum=0
      read(header,'(200a1)') (H(i),i=1,200)
      do i=1,200
       if(H(i).eq.' '.or.H(i).eq.'#') then
        if(H(i+1).ne.' ') then
         write(coll,'(3a1)') H(i+1),H(i+2),H(i+3)
         if(coll.eq.'Sta'.or.coll.eq.'sta') then
          nc_star=ncollum
          elseif(coll.eq.'Mag'.or.coll.eq.'mag') then
          ncollum=ncollum+1
          nc_mag=ncollum
          elseif(coll.eq.'Mer'.or.coll.eq.'mer') then
          ncollum=ncollum+1
          nc_merr=ncollum
          elseif(coll.eq.'V  '.or.coll.eq.'v  ') then
          ncollum=ncollum+1
          nc_v=ncollum
          elseif(coll.eq.'B-V'.or.coll.eq.'b-v') then
          ncollum=ncollum+1
          nc_bv=ncollum
          elseif(coll.eq.'U-B'.or.coll.eq.'u-b') then
          ncollum=ncollum+1
          nc_ub=ncollum
          elseif(coll.eq.'V-R'.or.coll.eq.'v-r') then
          ncollum=ncollum+1
          nc_vr=ncollum
          elseif(coll.eq.'V-I'.or.coll.eq.'v-i') then
          ncollum=ncollum+1
          nc_vi=ncollum
          elseif(coll.eq.'UT '.or.coll.eq.'ut '.or.coll.eq.'UT('.or.
     &    coll.eq.'ut(') then
          ncollum=ncollum+1
          nc_ut=ncollum
          elseif(coll.eq.'Air'.or.coll.eq.'air'.or.coll.eq.'X  '.or.
     &    coll.eq.'x  ') then
          ncollum=ncollum+1
          nc_air=ncollum
          elseif(coll.eq.'wei'.or.coll.eq.'Wei') then
          ncollum=ncollum+1
          nc_wei=ncollum
          elseif(coll.eq.'dat'.or.coll.eq.'Dat') then
          ncollum=ncollum+1
          nc_date=ncollum
          elseif(coll.eq.'Nob'.or.coll.eq.'nob') then
          ncollum=ncollum+1
          nc_nobs=ncollum
         endif
        endif
       endif
      enddo

      do i=1,ns
       read(31,'(T1,a16,T17,a200)',end=100) star(i),line
       read(line,*) (dat(j,i),j=1,ncollum)    
       mag(i)=dat(nc_mag,i)
       merr(i)=dat(nc_merr,i)
       bv(i)=dat(nc_bv,i)
       ub(i)=dat(nc_ub,i)
       vr(i)=dat(nc_vr,i)
       vi(i)=dat(nc_vi,i)
       v(i)=dat(nc_v,i)
       ut(i)=dat(nc_ut,i)
       wei(i)=dat(nc_wei,i)
       air(i)=dat(nc_air,i)
       date(i)=dat(nc_date,i)
       nobs(i)=dat(nc_nobs,i)
      enddo
100   nstar=i-1      
      close(31)


      open(32,file='output.f_u')
      write(32,*) 'Star             Mag      Merr   V      B-V    U-B  
     & V-R    V-I   f[(B-V)0]  Nobs  UT          Airmass  Weight date'
30    format(A16,f8.3,f7.3,f8.3,5f7.3,i8,2f12.7,f5.1,i5)
      do i=1,nstar
       bv_col=bv(i)
       call lint(bv_col,c_bv,c_u,f_col,21,21)
       fbv(i)=f_col
       write(32,30) star(i),mag(i),merr(i),v(i),bv(i),ub(i),vr(i),vi(i)
     & ,fbv(i),nobs(i),ut(i),air(i),wei(i),date(i)      
      enddo
      

      stop
      end

! linear interpolation program for fortran 95.
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
12    y=yt1(i)
      goto 14
10    continue
11    y=yt1(k)+(yt1(k+1)-yt1(k))*((x-xt1(k))/(xt1(k+1)-xt1(k)))
14    return
      end





