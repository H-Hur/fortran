c     for fortran95.
c     The firstrun at 06.30.2008 
c     This progran need wcout2.fs,coeff.302.figSung2001.ZAMS.cc as
c     iniput files.
c     06,30,2008 last updated by Hur Hyeon-Oh
      parameter (nf=30,ns=100000)
      character input*16,coeff*16,filter(nf)*10
      real*8 airmass(nf),UT(nf)
      integer num_f(nf),num_c(nf),ntransform

c      print *,'Input data file name=?'
c      read(*,'(a16)') input
      write(input,'(a16)') 'wcout2.fs       '
      call readwchead(input,airmass,UT,filter)

c      write(coeff,'(a16)') 'coeff.228.fig_vi'
      print *,'Input coeff data file name=?'
      read(*,'(a16)') coeff
      call readcoeff(coeff,num_f,num_c,ntransform)

      write(input,'(a16)') 'Sung2001.ZAMS.cc'
      call read_ZAMS(input)

      call cal_outmag(ntransform,airmass,UT,filter)




      stop
      end

      subroutine cal_outmag(ntrans,airmass,UT,filter)
      parameter (nf=30,ns=100000)
      integer ntrans,nreiter,niter,infor_iter(nf,ns),set_f,nfail,auto
      character filter(nf)*10,cname(nf)*10,H*1
      real*8 airmass(nf),UT(nf),xwt,ywt
      real*8 color(nf),Q,slope,red(ns,2),wt,wtsum,ebvwsum
      real*8 ebv_mean,ncorred
      real*8 dmag(nf),omag1(nf)
      real bv0,ub0,bv1,ub1,c_bv(21),c_u(21),f,ebv(ns),x1,y1,ebv_in
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),coeff1(nf),
     &coeff2(nf),zeropoint(nf),delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf)
      common /ncw/nc_U,nc_B,nc_V,nc_R,nc_I,nc_Ha,nc_K
     &,nfilter,nc_H,nc_J,nstars
      common /wcdata/ID(nf,ns),xc(nf,ns),yc(nf,ns),
     &cmag(nf,ns),cmerr(nf,ns),omag(nf,ns)
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),nzams
      data c_bv/-.35d0,-.33d0,-0.3d0,-0.2d0,-0.1d0, 0.0d0, 0.1d0, 0.2d0,
     &0.3d0, 0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,1.2d0,1.4d0,1.6d0
     &,1.8d0,2.0d0/
      data c_u/-.132d0,-.132d0,-.126d0,-.074d0,-.030d0,0.002d0,0.024d0,
     &0.026d0,0.010d0,-.017d0,-.06d0,-.104d0,-.132d0,-.132d0,-.132d0,
     &-.132d0,-.132d0,-.132d0,-.132d0,-.132d0,-.132d0/

      do i=1,nfilter
       do j=1,ntrans
        if(num_fil(j).eq.0) then
         write(cname(j),'(a10)') '   none   '   
        elseif(num_fil(j).eq.1) then
         write(cname(j),'(a10)') '     U-B  '
        elseif(num_fil(j).eq.2) then
         write(cname(j),'(a10)') '     B-V  '
        elseif(num_fil(j).eq.3) then
         write(cname(j),'(a10)') '     V-R  '
        elseif(num_fil(j).eq.4) then
         write(cname(j),'(a10)') '     R-I  '
        elseif(num_fil(j).eq.5) then
         write(cname(j),'(a10)') '     V-I  '      
        endif
       enddo
      enddo

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c     write head
      open(41,file='out1.aex')
      write(41,*) ' '
      write(41,*) " F : Filrer"
      write(41,*) " C : Color"
      write(41,*) " T : UT"
      write(41,*) " X : Airmass"
      write(41,*) " *"
      write(41,*) ' ID(of all detected stars)'
      write(41,*) ' Obtained the outside magnitude or not (0 = not)'
      write(41,*) ' Outside magnitude(aperture corrected)'
      write(41,*) ' Instrument magnitude(aperture corrected)'
      write(41,*) ' Instrument magnitude error'
      write(41,*) ' X center(on the Xid master image)'
      write(41,*) ' Y center(on the Xid master image)'
      write(41,*) ' Mean e(B-V) from early type stars'
      write(41,*) ""
      write(41,*) ""
      write(41,'(t1,a1,30a10)')'F',(filter(i),i=1,nfilter)
      write(41,'(t1,a1,30a10)') 'C',(cname(i),i=1,nfilter)
      write(41,'(t1,a1,30f10.6)')'T',(UT(i),i=1,nfilter)
      write(41,'(t1,a1,30f10.7)')'X',(airmass(i),i=1,nfilter)
      write(41,'(t1,a1)') "*"
      
      niter=0

      print *,'Do you want to correct the non-linear transform of Landol
     &t U filter? (yes = 1)'
      read(*,*) set_f
      print *,'Determine e(B-V) from OB stars? (yes= 1 ,otherinput= 0)'
      read(*,*) auto
      if(auto.ne.1) then
       print *,'What is mean e(B-V)=?'
       read(*,*) ebv_in
      endif
c      if(set_f.eq.1) then
c       print *,'e(U-B)/e(B-V)=?'
c       read(*,*) slope
c      endif
      ncorred=0
 
      do i=1,nstars
       nreiter=0
       ebv(i)=0.
100    niter=niter+1
       do j=1,nfilter
        dmag(j)=0.
       enddo
       if(niter.le.1) then
        do j=1,nfilter
         infor_iter(j,i)=0
         omag(j,i)=cmag(j,i)
        enddo
       endif
       do k=1,ntrans
        color(nc_U)=omag(nc_U,i)-omag(nc_B,i)
        color(nc_B)=omag(nc_B,i)-omag(nc_V,i)
        if     (num_fil(k).eq.3.and.num_col(k).eq.2) then
         color(nc_V)=omag(nc_B,i)-omag(nc_V,i)
         elseif(num_fil(k).eq.3.and.num_col(k).eq.5) then
         color(nc_V)=omag(nc_V,i)-omag(nc_I,i)
        endif 
        if     (num_fil(k).eq.4.and.num_col(k).eq.3) then
         color(nc_R)=omag(nc_V,i)-omag(nc_R,i)
         elseif(num_fil(k).eq.4.and.num_col(k).eq.4) then
         color(nc_R)=omag(nc_R,i)-omag(nc_I,i)
        endif
        if     (num_fil(k).eq.5.and.num_col(k).eq.5) then
         color(nc_I)=omag(nc_V,i)-omag(nc_I,i)
         elseif(num_fil(k).eq.5.and.num_col(k).eq.4) then
         color(nc_I)=omag(nc_R,i)-omag(nc_I,i)
        endif
        color(nc_Ha)=0.d0
       enddo

        do jf=1,nfilter
         if(jf.eq.nc_V) then !V filter
          do k=1,ntrans 
           if(UT(nc_V).ge.timeran(k,1).and.UT(nc_V).lt.timeran(k,2).and.
     &     color(nc_V).ge.colorran(k,1).and.color(nc_V).lt.colorran(k,2)
     &     .and.num_fil(k).eq.3.and.num_col(k).eq.2.and.cmag(nc_V,i)
     &     .ne.0.and.cmag(nc_B,i).ne.0) then ! V : B-V
             infor_iter(nc_V,i)=1
             omag1(nc_V)=omag(nc_V,i)
             omag(nc_V,i)=cmag(nc_V,i)
     &       -(coeff1(k)-coeff2(k)*color(nc_V))*airmass(nc_V)
     &       +color(nc_V)*transcoeff(k)+transzero(k)
     &       +(UT(nc_V)-delta_ut(k))*timevar(k)+timezero(k)
     &       +zeropoint(k)
             dmag(nc_V)=sqrt((omag1(nc_V)-omag(nc_V,i))**2.d0)
           endif
           if(UT(nc_V).ge.timeran(k,1).and.UT(nc_V).lt.timeran(k,2).and.
     &     color(nc_V).ge.colorran(k,1).and.color(nc_V).lt.colorran(k,2)
     &     .and.num_fil(k).eq.3.and.num_col(k).eq.5.and.cmag(nc_V,i)
     &     .ne.0.and.cmag(nc_I,i).ne.0) then ! V : V-I
             infor_iter(nc_V,i)=1
             omag1(nc_V)=omag(nc_V,i)
             omag(nc_V,i)=cmag(nc_V,i)
     &       -(coeff1(k)-coeff2(k)*color(nc_V))*airmass(nc_V)
     &       +color(nc_V)*transcoeff(k)+transzero(k)
     &       +(UT(nc_V)-delta_ut(k))*timevar(k)+timezero(k)
     &       +zeropoint(k)
             dmag(nc_V)=sqrt((omag1(nc_V)-omag(nc_V,i))**2.d0)
           endif
          enddo
         endif

         if(jf.eq.nc_B) then !B filter
          do k=1,ntrans
           if(UT(nc_B).ge.timeran(k,1).and.UT(nc_B).lt.timeran(k,2).and.
     &     color(nc_B).ge.colorran(k,1).and.color(nc_B).lt.colorran(k,2)
     &     .and.num_fil(k).eq.2.and.num_col(k).eq.2.and.cmag(nc_B,i)
     &     .ne.0.and.infor_iter(nc_V,i).eq.1) then ! B : B-V
            infor_iter(nc_B,i)=1
            omag1(nc_B)=omag(nc_B,i)
            omag(nc_B,i)=cmag(nc_B,i)
     &      -(coeff1(k)-coeff2(k)*color(nc_B))*airmass(nc_B)
     &      +color(nc_B)*transcoeff(k)+transzero(k)
     &      +(UT(nc_B)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)
            dmag(nc_B)=sqrt((omag1(nc_B)-omag(nc_B,i))**2.d0)
           endif
          enddo
         endif

         if(jf.eq.nc_I) then !I filter
          do k=1,ntrans
           if(UT(nc_I).ge.timeran(k,1).and.UT(nc_I).lt.timeran(k,2).and.
     &     color(nc_I).ge.colorran(k,1).and.color(nc_I).lt.colorran(k,2)
     &     .and.num_fil(k).eq.5.and.num_col(k).eq.5.and.cmag(nc_I,i)
     &     .ne.0.and.infor_iter(nc_V,i).eq.1) then ! I : V-I
            infor_iter(nc_I,i)=1
            omag1(nc_I)=omag(nc_I,i)
            omag(nc_I,i)=cmag(nc_I,i)
     &      -(coeff1(k)-coeff2(k)*color(nc_I))*airmass(nc_I)
     &      +color(nc_I)*transcoeff(k)+transzero(k)
     &      +(UT(nc_I)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)
            dmag(nc_I)=sqrt((omag1(nc_I)-omag(nc_I,i))**2.d0)
           endif
           if(UT(nc_I).ge.timeran(k,1).and.UT(nc_I).lt.timeran(k,2).and.
     &     color(nc_I).ge.colorran(k,1).and.color(nc_I).lt.colorran(k,2)
     &     .and.num_fil(k).eq.5.and.num_col(k).eq.4.and.cmag(nc_I,i)
     &     .ne.0.and.infor_iter(nc_R,i).eq.1) then ! I : R-I
            infor_iter(nc_I,i)=1
            omag1(nc_I)=omag(nc_I,i)
            omag(nc_I,i)=cmag(nc_I,i)
     &      -(coeff1(k)-coeff2(k)*color(nc_I))*airmass(nc_I)
     &      +color(nc_I)*transcoeff(k)+transzero(k)
     &      +(UT(nc_I)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)
            dmag(nc_I)=sqrt((omag1(nc_I)-omag(nc_I,i))**2.d0)
           endif
          enddo
         endif

         if(jf.eq.nc_R) then !R filter
          do k=1,ntrans
           if(UT(nc_R).ge.timeran(k,1).and.UT(nc_R).lt.timeran(k,2).and.
     &     color(nc_R).ge.colorran(k,1).and.color(nc_R).lt.colorran(k,2)
     &     .and.num_fil(k).eq.4.and.num_col(k).eq.3.and.cmag(nc_R,i)
     &     .ne.0.and.infor_iter(nc_V,i).eq.1) then ! R : V-R
            infor_iter(nc_R,i)=1
            omag1(nc_R)=omag(nc_R,i)
            omag(nc_R,i)=cmag(nc_R,i)
     &      -(coeff1(k)-coeff2(k)*color(nc_R))*airmass(nc_R)
     &      +color(nc_R)*transcoeff(k)+transzero(k)
     &      +(UT(nc_R)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)
            dmag(nc_R)=sqrt((omag1(nc_R)-omag(nc_R,i))**2.d0)
           endif
           if(UT(nc_R).ge.timeran(k,1).and.UT(nc_R).lt.timeran(k,2).and.
     &     color(nc_R).ge.colorran(k,1).and.color(nc_R).lt.colorran(k,2)
     &     .and.num_fil(k).eq.5.and.num_col(k).eq.4.and.cmag(nc_R,i)
     &     .ne.0.and.infor_iter(nc_I,i).eq.1) then ! R : R-I
            infor_iter(nc_R,i)=1
            omag1(nc_R)=omag(nc_R,i)
            omag(nc_R,i)=cmag(nc_R,i)
     &      -(coeff1(k)-coeff2(k)*color(nc_R))*airmass(nc_R)
     &      +color(nc_R)*transcoeff(k)+transzero(k)
     &      +(UT(nc_R)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)
            dmag(nc_R)=sqrt((omag1(nc_R)-omag(nc_R,i))**2.d0)
           endif
          enddo
         endif

         if(jf.eq.nc_Ha) then !Ha filter
          do k=1,ntrans
           if(UT(nc_Ha).ge.timeran(k,1).and.UT(nc_Ha).lt.timeran(k,2)
     &     .and.color(nc_Ha).ge.colorran(k,1).and.color(nc_Ha).lt.
     &     colorran(k,2).and.num_fil(k).eq.6.and.num_col(k).eq.0
     &     .and.cmag(nc_Ha,i).ne.0) then ! Ha
            infor_iter(nc_Ha,i)=1
            omag1(nc_Ha)=omag(nc_Ha,i)
            omag(nc_Ha,i)=cmag(nc_Ha,i)
     &      -(coeff1(k)-coeff2(k)*color(nc_Ha))*airmass(nc_Ha)
     &      +color(nc_Ha)*transcoeff(k)+transzero(k)
     &      +(UT(nc_Ha)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)
            dmag(nc_Ha)=sqrt((omag1(nc_Ha)-omag(nc_Ha,i))**2.d0)
           endif
          enddo
         endif

         if(jf.eq.nc_U) then !U filter
          do k=1,ntrans
           if(UT(nc_U).ge.timeran(k,1).and.UT(nc_U).lt.timeran(k,2)
     &     .and.color(nc_U).ge.colorran(k,1).and.color(nc_U).lt.
     &     colorran(k,2).and.num_fil(k).eq.1.and.num_col(k).eq.1
     &     .and.cmag(nc_U,i).ne.0.and.infor_iter(nc_B,i).eq.1) then 
            infor_iter(nc_U,i)=1   ! U : U-B
            omag1(nc_U)=omag(nc_U,i) 
            omag(nc_U,i)=cmag(nc_U,i)
     &      -(coeff1(k)-coeff2(k)*color(nc_U))*airmass(nc_U)
     &      +color(nc_U)*transcoeff(k)+transzero(k)
     &      +(UT(nc_U)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)
            dmag(nc_U)=sqrt((omag1(nc_U)-omag(nc_U,i))**2.d0)
           endif
          enddo
         endif
        enddo
        nreiter=0
        do j=1,nfilter
         if(dmag(j).gt.0.0001d0) then
          nreiter=nreiter+1
         endif
        enddo
        if(nreiter.eq.0) then
         goto 190
         elseif(nreiter.ne.0) then
         goto 100
        endif
190   niter=0
       enddo
      if(set_f.ne.1) then
       goto 400
      endif

      wtsum=0.d0
      ebvwsum=0.d0
      niter=0
      n=0
      do i=1,nstars
       color(nc_U)=omag(nc_U,i)-omag(nc_B,i)
       color(nc_B)=omag(nc_B,i)-omag(nc_V,i)

200    niter=niter+1
       do k=1,ntrans
        if(UT(nc_U).ge.timeran(k,1).and.UT(nc_U).lt.timeran(k,2)
     &  .and.color(nc_U).ge.colorran(k,1).and.color(nc_U).lt.
     &  colorran(k,2).and.num_fil(k).eq.1.and.num_col(k).eq.1
     &  .and.cmag(nc_U,i).ne.0.and.infor_iter(nc_B,i).eq.1) then
         if(color(nc_U).lt.-0.05.and.color(nc_B).lt.0.5) then
          y1=color(nc_U) 
          x1=color(nc_B)
          call f_iter(y1,x1,bv0,ub0,nfail)
          if(nfail.eq.1) then
           ebv(i)=0.
           goto 290
          endif
          Q=ub0 -0.72d0*bv0
          if(Q.lt.-0.5) then
           n=n+1
           ebv(i)=x1-bv0
           xwt=(0.01/(cmerr(nc_B,i)+cmerr(nc_V,i)))**2.d0
           if(cmerr(nc_B,i)+cmerr(nc_V,i).lt.0.01) then
            xwt=1.d0
           endif
           ywt=(0.01/(cmerr(nc_U,i)+cmerr(nc_B,i)))**2.d0
           if(cmerr(nc_U,i)+cmerr(nc_B,i).lt.0.01) then
            ywt=1.d0
           endif
           wt=xwt*ywt
           wtsum=wtsum+wt
           ebvwsum=ebvwsum+ebv(i)*wt     
          endif 
c          print *,wt,ebv(i)
         endif  
        endif
       enddo
290   enddo
      ebv_mean=ebvwsum/wtsum
      print *,ebv_mean,ebvwsum,wtsum,n
      do i=1,nstars
       if(infor_iter(nc_U,i).eq.1) then
        if(ebv(i).ne.0.d0) then
         if(auto.ne.1) then
          ebv(i)=ebv_in
         endif
         bv0=omag(nc_B,i)-omag(nc_V,i)-ebv(i)
         call lint(bv0,c_bv,c_u,f,21,21)
         color(nc_U)=omag(nc_U,i)-omag(nc_B,i)+f
         color(nc_B)=omag(nc_B,i)-omag(nc_V,i)
         omag(nc_U,i)=omag(nc_U,i)+f
        else
         if(auto.ne.1) then
          ebv_mean=ebv_in
         endif
         bv0=omag(nc_B,i)-omag(nc_V,i)-ebv_mean
         call lint(bv0,c_bv,c_u,f,21,21)
         color(nc_U)=omag(nc_U,i)-omag(nc_B,i)+f
         color(nc_B)=omag(nc_B,i)-omag(nc_V,i)
         omag(nc_U,i)=omag(nc_U,i)+f
         ebv(i)=ebv_mean
        endif
       endif
      enddo  

400   do i=1,nstars 
        write(41,'(30i10)') (ID(j,i),j=1,nfilter)
        write(41,'(30i10)') (infor_iter(j,i),j=1,nfilter)
        write(41,'(30f10.3)') (omag(j,i),j=1,nfilter)
        write(41,'(30f10.3)') (cmag(j,i),j=1,nfilter)
        write(41,'(30f10.3)') (cmerr(j,i),j=1,nfilter)
        write(41,'(30f10.3)') (xc(j,i),j=1,nfilter)
        write(41,'(30f10.3)') (yc(j,i),j=1,nfilter)
        write(41,'(30f10.3)') (ebv(i),j=1,nfilter)
        write(41,*) ""
        write(41,*) ""
      enddo
      close(41)
      return
      end

      subroutine f_iter(ub_s0,bv_s,x,y,fail)
      parameter (ns=100000)
      real bv_s,f0,f,c_bv(21),c_u(21)
      real x,x0,y,ub_s0
      integer niter,fail
      common /zams/zMv(34),zVI(34),zBV(34),zUB(34)
      data c_bv/-.35d0,-.33d0,-0.3d0,-0.2d0,-0.1d0, 0.0d0, 0.1d0, 0.2d0,
     &0.3d0, 0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,1.2d0,1.4d0,1.6d0
     &,1.8d0,2.0d0/
      data c_u/-.132d0,-.132d0,-.126d0,-.074d0,-.030d0,0.002d0,0.024d0,
     &0.026d0,0.010d0,-.017d0,-.06d0,-.104d0,-.132d0,-.132d0,-.132d0,
     &-.132d0,-.132d0,-.132d0,-.132d0,-.132d0,-.132d0/
      call lint(bv_s,c_bv,c_u,f,21,21)
      x=bv_s
      fail=0
      do i=1,100
       f0=f
       x0=x
       ub_s=ub_s0+f
       call cor_red(bv_s,ub_s,x,y)
       call lint(x,c_bv,c_u,f,21,21)
       if(sqrt((f-f0)**2.d0).le.0.001.and.sqrt((x-x0)*2.d0).le.
     & 0.001) then
        goto 500
       endif
c       print *,i,x,x0,f,f0
      enddo
      fail=1
500   ub_int=y
      bv_int=x
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
      return
      end


      subroutine readcoeff(coeff,num_f,num_c,ntransform)
      parameter (nf=30,ns=100000)
      character coeff*16,H*1
      integer ntrans,num_f(nf),num_c(nf),ntransform,numf
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),coeff1(nf),
     &coeff2(nf),zeropoint(nf),delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf)
      open(31,file=coeff)
      ntrans=0
      numf=1
      do i=1,50
       read(31,'(a1)',end=500) H
       if(H.eq.'@'.or.H.eq.'#') then
        goto 100
        else
        backspace(31)
        ntrans=ntrans+1
        read(31,*) num_fil(ntrans),num_col(ntrans),colorran(ntrans,1),
     &  colorran(ntrans,2),transcoeff(ntrans),transzero(ntrans),
     &  coeff1(ntrans),coeff2(ntrans),timevar(ntrans),timeran(ntrans,1),
     &  timeran(ntrans,2),timezero(ntrans),delta_ut(ntrans),
     &  zeropoint(ntrans)
        num_f(ntrans)=num_fil(ntrans)
        num_c(ntrans)=num_col(ntrans)
        num_f(ntrans)=num_fil(ntrans)
        num_c(ntrans)=num_col(ntrans)
       endif
100   enddo 
500   close(31)
      ntransform=ntrans
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


      subroutine readwchead(input,airmass,UT,filter)
      parameter (nf=30,ns=100000) 
      character header*200,fil(nf)*10,input*16,H*1,coll(10)*1
      character filter(nf)*10
      integer nH
      real*8 airmass(nf),UT(nf)
      common /ncw/nc_U,nc_B,nc_V,nc_R,nc_I,nc_Ha,nc_K
     &,nfilter,nc_H,nc_J,nstars
      common /wcdata/ID(nf,ns),xc(nf,ns),yc(nf,ns),
     &cmag(nf,ns),cmerr(nf,ns)
      open(31,file=input)
      nfilter=0
      nstars=0

      do i=1,ns
       read(31,'(a1)') H
       if(H.eq.'F') then
        backspace(31)
        read(31,'(t2,30a10)') (filter(m),m=1,nf)
        do j=1,nf
         read(filter(j),'(10a1)') (coll(m),m=1,10)
         do k=1,10
          if(coll(k).eq.'U') then
           nfilter=nfilter+1
           nc_U=nfilter
           elseif(coll(k).eq.'B') then
           nfilter=nfilter+1
           nc_B=nfilter
           elseif(coll(k).eq.'V') then
           nfilter=nfilter+1
           nc_V=nfilter
           elseif(coll(k).eq.'R') then
           nfilter=nfilter+1
           nc_R=nfilter
           elseif(coll(k).eq.'I') then
           nfilter=nfilter+1
           nc_I=nfilter
           elseif(coll(k).eq.'H'.and.coll(k+1).eq.'a') then
           nfilter=nfilter+1
           nc_Ha=nfilter
          endif
         enddo
        enddo     
        elseif(H.eq.'X') then
        backspace(31)
        read(31,'(t2,30a10)') (fil(m),m=1,nf)
        do j=1,nfilter
         read(fil(j),*) airmass(j)
        enddo
        elseif(H.eq.'T') then
        backspace(31)
        read(31,'(t2,30a10)') (fil(m),m=1,nf)
        do j=1,nfilter
         read(fil(j),*) UT(j)
        enddo
        elseif(H.eq.'*') then
        goto 100
       endif
      enddo

      nstars=0
100   do i=1,ns
       read(31,*,end=300) (ID(j,i),j=1,nfilter)
       read(31,*,end=300) nH
       read(31,*,end=300) (xc(j,i),j=1,nfilter)
       read(31,*,end=300) (yc(j,i),j=1,nfilter)
       read(31,*,end=300) (cmag(j,i),j=1,nfilter)
c       print *,i,j,cmag(j,i)
       read(31,*,end=300) (cmerr(j,i),j=1,nfilter)
       read(31,'(a1)',end=300) H
       read(31,'(a1)',end=300) H
       read(31,'(a1)',end=300) H
       read(31,'(a1)',end=300) H
       nstars=nstars+1
      enddo
300   close(31)
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
       if(x-xt1(i).lt.0.) then
        goto 11
        elseif(x-xt1(i).eq.0.) then
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

