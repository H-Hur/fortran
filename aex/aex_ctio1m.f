c     for fortran95.
c     The firstrun at 06.30.2008 
c     This progran need wcout2.wc,coeff.ctio.fig                    
c     iniput files.
c     06,30,2008 last updated by Hur Hyeon-Oh
c     f95 -o aex_ctio1m ~/fortran/subroutine/aex/aex_ctio1m.f 
c     f95 -o aex_ctio1m  aex_ctio1m.f 
      parameter (nf=30,ns=100000)
      character input*16,coeff*16,filter(nf)*10
      real*8 airmass(nf),UT(nf)
      integer num_f(nf),num_c(nf),ntransform
      call read_ZAMS()
      call reddening_input()
      write(input,'(a16)') 'wcout2.wc                 '
      call readwchead(input,airmass,UT,filter)

      print *,'Input coeff data file name=?'
      read(*,'(a16)') coeff
 
      call readcoeff(coeff,num_f,num_c,ntransform)
      call cal_outmag(input,coeff,ntransform,airmass,UT,filter)


      stop
      end

      subroutine reddening_input()
      parameter (ns=100000)
      common/red_inp/ebv_fg,ebv_cl,ebv_bg,rv_fg,rv_cl,rv_bg,id_wc(ns),
     &ebv_wc(ns),nredip
      character H*1
      ebv_fg=0.
      ebv_cl=0.
      ebv_bg=0.
      rv_fg=3.1
      rv_cl=3.1
      rv_bg=3.1
      open(21,file='reddening.aex')
      read(21,*,end=99) ebv_fg,ebv_cl,ebv_bg
      read(21,*,end=99) rv_fg,rv_cl,rv_bg
      nredip=0
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.ne.'#') then
        backspace(21)
        nredip=nredip+1
        read(21,*) id_wc(nredip),ebv_wc(nredip)
       endif
      enddo
99    close(21)
      do i=1,nredip
       print *,id_wc(i),ebv_wc(i)
      enddo
      return
      end

      subroutine cal_outmag(input,coeff,ntrans,airmass,UT,filter)
      parameter (nf=30,ns=100000)
      integer ntrans,nreiter,niter,infor_iter(nf,ns),set_f,nfail,auto
      integer ncolor,nvcal(ns,2),nvvi,niri
      integer nuv,nbv,nvi,nri,nrh,nid
      character filter(nf)*10,cname(nf)*10,H*1,input*16,coeff*16
      real*8 airmass(nf),UT(nf),xwt,ywt
      real*8 color(nf),Q,slope,red(ns,2),wt,wtsum,ebvwsum,irie,come
      real*8 ebv_mean,ncorred,omag(nf,ns),ovvi,oiri,ovvi2,oiri2
      real*8 dmag(nf),omag1(nf),dvvi,diri,cvvi,ciri,vbve,vvie,ivie
      real*8 wub,wbv,wvi,wri,wrh,ovbv,weub,webv,wevi,weri,werh,rad
      real bv0,ub0,bv1,ub1,c_bv(16),c_u(16),f,ebv(ns),x1,y1,ebv_in
      real commag(5),rrr
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),coeff1(nf),
     &coeff2(nf),zeropoint(nf),delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf),r1coeff(nf),r2coeff(nf) 
      common /ncw/nc_U,nc_B,nc_V,nc_R,nc_I,nc_Ha,nc_K
     &,nfilter,nc_H,nc_J,nstars
      common /wcdata/ID(ns),xc(ns),yc(ns),xe(ns),ye(ns),nobs(nf,ns),
     &cmag(nf,ns),cmerr(nf,ns),id_f(nf,ns)
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),nzams
      data c_bv/-1.25,0.1,0.19783,0.22391,0.31087,0.39783,0.48478,
     &0.57174,0.65870,0.74565,0.83261,0.91957,1.00652,1.09348,1.2674,4./
      data c_u /0.   ,0. ,0.     ,0.00095,0.00918,0.01614,0.02120,
     &0.02468,0.02500,0.02247,0.01820,0.01139,0.00443,0.00000,0.0000,0./
      common/red_inp/ebv_fg,ebv_cl,ebv_bg,rv_fg,rv_cl,rv_bg,id_wc(ns),
     &ebv_wc(ns),nredip

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
      write(41,*) " T : UT"
      write(41,*) " X : Airmass"
      write(41,*) " *"
      write(41,*) ' ID,ID(of all detected stars)'
      write(41,*) ' Xc, Obtained the outside magnitude or not (0 = not)'
      write(41,*) ' Yc, Outside magnitude(aperture corrected)'
      write(41,*) ' Xc error,Instrument magnitude(aperture corrected)'
      write(41,*) ' Yc error, Instrument magnitude error'
      write(41,*) ' Distance from 2032,2032                '
      write(41,*) ""
      write(41,*) ""
      write(41,'(t1,a1,t11,30a10)')'F',(filter(i),i=1,nfilter)
      write(41,'(t1,a1,t11,30f10.6)')'T',(UT(i),i=1,nfilter)
      write(41,'(t1,a1,t11,30f10.7)')'X',(airmass(i),i=1,nfilter)
      write(41,'(t1,a1)') "*"

      open(42,file='out2.aex')
      write(42,'(t1,a10,a16)') '# input = ',input
      write(42,'(t1,a10,a16)') '# coeff = ',coeff
      write(42,'(t1,a7,t11,30a10)')'#Filter',(filter(i),i=1,nfilter)
      write(42,'(t1,a5,t11,30f10.6)')'#Time',(UT(i),i=1,nfilter)
      write(42,'(t1,a8,t11,30f10.7)')'#Airmass',(airmass(i),i=1,nfilter)



      write(42,'(t1,a130)')  '#rawID    ID    Xc       Yc      V       U
     &-B    B-V    V-I    R-I    Ha    eV    e(U-B) e(B-V) e(V-I) e(R-I)
     & e(Ha)   nobs-------------                                       '
50    format(t1,2i6,2f9.3,f8.3,11f7.3,6i3)
      niter=0

      print *,'Do you want to correct the non-linear transform of 
     &U filter? (yes = 1)'
      read(*,*) set_f
c      print *,'Determine e(B-V) from OB stars? (yes= 1 ,otherinput= 0)'
c      read(*,*) auto
      auto=0
c      if(auto.ne.1) then
c       print *,'What is mean e(B-V)=?'
c       read(*,*) ebv_in
c      endif
      ebv_in=ebv_cl
c      if(set_f.eq.1) then
c       print *,'e(U-B)/e(B-V)=?'
c       read(*,*) slope
c      endif
      slope=0.72
      ncorred=0
      nid=0
      do i=1,nstars
       rrr=sqrt((xc(i)-2032.)**2+(yc(i)-2032.)**2)/1000.
       nvvi=0
       niri=0
       nvbv=0
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
           ovvi=omag(nc_V,i)
           cvvi=omag(nc_V,i)-omag(nc_I,i)
           oiri=omag(nc_I,i)
           ciri=omag(nc_R,i)-omag(nc_I,i)
        enddo
       endif
       do k=1,ntrans
        color(nc_U)=omag(nc_U,i)-omag(nc_B,i)
        color(nc_B)=omag(nc_B,i)-omag(nc_V,i)
        color(nc_Ha)=0.d0
       enddo


       do jf=1,nfilter
        if(jf.eq.nc_V) then !V filter
         do k=1,ntrans 
          color(nc_V)=omag(nc_B,i)-omag(nc_V,i)
          if(UT(nc_V).ge.timeran(k,1).and.UT(nc_V).lt.timeran(k,2).and.
     &    color(nc_V).ge.colorran(k,1).and.color(nc_V).lt.colorran(k,2)
     &    .and.num_fil(k).eq.3.and.num_col(k).eq.2.and.cmag(nc_V,i)
     &    .ne.0.and.cmag(nc_B,i).ne.0) then ! V : B-V
            nvbv=1
            infor_iter(nc_V,i)=1
            omag1(nc_V)=omag(nc_V,i)
            omag(nc_V,i)=cmag(nc_V,i)
     &      -(coeff1(k)-coeff2(k)*color(nc_V))*airmass(nc_V)
     &      +color(nc_V)*transcoeff(k)+transzero(k)
     &      +(UT(nc_V)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr

            dmag(nc_V)=sqrt((omag1(nc_V)-omag(nc_V,i))**2.d0)
          endif
          color(nc_V)=omag(nc_V,i)-omag(nc_I,i)
          cvvi=ovvi-omag(nc_i,i)

          if(UT(nc_V).ge.timeran(k,1).and.UT(nc_V).lt.timeran(k,2).and.
     &    cvvi.ge.colorran(k,1).and.cvvi.lt.colorran(k,2)
     &    .and.num_fil(k).eq.3.and.num_col(k).eq.5.and.cmag(nc_V,i)
     &    .ne.0.and.cmag(nc_I,i).ne.0.and.infor_iter(nc_V,i).eq.1
     &    .and.infor_iter(nc_I,i).eq.1.and.cmag(nc_I,i).ne.0)then!
           nvvi=1
           cvvi=ovvi-omag(nc_i,i)
           ovvi2=ovvi
           ovvi=cmag(nc_V,i)
     &      -(coeff1(k)-coeff2(k)*cvvi)*airmass(nc_V)
     &      +cvvi*transcoeff(k)+transzero(k)
     &      +(UT(nc_V)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr

           dvvi=sqrt((ovvi2-ovvi)**2.d0)
          endif
          if(UT(nc_V).ge.timeran(k,1).and.UT(nc_V).lt.timeran(k,2).and.
     &    color(nc_V).ge.colorran(k,1).and.color(nc_V).lt.colorran(k,2)
     &    .and.num_fil(k).eq.3.and.num_col(k).eq.5.and.cmag(nc_V,i)
     &    .ne.0.and.cmag(nc_I,i).ne.0.and.infor_iter(nc_V,i).eq.0)then! V : V-I
            infor_iter(nc_V,i)=1
            omag1(nc_V)=omag(nc_V,i)
            omag(nc_V,i)=cmag(nc_V,i)
     &      -(coeff1(k)-coeff2(k)*color(nc_V))*airmass(nc_V)
     &      +color(nc_V)*transcoeff(k)+transzero(k)
     &      +(UT(nc_V)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr
            dmag(nc_V)=sqrt((omag1(nc_V)-omag(nc_V,i))**2.d0)
          endif
         enddo
        endif

        if(jf.eq.nc_B) then !B filter
         do k=1,ntrans
          if(UT(nc_B).ge.timeran(k,1).and.UT(nc_B).lt.timeran(k,2).and.
     &    color(nc_B).ge.colorran(k,1).and.color(nc_B).lt.colorran(k,2)
     &    .and.num_fil(k).eq.2.and.num_col(k).eq.2.and.cmag(nc_B,i)
     &    .ne.0.and.infor_iter(nc_V,i).eq.1) then ! B : B-V
           infor_iter(nc_B,i)=1
           omag1(nc_B)=omag(nc_B,i)
           omag(nc_B,i)=cmag(nc_B,i)
     &     -(coeff1(k)-coeff2(k)*color(nc_B))*airmass(nc_B)
     &     +color(nc_B)*transcoeff(k)+transzero(k)
     &     +(UT(nc_B)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr

           dmag(nc_B)=sqrt((omag1(nc_B)-omag(nc_B,i))**2.d0)
          endif
         enddo
        endif

        if(jf.eq.nc_I) then !I filter
         do k=1,ntrans
          color(nc_I)=omag(nc_V,i)-omag(nc_I,i)
          if(UT(nc_I).ge.timeran(k,1).and.UT(nc_I).lt.timeran(k,2).and.
     &    color(nc_I).ge.colorran(k,1).and.color(nc_I).lt.colorran(k,2)
     &    .and.num_fil(k).eq.5.and.num_col(k).eq.5.and.cmag(nc_I,i)
     &    .ne.0.and.infor_iter(nc_V,i).eq.1) then ! I : V-I
           infor_iter(nc_I,i)=1
           omag1(nc_I)=omag(nc_I,i)
           omag(nc_I,i)=cmag(nc_I,i)
     &     -(coeff1(k)-coeff2(k)*color(nc_I))*airmass(nc_I)
     &     +color(nc_I)*transcoeff(k)+transzero(k)
     &     +(UT(nc_I)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr

           dmag(nc_I)=sqrt((omag1(nc_I)-omag(nc_I,i))**2.d0)

          endif
          color(nc_I)=omag(nc_R,i)-omag(nc_I,i)
          ciri=omag(nc_R,i)-oiri
          if(UT(nc_I).ge.timeran(k,1).and.UT(nc_I).lt.timeran(k,2).and.
     &    ciri.ge.colorran(k,1).and.ciri.lt.colorran(k,2)
     &    .and.num_fil(k).eq.5.and.num_col(k).eq.4.and.cmag(nc_I,i)
     &    .ne.0.and.infor_iter(nc_R,i).eq.1.and.infor_iter(nc_I,i).eq.
     &    1.and.cmag(nc_I,i).ne.0.) then ! I : R-I
           niri=1
           infor_iter(nc_I,i)=1
           ciri=omag(nc_R,i)-oiri
           oiri2=oiri
           oiri=cmag(nc_I,i)
     &     -(coeff1(k)-coeff2(k)*ciri)*airmass(nc_I)
     &     +ciri*transcoeff(k)+transzero(k)
     &     +(UT(nc_I)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr

           diri=sqrt((oiri2-oiri)**2.d0)
          endif

          if(UT(nc_I).ge.timeran(k,1).and.UT(nc_I).lt.timeran(k,2).and.
     &    color(nc_I).ge.colorran(k,1).and.color(nc_I).lt.colorran(k,2)
     &    .and.num_fil(k).eq.5.and.num_col(k).eq.4.and.cmag(nc_I,i)
     &    .ne.0.and.infor_iter(nc_R,i).eq.1.and.infor_iter(nc_I,i).eq.
     &    0) then ! I : R-I
           infor_iter(nc_I,i)=1
           omag1(nc_I)=omag(nc_I,i)
           omag(nc_I,i)=cmag(nc_I,i)
     &     -(coeff1(k)-coeff2(k)*color(nc_I))*airmass(nc_I)
     &     +color(nc_I)*transcoeff(k)+transzero(k)
     &     +(UT(nc_I)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr

           dmag(nc_I)=sqrt((omag1(nc_I)-omag(nc_I,i))**2.d0)
          endif
         enddo
        endif

        if(jf.eq.nc_R) then !R filter
         do k=1,ntrans
          color(nc_R)=omag(nc_R,i)-omag(nc_I,i)
          if(UT(nc_R).ge.timeran(k,1).and.UT(nc_R).lt.timeran(k,2).and.
     &    color(nc_R).ge.colorran(k,1).and.color(nc_R).lt.colorran(k,2)
     &    .and.num_fil(k).eq.4.and.num_col(k).eq.4.and.cmag(nc_R,i)
     &    .ne.0.and.infor_iter(nc_I,i).eq.1) then ! R : R-I
           infor_iter(nc_R,i)=1
           omag1(nc_R)=omag(nc_R,i)
           omag(nc_R,i)=cmag(nc_R,i)
     &     -(coeff1(k)-coeff2(k)*color(nc_R))*airmass(nc_R)
     &     +color(nc_R)*transcoeff(k)+transzero(k)
     &     +(UT(nc_R)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr

           dmag(nc_R)=sqrt((omag1(nc_R)-omag(nc_R,i))**2.d0)
          endif
         enddo
        endif

        if(jf.eq.nc_Ha) then !Ha filter
         do k=1,ntrans
          if(UT(nc_Ha).ge.timeran(k,1).and.UT(nc_Ha).lt.timeran(k,2)
     &    .and.color(nc_Ha).ge.colorran(k,1).and.color(nc_Ha).lt.
     &    colorran(k,2).and.num_fil(k).eq.6.and.num_col(k).eq.0
     &    .and.cmag(nc_Ha,i).ne.0) then ! Ha
           infor_iter(nc_Ha,i)=1
           omag1(nc_Ha)=omag(nc_Ha,i)
           omag(nc_Ha,i)=cmag(nc_Ha,i)
     &     -(coeff1(k)-coeff2(k)*color(nc_Ha))*airmass(nc_Ha)
     &     +color(nc_Ha)*transcoeff(k)+transzero(k)
     &     +(UT(nc_Ha)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr

           dmag(nc_Ha)=sqrt((omag1(nc_Ha)-omag(nc_Ha,i))**2.d0)
          endif
         enddo
        endif

        if(jf.eq.nc_U) then !U filter
         do k=1,ntrans
          if(UT(nc_U).ge.timeran(k,1).and.UT(nc_U).lt.timeran(k,2)
     &    .and.color(nc_U).ge.colorran(k,1).and.color(nc_U).lt.
     &    colorran(k,2).and.num_fil(k).eq.1.and.num_col(k).eq.1
     &    .and.cmag(nc_U,i).ne.0.and.infor_iter(nc_B,i).eq.1) then 
           infor_iter(nc_U,i)=1   ! U : U-B
           omag1(nc_U)=omag(nc_U,i) 
           omag(nc_U,i)=cmag(nc_U,i)
     &     -(coeff1(k)-coeff2(k)*color(nc_U))*airmass(nc_U)
     &     +color(nc_U)*transcoeff(k)+transzero(k)
     &     +(UT(nc_U)-delta_ut(k))*timevar(k)+timezero(k)
     &      +zeropoint(k)+r1coeff(k)*rrr+r2coeff(k)*rrr*rrr

           dmag(nc_U)=sqrt((omag1(nc_U)-omag(nc_U,i))**2.d0)
          endif
         enddo
        endif
       enddo
c      if(id(i).eq.206) then
c       write(*,'(i6,4f8.3,3i2)') ID(i),omag(nc_V,i),ovvi,omag(nc_I,i),
c     &  oiri,nvvi,niri,niter
c      endif

       nreiter=0
       do j=1,nfilter
        if(dmag(j).gt.0.001d0) then
         nreiter=nreiter+1
        endif
        if(niri.ne.0.and.diri.gt.0.001d0) then
         nreiter=nreiter+1
        endif
        if(nvvi.ne.0.and.dvvi.gt.0.001d0) then
         nreiter=nreiter+1
        endif
       enddo
       ovbv=omag(nc_v,i)
       if(nreiter.eq.0.or.niter.ge.100) then
        if(nvvi.eq.1.and.nvbv.eq.1) then
         call cal_error(cmerr(nc_V,i),cmerr(nc_B,i),vbve)
         call cal_error(cmerr(nc_V,i),cmerr(nc_I,i),vvie)
         if(vbve.lt.0.1) then
          vbve=1.d0
         else
          vbve=0.1d0
         endif
         if(vvie.lt.0.1) then
          vvie=1.d0
         else
          vvie=0.1d0
         endif
         omag(nc_V,i)=(omag(nc_V,i)*vbve+ovvi*vvie)/(vbve+vvie)
        elseif(nvvi.eq.1.and.nvbv.eq.0) then
         omag(nc_V,i)=ovvi
        endif
        if(niri.eq.1) then
         call cal_error(cmerr(nc_R,i),cmerr(nc_I,i),irie)
         call cal_error(cmerr(nc_V,i),cmerr(nc_I,i),ivie)
         if(irie.lt.0.1) then
          irie=1.d0
         else
          irie=0.1d0
         endif
         if(vvie.lt.0.1) then
          ivie=1.d0
         else
          ivie=0.1d0
         endif
         omag(nc_I,i)=(omag(nc_I,i)*ivie+oiri*irie)/(irie+ivie)
        endif
        goto 190
        elseif(nreiter.ne.0.or.niter.lt.100) then
        goto 100
       endif
190    niter=0
c       if(id(i).eq.206) then
c        write(*,'(i6,4f8.3,2i2)') ID(i),omag(nc_V,i),ovvi,omag(nc_I,i),
c     &  oiri,nvvi,niri
c       endif
       if(infor_iter(nc_v,i).ne.0.and.infor_iter(nc_b,i).ne.0.and.
     & infor_iter(nc_i,i).ne.0.) then!and.cvvi.ge.1.75) then
c        if(cvvi.gt.2.5.and.ovbv-ovvi.gt.0.005) then
c         write(*,'(I6,5f9.3)') id(i),cvvi,ovvi,ovbv,ovbv-ovvi
c     &   ,sqrt(2*cmerr(nc_V,i)**2+cmerr(nc_I,i)**2+cmerr(nc_B,i)**2)
c        endif
       endif
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
          if((Q.lt.-0.5.and.omag(nc_V,i).lt.10.d0).or.(y1.lt.0.1d0.and.
     &    omag(nc_V,i).lt.10.0d0)) then
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
          print *,id(i),wt,ebv(i),bv0
         endif  
        endif
       enddo
290    i290=i290
      enddo
      ebv_mean=ebvwsum/wtsum
      print *,ebv_mean,ebvwsum,wtsum,n
      do i=1,nstars
       if(infor_iter(nc_U,i).eq.1) then
        if(ebv(i).ne.0.d0) then
         bv0=omag(nc_B,i)-omag(nc_V,i)-ebv(i)
         call lint(bv0,c_bv,c_u,f,16,16)
         color(nc_U)=omag(nc_U,i)-omag(nc_B,i)+f
         color(nc_B)=omag(nc_B,i)-omag(nc_V,i)
         omag(nc_U,i)=omag(nc_U,i)+f
        else
         if(auto.ne.1) then
          ebv_mean=ebv_in
         endif
         bv0=omag(nc_B,i)-omag(nc_V,i)-ebv_mean
         if(omag(nc_V,i).lt.10.0d0) then
          bv0=omag(nc_B,i)-omag(nc_V,i)
         else
          if(auto.ne.1) then
           ebv_mean=ebv_in
          endif
          bv0=omag(nc_B,i)-omag(nc_V,i)-ebv_mean
         endif
         call lint(bv0,c_bv,c_u,f,21,21)
c         print *,id(i),color(nc_B),color(nc_U),bv0,f
         color(nc_U)=omag(nc_U,i)-omag(nc_B,i)+f
         color(nc_B)=omag(nc_B,i)-omag(nc_V,i)
         omag(nc_U,i)=omag(nc_U,i)+f
         ebv(i)=ebv_mean
        endif
       endif
      enddo
400   do i=1,nstars
        rad=sqrt((xc(i)-2032.d0)**2+(yc(i)-2032.d0)**2)
        write(41,'(t1,31i10)') i,(ID(i),j=1,nfilter)
        write(41,'(t1,f10.3,30i10)') xc(i),(infor_iter(j,i),j=1,nfilter)
        write(41,'(t1,31f10.3)') yc(i),(omag(j,i),j=1,nfilter)
        write(41,'(t1,31f10.3)') xe(i),(cmag(j,i),j=1,nfilter)
        write(41,'(t1,31f10.3)') ye(i),(cmerr(j,i),j=1,nfilter)
        write(41,'(t1,f10.3,30i10)') rad,(nobs(j,i),j=1,nfilter)
        write(41,*) ""
        write(41,*) ""
        write(41,*) ""
        write(41,*) ""
        do j=1,nfilter
         if(nobs(j,i).eq.0) omag(j,i)=0.d0
        enddo
        nbv=min(nobs(nc_v,i),nobs(nc_b,i))
        nub=min(nobs(nc_u,i),nobs(nc_b,i))
        nvi=min(nobs(nc_v,i),nobs(nc_i,i))
        nri=min(nobs(nc_r,i),nobs(nc_i,i))
        nrh=min(nobs(nc_r,i),nobs(nc_ha,i))
        if(infor_iter(nc_v,i).eq.0) nobs(nc_v,i)=0
        if(infor_iter(nc_u,i).eq.0) nub=0
        if(infor_iter(nc_b,i).eq.0) nub=0
        if(infor_iter(nc_b,i).eq.0) nbv=0
        if(infor_iter(nc_i,i).eq.0) nvi=0
        if(infor_iter(nc_i,i).eq.0) nri=0
        if(infor_iter(nc_r,i).eq.0) nri=0
        if(infor_iter(nc_r,i).eq.0) nrh=0

        wub=omag(nc_u,i)-omag(nc_b,i)
        wbv=omag(nc_b,i)-omag(nc_v,i)
        wvi=omag(nc_v,i)-omag(nc_i,i)
        wri=omag(nc_r,i)-omag(nc_i,i)
        wrh=omag(nc_r,i)-omag(nc_Ha,i)

        call cal_error2(cmerr(nc_u,i),cmerr(nc_b,i),weub)
        call cal_error2(cmerr(nc_b,i),cmerr(nc_v,i),webv)
        call cal_error2(cmerr(nc_v,i),cmerr(nc_i,i),wevi)
        call cal_error2(cmerr(nc_r,i),cmerr(nc_i,i),weri)
        call cal_error2(cmerr(nc_r,i),cmerr(nc_Ha,i),werh)

        if(nub.eq.0) wub=0.d0
        if(nbv.eq.0) wbv=0.d0
        if(nvi.eq.0) wvi=0.d0
        if(nri.eq.0) wri=0.d0
        if(nrh.eq.0) wrh=0.d0
        if(infor_iter(nc_v,i).eq.0) omag(nc_v,i)=0.d0
        if(infor_iter(nc_v,i).eq.0) cmerr(nc_v,i)=0.d0
        if(nobs(nc_v,i).ne.0) then 
         nid=nid+1
         write(42,50)i,nid,xc(i),yc(i),omag(nc_v,i),wub,wbv,wvi,wri,
     &   omag(nc_Ha,i),cmerr(nc_v,i),weub,webv,wevi,weri,cmerr(nc_ha,i),
     &   nobs(nc_v,i),nub,nbv,nvi,nri,nobs(nc_ha,i)
c         if(omag(nc_v,i).lt.15.and.wvi.gt.2.5) print*,i,omag(nc_v,i),wvi
        endif
      enddo
      close(41)
      close(42)
      return
      end
!     combine two errors.     
      subroutine cal_error2(e1,e2,e3)
      real e1,e2
      real*8 e3
      if(e1.le.0.d0.or.e2.le.0.d0) then
       e3= 0.d0
      else
       e3=sqrt(e1**2.d0+e2**2.d0)
      endif
      return
      end

!     combine two errors.     
      subroutine cal_error(e1,e2,e3)
      real*8 e1,e2,e3
      if(e1.lt.0.d0.or.e2.lt.0.d0) then
       e3=-1.d0
      else
       e3=sqrt(e1**2.d0+e2**2.d0)
      endif
      return
      end

      subroutine f_iter(ub_s0,bv_s,x,y,fail)
      real bv_s,f0,f,c_bv(16),c_u(16)
      real x,x0,y,ub_s0
      integer niter,fail
      common /zams/zMv(34),zVI(34),zBV(34),zUB(34)
      data c_bv/-1.25,0.1,0.19783,0.22391,0.31087,0.39783,0.48478,
     &0.57174,0.65870,0.74565,0.83261,0.91957,1.00652,1.09348,1.2674,4./
      data c_u /0.   ,0. ,0.     ,0.00095,0.00918,0.01614,0.02120,
     &0.02468,0.02500,0.02247,0.01820,0.01139,0.00443,0.00000,0.0000,0./
      call lint(bv_s,c_bv,c_u,f,16,16)
      x=bv_s
      fail=0
      do i=1,100
       f0=f
       x0=x
       ub_s=ub_s0+f
       call cor_red(bv_s,ub_s,x,y)
       call lint(x,c_bv,c_u,f,16,16)
       if(sqrt((f-f0)**2.d0).le.0.001.and.sqrt((x-x0)*2.d0).le.
     & 0.001) then
        goto 500
       endif
      enddo
      fail=1
500   ub_int=y
      bv_int=x
      return
      end

      subroutine cor_red(sBV,sUB,unredx,unredy)
      parameter (nf=10,ns=1000,maxiter=1000)
      common /zams/zMv(34),zVI(34),zBV(34),zUB(34)
      real fx,fy,dy,unredx,unredy,sBV,sUB,slope,eBV,eUB
      slope=0.72
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
     &,timevar(nf),transzero(nf),timezero(nf),r1coeff(nf),r2coeff(nf) 
      open(31,file=coeff)
      ntrans=0
      numf=1
      do i=1,50
       read(31,'(t1,a1)',end=500) H
       if(H.eq.'@'.or.H.eq.'#') then
        goto 100
        else
        backspace(31)
        ntrans=ntrans+1
        read(31,*) num_fil(ntrans),num_col(ntrans),colorran(ntrans,1),
     &  colorran(ntrans,2),transcoeff(ntrans),transzero(ntrans),
     &  coeff1(ntrans),coeff2(ntrans),timevar(ntrans),timeran(ntrans,1),
     &  timeran(ntrans,2),timezero(ntrans),delta_ut(ntrans),
     &  zeropoint(ntrans),r1coeff(ntrans),r2coeff(ntrans)
c        print *,num_fil(ntrans),num_col(ntrans),colorran(ntrans,1),
c     &  colorran(ntrans,2),transcoeff(ntrans),transzero(ntrans),
c     &  coeff1(ntrans),coeff2(ntrans),timevar(ntrans),timeran(ntrans,1),
c     &  timeran(ntrans,2),timezero(ntrans),delta_ut(ntrans),
c     &  zeropoint(ntrans),r1coeff(ntrans),r2coeff(ntrans)

        num_f(ntrans)=num_fil(ntrans)
        num_c(ntrans)=num_col(ntrans)
        num_f(ntrans)=num_fil(ntrans)
        num_c(ntrans)=num_col(ntrans)
       endif
100    i100=i100
      enddo 
500   close(31)
      ntransform=ntrans
      print *,ntransform,i
      return
      end


      subroutine read_ZAMS()
      character H*1,input*16
      common /zams/zMv(34),zVI(34),zBV(34),zUB(34)
      open(31,file='Sung2001.ZAMS.cc')
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
      common /wcdata/ID(ns),xc(ns),yc(ns),xe(ns),ye(ns),nobs(nf,ns),
     &cmag(nf,ns),cmerr(nf,ns),id_f(nf,ns)
      open(31,file=input)
      nfilter=0
      nstars=0

      do i=1,ns
       read(31,'(a1)') H
       if(H.eq.'F') then
        backspace(31)
        read(31,'(t11,30a10)') (filter(m),m=1,nf)
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
           elseif(coll(k).eq.'H'.and.coll(k+1).eq.'a'.and.coll(k-1).eq.
     &     ' ') then
           nfilter=nfilter+1
           nc_Ha=nfilter
          endif
         enddo
        enddo     
        elseif(H.eq.'X') then
        backspace(31)
        read(31,'(t11,30a10)') (fil(m),m=1,nf)
        do j=1,nfilter
         read(fil(j),*) airmass(j)
        enddo
        elseif(H.eq.'T') then
        backspace(31)
        read(31,'(t11,30a10)') (fil(m),m=1,nf)
        do j=1,nfilter
         read(fil(j),*) UT(j)
        enddo
        elseif(H.eq.'*') then
        goto 100
       endif
      enddo

      nstars=0
100   do i=1,ns
       read(31,*,end=300) ID(i),(id_f(j,i),j=1,nfilter)
       read(31,*,end=300) xc(i),(cmag(j,i),j=1,nfilter)
       read(31,*,end=300) yc(i),(cmerr(j,i),j=1,nfilter)
       read(31,*,end=300) xe(i),(nobs(j,i),j=1,nfilter)
       read(31,*,end=300) ye(i)
       read(31,'(a1)',end=300) H
       read(31,'(a1)',end=300) H
       read(31,'(a1)',end=300) H
       read(31,'(a1)',end=300) H
       read(31,'(a1)',end=300) H
       nstars=nstars+1
c       print *,ID(i),xc(i),yc(i),xe(i),ye(i),(cmag(j,i),j=1,nfilter)

      enddo
300   close(31)
      print *,nstars
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

