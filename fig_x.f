c for only fortran77.
c f77 -o 'fig.exe' 'fig_x.f' -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11

      parameter (nf=10,ns=1000)
      character filter(nf)*2,star(nf,ns)*16,input*16,coeff*16,output*22
      real*8 transy(nf,ns),airy(nf,ns),timey(nf,ns),colx(nf,ns)
      real*8 timex(nf,ns)
      integer num_f(nf),num_c(nf),ntransform,lastnum_f,ne,nne

      print *,'Input data file name=?'
      read(*,'(a16)') input
c      write(input,'(a16)') '228_9703.fig    '
      call readhead(input)
      call readmag(input,filter,star)

c      write(coeff,'(a16)') 'coeff.228.fig   '
      print *,'Input coeff data file name=?'
      read(*,'(a16)') coeff
      call readcoeff(coeff,num_f,num_c,ntransform)

      write(input,'(a16)') 'Sung2001.ZAMS.cc'
      call readhead(input)
      call read_ZAMS(input)




      call caly(ntransform,filter,star,transy,airy,timey,colx,timex)
      call dataran()


      do i=1,ntransform

       if(num_f(i).eq.1) then
        write(output,'(a21)') 'postportfile fig_U.ps'
        elseif(num_f(i).eq.2) then
        write(output,'(a21)') 'postportfile fig_B.ps'
        elseif(num_f(i).eq.3) then
        write(output,'(a21)') 'postportfile fig_V.ps'
        elseif(num_f(i).eq.4) then
        write(output,'(a21)') 'postportfile fig_R.ps'
        elseif(num_f(i).eq.5) then
        write(output,'(a21)') 'postportfile fig_I.ps'
        elseif(num_f(i).eq.6) then
        write(output,'(a22)') 'postportfile fig_Ha.ps'
       endif

       if(i.le.ntransform-1) then
        if(num_f(i).eq.num_f(i+1)) then
         nne=1
         else
         nne=0
        endif
       endif

       if(num_f(i).ne.lastnum_f) then
        call sm_device (output)
        call sm_graphics
        call sm_erase()
        ne=0
        else
        ne=1
       endif

       if(i.eq.1) then
        lastnum_f=num_f(i)
        call sm_device (output)
        call sm_graphics
        call sm_erase()
        ne=0
        elseif(ne.eq.0) then
        lastnum_f=num_f(i)
        call sm_device (output)
        call sm_graphics
        call sm_erase()
       endif
       lastnum_f=num_f(i)

       call sm_location(5000,31000,26000,31000)    
       call Xdraw(filter,airy,num_f(i),num_c(i),i,ne)

       call sm_location(5000,31000,18000,23000)
       call timedraw(star,filter,timey,timex,num_f(i),num_c(i),i,ne)

       call sm_location(5000,31000,10000,15000)
       call transdraw(star,filter,transy,colx,num_f(i),num_c(i),i,ne)

       if(nne.ne.1) then
        call sm_gflush()
        call sm_hardcopy
        call sm_alpha
        elseif(i.eq.ntransform) then
        call sm_gflush()
        call sm_hardcopy
        call sm_alpha
       endif        

      enddo

      stop
      end



      subroutine dataran()
      parameter(nf=10,ns=1000)
      character input*16
      integer nstars
      common /ran/BVran(2),UBran(2),VRran(2),VIran(2),UTran(2),Airran(2)
      common/dat1/dmag(nf,ns),dmerr(nf,ns),ut(nf,ns),air(nf,ns),
     &weight(nf,ns),nall
      common/dat2/V(nf,ns),BV(nf,ns),VR(nf,ns),VI(nf,ns),UB(nf,ns),
     &nstar(nf),nfilter,mfilter(nf)
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),coeff1(nf),
     &coeff2(nf),zeropoint(nf),delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf),Lfilter(nf)
      nstars=0
      do i=1,nfilter
       nstars=nstars+nstar(i)
      enddo

      BVran(1)=1.
      UBran(1)=1.
      VRran(1)=1.
      VIran(1)=1.
      UTran(1)=10.
      Airran(1)=0.95
      BVran(2)=0.
      UBran(2)=0.
      VRran(2)=0.
      VIran(2)=0.
      UTran(2)=-10.
      Airran(2)=2.3

      do j=1,nfilter
       do i=1,nstar(j)
        if(BV(j,i).lt.BVran(1)) then
         BVran(1)=BV(j,i)
         elseif(BV(j,i).gt.BVran(2)) then
         BVran(2)=BV(j,i)
        endif
        if(UB(j,i).lt.UBran(1)) then
         UBran(1)=UB(j,i)
         elseif(UB(j,i).gt.UBran(2)) then
         UBran(2)=UB(j,i)
        endif
        if(VR(j,i).lt.VRran(1)) then
         VRran(1)=VR(j,i)
         elseif(VR(j,i).gt.VRran(2)) then
         VRran(2)=VR(j,i)
        endif
        if(VI(j,i).lt.VIran(1)) then
         VIran(1)=VI(j,i)
         elseif(VI(j,i).gt.VIran(2)) then
         VIran(2)=VI(j,i)
        endif
        if(UT(j,i)-delta_ut(j).lt.UTran(1)) then
         UTran(1)=UT(j,i)-delta_ut(j)
         elseif(UT(j,i)-delta_ut(j).gt.UTran(2)) then
         UTran(2)=UT(j,i)-delta_ut(j)
        endif
       enddo
      enddo

c      BVran(1)=BVran(1)-0.1
c      UBran(1)=UBran(1)-0.1
c      VRran(1)=VRran(1)-0.1
c      VIran(1)=VIran(1)-0.1
      UTran(1)=UTran(1)-1.0
c      BVran(2)=BVran(2)+0.1
c      UBran(2)=UBran(2)+0.1
c      VRran(2)=VRran(2)+0.1
c      VIran(2)=VIran(2)+0.1
      UTran(2)=UTran(2)+1.0

      BVran(1)=-0.5
      UBran(1)=-1.2
      VRran(1)=-0.5
      VIran(1)=-0.5
      BVran(2)=2.
      UBran(2)=1.5
      VRran(2)=1.2
      VIran(2)=2.5

      return
      end 

      



      subroutine transdraw(star,filter,transy,colx,numfil,numcol,numtr,
     &ne)
      parameter (nf=10,ns=1000)
      real*8 transy(nf,ns),colx(nf,ns),wei_mean
      real c_bv(21),c_u(21),fbv,xmin,xmax,xran,ymin,ymax,yran
      common /ran/BVran(2),UBran(2),VRran(2),VIran(2),UTran(2),Airran(2)
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),coeff1(nf),
     &coeff2(nf),zeropoint(nf),delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf),Lfilter(nf)
      common/dat1/dmag(nf,ns),dmerr(nf,ns),ut(nf,ns),air(nf,ns),
     &weight(nf,ns),nall
      common/dat2/V(nf,ns),BV(nf,ns),VR(nf,ns),VI(nf,ns),UB(nf,ns),
     &nstar(nf),nfilter,mfilter(nf)
      real x(ns),y(ns),ylocate,xlocate,pt1,pt2,pt3,pp,w(ns)
      integer npoint,numfil,numcol,numtr,ne
      character filter(nf)*2,label*20,star(nf,ns)*16
      npoint=0
      wei_mean=0.d0
      do i=1,nstar(mfilter(numfil))
       npoint=npoint+1
       y(i)=transy(mfilter(numfil),i)
       x(i)=colx(mfilter(numfil),i)
       w(i)=weight(mfilter(numfil),i)
       wei_mean=wei_mean+w(i)
      enddo
      wei_mean=wei_mean/real(npoint)
      
      if(numcol.eq.1) then
       xmin=UBran(1)
       xmax=UBran(2)
       elseif(numcol.eq.2) then
       xmin=BVran(1)
       xmax=BVran(2)
       elseif(numcol.eq.3) then
       xmin=VRran(1)
       xmax=VRran(2)
       elseif(numcol.eq.5) then
       xmin=VIran(1)
       xmax=VIran(2)
      endif
      xran=xmax-xmin
      ymin=-0.26
      ymax=0.26
      yran=ymax-ymin 

      call sm_limits(xmin,xmax,-0.26,0.26)
      call sm_ctype('black')
      call sm_expand(0.5)
      if(ne.eq.0) then
       call sm_box(1,2,1,2)
      endif
      call sm_expand(1.)

      if(numcol.eq.1) then
       write(label,'(a4)') 'U-B '
       elseif(numcol.eq.2) then
       write(label,'(a4)') 'B-V '
       elseif(numcol.eq.3) then
       write(label,'(a4)') 'V-R '
       elseif(numcol.eq.4) then
       write(label,'(a4)') 'R-I '
       elseif(numcol.eq.5) then
       write(label,'(a4)') 'V-I '
       elseif(numcol.eq.6) then
       write(label,'(a4)') 'R-Ha'
      endif
      if(ne.eq.0) then
       call sm_xlabel(label)
      endif

      call sm_expand(0.5)
      if(numfil.eq.1.and.ne.eq.0) then
       call sm_ylabel('M-m+k\\d1X-k\\d2CX-\\ga(UT)-\\gx-f[(B-V)\\d0]')
       elseif(ne.eq.0) then
       call sm_ylabel('M-m+k\\d1X-k\\d2CX-\\ga(UT)-\\gx')
      endif

      ylocate=ymin+yran*0.85
      xlocate=xmin
      call sm_relocate(xlocate,ylocate)
10    format(a7,a2,a1,f7.4)
      call sm_expand(0.7)
      write(label,10) '\\2\\gh',filter(mfilter(numfil)),'=',
     &transcoeff(numtr)
      if(colorran(numtr,2)-colorran(numtr,1) .gt.5.) then
       call sm_label(label)
      endif
      call sm_expand(1.)

      if(xmin.lt.colorran(numtr,1)) then
       xmin=colorran(numtr,1)
      endif
      if(xmax.gt.colorran(numtr,2)) then
       xmax=colorran(numtr,2)
      endif
      xlocate=xmin
      ylocate=xlocate*transcoeff(numtr)+transzero(numtr)
      call sm_relocate(xlocate,ylocate)
      xlocate=xmax
      ylocate=xlocate*transcoeff(numtr)+transzero(numtr)
      call sm_draw(xlocate,ylocate)

      pt1=400.
      pt2=40.
      pt3=1.01
      do i=1,npoint
       if(x(i).gt.xmin.and.x(i).le.xmax) then
        call sm_expand(2.5)
        pp=pt1+w(i)/pt3
        call sm_ptype(pp,1)
        call sm_points(x(i),y(i),1)
        if(w(i).ge.wei_mean*1.3d0) then
         write(label,'(a20)') star(mfilter(numfil),i)
         call sm_expand(0.4)
         call sm_label(label)
        endif 
       endif
      enddo

      call sm_expand(1.)

      return
      end


      subroutine timedraw(star,filter,timey,timex,numfil,numcol,numtr,
     &ne)
      parameter (nf=10,ns=1000)
      real c_bv(21),c_u(21),fbv,xran,yran,ymin,ymax,xmin,xmax
      real*8 timey(nf,ns),timex(nf,ns),wei_mean
      common /ran/BVran(2),UBran(2),VRran(2),VIran(2),UTran(2),Airran(2)
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),coeff1(nf),
     &coeff2(nf),zeropoint(nf),delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf),Lfilter(nf)
      common/dat1/dmag(nf,ns),dmerr(nf,ns),ut(nf,ns),air(nf,ns),
     &weight(nf,ns),nall
      common/dat2/V(nf,ns),BV(nf,ns),VR(nf,ns),VI(nf,ns),UB(nf,ns),
     &nstar(nf),nfilter,mfilter(nf)
      real x(ns),y(ns),ylocate,xlocate,pt1,pt2,pt3,pp,w(ns)
      integer npoint,numfil,numcol,numtr,ne
      character filter(nf)*2,label*20,star(nf,ns)*16
      npoint=0
      wei_mean=0.d0
      do i=1,nstar(mfilter(numfil))
       npoint=npoint+1
       y(i)=timey(mfilter(numfil),i)
       x(i)=timex(mfilter(numfil),i)
       w(i)=weight(mfilter(numfil),i)
       wei_mean=wei_mean+w(i)
      enddo
      wei_mean=wei_mean/real(npoint)
      ymin=-0.11
      ymax=0.11
      xmin=UTran(1)
      xmax=UTran(2)
      call sm_limits(xmin,xmax,ymin,ymax)
      call sm_ctype('black')
      call sm_expand(0.5)
      call sm_box(1,2,1,2)
      call sm_expand(1.)

      write(label,*) 'UT -',delta_ut(numtr)
      call sm_xlabel(label)
      call sm_expand(0.5)
      if(numfil.eq.1) then
       call sm_ylabel('M-m+k\\d1X-k\\d2CX-\\ghC-\\gx-f[(B-V)\\d0]')
       else
       call sm_ylabel('M-m+k\\d1X-k\\d2CX-\\ghC-\\gx')
      endif
      yran=ymax-ymin
      ylocate=ymin+yran*0.85
      xran=xmax-xmin
      xlocate=xmin
      call sm_relocate(xlocate,ylocate)
10    format(a10,a2,a1,f6.4)
      call sm_expand(0.7)
      write(label,10) '\\2\\2\\ga',filter(mfilter(numfil)),'=',
     &timevar(numtr)
      if(timeran(numtr,2)-timeran(numtr,1) .gt.12.) then
       call sm_label(label)
      endif

      xlocate=UTran(1)
      ylocate=UTran(1)*timevar(numtr)
      call sm_relocate(xlocate,ylocate)
      xlocate=UTran(2)
      ylocate=UTran(2)*timevar(numtr)
      call sm_draw(xlocate,ylocate)
      call sm_expand(2.5)
      pt1=400.
      pt2=40.
      pt3=1.01
      do i=1,npoint
       call sm_expand(2.5)
       pp=pt1+w(i)/pt3
       call sm_ptype(pp,1)
       call sm_points(x(i),y(i),1)
       if(w(i).ge.wei_mean*1.3d0) then
        write(label,'(a20)') star(mfilter(numfil),i)
        call sm_expand(0.4)
        call sm_label(label)
       endif
      enddo
      call sm_expand(1.)
      return
      end




      subroutine Xdraw(filter,airy,numfil,numcol,numtr,ne)
      parameter (nf=10,ns=1000)
      real*8 airy(nf,ns)
      real c_bv(21),c_u(21),fbv,xmin,xmax,xran,ymin,ymax,yran
      common /ran/BVran(2),UBran(2),VRran(2),VIran(2),UTran(2),Airran(2)
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),coeff1(nf),
     &coeff2(nf),zeropoint(nf),delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf),Lfilter(nf)
      common/dat1/dmag(nf,ns),dmerr(nf,ns),ut(nf,ns),air(nf,ns),
     &weight(nf,ns),nall
      common/dat2/V(nf,ns),BV(nf,ns),VR(nf,ns),VI(nf,ns),UB(nf,ns),
     &nstar(nf),nfilter,mfilter(nf)
      real x(ns),y(ns),ylocate,xlocate,pt1,pt2,pt3,pp,w(ns)
      integer npoint,numfil,numcol,numtr,ne
      character filter(nf)*2,label*20
      npoint=0
      do i=1,nstar(mfilter(numfil))
       npoint=npoint+1
       y(i)=airy(mfilter(numfil),i)
       x(i)=air(mfilter(numfil),i)
       w(i)=weight(mfilter(numfil),i)
       if(numfil.eq.5) then
       endif
      enddo

      if(numfil.eq.1) then
       ymin=-1.3
       elseif(numfil.eq.2) then
       ymin=-0.8
       elseif(numfil.eq.3) then
       ymin=-0.45
       elseif(numfil.eq.4) then
       ymin=-0.4
       elseif(numfil.eq.5) then
       ymin=-0.25
      endif
      ymax=0.1
      yran=ymax-ymin
      xmin=Airran(1)
      xmax=Airran(2)
      xran=xmax-xmin
      call sm_limits(xmin,xmax,ymin,ymax)
      call sm_ctype('black')
      call sm_expand(0.5)
      call sm_box(1,2,1,2)
      call sm_expand(1.)
      call sm_xlabel('Airmass')
      call sm_expand(0.5)
      if(numfil.eq.1) then
       call sm_ylabel('M-m-\\ghC-k\\d2CX-\\gx-\\ga(UT)-f[(B-V)\\d0]')
       else
       call sm_ylabel('M-m-\\ghC-k\\d2CX-\\gx-\\ga(UT)')
      endif

      ylocate=ymin+yran*0.85
      xlocate=xmin
      call sm_relocate(xlocate,ylocate)
10    format(a11,a2,a1,f6.4)
      call sm_expand(0.7)
      write(label,10) '\\2k\\2\\d1',filter(mfilter(numfil)),'=',
     &coeff1(numtr)
      call sm_label(label)
      call sm_expand(1.)
      ylocate=-2.3*coeff1(numtr)
      call sm_relocate(2.3,ylocate)
      call sm_draw(0.,0.)
      call sm_expand(2.5)

      pt1=400.
      pt2=40.
      pt3=1.01
      do i=1,npoint
       pp=pt1+w(i)/pt3
       call sm_ptype(pp,1)
       call sm_points(x(i),y(i),1)
      enddo
      call sm_expand(1.)
      return
      end



      subroutine caly(ntrans,filter,star,transy,airy,timey,colx,timex)
      parameter (nf=10,ns=1000)
      integer i1,ntrans
      character filter(nf)*2,star(nf,ns)*16
      real*8 transy(nf,ns),airy(nf,ns),timey(nf,ns),colx(nf,ns)
      real*8 timex(nf,ns)
      real c_bv(21),c_u(21),fbv,bv0,zamsUB,ub0
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),coeff1(nf),
     &coeff2(nf),zeropoint(nf),delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf),Lfilter(nf)
      common/dat1/dmag(nf,ns),dmerr(nf,ns),ut(nf,ns),air(nf,ns),
     &weight(nf,ns),nall
      common/dat2/V(nf,ns),BV(nf,ns),VR(nf,ns),VI(nf,ns),UB(nf,ns),
     &nstar(nf),nfilter,mfilter(nf)
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),zVR(ns),zRI(ns),nzams
      data c_bv/-.35,-.33,-0.3,-0.2,-0.1, 0.0, 0.1, 0.2,
     &0.3, 0.4,0.5,0.6,0.7,0.8,0.9,1.0,1.2,1.4,1.6
     &,1.8,2.0/
      data c_u/-.132,-.132,-.126,-.074,-.030,0.002,0.024,
     &0.026,0.010,-.017,-.06,-.104,-.132,-.132,-.132,
     &-.132,-.132,-.132,-.132,-.132,-.132/

      do i=1,nfilter
       if(i.eq.mfilter(1)) then!U
        do j=1,nstar(i)
         do k=1,ntrans
          if(ut(i,j).ge.timeran(k,1).and.ut(i,j).lt.timeran(k,2).and.
     &    num_fil(k).eq.1)then
           if(UB(i,j).ge.colorran(k,1).and.UB(i,j).lt.colorran(k,2))
     &     then
            call lint(UB(i,j),zUB,zBV,bv0,nzams,nzams)
            if(BV(i,j).le.0.and.bv0.lt.BV(i,j)) then
             call cor_red(BV(i,j),UB(i,j),bv0,ub0)
             print *,star(i,j),'E(B-V)=',BV(i,j)-bv0
            else
             bv0=BV(i,j)
            endif
            call lint(bv0,c_bv,c_u,fbv,21,21)
            airy(i,j)=V(i,j)+BV(i,j)+UB(i,j)-dmag(i,j)-transcoeff(k)
     &      *UB(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)-
     &      coeff2(k)*UB(i,j)*air(i,j)-fbv-transzero(k)-timezero(k)
            timey(i,j)=V(i,j)+BV(i,j)+UB(i,j)-dmag(i,j)
     &      -transcoeff(k)*UB(i,j)+coeff1(k)*air(i,j)-zeropoint(k)-
     &      coeff2(k)*UB(i,j)*air(i,j)-fbv-transzero(k)
            transy(i,j)=V(i,j)+BV(i,j)+UB(i,j)-dmag(i,j)+coeff1(k)*
     &      air(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)-
     &      coeff2(k)*UB(i,j)*air(i,j)-fbv-timezero(k)
            colx(i,j)=UB(i,j)
            timex(i,j)=ut(i,j)-delta_ut(k)
c            print *,airy(i,j)
           endif
          endif
         enddo
        enddo

       elseif(i.eq.mfilter(2)) then!B
        do j=1,nstar(i)
         do k=1,ntrans
          if(ut(i,j).ge.timeran(k,1).and.ut(i,j).lt.timeran(k,2).and.
     &    num_fil(k).eq.2) then
           if(BV(i,j).ge.colorran(k,1).and.BV(i,j).lt.colorran(k,2))
     &     then
            airy(i,j)=V(i,j)+BV(i,j)-dmag(i,j)-transcoeff(k)*
     &      BV(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &      -coeff2(k)*BV(i,j)*air(i,j)-transzero(k)-timezero(k)
            timey(i,j)=V(i,j)+BV(i,j)-dmag(i,j)-transcoeff(k)*
     &      BV(i,j)+coeff1(k)*air(i,j)-zeropoint(k)-transzero(k)-
     &      coeff2(k)*BV(i,j)*air(i,j)
            transy(i,j)=V(i,j)+BV(i,j)-dmag(i,j)+coeff1(k)*
     &      air(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &      -coeff2(k)*BV(i,j)*air(i,j)-timezero(k)
            colx(i,j)=BV(i,j)
            timex(i,j)=ut(i,j)-delta_ut(k)
           endif
          endif
         enddo
        enddo

       elseif(i.eq.mfilter(3)) then!V
        do j=1,nstar(i)
         do k=1,ntrans
          if(ut(i,j).ge.timeran(k,1).and.ut(i,j).lt.timeran(k,2).and.
     &    num_fil(k).eq.3) then
           if(num_col(k).eq.2) then!B-V
            if(BV(i,j).ge.colorran(k,1).and.BV(i,j).lt.colorran(k,2))
     &      then
             airy(i,j)=V(i,j)-dmag(i,j)-transcoeff(k)*
     &       BV(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &       -coeff2(k)*BV(i,j)*air(i,j)-transzero(k)-timezero(k)
             timey(i,j)=V(i,j)-dmag(i,j)-transcoeff(k)*
     &       BV(i,j)+coeff1(k)*air(i,j)-zeropoint(k)
     &       -coeff2(k)*BV(i,j)*air(i,j)-transzero(k)
             transy(i,j)=V(i,j)-dmag(i,j)+coeff1(k)*
     &       air(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &       -coeff2(k)*BV(i,j)*air(i,j)-timezero(k)
             colx(i,j)=BV(i,j)
             timex(i,j)=ut(i,j)-delta_ut(k)
            endif
           elseif(num_col(k).eq.5) then!V-I
            if(VI(i,j).ge.colorran(k,1).and.VI(i,j).lt.colorran(k,2))
     &      then
             airy(i,j)=V(i,j)-dmag(i,j)-transcoeff(k)*
     &       VI(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &       -coeff2(k)*VI(i,j)*air(i,j)-transzero(k)-timezero(k)
             timey(i,j)=V(i,j)-dmag(i,j)-transcoeff(k)*
     &       VI(i,j)+coeff1(k)*air(i,j)-zeropoint(k)
     &       -coeff2(k)*VI(i,j)*air(i,j)-transzero(k)
             transy(i,j)=V(i,j)-dmag(i,j)+coeff1(k)*
     &       air(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &       -coeff2(k)*VI(i,j)*air(i,j)-timezero(k)
             colx(i,j)=VI(i,j)
             timex(i,j)=ut(i,j)-delta_ut(k)
            endif
           endif
          endif
         enddo
        enddo

       elseif(i.eq.mfilter(4)) then!R
        do j=1,nstar(i)
         do k=1,ntrans
          if(ut(i,j).ge.timeran(k,1).and.ut(i,j).lt.timeran(k,2).and.
     &    num_fil(k).eq.4) then
           if(VR(i,j).ge.colorran(k,1).and.VR(i,j).lt.colorran(k,2))
     &     then
            airy(i,j)=V(i,j)-VR(i,j)-dmag(i,j)-transcoeff(k)*
     &      VR(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &      -coeff2(k)*VR(i,j)*air(i,j)-transzero(k)-timezero(k)
            timey(i,j)=V(i,j)-VR(i,j)-dmag(i,j)-transcoeff(k)*
     &      VR(i,j)+coeff1(k)*air(i,j)-zeropoint(k)
     &      -coeff2(k)*VR(i,j)*air(i,j)-transzero(k)
            transy(i,j)=V(i,j)-VR(i,j)-dmag(i,j)+coeff1(k)*
     &      air(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &      -coeff2(k)*VR(i,j)*air(i,j)-timezero(k)
            colx(i,j)=VR(i,j)
            timex(i,j)=ut(i,j)-delta_ut(k)
           endif
          endif
         enddo
        enddo

       elseif(i.eq.mfilter(5)) then!I
        do j=1,nstar(i)
         do k=1,ntrans
          if(ut(i,j).ge.timeran(k,1).and.ut(i,j).lt.timeran(k,2).and.
     &    num_fil(k).eq.5) then
           if(VI(i,j).ge.colorran(k,1).and.VI(i,j).lt.colorran(k,2))
     &     then
            airy(i,j)=V(i,j)-VI(i,j)-dmag(i,j)-transcoeff(k)*
     &      VI(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &      -coeff2(k)*VI(i,j)*air(i,j)-transzero(k)-timezero(k)
            timey(i,j)=V(i,j)-VI(i,j)-dmag(i,j)-transcoeff(k)*
     &      VI(i,j)+coeff1(k)*air(i,j)-zeropoint(k)
     &      -coeff2(k)*VI(i,j)*air(i,j)-transzero(k)
            transy(i,j)=V(i,j)-VI(i,j)-dmag(i,j)+coeff1(k)*
     &      air(i,j)-timevar(k)*(ut(i,j)-delta_ut(k))-zeropoint(k)
     &      -coeff2(k)*VI(i,j)*air(i,j)-timezero(k)
            colx(i,j)=VI(i,j)
            timex(i,j)=ut(i,j)-delta_ut(k)
           endif
          endif
         enddo
        enddo

       elseif(i.eq.mfilter(6)) then!Ha
        do j=1,nstar(i)
         do k=1,ntrans
          if(ut(i,j).ge.timeran(k,1).and.ut(i,j).lt.timeran(k,2).and.
     &    num_fil(k).eq.6) then
           airy(i,j)=dmag(i,j)-zeropoint(k)
           timey(i,j)=0.
           transy(i,j)=0.
           colx(i,j)=0.
           timex(i,j)=ut(i,j)-delta_ut(k)
          endif
         enddo
        enddo
       endif

      enddo

      return
      end
 

      subroutine readcoeff(coeff,num_f,num_c,ntransform)
      parameter (nf=10,ns=1000)
      character coeff*16,H*1
      integer ntrans,num_f(nf),num_c(nf),ntransform,numf
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),coeff1(nf),
     &coeff2(nf),zeropoint(nf),delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf),Lfilter(nf)
      common/dat1/dmag(nf,ns),dmerr(nf,ns),ut(nf,ns),air(nf,ns),
     &weight(nf,ns),nall
      common/dat2/V(nf,ns),BV(nf,ns),VR(nf,ns),VI(nf,ns),UB(nf,ns),
     &nstar(nf),nfilter,mfilter(nf)
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
        Lfilter(mfilter(ntrans))=ntrans
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


      subroutine readmag(input,filter,star)
      parameter (nf=10,ns=1000,nmcol=20)
      integer numcol,nc1,nc2,ncol,nc10,nc20
      common /nc/nc_filter,nc_star,nc_mag,nc_merr,nc_UT,nc_air,nc_weight
     &,ncollum,nc_V,nc_BV,nc_VR,nc_RI,nc_VI,nc_UB,nc_Mv
      common/dat1/dmag(nf,ns),dmerr(nf,ns),ut(nf,ns),air(nf,ns),
     &weight(nf,ns),nall
      common/dat2/V(nf,ns),BV(nf,ns),VR(nf,ns),VI(nf,ns),UB(nf,ns),
     &nstar(nf),nfilter,mfilter(nf)
      character line*200,H(200)*1,input*16,coll(nmcol)*32,filter(nf)*2
      character star(nf,ns)*16,Hv(32)*1,lastcoll1*32
      do i=1,32
       write(Hv(i),'(a1)') ' '
      enddo
      open(31,file=input)
      nfilter=0
      write(lastcoll1,'(32a1)') (Hv(k),k=1,32)

      do i=1,ns

       read(31,'(a200)',end=500) line
       read(line,'(200a1)') (H(j),j=1,200)
       if(H(1).eq.'#') then    
        goto 100
       endif
       do j=1,16
        if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &  '^')then
         nc1=1
         nc2=j
         goto 10
        endif
       enddo

10     write(coll(1),'(32a1)') (Hv(k),k=1,32)
       write(coll(1),'(32a1)') (H(k),k=nc1,nc2)
       if(coll(1).ne.lastcoll1) then
        nfilter=nfilter+1
        read(coll(1),'(a2)') filter(nfilter)
        if(filter(nfilter).eq.'U ') then
         mfilter(1)=nfilter
         elseif(filter(nfilter).eq.'B ') then
         mfilter(2)=nfilter
         elseif(filter(nfilter).eq.'V ') then
         mfilter(3)=nfilter
         elseif(filter(nfilter).eq.'R ') then
         mfilter(4)=nfilter
         elseif(filter(nfilter).eq.'I ') then
         mfilter(5)=nfilter
         elseif(filter(nfilter).eq.'Ha') then
         mfilter(6)=nfilter
        endif
        nstar(nfilter)=1
        else
        nstar(nfilter)=nstar(nfilter)+1
       endif
       write(lastcoll1,'(a32)') coll(1)

       do numcol=2,ncollum
        do j=nc2+1,nc2+1+32
         if(H(j).eq.' '.and.H(j+1).ne.' ') then
          nc1=j+1
         endif
         if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &   '^')  then
          nc2=j
          goto 20
         endif        
        enddo
20      write(coll(numcol),'(32a1)') (Hv(k),k=1,32)
        write(coll(numcol),'(32a1)') (H(k),k=nc1,nc2)
       enddo

       read(coll(nc_star),*) star(nfilter,nstar(nfilter))
       read(coll(nc_mag),*) dmag(nfilter,nstar(nfilter))
       read(coll(nc_merr),*) dmerr(nfilter,nstar(nfilter))
       read(coll(nc_UT),*) ut(nfilter,nstar(nfilter))
       read(coll(nc_air),*) air(nfilter,nstar(nfilter))
       read(coll(nc_V),*) V(nfilter,nstar(nfilter))
       read(coll(nc_BV),*) BV(nfilter,nstar(nfilter))
       read(coll(nc_VR),*) VR(nfilter,nstar(nfilter))
       read(coll(nc_VI),*) VI(nfilter,nstar(nfilter))
       read(coll(nc_UB),*) UB(nfilter,nstar(nfilter))
       read(coll(nc_weight),*) weight(nfilter,nstar(nfilter))

c       print *,ut(nfilter,nstar(nfilter))
100   enddo 
500   close(31)

      nall=i-2
      return
      end
    

      subroutine readhead(input)
      parameter (nf=10,ns=1000)
      character header*200,H(200)*1,input*16,coll*3
      common /nc/nc_filter,nc_star,nc_mag,nc_merr,nc_UT,nc_air,nc_weight
     &,ncollum,nc_V,nc_BV,nc_VR,nc_RI,nc_VI,nc_UB,nc_Mv
      open(31,file=input)
      read(31,'(a200)') header
      ncollum=0
      read(header,'(200a1)') (H(i),i=1,200)

      do i=1,200
       if(H(i).eq.' '.or.H(i).eq.'#') then
        if(H(i+1).ne.' ') then
         write(coll,'(3a1)') H(i+1),H(i+2),H(i+3)
         if(coll.eq.'Sta'.or.coll.eq.'sta') then
          ncollum=ncollum+1
          nc_star=ncollum
          elseif(coll.eq.'Mag'.or.coll.eq.'mag') then
          ncollum=ncollum+1
          nc_mag=ncollum
          elseif(coll.eq.'Mer'.or.coll.eq.'mer') then
          ncollum=ncollum+1
          nc_merr=ncollum
          elseif(coll.eq.'V  '.or.coll.eq.'v  ') then
          ncollum=ncollum+1
          nc_V=ncollum
          elseif(coll.eq.'B-V'.or.coll.eq.'b-v') then
          ncollum=ncollum+1
          nc_BV=ncollum
          elseif(coll.eq.'U-B'.or.coll.eq.'u-b') then
          ncollum=ncollum+1
          nc_UB=ncollum
          elseif(coll.eq.'V-R'.or.coll.eq.'v-r') then
          ncollum=ncollum+1
          nc_VR=ncollum
          elseif(coll.eq.'R-I'.or.coll.eq.'r-i') then
          ncollum=ncollum+1
          nc_RI=ncollum
          elseif(coll.eq.'V-I'.or.coll.eq.'v-i') then
          ncollum=ncollum+1
          nc_VI=ncollum
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
          nc_weight=ncollum
          elseif(coll.eq.'fil'.or.coll.eq.'Fil') then
          ncollum=ncollum+1
          nc_filter=ncollum
          elseif(coll.eq.'Mv') then
          ncollum=ncollum+1
          nc_Mv=ncollum
         endif
        endif
       endif
      enddo
      close(31)
      return
      end

      subroutine read_ZAMS(input)
      parameter (nf=10,ns=1000,nmcol=20)
      character line*200,H(200)*1,input*16,coll(nmcol)*32
      character Hv(32)*1,lastcoll1*32
      integer numcol,nc1,nc2,ncol,nc10,nc20
      common /nc/nc_filter,nc_star,nc_mag,nc_merr,nc_UT,nc_air,nc_weight
     &,ncollum,nc_V,nc_BV,nc_VR,nc_RI,nc_VI,nc_UB,nc_Mv
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),zVR(ns),zRI(ns),nzams
      do i=1,32
       write(Hv(i),'(a1)') ' '
      enddo
      open(31,file=input)
      write(lastcoll1,'(32a1)') (Hv(k),k=1,32)
      nzams=0
      do i=1,ns
       read(31,'(a200)',end=500) line
       nzams=nzams+1
       read(line,'(200a1)') (H(j),j=1,200)
       if(H(1).eq.'#') then
        nzams=nzams-1
        goto 100
       endif
       do j=1,16
        if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &  '^')then
         nc2=j
         goto 10
        endif
       enddo

10     do numcol=1,ncollum
        if(numcol.eq.1) then
         nc1=1
         do j=1,nc2+1+32
          if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &    '^')  then
           nc2=j
           goto 20
          endif
         enddo
        else
         do j=nc2+1,nc2+1+32
          if(H(j).eq.' '.and.H(j+1).ne.' ') then
           nc1=j+1
          endif
          if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &    '^')  then
           nc2=j
           goto 20
          endif
         enddo
        endif
20      write(coll(numcol),'(32a1)') (Hv(k),k=1,32)
        write(coll(numcol),'(32a1)') (H(k),k=nc1,nc2)
       enddo

       read(coll(nc_Mv),*) zMv(nzams)
       read(coll(nc_BV),*) zBV(nzams)
c       read(coll(nc_VR),*) zVR(nzams)
c       read(coll(nc_RI),*) zRI(nzams)
       read(coll(nc_VI),*) zVI(nzams)
       read(coll(nc_UB),*) zUB(nzams)

100   enddo
500   close(31)
      
      return
      end


      subroutine cor_red(sBV,sUB,unredx,unredy)
      parameter (nf=10,ns=1000,maxiter=1000)
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),zVR(ns),zRI(ns),nzams
      real fx,fy,dy,unredx,unredy,sBV,sUB,slope,eBV,eUB,Rv
      slope=0.72
      Rv=3.1
      call lint(sUB,zUB,zBV,fx,nzams,nzams)
      dx=fx-sBV
      call lint(sBV,zBV,zUB,fy,nzams,nzams)
      dy=sUB-fy
      unredx=sBV
      unredy=sUB
      do j=1,1000
       call lint(unredy,zUB,zBV,fx,ns,ns)
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

