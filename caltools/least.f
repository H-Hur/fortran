      subroutine cal_least_xy_dxdy_exc(ntype,
     &x,x1,x2,
     &y,y1,y2,
     &w,n,  
     &thre,nmin,maxiter,num_cal,ncal,
     &s_in,s1,s2,ds,z_in,z1,z2,dz,
     &s_out,z_out,wsum,std,sstd,zstd)
      integer ntype,n,nc
      integer num_cal0(n),ncal,nmin,maxiter,niter,ncal0,num_cal(n)
      real x(n),y(n),w(n),y1,y2,s_in,s1,s2,ds,z_in,z1,z2,dz,d_s,d_z
      real s_out,z_out,s,z,thre,wsum,sstd,zstd,xcal(n),ycal(n),wcal(n)
      print *,thre,nmin,maxiter,ntype,s_in,z_in
      niter=0
      ncal=0
100   niter=niter+1
      if(niter.eq.1) then
       do i=1,n
        num_cal(i)=0
        if(x(i).ge.x1.and.x(i).le.x2) then
         num_cal0(i)=1
         num_cal(i)=1
         ncal=ncal+1
         xcal(ncal)=x(i)
         ycal(ncal)=y(i)
         wcal(ncal)=w(i)
        endif
       enddo
      endif
      call cal_least_xy_dxdy(ntype,xcal,x1,x2,ycal,y1,y2,wcal,ncal,
     &s_in,s1,s2,ds,z_in,z1,z2,dz,
     &s_out,z_out)
      std=0.
      sstd=0.
      zstd=0.
      wsum=0.
      do i=1,n
       if(x(i).ge.x1.and.x(i).le.x2.and.num_cal(i).ne.0) then
        d_s=(y(i)-z_out)/x(i)-s_out
        d_z=y(i)-s_out*x(i)-z_out
        ddd=abs(s_out*x(i)-1.*y(i)+z_out)/sqrt(s_out*s_out+1.*1.)
        sstd=sstd+w(i)*d_s*d_s
        zstd=zstd+w(i)*d_z*d_z
        std=std+w(i)*ddd*ddd
        wsum=wsum+w(i)
       endif
      enddo
      std=sqrt(std/(real(ncal-1)*wsum/real(ncal)))
      sstd=sqrt(sstd/(real(ncal-1)*wsum/real(ncal)))
      zstd=sqrt(zstd/(real(ncal-1)*wsum/real(ncal)))
      if(niter.eq.maxiter) goto 900
      ncal0=ncal
      ncal=0
      nc=0
      do i=1,n
       num_cal0(i)=num_cal(i)
       num_cal(i)=0
       if(x(i).ge.x1.and.x(i).le.x2) then
        ddd=abs(s_out*x(i)-1.*y(i)+z_out)/sqrt(s_out*s_out+1.*1.)
        if(ddd.le.thre*std) then
         num_cal(i)=1
         ncal=ncal+1
         xcal(ncal)=x(i)
         ycal(ncal)=y(i)
         wcal(ncal)=w(i)
        endif
       endif
       if(num_cal(i).eq.num_cal0(i)) nc=nc+1
      enddo
      if(nc.ne.n.and.niter.lt.maxiter) goto 100
    
900   print *,'N iter= ',niter
      return
      end
!     ntype=1 -> compute slope,y
!     2 -> compute slope with coordinate (s_in,y_in)
cccc!     3 -> compute y with slope input
      subroutine cal_least_xy_dxdy(ntype,
     &x,x1,x2,
     &y,y1,y2,
     &w,n,
     &s_in,s1,s2,ds,z_in,z1,z2,dz,
     &s_out,z_out)
      integer ntype,n
      real x(n),y(n),w(n),y1,y2,s_in,s1,s2,ds,z_in,z1,z2,dz
      real s_out,z_out,wsum,s,z,dmin,d,dsum
      wsum=0.
      s_out=0.
      z_out=0. 
      if(ntype.eq.1) then
       do s=s1,s2,ds
        do z=z1,z2,dz
         dsum=0.
         wsum=0.
         do i=1,n
          if(x(i).ge.x1.and.x(i).le.x2.and.y(i).ge.y1.and.y(i).le.y2) 
     &    then
           d=w(i)*abs(s*x(i)-1.*y(i)+z)/sqrt(s*s+1.*1.)
           dsum=dsum+d
           wsum=wsumn+w(i)
          endif
         enddo
         if(wsum.eq.0.) then
          print *,'cal_least_xy_dxdy : error, weight sum = 0'
         endif
         if((s.eq.s1.and.z.eq.z1).or.(dsum/wsum.lt.dmin)) then
          dmin=dsum/wsum
          s_out=s
          z_out=z
         endif
        enddo
       enddo  
      elseif(ntype.eq.2) then
       do s=s1,s2,ds
        z=z_in-s*s_in
        dsum=0.
        wsum=0.
        do i=1,n
         if(x(i).ge.x1.and.x(i).le.x2.and.y(i).ge.y1.and.y(i).le.y2)
     &   then
          d=w(i)*abs(s*x(i)-1.*y(i)+z)/sqrt(s*s+1.*1.)
          dsum=dsum+d
          wsum=wsumn+w(i)
         endif
        enddo
        if(wsum.eq.0.) then
         print *,'cal_least_xy_dxdy : error, weight sum = 0'
        endif
        if((s.eq.s1).or.(dsum/wsum.lt.dmin)) then
         dmin=dsum/wsum
         s_out=s
         z_out=z
        endif
       enddo
      endif
c      print *,s_in,z_in,s,z,s_out,z_out
      return
      end

!     ntype=1 -> compute slope,y
!     2 -> compute slope with (0,0)
!     3 -> compute slope with y input
!     4 -> compute slope with coordinate (slopein,yin)
!     5 -> compute y with slope input
      subroutine cal_least_dp(x,y,n,ntype,slopein,yin,slopeout,yout)
      real*8 x(n),y(n),sx,sy,sxx,sxy,slopein,yin,slopeout,yout
      real*8 sxslopein
      integer n,ntype
      sx=0.d0
      sy=0.d0
      sxx=0.d0
      sxy=0.d0
      sxslopein=0.d0
      do i=1,n
       sx=sx+x(i)
       sxx=sxx+x(i)**2.d0
       sxy=sxy+x(i)*y(i)
       sy=sy+y(i)
       sxslopein=sxslopein+(x(i)-slopein)**2.d0
      enddo
      if(ntype.eq.1) then
       slopeout=(real(n)*sxy-sx*sy)/(real(n)*sxx-sx**2.d0)
       yout=(sxx*sy-sx*sxy)/(real(n)*sxx-sx**2.d0)
      elseif(ntype.eq.2) then
       slopeout=sy/sx
       yout=0.d0
      elseif(ntype.eq.3) then
       slopeout=(sxy-yin*sx)/sxx
       yout=yin
      elseif(ntype.eq.4) then
       slopeout=(sxy-slopein*sy-yin*sx+real(n)*slopein*yin)/sxslopein
       yout=yin-slopeout*slopein
      elseif(ntype.eq.5) then
       yout=(sy-slopein*sx)/real(n)
       slopeout=slopein
      endif
      return
      end

!     ntype=1 -> compute slope,y
!     2 -> compute slope with coordinate (x,y)=(s_in,y_in)
!     3 -> compute y(yo) with slope(s_in) input
      subroutine cal_wleast(x,y,w,n,ntype,s_in,y_in,so,yo,es,ey)
      real x(n),y(n),w(n),wx,wy,wxx,wxy,ws,s_in,y_in,so,yo,wdd,es,ey
      integer n,i,ntype
      wx=0.
      wy=0.
      wxx=0.
      wxy=0.
      ws=0.
      wdd=0.
      do i=1,n
       ws    =  ws + w(i)
       wx    =  wx + w(i) * x(i)
       wy    =  wy + w(i) * y(i)
       wxx   =  wxx + w(i) * x(i) * x(i)
       wxy   =  wxy + w(i) * x(i) * y(i)
      enddo
      if(ntype.eq.1) then
       so = (ws*wxy-wx*wy)/(ws*wxx-wx*wx)
       yo = (wxx*wy-wx*wxy)/(ws*wxx-wx*wx)
      elseif(ntype.eq.2) then
       so=0.
       do i=1,n
        so=so+w(i)*(y(i)-y_in)/(x(i)-s_in)
       enddo
       so=so/ws
       yo=0.
       do i=1,n
        yo=yo+w(i)*(-so*x(i)+y(i))
       enddo
       yo=yo/ws
      elseif(ntype.eq.3) then
       yo=0.
       do i=1,n
        yo=yo+w(i)*(-s_in*x(i)+y(i))
       enddo
       yo=yo/ws
       so=s_in
      endif
      do i=1,n
       wdd =  wdd+w(i)*(y(i)-yo-x(i)*so)*(y(i)-yo-x(i)*so)
      enddo
      es=sqrt(ws*wdd/((ws*wxx-wx*wx)*real(n-2)))
      ey=sqrt(wxx*wdd/((ws*wxx-wx*wx)*real(n-2)))
      return
      end


!     ntype=1 -> compute slope,y
!     2 -> compute slope with coordinate (slopein,yin)
!     3 -> compute y with slope input
      subroutine cal_wleast_eq(x,y,w,n,x1,x2,ntype,
     &sl_in,y_in,sl_out,yout,esl,ey,n_ex)
      parameter(ns=100000)
      real x(n),y(n),w(n),sl_in,y_in,xx(ns),yy(ns),ww(ns),x1,x2
      real*8 wx,wy,wxx,wxy,wsum,sl_out,yout
      real*8  wdydy,esl,ey,wxslin,st
      integer n,i,nn,n_ex_last,n_ex
      do i=1,n
       xx(i)=x(i)
       yy(i)=y(i)
       ww(i)=w(i)
      enddo
      nn=n
      n_ex=0
100   wx=0.d0
      wy=0.d0
      wxx=0.d0
      wxy=0.d0
      wsum=0.d0
      wdydy=0.d0
      wxslin=0.d0
      n_ex_last=n_ex
      do i=1,nn
       if(xx(i).ge.x1.and.xx(i).lt.x2) then
        wsum  =  wsum + ww(i)
        wx    =  wx + ww(i) * xx(i)
        wy    =  wy + ww(i) * yy(i)
        wxx   =  wxx + ww(i) * xx(i) * xx(i)
        wxy   =  wxy + ww(i) * xx(i) * yy(i)
       endif
      enddo
      if(ntype.eq.1) then
       sl_out = (wsum*wxy-wx*wy)/(wsum*wxx-wx*wx)
       yout    = (wxx*wy-wx*wxy)/(wsum*wxx-wx*wx)
      elseif(ntype.eq.2) then
       sl_out=0.d0
       do i=1,nn
        if(xx(i).ge.x1.and.xx(i).lt.x2)
     &  sl_out=ww(i)*(yy(i)-yin)/(xx(i)-sl_in)
       enddo
       sl_out=sl_out/wsum
       yout=0.d0
       do i=1,nn
        if(xx(i).ge.x1.and.xx(i).lt.x2) yout=ww(i)*(sl_out*xx(i)+yy(i))
       enddo
       yout=yout/wsum
      elseif(ntype.eq.3) then
       yout=0.d0
       do i=1,nn
        if(xx(i).ge.x1.and.xx(i).lt.x2) yout=ww(i)*(sl_in*xx(i)+yy(i))
       enddo
       yout=yout/wsum
       sl_out=sl_in
      endif
      st=0.d0
      do i=1,nn
       if(xx(i).ge.x1.and.xx(i).lt.x2) then
       st=st+(yy(i)-yout-xx(i)*sl_out)**2
       wdydy =  (yy(i)-yout-xx(i)*sl_out)*(yy(i)-yout-xx(i)*sl_out)
     & /(wsum-2.d0*wsum/(real(nn)))
       endif
      enddo
      st=sqrt(st/real(nn-1))
      nn=0
      do i=1,n
       if(abs(y(i)-yout-x(i)*sl_out).le.st*2.5d0) then
        nn=nn+1
        xx(nn)=x(i)
        yy(nn)=y(i)
        ww(nn)=w(i)
       endif
      enddo
      n_ex=n-nn
      if(n_ex.ne.n_ex_last) goto 100
      esl=sqrt(wsum*wdydy/(wsum*wxx-wx*wx))
      ey=sqrt((wxx*wdydy)/(wsum*wxx-wx*wx))
      return
      end

      subroutine least_xy(xl,yl,n,slope,yp)
      real xl(n),yl(n),slope,yp
      integer n
      real xx,yy,xy,x,y,a,b,c,d,xc,yc,xyc,dslope,sq,rt
      xx=0.
      yy=0.
      xy=0.
      x=0.
      y=0.
      do i=1,n
       xx=xx+xl(i)*xl(i)
       yy=yy+yl(i)*yl(i)
       xy=xy+xl(i)*yl(i)
       x=x+xl(i)
       y=y+yl(i)
      enddo
      xc=real(n)*xx - x*x
      yc=real(n)*yy - y*y
      xyc = real(n)*xy - x*y  
      a = 2.*xc
      b = 2.*xc*yc - xc
      c = x - 2.*xc*yc
      d = 2.*xc*yc
  
      rt = 2.*b*b*b - 9.*a*b*c + 27.*a*a*d  
      sq = sqrt(rt*rt - 4.*(b*b-3.*a*c)**3) 
      dslope = (-1./(3.*a))*
     &(b + ((rt+sq)/2.)**(1./3.) + ((rt-sq)/2.)**(1./3.))

      slope = real(dslope)
      yp = - (dslope*x + y)/real(n)
      return
      end
