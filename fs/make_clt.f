c     gfortran -o make_clt ~/work/fortran/subroutine/fs/make_clt.f
c     21,22,31

      parameter(ns=100000)
      real xref1,yref1,xref2,yref2,xtrim1,ytrim1,xtrim2,ytrim2,cidrad
      real distcut,errorcut,chicut,sharpcut,skycut
      real sky1(ns),sky2(ns)
      real min_1,max_1,min_2,max_2,xx1(ns),yy1(ns),xx2(ns),yy2(ns),dist
      real x1(ns),y1(ns),mag1(ns),merr1(ns),sharp1(ns),chi1(ns),fwhm1
      real x2(ns),y2(ns),mag2(ns),merr2(ns),sharp2(ns),chi2(ns),fwhm2
      integer id1(ns),id2(ns),nstar1,nstar2,nn1,nn2,id11(ns),id22(ns)
      integer num_im,num_ref,iter1(ns),iter2(ns),itercut
      character als1*32,als2*32,H*1,line*100,als*32,output*32,c(32)
      open(22,file='input.make_clt')
      do i=1,ns
       read(22,'(t1,a1,t3,a32,t35,a100)',end=900) H,als,line
c       print *,als
       if(H.eq.'R') then
        print *,H,als,line
        read(als,'(32a1)') (c(j),j=1,32)
        do j=1,32
         if(c(j).eq.'.') goto 30
        enddo
30      write(c(j+1),'(a1)') 'c'
        write(c(j+2),'(a1)') 'l'
        write(c(j+3),'(a1)') 't'
        write(output,'(32a1)') (c(j),j=1,32)
        read(line,*) num_im,distcut,errorcut,chicut,sharpcut
     &  ,skycut,itercut
        num_ref=num_im
        call read_als(als,id1,x1,y1,mag1,merr1,sharp1,chi1,sky1,iter1,
     &  nstar1,fwhm1,min_1,max_1)
        call make_candidate(id1,x1,y1,mag1,merr1,sharp1,chi1,sky1,iter1
     &  ,nstar1,min_1,max_1,fwhm1,distcut,errorcut,chicut,sharpcut,
     &  skycut,itercut,0.05,1.5,id11,xx1,yy1,nn1)
        call write_clt(output,num_ref,num_im,fwhm1,cidrad,
     &  xref1,yref1,xref2,yref2,xtrim1,ytrim1,xtrim2,ytrim2,
     &  xx1,yy1,nn1,xx2,yy2,nn2)
       elseif(H.eq.'T') then
        print *,H,als,line
        read(als,'(32a1)') (c(j),j=1,32)
        do j=1,32
         if(c(j).eq.'.') goto 50
        enddo
50     write(c(j+1),'(a1)') 'c'
        write(c(j+2),'(a1)') 'l'
        write(c(j+3),'(a1)') 't'
        write(output,'(32a1)') (c(j),j=1,32)
        read(line,*) num_im,xref1,yref1,xref2,yref2,
     &  xtrim1,ytrim1,xtrim2,ytrim2,
     &  distcut2,errorcut2,chicut2,sharpcut2,cidrad
        call read_als(als,id2,x2,y2,mag2,merr2,sharp2,chi2,sky2,iter2,
     &  nstar2,fwhm2,min_2,max_2)
        call search_star(xref1,yref1,id1,x1,y1,mag1,merr1,sharp1,chi1,
     &  nstar1,fwhm1,xref1,yref1)
        call search_star(xref2,yref2,id1,x1,y1,mag1,merr1,sharp1,chi1,
     &  nstar1,fwhm1,xref2,yref2)
        call search_star(xtrim1,ytrim1,id2,x2,y2,mag2,merr2,sharp2,chi2,
     &  nstar2,fwhm1,xtrim1,ytrim1)
        call search_star(xtrim2,ytrim2,id2,x2,y2,mag2,merr2,sharp2,chi2,
     &  nstar2,fwhm2,xtrim2,ytrim2)
        call make_candidate(id2,x2,y2,mag2,merr2,sharp2,chi2,sky2,iter2,
     &  nstar2,min_2,max_2,fwhm2,distcut,errorcut,chicut,sharpcut,
     &  skycut,itercut,0.05,1.5,id22,xx2,yy2,nn2)
        call write_clt(output,num_ref,num_im,fwhm1,cidrad,
     &  xref1,yref1,xref2,yref2,xtrim1,ytrim1,xtrim2,ytrim2,
     &  xx1,yy1,nn1,xx2,yy2,nn2)
       endif
      enddo
900   close(22)
 
      stop
      end

      subroutine write_clt(output,num_refim,num_trim,fwhm1,cidrad,
     &xref1,yref1,xref2,yref2,xtrim1,ytrim1,xtrim2,ytrim2,
     &xx1,yy1,nn1,xx2,yy2,nn2)
      parameter(ns=100000)
      integer nn1,nn2,ncal
      real xref1,yref1,xref2,yref2,xtrim1,ytrim1,xtrim2,ytrim2,fwhm1
      real xx1(ns),yy1(ns),xx2(ns),yy2(ns)
      real scl,th,xcal,ycal,dist,cidrad
      character output*32
      open(31,file=output)
      write(31,'(t1,a1)') '#'
      write(31,*) num_refim,num_trim
      if(num_refim.eq.num_trim) goto 100
      call cal_2star(xtrim1,ytrim1,xref1,yref1,xtrim2,ytrim2,xref2,yref2
     &,scl,th)
      ncal=0
      do i=1,nn1
       call cal_2star2(xref1,yref1,xtrim1,ytrim1,xx1(i),yy1(i),scl,th,
     & xcal,ycal)
       do j=1,nn2
        dist=sqrt((xx2(j)-xcal)**2+(yy2(j)-ycal)**2)
        if(dist.le.fwhm1*cidrad) then
         write(31,'(2f9.3,1x,2f9.3,1x,f10.3)') 
     &   xx1(i),yy1(i),xx2(j),yy2(j),dist
         ncal=ncal+1
        endif
       enddo
      enddo      
      print *,ncal
100   close(31)
      return
      end

      subroutine make_candidate(id,x,y,amag,amerr,sharp,chi,sky,iter,
     &nstar,amin,amax,fwhm,th_dist,th_err,th_chi,th_sharp,th_sky,
     &iter_th,th_min,th_max,id_out,x_out,y_out,n_out)
      parameter(ns=100000)
      integer nstar,id(ns),id_out(ns),n_out,nclose(ns),iter(ns),iter_th
      real x(ns),y(ns),amag(ns),amerr(ns),sharp(ns),chi(ns),amin,amax
      real fwhm,th_err,th_chi,th_sharp,th_min,th_max,x_out(ns),y_out(ns)
      real th_dist,sky(ns),th_sky,sky_min,sky_max
      n_out=0
      do i=1,nstar
       nclose(i)=0
       if(i.eq.1) sky_min=sky(i)
       if(i.eq.1) sky_max=sky(i)
       if(i.ne.1) sky_min=min(sky_min,sky(i))
       if(i.ne.1) sky_max=max(sky_max,sky(i))
      enddo
      do i=1,nstar
       do j=1,nstar
        dist=sqrt((x(i)-x(j))**2+(y(i)-y(j))**2)
        if(i.ne.j.and.dist.le.fwhm*th_dist) then
         nclose(i)=nclose(i)+1
         nclose(j)=nclose(j)+1
        endif
       enddo
      enddo
      do i=1,nstar
       if(amag(i).gt.amin+th_min.and.amag(i).lt.amax-th_max.and.
     & amerr(i).le.th_err.and.chi(i).le.th_chi.and.nclose(i).eq.0.and.
     & sharp(i).ge.-th_sharp.and.sharp(i).le.th_sharp.and.
     & iter(i).le.iter_th.and.sky(i).le.sky_min+th_sky*(sky_max-sky_min)
     & ) then
        n_out=n_out+1
        id_out(n_out)=id(i)
        x_out(n_out)=x(i)
        y_out(n_out)=y(i)
       endif 
      enddo
      return
      end


      subroutine search_star(xc_in,yc_in,id,x,y,amag,amerr,sharp,chi,
     &nstar,fwhm,xc_out,yc_out)
      parameter(ns=100000)
      integer nstar,id(ns)
      real x(ns),y(ns),amag(ns),amerr(ns),sharp(ns),chi(ns)
      real xc_in,yc_in,xc_out,yc_out,dist,fwhm
      do i=1,nstar
       dist=sqrt((x(i)-xc_in)**2+(y(i)-yc_in)**2)
       if(dist.le.fwhm*0.5) then
        xc_out=x(i)
        yc_out=y(i)
        goto 100
       endif
      enddo
100   return
      end


      subroutine read_als(input,id,x,y,amag,amerr,sharp,chi,sky,iter,
     &nstar,fwhm,amag1,amag2)
      parameter(ns=100000)
      integer nstar,id(ns),iter(ns)
      real x(ns),y(ns),amag(ns),amerr(ns),sharp(ns),chi(ns),amag1,amag2
      real fwhm,sky(ns)
      character input*32
      open(21,file=input)
      do i=1,44
       read(21,'(t17,a32)') input
       if(i.eq.26) read(input,*) fwhm
      enddo
      do i=1,ns
       read(21,*,end=99) id(i),x(i),y(i),amag(i),amerr(i),sky(i),iter(i)
       read(21,*,end=99) sharp(i),chi(i)
       if(i.eq.1) amag1=amag(i)
       if(i.eq.1) amag2=amag(i)
       if(i.ne.1) amag1=min(amag1,amag(i))
       if(i.ne.1) amag2=max(amag2,amag(i))
      enddo
99    nstar=i-1
      close(21)
      return
      end
      subroutine cal_2star(x1,y1,xx1,yy1,x2,y2,xx2,yy2,scl,th)
      real x1,y1,x2,y2,xx1,yy1,xx2,yy2
      real dx,dy,dxx,dyy,dist_xy,dist_xxyy,th_xy,th_xxyy,th
      real pi,scl,xx22,yy22,ra33
      data pi/3.14159265358979d0/
      dx=x2-x1
      dy=y2-y1
      dxx=xx2-xx1
      dyy=yy2-yy1
      dist_xy=sqrt(dx*dx+dy*dy)
      dist_xxyy=sqrt(dxx*dxx+dyy*dyy)
      scl=dist_xxyy/dist_xy
      if(dx.gt.0.d0.and.dy.gt.0.d0) then!1st quater
       th_xy=atan(dy/dx)*180.d0/pi
      elseif(dx.lt.0.d0.and.dy.gt.0.d0) then!2nd quater
       th_xy=180.d0-atan(dy/(-dx))*180.d0/pi
      elseif(dx.lt.0.d0.and.dy.lt.0.d0) then!3rd quater
       th_xy=180.d0+atan(-dy/(-dx))*180.d0/pi
      elseif(dx.gt.0.d0.and.dy.lt.0.d0) then!4th quater
       th_xy=360.d0-atan(-dy/dx)*180.d0/pi
      endif
      if(dxx.gt.0.d0.and.dyy.gt.0.d0) then!1st quater
       th_xxyy=atan(dyy/(dxx))*180.d0/pi
      elseif(dxx.lt.0.d0.and.dyy.gt.0.d0) then!2nd quater
       th_xxyy=180.d0-atan(dyy/(-dxx))*180.d0/pi
      elseif(dxx.lt.0.d0.and.dyy.lt.0.d0) then!3rd quater
       th_xxyy=180.d0+atan(-dyy/(-dxx))*180.d0/pi
      elseif(dxx.gt.0.d0.and.dyy.lt.0.d0) then!4th quater
       th_xxyy=360.d0-atan(-dyy/(dxx))*180.d0/pi
      endif
      th=th_xxyy-th_xy
180   if(th.ge.0.d0.and.th.le.360.d0) goto 200
      if(th.lt.0.d0) th=th+360.d0
      if(th.gt.360.d0) th=th-360.d0
      goto 180
200   if(th.gt.0.d0.and.th.le.90.d0) then
       xx22=xx1+dx*scl*cos(pi*th/180.d0)-dy*scl*sin(pi*th/180.d0)
       yy22=yy1+dx*scl*sin(pi*th/180.d0)+dy*scl*cos(pi*th/180.d0)
      elseif(th.gt.90.d0.and.th.le.180.d0) then
       xx22=xx1-dx*scl*cos(pi*th/180.d0)-dy*scl*sin(pi*th/180.d0)
       yy22=yy1+dx*scl*sin(pi*th/180.d0)-dy*scl*cos(pi*th/180.d0)
      elseif(th.gt.180.d0.and.th.le.270.d0) then
       xx22=xx1-dx*scl*cos(pi*th/180.d0)+dy*scl*sin(pi*th/180.d0)
       yy22=yy1-dx*scl*sin(pi*th/180.d0)-dy*scl*cos(pi*th/180.d0)
      elseif(th.gt.270.d0.and.th.le.360.d0) then
       xx22=xx1+dx*scl*cos(pi*th/180.d0)+dy*scl*sin(pi*th/180.d0)
       yy22=yy1-dx*scl*sin(pi*th/180.d0)+dy*scl*cos(pi*th/180.d0)
      endif
      return
      end
      subroutine cal_2star2(x1,y1,xx1,yy1,x2,y2,scl,th,xx,yy)
      real scl,th,dx,dy,pi,x1,y1,xx1,yy1,x2,y2,xx,yy
      data pi/3.14159265358979d0/
      dx=x2-x1
      dy=y2-y1
      if(th.gt.0.and.th.le.90.d0) then
       xx=xx1+dx*scl*cos(pi*th/180.d0)+dy*scl*sin(pi*th/180.d0)
       yy=yy1-dx*scl*sin(pi*th/180.d0)+dy*scl*cos(pi*th/180.d0)
      elseif(th.gt.90.d0.and.th.le.180.d0) then
       xx=xx1-dx*scl*cos(pi*th/180.d0)+dy*scl*sin(pi*th/180.d0)
       yy=yy1-dx*scl*sin(pi*th/180.d0)-dy*scl*cos(pi*th/180.d0)
      elseif(th.gt.180.d0.and.th.le.270.d0) then
       xx=xx1-dx*scl*cos(pi*th/180.d0)-dy*scl*sin(pi*th/180.d0)
       yy=yy1+dx*scl*sin(pi*th/180.d0)-dy*scl*cos(pi*th/180.d0)
      elseif(th.gt.270.d0.and.th.le.360.d0) then
       xx=xx1+dx*scl*cos(pi*th/180.d0)-dy*scl*sin(pi*th/180.d0)
       yy=yy1+dx*scl*sin(pi*th/180.d0)+dy*scl*cos(pi*th/180.d0)
      endif
      return
      end

