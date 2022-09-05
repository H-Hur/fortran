c     gfortran -o merge_cid ~/work/fortran/subroutine/daotools/ctio4m_mosaic/merge_cid.f ~/work/fortran/subroutine/caltools/mat.f
c     21,22,31-42,81
      parameter (ns=500000,ni=20)
      integer nimage,nstars,IDF(ns,ni),mas(ni),nstar(ni),nchip_out
      character input(ni)*10,filter(ni)*10,inwc*32,fs*16
      real*8 MAGF(ns,ni,8),air(ni),ut(ni),psf(ni)
      real fitrad
      data fitrad/1.5/
c      write(inwc,'(a32)') 'i3s_1_to_2.dat                              '
      print *,' input =?'
      read(*,'(t1,a32)') inwc 
      call read_coo(inwc,nchip_out)
      call merge_fs(0,fitrad)

      write(inwc,'(a32)') 'out2_c1.cid_3rd                             '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstar,nstars)
      call comb_nobs_fs(1,nimage,input,filter,IDF,MAGF,mas,
     &air,psf,ut,nstar,nstars)
      call merge_fs(1,fitrad)
      write(inwc,'(a32)') 'out2_c2.cid_3rd                             '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstar,nstars)
      call comb_nobs_fs(2,nimage,input,filter,IDF,MAGF,mas,
     &air,psf,ut,nstar,nstars)
      call merge_fs(2,fitrad)
      write(inwc,'(a32)') 'out2_c3.cid_3rd                             '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstar,nstars)
      call comb_nobs_fs(3,nimage,input,filter,IDF,MAGF,mas,
     &air,psf,ut,nstar,nstars)
      call merge_fs(3,fitrad)
      write(inwc,'(a32)') 'out2_c4.cid_3rd                             '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstar,nstars)
      call comb_nobs_fs(4,nimage,input,filter,IDF,MAGF,mas,
     &air,psf,ut,nstar,nstars)
      call merge_fs(4,fitrad)
      write(inwc,'(a32)') 'out2_c5.cid_3rd                             '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstar,nstars)
      call comb_nobs_fs(5,nimage,input,filter,IDF,MAGF,mas,
     &air,psf,ut,nstar,nstars)
      call merge_fs(5,fitrad)
      write(inwc,'(a32)') 'out2_c6.cid_3rd                             '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstar,nstars)
      call comb_nobs_fs(6,nimage,input,filter,IDF,MAGF,mas,
     &air,psf,ut,nstar,nstars)
      call merge_fs(6,fitrad)
      write(inwc,'(a32)') 'out2_c7.cid_3rd                             '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstar,nstars)
      call comb_nobs_fs(7,nimage,input,filter,IDF,MAGF,mas,
     &air,psf,ut,nstar,nstars)
      call merge_fs(7,fitrad)
      write(inwc,'(a32)') 'out2_c8.cid_3rd                             '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstar,nstars)
      call comb_nobs_fs(8,nimage,input,filter,IDF,MAGF,mas,
     &air,psf,ut,nstar,nstars)
      call merge_fs(8,fitrad)


          
      call write_trs(nchip_out,fitrad)

      stop
      end

      subroutine write_trs(nc,fitrad)
      parameter (ns=500000,ni=20,nx=1000)
      common/coo_in/xc1,yc1,xc2,yc2,ncoo,filter1,filter2,exps1,exps2
      real*8 xc1(ns),yc1(ns),xc2(ns),yc2(ns),xq(nx),yq(nx),x,y,xx,yy
      real*8 xeq(nx),yeq(nx),xxx,yyy 
      character filter1*2,filter2*2,exps1*5,exps2*5
      common/chips/id_fs(ns,ni),nchips(ns),n_stars,
     &x_im(ns,ni),y_im(ns,ni),xc(ns,ni),yc(ns,ni),
     &rmag(ns,ni),rmerr(ns,ni),
     &x_f_im(ns,ni),y_f_im(ns,ni),xc_f(ns,ni),yc_f(ns,ni),
     &nul,nus,nbl,nbs,nvl,nvs,nrl,nrs,nil,nis,nhal,nhas
      common/table_merge/xus(ns),yus(ns),n_us,xul(ns),yul(ns),n_ul,
     &xbs(ns),ybs(ns),n_bs,xbl(ns),ybl(ns),n_bl,
     &xvs(ns),yvs(ns),n_vs,xvl(ns),yvl(ns),n_vl,
     &xrs(ns),yrs(ns),n_rs,xrl(ns),yrl(ns),n_rl,
     &xis(ns),yis(ns),n_is,xil(ns),yil(ns),n_il,
     &xhs(ns),yhs(ns),n_hs,xhl(ns),yhl(ns),n_hl
      common/chip/x0ccd(8),y0ccd(8),rot(8)
      integer nc
      real fitrad,d
      call make_3xy(10,xc1,yc1,xc2,ncoo,xq,1.d0)
      call make_3xy(10,xc1,yc1,yc2,ncoo,yq,1.d0)
c      if(n_ul.ne.0) open(31,file='u_long.coo')
c      if(n_us.ne.0) open(32,file='u_short.coo')
c      if(n_bl.ne.0) open(33,file='b_long.coo')
c      if(n_bs.ne.0) open(34,file='b_short.coo')
      if(n_vl.ne.0) then
       open(35,file='v_long.coo')
       call make_transform_filter(nil,nvl,xeq,yeq) 
       do i=1,n_vl
        x=xvl(i)
        y=yvl(i)
        xx=xq(1)+xq(2)*x + xq(3)*y + xq(4)*x*x + xq(5)*x*y + xq(6)*y*y +
     &  xq(7)*x*x*x + xq(8)*x*x*y + xq(9)*x*y*y + xq(10)*y*y*y
        yy=yq(1)+yq(2)*x + yq(3)*y + yq(4)*x*x + yq(5)*x*y + yq(6)*y*y +
     &  yq(7)*x*x*x + yq(8)*x*x*y + yq(9)*x*y*y + yq(10)*y*y*y
        x=xx
        y=yy
        xx=xeq(1)+xeq(2)*x+xeq(3)*y+xeq(4)*x*x+xeq(5)*x*y+xeq(6)*y*y +
     &  xeq(7)*x*x*x+xeq(8)*x*x*y+xeq(9)*x*y*y+xeq(10)*y*y*y
        yy=yeq(1)+yeq(2)*x+yeq(3)*y+yeq(4)*x*x+yeq(5)*x*y+yeq(6)*y*y +
     &  yeq(7)*x*x*x +yeq(8)*x*x*y +yeq(9)*x*y*y +yeq(10)*y*y*y
c        write(*,'(10f10.3)') xvl(i),yvl(i),x,y,xx,yy
        if(xx-x0ccd(nc).ge.fitrad.and.xx-x0ccd(nc).le.2048.-fitrad.and.
     &  yy-y0ccd(nc).ge.fitrad.and.yy-y0ccd(nc).le.4096.-fitrad)
     &  write(35,'(2f10.3)') xx-x0ccd(nc),yy-y0ccd(nc)
       enddo
      endif
      if(n_vs.ne.0) then
       open(36,file='v_short.coo')
       call make_transform_filter(nil,nvs,xeq,yeq)
       do i=1,n_vs
        x=xvs(i)
        y=yvs(i)
        xx=xq(1)+xq(2)*x + xq(3)*y + xq(4)*x*x + xq(5)*x*y + xq(6)*y*y +
     &  xq(7)*x*x*x + xq(8)*x*x*y + xq(9)*x*y*y + xq(10)*y*y*y
        yy=yq(1)+yq(2)*x + yq(3)*y + yq(4)*x*x + yq(5)*x*y + yq(6)*y*y +
     &  yq(7)*x*x*x + yq(8)*x*x*y + yq(9)*x*y*y + yq(10)*y*y*y
        x=xx
        y=yy
        xx=xeq(1)+xeq(2)*x+xeq(3)*y+xeq(4)*x*x+xeq(5)*x*y+xeq(6)*y*y +
     &  xeq(7)*x*x*x+xeq(8)*x*x*y+xeq(9)*x*y*y+xeq(10)*y*y*y
        yy=yeq(1)+yeq(2)*x+yeq(3)*y+yeq(4)*x*x+yeq(5)*x*y+yeq(6)*y*y +
     &  yeq(7)*x*x*x +yeq(8)*x*x*y +yeq(9)*x*y*y +yeq(10)*y*y*y
        if(xx-x0ccd(nc).ge.fitrad.and.xx-x0ccd(nc).le.2048.-fitrad.and.
     &  yy-y0ccd(nc).ge.fitrad.and.yy-y0ccd(nc).le.4096.-fitrad)
     &  write(36,'(2f10.3)') xx-x0ccd(nc),yy-y0ccd(nc)
       enddo
      endif 
      if(n_rl.ne.0) then
       open(37,file='r_long.coo')
       call make_transform_filter(nil,nrl,xeq,yeq)
       do i=1,n_rl
        x=xrl(i)
        y=yrl(i)
        xx=xq(1)+xq(2)*x + xq(3)*y + xq(4)*x*x + xq(5)*x*y + xq(6)*y*y +
     &  xq(7)*x*x*x + xq(8)*x*x*y + xq(9)*x*y*y + xq(10)*y*y*y
        yy=yq(1)+yq(2)*x + yq(3)*y + yq(4)*x*x + yq(5)*x*y + yq(6)*y*y +
     &  yq(7)*x*x*x + yq(8)*x*x*y + yq(9)*x*y*y + yq(10)*y*y*y
        x=xx
        y=yy
        xx=xeq(1)+xeq(2)*x+xeq(3)*y+xeq(4)*x*x+xeq(5)*x*y+xeq(6)*y*y +
     &  xeq(7)*x*x*x+xeq(8)*x*x*y+xeq(9)*x*y*y+xeq(10)*y*y*y
        yy=yeq(1)+yeq(2)*x+yeq(3)*y+yeq(4)*x*x+yeq(5)*x*y+yeq(6)*y*y +
     &  yeq(7)*x*x*x +yeq(8)*x*x*y +yeq(9)*x*y*y +yeq(10)*y*y*y
        if(xx-x0ccd(nc).ge.fitrad.and.xx-x0ccd(nc).le.2048.-fitrad.and.
     &  yy-y0ccd(nc).ge.fitrad.and.yy-y0ccd(nc).le.4096.-fitrad)
     &  write(37,'(2f10.3)') xx-x0ccd(nc),yy-y0ccd(nc)
       enddo
      endif
      if(n_rs.ne.0) then
       open(38,file='r_short.coo')
       call make_transform_filter(nil,nrs,xeq,yeq)
       do i=1,n_rs
        x=xrs(i)
        y=yrs(i)
        xx=xq(1)+xq(2)*x + xq(3)*y + xq(4)*x*x + xq(5)*x*y + xq(6)*y*y +
     &  xq(7)*x*x*x + xq(8)*x*x*y + xq(9)*x*y*y + xq(10)*y*y*y
        yy=yq(1)+yq(2)*x + yq(3)*y + yq(4)*x*x + yq(5)*x*y + yq(6)*y*y +
     &  yq(7)*x*x*x + yq(8)*x*x*y + yq(9)*x*y*y + yq(10)*y*y*y
        x=xx
        y=yy
        xx=xeq(1)+xeq(2)*x+xeq(3)*y+xeq(4)*x*x+xeq(5)*x*y+xeq(6)*y*y +
     &  xeq(7)*x*x*x+xeq(8)*x*x*y+xeq(9)*x*y*y+xeq(10)*y*y*y
        yy=yeq(1)+yeq(2)*x+yeq(3)*y+yeq(4)*x*x+yeq(5)*x*y+yeq(6)*y*y +
     &  yeq(7)*x*x*x +yeq(8)*x*x*y +yeq(9)*x*y*y +yeq(10)*y*y*y
        if(xx-x0ccd(nc).ge.fitrad.and.xx-x0ccd(nc).le.2048.-fitrad.and.
     &  yy-y0ccd(nc).ge.fitrad.and.yy-y0ccd(nc).le.4096.-fitrad)
     &  write(38,'(2f10.3)') xx-x0ccd(nc),yy-y0ccd(nc)
       enddo
      endif
      if(n_il.ne.0) then
       open(39,file='i_long.coo')
       do i=1,n_il
        x=xil(i)
        y=yil(i)
        xx=xq(1)+xq(2)*x + xq(3)*y + xq(4)*x*x + xq(5)*x*y + xq(6)*y*y +
     &  xq(7)*x*x*x + xq(8)*x*x*y + xq(9)*x*y*y + xq(10)*y*y*y
        yy=yq(1)+yq(2)*x + yq(3)*y + yq(4)*x*x + yq(5)*x*y + yq(6)*y*y +
     &  yq(7)*x*x*x + yq(8)*x*x*y + yq(9)*x*y*y + yq(10)*y*y*y
c        x=xx
c        y=yy
c        d=sqrt((xx-x0ccd(nc)-1024.)**2+(yy-y0ccd(nc)-2048)**2.)
c        xxx=xx-x0ccd(nc)-d*cos(rot(nc))
c        yyy=yy-y0ccd(nc)-d*sin(rot(nc))
c        if(xxx.ge.fitrad.and.xxx.le.2048.-fitrad.and.
c     &  yyy.ge.fitrad.and.yyy.le.4096.-fitrad)
c     &  write(39,'(2f10.3)') xxx,yyy
        if(xx-x0ccd(nc).ge.fitrad.and.xx-x0ccd(nc).le.2048.-fitrad.and.
     &  yy-y0ccd(nc).ge.fitrad.and.yy-y0ccd(nc).le.4096.-fitrad)
     &  write(39,'(2f10.3)') xx-x0ccd(nc),yy-y0ccd(nc)

       enddo
      endif
      if(n_is.ne.0) then
       open(40,file='i_short.coo')
       call make_transform_filter(nil,nis,xeq,yeq)
       do i=1,n_is
        x=xis(i)
        y=yis(i)
        xx=xq(1)+xq(2)*x + xq(3)*y + xq(4)*x*x + xq(5)*x*y + xq(6)*y*y +
     &  xq(7)*x*x*x + xq(8)*x*x*y + xq(9)*x*y*y + xq(10)*y*y*y
        yy=yq(1)+yq(2)*x + yq(3)*y + yq(4)*x*x + yq(5)*x*y + yq(6)*y*y +
     &  yq(7)*x*x*x + yq(8)*x*x*y + yq(9)*x*y*y + yq(10)*y*y*y
        x=xx
        y=yy
        xx=xeq(1)+xeq(2)*x+xeq(3)*y+xeq(4)*x*x+xeq(5)*x*y+xeq(6)*y*y +
     &  xeq(7)*x*x*x+xeq(8)*x*x*y+xeq(9)*x*y*y+xeq(10)*y*y*y
        yy=yeq(1)+yeq(2)*x+yeq(3)*y+yeq(4)*x*x+yeq(5)*x*y+yeq(6)*y*y +
     &  yeq(7)*x*x*x +yeq(8)*x*x*y +yeq(9)*x*y*y +yeq(10)*y*y*y
        if(xx-x0ccd(nc).ge.fitrad.and.xx-x0ccd(nc).le.2048.-fitrad.and.
     &  yy-y0ccd(nc).ge.fitrad.and.yy-y0ccd(nc).le.4096.-fitrad) !.and.
c     &  nchips(i).eq.nc)
     &  write(40,'(2f10.3)') xx-x0ccd(nc),yy-y0ccd(nc)
       enddo
      endif
      if(n_hl.ne.0) then
       open(41,file='ha_long.coo')
       call make_transform_filter(nil,nhal,xeq,yeq)
       do i=1,n_hl
        x=xhl(i)
        y=yhl(i)
        xx=xq(1)+xq(2)*x + xq(3)*y + xq(4)*x*x + xq(5)*x*y + xq(6)*y*y +
     &  xq(7)*x*x*x + xq(8)*x*x*y + xq(9)*x*y*y + xq(10)*y*y*y
        yy=yq(1)+yq(2)*x + yq(3)*y + yq(4)*x*x + yq(5)*x*y + yq(6)*y*y +
     &  yq(7)*x*x*x + yq(8)*x*x*y + yq(9)*x*y*y + yq(10)*y*y*y
        x=xx
        y=yy
c        xx=xeq(1)+xeq(2)*x+xeq(3)*y+xeq(4)*x*x+xeq(5)*x*y+xeq(6)*y*y +
c     &  xeq(7)*x*x*x+xeq(8)*x*x*y+xeq(9)*x*y*y+xeq(10)*y*y*y
c        yy=yeq(1)+yeq(2)*x+yeq(3)*y+yeq(4)*x*x+yeq(5)*x*y+yeq(6)*y*y +
c     &  yeq(7)*x*x*x +yeq(8)*x*x*y +yeq(9)*x*y*y +yeq(10)*y*y*y
        if(xx-x0ccd(nc).ge.fitrad.and.xx-x0ccd(nc).le.2048.-fitrad.and.
     &  yy-y0ccd(nc).ge.fitrad.and.yy-y0ccd(nc).le.4096.-fitrad) !.and.
c     &  nchips(i).eq.nc)
     &  write(41,'(2f10.3)') xx-x0ccd(nc),yy-y0ccd(nc)
       enddo
      endif
      if(n_hs.ne.0) then
       open(42,file='ha_short.coo')
       call make_transform_filter(nil,nhas,xeq,yeq)
       do i=1,n_hs
        x=xhs(i)
        y=yhs(i)
        xx=xq(1)+xq(2)*x + xq(3)*y + xq(4)*x*x + xq(5)*x*y + xq(6)*y*y +
     &  xq(7)*x*x*x + xq(8)*x*x*y + xq(9)*x*y*y + xq(10)*y*y*y
        yy=yq(1)+yq(2)*x + yq(3)*y + yq(4)*x*x + yq(5)*x*y + yq(6)*y*y +
     &  yq(7)*x*x*x + yq(8)*x*x*y + yq(9)*x*y*y + yq(10)*y*y*y
        x=xx
        y=yy
c        xx=xeq(1)+xeq(2)*x+xeq(3)*y+xeq(4)*x*x+xeq(5)*x*y+xeq(6)*y*y +
c     &  xeq(7)*x*x*x+xeq(8)*x*x*y+xeq(9)*x*y*y+xeq(10)*y*y*y
c        yy=yeq(1)+yeq(2)*x+yeq(3)*y+yeq(4)*x*x+yeq(5)*x*y+yeq(6)*y*y +
c     &  yeq(7)*x*x*x +yeq(8)*x*x*y +yeq(9)*x*y*y +yeq(10)*y*y*y
        if(xx-x0ccd(nc).ge.fitrad.and.xx-x0ccd(nc).le.2048.-fitrad.and.
     &  yy-y0ccd(nc).ge.fitrad.and.yy-y0ccd(nc).le.4096.-fitrad) !.and.
c     &  nchips(i).eq.nc)
     &  write(42,'(2f10.3)') xx-x0ccd(nc),yy-y0ccd(nc)
       enddo
      endif 
      do i=31,42
       close(i)
      enddo
      return
      end

      subroutine make_transform_filter(im1,im2,xeq,yeq)
      parameter (ns=500000,ni=20,nx=1000)
      real*8 xeq(nx),yeq(nx),x1(ns),y1(ns),x2(ns),y2(ns)
      common/chips/id_fs(ns,ni),nchips(ns),n_stars,
     &x_im(ns,ni),y_im(ns,ni),xc(ns,ni),yc(ns,ni),
     &rmag(ns,ni),rmerr(ns,ni),
     &x_f_im(ns,ni),y_f_im(ns,ni),xc_f(ns,ni),yc_f(ns,ni),
     &nul,nus,nbl,nbs,nvl,nvs,nrl,nrs,nil,nis,nhal,nhas
      integer n,n1,n2,n3,n4,n5,n6,n7,n8,n9
      real cut
      n1=0
      n2=0
      n3=0
      n4=0
      n5=0
      do i=1,n_stars
       if(id_fs(i,im1).ne.0.and.rmerr(i,im1).le.0.01.and.
     & id_fs(i,im2).ne.0.and.rmerr(i,im2).le.0.01) n1=n1+1
       if(id_fs(i,im1).ne.0.and.rmerr(i,im1).le.0.02.and.
     & id_fs(i,im2).ne.0.and.rmerr(i,im2).le.0.02) n2=n2+1
       if(id_fs(i,im1).ne.0.and.rmerr(i,im1).le.0.03.and.
     & id_fs(i,im2).ne.0.and.rmerr(i,im2).le.0.03) n3=n3+1
       if(id_fs(i,im1).ne.0.and.rmerr(i,im1).le.0.04.and.
     & id_fs(i,im2).ne.0.and.rmerr(i,im2).le.0.04) n4=n4+1
       if(id_fs(i,im1).ne.0.and.rmerr(i,im1).le.0.05.and.
     & id_fs(i,im2).ne.0.and.rmerr(i,im2).le.0.05) n5=n5+1
      enddo
      cut=0.05
      if(n5.gt.1000) cut=0.04
      if(n4.gt.1000) cut=0.03
      if(n3.gt.1000) cut=0.02
      if(n2.gt.1000) cut=0.01
      n=0
      do i=1,n_stars
       if(id_fs(i,im1).ne.0.and.rmerr(i,im1).le.cut.and.
     & id_fs(i,im2).ne.0.and.rmerr(i,im2).le.cut) then
        n=n+1
        x1(n)=xc(i,im1)   
        y1(n)=yc(i,im1)
        x2(n)=xc(i,im2)
        y2(n)=yc(i,im2)
       endif
      enddo
      call make_3xy(10,x1,y1,x2,n,xeq,1.d0)
      call make_3xy(10,x1,y1,y2,n,yeq,1.d0)
      return
      end
 
      subroutine read_coo(input,nchip_out)
      parameter (ns=500000)
      common/coo_in/xc1,yc1,xc2,yc2,ncoo,filter1,filter2,exps1,exps2
      common/chip/x0ccd(8),y0ccd(8),rot(8)
      real*8 xc1(ns),yc1(ns),xc2(ns),yc2(ns),d
      character input*32,H*1,filter1*2,filter2*2,exps1*5,exps2*5
      integer nchip_out
      open(22,file=input)
      read(22,'(t1,a1)') H
      read(22,'(t1,a2)') filter1
      read(22,'(t1,a5)') exps1
      read(22,'(t1,a2)') filter2
      read(22,'(t1,a5)') exps2
      read(22,*) nchip_out
      nc=nchip_out
      read(22,'(t1,a1)') H
      do i=1,ns
       read(22,*,end=99) xc1(i),yc1(i),xc2(i),yc2(i)
       xc1(i)=xc1(i)+x0ccd(nc)
       xc2(i)=xc2(i)+x0ccd(nc)
       yc1(i)=yc1(i)+y0ccd(nc)
       yc2(i)=yc2(i)+y0ccd(nc)
      enddo
99    ncoo=i-1
      return
      end
      subroutine merge_fs(nchip,fit)
      parameter (ns=500000,ni=20)
      common/chips/id_fs(ns,ni),nchips(ns),n_stars,
     &x_im(ns,ni),y_im(ns,ni),xc(ns,ni),yc(ns,ni),
     &rmag(ns,ni),rmerr(ns,ni),
     &x_f_im(ns,ni),y_f_im(ns,ni),xc_f(ns,ni),yc_f(ns,ni),
     &nul,nus,nbl,nbs,nvl,nvs,nrl,nrs,nil,nis,nhal,nhas
      common/table_merge/xus(ns),yus(ns),n_us,xul(ns),yul(ns),n_ul,
     &xbs(ns),ybs(ns),n_bs,xbl(ns),ybl(ns),n_bl,
     &xvs(ns),yvs(ns),n_vs,xvl(ns),yvl(ns),n_vl,
     &xrs(ns),yrs(ns),n_rs,xrl(ns),yrl(ns),n_rl,
     &xis(ns),yis(ns),n_is,xil(ns),yil(ns),n_il,
     &xhs(ns),yhs(ns),n_hs,xhl(ns),yhl(ns),n_hl
      real fit
      integer n_c
      open(81,file='dummy.dat')
      if(nchip.eq.0) then!initialyzing
       n_ul=0
       n_us=0
       n_bl=0
       n_bs=0
       n_vl=0
       n_vs=0
       n_rl=0
       n_rs=0
       n_il=0
       n_is=0
       n_hl=0
       n_hs=0
       goto 900
      endif
      do i=1,nstars
       if(id_fs(i,nul).ne.0) then
        n_c=0
        if(n_ul.ne.0) 
     &  call search_close(fit,xc_f(i,nul),yc_f(i,nul),xul,yul,n_ul,n_c)
        if(n_ul.eq.0.or.n_c.eq.0) then
         n_ul=n_ul+1
         xul(n_ul)=xc_f(i,nul)
         yul(n_ul)=yc_f(i,nul)
        endif
       endif
       if(id_fs(i,nus).ne.0) then
        n_c=0
        if(n_us.ne.0) 
     &  call search_close(fit,xc_f(i,nus),yc_f(i,nus),xus,yus,n_us,n_c)
        if(n_us.eq.0.or.n_c.eq.0) then
         n_us=n_us+1
         xus(n_us)=xc_f(i,nus)
         yus(n_us)=yc_f(i,nus)
        endif
       endif

       if(id_fs(i,nbl).ne.0) then
        n_c=0
        if(n_bl.ne.0)
     &  call search_close(fit,xc_f(i,nbl),yc_f(i,nbl),xbl,ybl,n_bl,n_c)
        if(n_bl.eq.0.or.n_c.eq.0) then
         n_bl=n_bl+1
         xbl(n_bl)=xc_f(i,nbl)
         ybl(n_bl)=yc_f(i,nbl)
        endif
       endif
       if(id_fs(i,nbs).ne.0) then
        n_c=0
        if(n_bs.ne.0)
     &  call search_close(fit,xc_f(i,nbs),yc_f(i,nbs),xbs,ybs,n_bs,n_c)
        if(n_bs.eq.0.or.n_c.eq.0) then
         n_bs=n_bs+1
         xbs(n_bs)=xc_f(i,nbs)
         ybs(n_bs)=yc_f(i,nbs)
        endif
       endif

       if(id_fs(i,nvl).ne.0) then
        n_c=0
        if(n_vl.ne.0)
     &  call search_close(fit,xc_f(i,nvl),yc_f(i,nvl),xvl,yvl,n_vl,n_c)
c        print *,i,n_vl
        if(n_vl.eq.0.or.n_c.eq.0) then
         n_vl=n_vl+1
         xvl(n_vl)=xc_f(i,nvl)
         yvl(n_vl)=yc_f(i,nvl)
        endif
       endif
       if(id_fs(i,nvs).ne.0) then
        n_c=0
        if(n_vs.ne.0)
     &  call search_close(fit,xc_f(i,nvs),yc_f(i,nvs),xvs,yvs,n_vs,n_c)
        if(n_vs.eq.0.or.n_c.eq.0) then
         n_vs=n_vs+1
         xvs(n_vs)=xc_f(i,nvs)
         yvs(n_vs)=yc_f(i,nvs)
        endif
       endif

       if(id_fs(i,nrl).ne.0) then
        n_c=0
        if(n_rl.ne.0)
     &  call search_close(fit,xc_f(i,nrl),yc_f(i,nrl),xrl,yrl,n_rl,n_c)
        if(n_rl.eq.0.or.n_c.eq.0) then
         n_rl=n_rl+1
         xrl(n_rl)=xc_f(i,nrl)
         yrl(n_rl)=yc_f(i,nrl)
        endif
       endif
       if(id_fs(i,nrs).ne.0) then
        n_c=0
        if(n_rs.ne.0)
     &  call search_close(fit,xc_f(i,nrs),yc_f(i,nrs),xrs,yrs,n_rs,n_c)
        if(n_rs.eq.0.or.n_c.eq.0) then
         n_rs=n_rs+1
         xrs(n_rs)=xc_f(i,nrs)
         yrs(n_rs)=yc_f(i,nrs)
        endif
       endif

       if(id_fs(i,nil).ne.0) then
        n_c=0
        if(n_il.ne.0)
     &  call search_close(fit,xc_f(i,nil),yc_f(i,nil),xil,yil,n_il,n_c)
        if(n_il.eq.0.or.n_c.eq.0) then
         n_il=n_il+1
         xil(n_il)=xc_f(i,nil)
         yil(n_il)=yc_f(i,nil)
        endif
       endif
       if(id_fs(i,nis).ne.0) then
        n_c=0
        if(n_is.ne.0) 
     &  call search_close(fit,xc_f(i,nis),yc_f(i,nis),xis,yis,n_is,n_c)
        if(n_is.eq.0.or.n_c.eq.0) then
         n_is=n_is+1
         xis(n_is)=xc_f(i,nis)
         yis(n_is)=yc_f(i,nis)
        endif
       endif

       if(id_fs(i,nhal).ne.0) then
        write(81,*) ,i,nhal,n_hl
        n_c=0
        if(n_hl.ne.0)
     &  call search_close(fit,xc_f(i,nhal),yc_f(i,nhal),
     &  xhl,yhl,n_hl,n_c)
        if(n_hl.eq.0.or.n_c.eq.0) then
         n_hl=n_hl+1
         xhl(n_hl)=xc_f(i,nhal)
         yhl(n_hl)=yc_f(i,nhal)
        endif
       endif
       if(id_fs(i,nhas).ne.0) then
        n_c=0
        if(n_hs.ne.0)
     &  call search_close(fit,xc_f(i,nhas),yc_f(i,nhas),
     &  xhs,yhs,n_hs,n_c)
        if(n_hs.eq.0.or.n_c.eq.0) then
         n_hs=n_hs+1
         xhs(n_hs)=xc_f(i,nhas)
         yhs(n_hs)=yc_f(i,nhas)
        endif     
       endif
      enddo
      print *,n_ul,n_us,n_bl,n_bs,n_vl,n_vs,n_rl,n_rs,n_il,n_is,
     &n_hl,n_hs
      close(81)
900   return
      end
      subroutine search_close(fitrad,xc,yc,xl,yl,nl,n_close)
      integer n_close,nl
      real fitrad,xc,yc,xl(nl),yl(nl)
      n_close=0
      do i=1,nl
       if(sqrt((xc-xl(i))**2+(yc-yl(i))**2).le.fitrad) n_close=n_close+1
      enddo 
      return
      end
      subroutine comb_nobs_fs(num_chip,nimage,input,filter,IDF,MAGF,mas,
     &air,psf,ut,nstar,nstars)
      parameter (ns=500000,ni=20)
      integer nimage,nstars,IDF(ns,ni),mas(ni),num_chip,nstar(ni)
      character input(ni)*10,filter(ni)*10,inwc*16,H*1,line*400,fh*2
      character ch(10)*1
      real*8 MAGF(ns,ni,8),air(ni),psf(ni),ut(ni)
      common/chip/x0ccd(8),y0ccd(8),rot(8)
      common/chips/id_fs(ns,ni),nchips(ns),n_stars,
     &x_im(ns,ni),y_im(ns,ni),xc(ns,ni),yc(ns,ni),
     &rmag(ns,ni),rmerr(ns,ni),
     &x_f_im(ns,ni),y_f_im(ns,ni),xc_f(ns,ni),yc_f(ns,ni),
     &nul,nus,nbl,nbs,nvl,nvs,nrl,nrs,nil,nis,nhal,nhas
      integer nu,nb,nv,nr,nii,nha,n_hl,n_hs
      integer n_ul,n_us,n_bl,n_bs,n_vl,n_vs,n_rl,n_rs,n_il,n_is
      nu=0
      nb=0
      nv=0
      nr=0 
      nii=0
      nha=0 
      n_ul=0
      n_us=0
      n_bl=0
      n_bs=0
      n_vl=0
      n_vs=0
      n_rl=0
      n_rs=0
      n_il=0
      n_is=0
      n_hl=0
      n_hs=0
      do i=1,nimage
       read(filter(i),'(10a1)') (ch(j),j=1,10)
       do j=1,10
        if(ch(j).eq.'U') write(fh,'(a2)') 'U '
        if(ch(j).eq.'U') goto 10
        if(ch(j).eq.'B') write(fh,'(a2)') 'B '
        if(ch(j).eq.'B') goto 10
        if(ch(j).eq.'V') write(fh,'(a2)') 'V '
        if(ch(j).eq.'V') goto 10
        if(ch(j).eq.'R') write(fh,'(a2)') 'R '
        if(ch(j).eq.'R') goto 10
        if(ch(j).eq.'I') write(fh,'(a2)') 'I '
        if(ch(j).eq.'I') goto 10
        if(ch(j).eq.'H'.and.ch(j+1).eq.'a'.and.j.le.9) 
     &  write(fh,'(a2)') 'Ha'
        if(ch(j).eq.'H'.and.ch(j+1).eq.'a'.and.j.le.9) goto 10
       enddo
       write(fh,'(a2)') '  '
10     if(fh.eq.'U ') then
        if(nu.eq.0) then
         nu=nu+1
         nul=i
         nus=i
         n_ul=nstar(i)
         n_us=nstar(i)
        elseif(nstar(i).gt.u_ul) then
         nu=nu+1
         nul=i
         n_ul=nstar(i)
        elseif(nstar(i).lt.n_us) then
         nu=nu+1
         nus=i
         n_us=nstar(i)
        endif
       endif
       if(fh.eq.'B ') then
        if(nb.eq.0) then
         nb=nb+1
         nbl=i
         nbs=i
         n_bl=nstar(i)
         n_bs=nstar(i)
        elseif(nstar(i).gt.n_bl) then
         nb=nb+1
         nbl=i
         n_bl=nstar(i)
        elseif(nstar(i).lt.n_bs) then
         nb=nb+1
         nbs=i
         n_bs=nstar(i)
        endif
       endif
       if(fh.eq.'V ') then
        if(nv.eq.0) then
         nv=nv+1
         nvl=i
         nvs=i
         n_vl=nstar(i)
         n_vs=nstar(i)
        elseif(nstar(i).gt.n_vl) then
         nv=nv+1
         vul=i
         n_vl=nstar(i)
        elseif(nstar(i).lt.n_vs) then
         nv=nv+1
         nvs=i
         n_vs=nstar(i)
        endif
       endif
       if(fh.eq.'R ') then
        if(nr.eq.0) then
         nr=nr+1
         nrl=i
         nrs=i
         n_rl=nstar(i)
         n_rs=nstar(i)
        elseif(nstar(i).gt.n_rl) then
         nr=nr+1
         nrl=i
         n_rl=nstar(i)
        elseif(nstar(i).lt.n_rs) then
         nr=nr+1
         nrs=i
         n_rs=nstar(i)
        endif
       endif
       if(fh.eq.'I ') then
        if(nii.eq.0) then
         nii=nii+1
         nil=i
         nis=i
         n_il=nstar(i)
         n_is=nstar(i)
        elseif(nstar(i).gt.n_il) then
         nii=nii+1
         nil=i
         n_il=nstar(i)
        elseif(nstar(i).lt.n_is) then
         nii=nii+1
         nis=i
         n_is=nstar(i)
        endif
       endif
       if(fh.eq.'Ha') then!###check it works well
        if(nha.eq.0) then
         nha=nha+1
         nhal=i
         nhas=i
         n_hl=nstar(i)
         n_hs=nstar(i)
        elseif(nstar(i).gt.n_hl) then
         nha=nha+1
         nhal=i
         n_hl=nstar(i)
        elseif(nstar(i).lt.n_hs) then
         nha=nha+1
         nhas=i
         n_hs=nstar(i)
        endif
       endif
      enddo

      do i=1,nstars
       nchips(i)=num_chip
       do j=1,nimage
        x_im(i,j)=magf(i,j,1)
        y_im(i,j)=magf(i,j,2)
        xc(i,j)=magf(i,j,1)+x0ccd(num_chip)
        yc(i,j)=magf(i,j,2)+y0ccd(num_chip)
        id_fs(i,j)=idf(i,j)
        x_f_im(i,j)=magf(i,j,3)
        y_f_im(i,j)=magf(i,j,4)
        xc_f(i,j)=magf(i,j,3)+x0ccd(num_chip)
        yc_f(i,j)=magf(i,j,4)+y0ccd(num_chip)
        rmag(i,j)=magf(i,j,5)
        rmerr(i,j)=magf(i,j,6)
       enddo
      enddo
      n_stars=nstars
      return
      end
      
      subroutine read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstar,nstars)
      parameter (ns=500000,ni=20)
      integer nimage,nstars,IDF(ns,ni),mas(ni),nstar(ni)
      character input(ni)*10,filter(ni)*10,inwc*32,H*1,line*400
      real*8 MAGF(ns,ni,8),air(ni),psf(ni),ut(ni)
      open(21,file=inwc)
      read(21,'(t1,30a10)') (input(j),j=1,ni)
      do i=1,ni
       psf(i)=0.d0
       nstar(i)=0
      enddo
      do i=1,100
       read(21,'(t1,a1,t2,a400)') H,line
       if(H.eq.'*') goto 99
       if(H.eq.'F') then
        read(line,'(t6,30a10)') (filter(j),j=1,ni)
       elseif(H.eq.'P') then
        read(line,'(t2,30f10.6)') (psf(j),j=1,ni)
       elseif(H.eq.'T') then
        read(line,'(t2,30f10.7)') (ut(j),j=1,ni)
       elseif(H.eq.'X') then
        read(line,'(t2,30f10.8)') (air(j),j=1,ni)
       elseif(H.eq.'M') then
        read(line,'(30I10)') (mas(j),j=1,ni)
       endif
      enddo
99    do i=1,ni
       if(psf(i).eq.0.d0) goto 109
      enddo
109   nimage=i-1
      do i=1,ns
       read(21,*,end=199) (IDF(i,j),j=1,nimage)
       read(21,*,end=199) (MAGF(i,j,1),j=1,nimage)
       read(21,*,end=199) (MAGF(i,j,2),j=1,nimage)
       read(21,*,end=199) (MAGF(i,j,3),j=1,nimage)
       read(21,*,end=199) (MAGF(i,j,4),j=1,nimage)
       read(21,*,end=199) (MAGF(i,j,5),j=1,nimage)
       read(21,*,end=199) (MAGF(i,j,6),j=1,nimage)
       read(21,*,end=199) (MAGF(i,j,7),j=1,nimage)
       read(21,*,end=199) (MAGF(i,j,8),j=1,nimage)
       do j=1,nimage
        if(idf(i,j).ne.0) then
         nstar(j)=nstar(j)+1 
        endif
       enddo
      enddo
199   nstars=i-1
      print *,nstars,' stars, ',nimage,' images read.'
      return
      end

      block data coo_chip
      common/chip/x0ccd(8),y0ccd(8),rot(8)
      data x0ccd/3.   , 2121., 4243., 6336.,    0., 2120., 4240., 6333./
      data y0ccd/2.   ,    7.,   4.5,    6., 4132., 4132., 4132., 4136./
      data rot  /0.    ,    0.,    0., .0015,    0.,    0.,    0.,   0./
c      data x0ccd/0.   ,   0.,   0.,      0.,   0.,   0.,   0.,   0./
c      data y0ccd/0.   ,   0.,   0.,      0.,   0.,   0.,   0.,   0./
      end
