c  Mad for fortran95.
c  Weight combine part only.
c  Last modified on 12 Sep. 2009 bu Hyeon-Oh Hur.
c  gfortran -o wc_cid_3rd_trans_to_master ~/work/fortran/subroutine/fs/wc_cid_3rd_trans_to_master.f ~/work/fortran/subroutine/caltools/std_excl.f
c  21,22,51,52,53,61,62,63,64,65,66,70---
      parameter (ns=100000,ni=30)
      integer nimage,nstars,IDF(ns,ni),mas(ni),n_apcor(ni)
      character input(ni)*10,filter(ni)*10,inwc*32,fs*32
      real*8 MAGF(ns,ni,8),air(ni),ut(ni),psf(ni)
      real minmax(ni,2)
      print *,'input output from fs file=  ? ( 1 = out2.cid_3rd )'
      read(*,'(a32)') inwc
      if(inwc.eq.'1') write(inwc,'(a32)') 'out2.cid_3rd                
     &                                                                 '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,n_apcor,
     & air,psf,ut,minmax,nstars)
      call read_dub()
      do i=1,nimage
       call readfilter(i,filter)     
      enddo
      open(65,file='dm.wc')
      write(65,'(t1,a45)')
     &'#ref  tar d_mag   Ncal calmag1 calmag2 Filter                   '
      close(65)
      call wccoo (nimage,nstars,psf,MAGF,IDF)
      call wcombpar(nimage,mas,IDF,MAGF,n_apcor,air,ut,filter,minmax,
     &nstars)
      call writewcomb(input,filter,ut,air,IDF,nstars)

      stop
      end

      subroutine wcombpar(ninput,mas,IDF,MAGF,n_apcor,
     &air,ut,filter,minmax,nstars)
      parameter (ns=100000,ni=30)
      character filter(ni)*10
      integer ninput,IDF(ns,ni),nfnum(20,ni),n,nf(ni),nc,nstars
      integer mas(ni),n_apcor(ni)
      real*8 MAGF(ns,ni,8),ut(ni),air(ni)
      real minmax(ni,2)             
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      common/wwc/nfilter
      common/dmm/master(ni),dm(ni,ni),num_im(ni),
     &xdm(ns,ni),ydm(ns,ni),adm(ns,ni),wdm(ns,ni),ncal(ni)
      do i=1,ni
       ims(i)=0
       nf(i)=0
       do j=1,20
        nfnum(j,i)=0
       enddo
      enddo    
      nfilter=0
      do j=1,20 
       n=0
       do i=1,ninput
        if(nfc(i,1).eq.j) then
         n=n+1
         nfnum(j,n)=i
        endif
       enddo
       if(nfnum(j,1).ne.0) then
        nfilter=nfilter+1
       endif
      enddo
      print*,nfilter,' filters found.'
      nc=0
      do i=1,10       
       do j=1,ni
        ims(j)=0
       enddo

       n=0
       do j=1,ninput
        if(nfnum(i,j).ne.0) then
         n=n+1
         ims(n)=nfnum(i,j)
         if(n.eq.1) then
          nc=nc+1
         endif
        endif
       enddo
       nf(i)=nc
       if(n.eq.0) then
        nf(i)=0
        goto 100
       endif
       call wcomb(nf(i),mas,IDF,MAGF,n_apcor,air,filter,minmax,nstars)
100    nc=nc
      enddo   
      return
      end 
 
      subroutine wcomb(ninput,mas,IDF,MAGF,n_apcor,air,filter,minmax,
     &nstars)
      parameter (ns=100000,ni=30)
      integer ninput,IDF(ns,ni),nstars,nfilt,nq
      integer npoint,mas(ni),n_apcor(ni),numb(ni)
      character filter(ni)*10,output*32
      real*8 MAGF(ns,ni,8),air(ni)
      real wsum,wmsum,wesum,wt,wer_er,wer_st,refmin,refmax,tarmax,tarmin
      real minmax(ni,2)
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      common/dmm/master(ni),dm(ni,ni),num_im(ni),
     &xdm(ns,ni),ydm(ns,ni),adm(ns,ni),wdm(ns,ni),ncal(ni)
      nnnn=70
      do i=1,ni
       ncal(i)=0
       if(ims(i).eq.0) then
        goto 100
       endif
      enddo
100   nfilt=i-1
      do i=1,ni
       do j=1,ni
        dm(i,j)=0.
       enddo
      enddo
      call selectmaster(ninput,mas,master,air,nstars)
      call cal_dm(ninput,nfilt,ims,IDF,MAGF,nstars)
      do i=1,nfilt 
       numb(i)=0
       write(output,'(a4,i1,a1,i1,a7)')
     & 'dmag_',master(ninput),'_',ims(i),'.wc_trs'
       if(master(ninput).ne.ims(i)) then
        nnnn=nnnn+1
        numb(i)=nnnn
        open(numb(i),file=output)
        write(numb(i),'(t1,a98)') 
     &  '#   xc_im     yc_im mag-mag_ref    wt #im #ref mag_im  mag_ref 
     & mag_trs   err_im  err_ref  tar-ref                              '
       endif
      enddo
      do i=1,nstars
       do j=1,nfilt 
        nq=ims(j) 
        if(MAGF(i,nq,6).gt.90.d0.or.MAGF(i,nq,5).gt.90.d0.or.
     &  MAGF(i,master(ninput),6).gt.90.d0.or.MAGF(i,master(ninput),5)
     &  .gt.90.d0.or.IDF(i,master(ninput)).eq.0.or.IDF(i,nq).eq.0
     &  .or.MAGF(i,nq,5).eq.0.or.MAGF(i,master(ninput),5).eq.0)
     &  then
        else
         wt=(0.005/MAGF(i,nq,6))**2.d0 + 
     &   (0.005/MAGF(i,master(ninput),6))**2.d0
         if(wt.gt.1.) wt=1.
         ncal(nq)=ncal(nq)+1
         wdm(ncal(nq),nq)=wt
         xdm(ncal(nq),nq)=MAGF(i,nq,1)
         ydm(ncal(nq),nq)=MAGF(i,nq,2)
         adm(ncal(nq),nq)=MAGF(i,master(ninput),5)-MAGF(i,nq,5) 
        endif
       enddo
      enddo

      do i=1,nstars
       do j=1,nfilt
        nq=ims(j)
        if(MAGF(i,nq,6).gt.90.d0.or.MAGF(i,nq,5).gt.90.d0.or.
     &  MAGF(i,master(ninput),6).gt.90.d0.or.MAGF(i,master(ninput),5)
     &  .gt.90.d0.or.IDF(i,master(ninput)).eq.0.or.IDF(i,nq).eq.0
     &  .or.MAGF(i,nq,5).eq.0.or.MAGF(i,master(ninput),5).eq.0)
     &  then
        else
         wt=(0.005/MAGF(i,nq,6))**2.d0 +
     &   (0.005/MAGF(i,master(ninput),6))**2.d0
         call cal_dm_wm(MAGF(i,ims(j),1),MAGF(i,ims(j),2),ims(j),minmax
     &   ,stdd,dmag)
         if(numb(j).ne.0) write(numb(j),'(t1,2f10.3,2f9.3,2i3,6f9.3)')
     &   MAGF(i,nq,1),MAGF(i,nq,2),
     &   MAGF(i,master(ninput),5)-MAGF(i,nq,5),wt
     &   ,nq,master(ninput)
     &   ,MAGF(i,nq,5),MAGF(i,master(ninput),5),MAGF(i,nq,5)+dmag
     &   ,MAGF(i,nq,6),MAGF(i,master(ninput),6),
     &    MAGF(i,nq,5)+dmag-MAGF(i,master(ninput),5)
        endif
       enddo
      enddo

      do i=71,nnnn
       close(nnnn)
      enddo

      do i=1,nstars
       wsum=0.
       wmsum=0.
       wesum=0.
       npoint=0
       wer_er=0.
       wer_st=0.
       do j=1,nfilt
        if(MAGF(i,ims(j),6).le.0.005d0) then
         wt=1.
         elseif(MAGF(i,ims(j),6).gt.90.d0.or.
     &   MAGF(i,ims(j),5).gt.90.d0) then
         wt=0.
         elseif(MAGF(i,ims(j),6).gt.0.005d0) then
         wt=(0.005/MAGF(i,ims(j),6))**2.d0
        endif
        if(MAGF(i,ims(j),5).eq.0) then
         wt=0.
        else
         if(n_apcor(ims(j)).eq.0) then
          call cal_dm_wm(MAGF(i,ims(j),1),MAGF(i,ims(j),2),ims(j),minmax
     &   ,stdd,dmag)
          wmsum=wmsum+(MAGF(i,ims(j),5)+dmag)*wt
         elseif(n_apcor(ims(j)).eq.1) then
          wmsum=wmsum+(MAGF(i,ims(j),5)+dm(ninput,j))*wt
         endif
         wesum=wesum+MAGF(i,ims(j),6)*wt
         npoint=npoint+1
         wsum=wsum+wt
         if(wt.ne.0.)
     &   wer_er=wer_er+1./(MAGF(i,ims(j),6)**2)
        endif
       enddo
       nobs(ninput,i)=npoint
       if(wsum.gt.0.) then
        IDW(ninput,i)=IDF(i,master(ninput))
        wcmag(ninput,i)=wmsum/wsum
        do j=1,nfilt
         if(MAGF(i,ims(j),6).le.0.005d0) then
          wt=1.
         elseif(MAGF(i,ims(j),6).gt.90.d0.or.
     &    MAGF(i,ims(j),5).gt.90.d0) then
          wt=0.
         elseif(MAGF(i,ims(j),6).gt.0.005d0) then
          wt=(0.005/MAGF(i,ims(j),6))**2.d0
         endif
         if(MAGF(i,ims(j),5).eq.0) then
          wt=0.
         endif
         if(n_apcor(ims(j)).eq.0) then
          call cal_dm_wm(MAGF(i,ims(j),1),MAGF(i,ims(j),2),ims(j),minmax
     &    ,stdd,dmag)
          wer_st=wer_st+
     &    wt*(MAGF(i,ims(j),5)+dmag-wcmag(ninput,i))**2
         elseif(n_apcor(ims(j)).eq.1) then
          wer_st=wer_st+
     &    wt*(MAGF(i,ims(j),5)+dm(ninput,j)-wcmag(ninput,i))**2
         endif
         if(wt.ne.0.)
     &   wer_er=wer_er+1./(MAGF(i,ims(j),6)**2)
        enddo
        wer_er=sqrt(1./wer_er)
        if(npoint.ge.2) then
         wer_st=sqrt(wer_st/real(npoint-1)/wsum)
        elseif(npoint.eq.1) then
         wer_st=wesum/wsum
        endif
        wcmerr(ninput,i)=max(wer_st,wer_er)
       elseif (wsum.eq.0.) then
        IDW(ninput,i)=0.
        wcmag(ninput,i)=0.
        wcmerr(ninput,i)=0.
       endif
      enddo
      return
      end

      subroutine cal_dm_wm(xc,yc,imnum,minmax,std,dmag_out)
      parameter (ns=100000,ni=30)
      common/dmm/master(ni),dm(ni,ni),num_im(ni),
     &xdm(ns,ni),ydm(ns,ni),adm(ns,ni),wdm(ns,ni),ncal(ni)
      real dmag_out,wsum,rad,std,th,minmax(ni,2)
      real*8 xc,yc,av,stdd,w(ns),dmag(ns),dist,w_dist,wconst
      integer imnum,n,nmin,nc
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      wsum=0.
      dmag_out=0.
      n=0
      do i=1,ncal(imnum)
       dist=sqrt((xc-xdm(i,imnum))**2+(yc-ydm(i,imnum))**2)
       n=n+1
       w_dist=exp(-dist/250.d0)/exp(-100.d0/250.d0)
       if(w_dist.gt.1.d0) w_dist=1.d0
       if(adm(i,imnum).lt.minmax(imnum,1)+0.5.and.
     & adm(i,imnum).gt.minmax(imnum,2)) w_dist=0.d0
       w(n)=wdm(i,imnum)*w_dist
       dmag(n)=adm(i,imnum)
      enddo
      call wstd_dp(n,dmag,w,av,stdd)
      std=real(stdd)
      dmag_out=real(av)
200   return
      end


      subroutine cal_dm(nf,nfimage,ims,IDF,MAGF,nstars)
      parameter (ns=100000,ni=30)
      common/dmm/master(ni),dm(ni,ni),num_im(ni),
     &xdm(ns,ni),ydm(ns,ni),adm(ns,ni),wdm(ns,ni),ncal(ni)
      real*8 wsum,dmsum,mmax1,mmax2,mmin1,mmin2
      real*8 MAGF(ns,ni,8),wt,ch1,ch2,wt1,wt2,err,err_st,err_psf
      real x1,x2,dm2(ns),mag_mas(ns),ww(ns),avg
      integer nfimage,IDF(ns,ni),nstars,ims(ni),nf,it1,it2,nit,nc,nt(ns)
      character H*1,ps*25,nn*5
      open(64,file='dm.wc')
      do i=1,ns
       read(64,'(a1)',end=9) H
      enddo
9     backspace(64)
      write(64,*) "   "
      do i=1,nstars
       if(IDF(i,master(nf)).ne.0) then
        mmax1=MAGF(i,master(nf),5)
        mmin1=MAGF(i,master(nf),5)
        goto 10
       endif
      enddo
10    do i=1,nstars
       if(MAGF(i,master(nf),5).gt.mmax1) then
        mmax1=MAGF(i,master(nf),5)
       endif
       if(MAGF(i,master(nf),5).lt.mmin1.and.
     & MAGF(i,master(nf),5).ne.0.d0)then
        mmin1=MAGF(i,master(nf),5)
       endif
      enddo

      do i=1,nfimage
       nit=0
       nc=0
       wt=0.d0
       wsum=0.d0
       dmsum=0.d0
       do j=1,nstars
        if(IDF(j,ims(i)).ne.0) then
         mmax2=MAGF(j,ims(i),5)
         mmin2=MAGF(j,ims(i),5)
         goto 50
        endif
       enddo
50     do j=1,nstars
        if(MAGF(j,ims(i),5).gt.mmax2) then
         mmax2=MAGF(j,ims(i),5)
        endif
        if(MAGF(j,ims(i),5).lt.mmin2.and.MAGF(j,ims(i),5).ne.0.d0) then
         mmin2=MAGF(j,ims(i),5)
        endif
       enddo
51    format(I2,I2,4f8.4)
       do j=1,nstars
        it1=0
        it2=0
        if(IDF(j,ims(i)).ne.0.and.IDF(j,master(nf)).ne.0.and.
     &   MAGF(j,ims(i),5).gt.mmin2+0.5d0.and.
     &   MAGF(j,ims(i),5).lt.mmax2-2.5d0.and.
     &   MAGF(j,master(nf),5).gt.mmin1+0.5d0.and.
     &   MAGF(j,master(nf),5).lt.mmax1-2.5d0) then
         nit=nit+1
         nc=nc+1
         nt(nc)=1
         dm2(nc)=MAGF(j,master(nf),5)-MAGF(j,ims(i),5)
         mag_mas(nc)=MAGF(j,master(nf),5)
         if(MAGF(j,ims(i),6).lt.0.005d0) then
          it1=1
          ch1=MAGF(j,ims(i),6)
          MAGF(j,ims(i),6)=0.005d0
         endif
         if(MAGF(j,master(nf),6).lt.0.005d0) then
          it2=1
          ch2=MAGF(j,master(nf),6)
          MAGF(j,master(nf),6)=0.005d0
         endif
c         wt=(15.d0*15.d0/(MAGF(j,ims(i),5)*MAGF(j,master(nf),5)))
c     &   *0.5*sqrt((0.005d0**2.d0+0.005d0**2.d0)/
c     &   (MAGF(j,ims(i),6)**2.d0+MAGF(j,master(nf),6)**2.d0))
         wt1=0.005d0*0.005d0/MAGF(j,ims(i),6)/MAGF(j,ims(i),6)
         wt2=0.005d0*0.005d0/MAGF(j,master(nf),6)/MAGF(j,master(nf),6)
         if(wt1.ge.1.d0) wt1=1.d0
         if(wt2.ge.1.d0) wt2=1.d0
         avg = (MAGF(j,ims(i),5)*wt1+MAGF(j,master(nf),5)*wt2)/(wt1+wt2)
         err_ st=sqrt((wt1*(MAGF(j,ims(i),5)-avg)**2+
     &   wt2*(MAGF(j,master(nf),5)-avg)**2)/(wt1+wt2))

         err_psf=sqrt(MAGF(j,master(nf),6)**2+MAGF(j,ims(i),6)**2)
         err=sqrt(MAGF(j,master(nf),6)**2+MAGF(j,ims(i),6)**2)
c         err=max(err_st,err_psf)
         wt=0.005d0*0.005d0/err/err
         if(wt.gt.1.d0) wt=1.d0 
         ww(nc)=wt
         wsum=wsum+wt
         dmsum=dmsum+wt*(MAGF(j,master(nf),5)-MAGF(j,ims(i),5))
         if(it1.eq.1) then
          MAGF(j,ims(i),6)=ch1
          it1=0
         endif
         if(it2.eq.1) then
          MAGF(j,master(nf),6)=ch2
          it2=0
         endif
        elseif(IDF(j,ims(i)).ne.0.and.IDF(j,master(nf)).ne.0) then
         nc=nc+1
         dm2(nc)=MAGF(j,master(nf),5)-MAGF(j,ims(i),5)
         nt(nc)=2
         mag_mas(nc)=MAGF(j,master(nf),5)
         ww(nc)=0.01
        endif
       enddo
       dm(nf,i)=dmsum/wsum
       write(64,'(i3,i3,f9.5,i5,2f8.3)') nf,i,dm(nf,i),nit,mmin2+0.5d0,
     & mmax2-2.5d0
      enddo
      write(64,'(a1)') ' '
      close(64)
c      if(dm(nf,i).ne.0.) then
c        print *,mag_mas(i),dm2(i)
c      endif
      return
      end



      subroutine wccoo (nimage,nstars,psf,MAGF,IDF)
      parameter (ns=100000,ni=30)
      integer nimage,nstars,npoint,IDF(ns,ni)
      real*8 wt,wtsum,MAGF(ns,ni,8),psf(ni),wt_f
      common/rcoo/wxc(ns),wyc(ns),wxe(ns),wye(ns),nwdouble(ns)
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      do i=1,nstars
       wxc(i)=0.d0
       wyc(i)=0.d0
       wxe(i)=0.d0
       wye(i)=0.d0
       wtsum=0.d0
       npoint=0
       nwdouble(i)=0
       do j=1,nimage
        wt_f=1.d0
        if(nfc(j,1).eq.6) wt_f=0.1d0 
        if(IDF(i,j).ne.0) then
         call cid_double(j,real(MAGF(i,j,1)),real(MAGF(i,j,2)),
     &   0.5*real(psf(j)),nwdouble(i))
         npoint=npoint+1
         wt=(0.005d0/MAGF(i,j,6))**2.d0
         if(wt.gt.1.d0) wt=1.d0
         wtsum=wtsum+wt*wt_f
         wxc(i)=wxc(i)+wt*wt_f*MAGF(i,j,3)
         wyc(i)=wyc(i)+wt*wt_f*MAGF(i,j,4)    
        endif    
       enddo
       wxc(i)=wxc(i)/wtsum
       wyc(i)=wyc(i)/wtsum
       if(npoint.eq.1) then
        wxe(i)=-1.d0
        wye(i)=-1.d0
       else
        do j=1,nimage
        wt_f=1.d0
        if(nfc(j,1).eq.6) wt_f=0.1d0
         if(IDF(i,j).ne.0) then
          wt=(0.005d0/MAGF(i,j,6))**2.d0
          if(wt.gt.1.d0) wt=1.d0
          wxe(i)=wxe(i)+wt*wt_f*(wxc(i)-MAGF(i,j,3))**2.d0
          wye(i)=wye(i)+wt*wt_f*(wyc(i)-MAGF(i,j,4))**2.d0
         endif      
        enddo 
        wxe(i)=sqrt(wxe(i)/(real(npoint-1)*wtsum/real(npoint)))
        wye(i)=sqrt(wye(i)/(real(npoint-1)*wtsum/real(npoint)))
       endif
      enddo
      return
      end

      subroutine cid_double(imnum,xc,yc,rad,ndouble)
      parameter (ns=100000,ni=30)
      common/double/ndub(ni),imnum_d(ni),xc_d(ni,ns),yc_d(ni,ns),
     &dub,nim_dub
      real xc,yc,rad,dist
      integer imnum,ndouble
      character dub(ni)*32
      do j=1,nim_dub
       if(imnum.eq.imnum_d(j)) then
        do i=1,ndub(j)
         dist=sqrt((xc-xc_d(j,i))**2+(yc-yc_d(j,i))**2)
         if(dist.le.rad) then
          ndouble=1
          goto 900
         endif
        enddo
       endif
      enddo
900   return
      end 

      subroutine writewcomb(input,filter,ut,air,IDF,nstars)
      parameter(ns=100000,ni=30)
      real*8 air(ni),ut(ni)
      character filter(ni)*10,input(ni)*10
      integer IDF(ns,ni),nid,nstars
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      common/rcoo/wxc(ns),wyc(ns),wxe(ns),wye(ns),nwdouble(ns)
      common/wwc/nfilter
      common/dmm/master(ni),dm(ni,ni),num_im(ni),
     &xdm(ns,ni),ydm(ns,ni),adm(ns,ni),wdm(ns,ni),ncal(ni)
      open(61,file='wcout0.wc')
      open(62,file='wcout1.wc')
      open(63,file='wcout2.wc')
c     0st file header##################
      write(61,'(t1,a13)') "# F : Filrer               "
      write(61,'(t1,a13)') "# T : UT                   "
      write(61,'(t1,a13)') "# X : Airmass              "
      write(61,'(t1,a1,a6,t26,30a14)')
     &'F','',(filter(master(i)),i=1,nfilter)
      write(61,'(t1,a1,t25,30f14.8)')'T',(ut(master(i)),i=1,nfilter)
      write(61,'(t1,a1,t25,30f14.7)')'X',(air(master(i)),i=1,nfilter)
      write(61,'(t1,a24,t25,30a14)')'#   ID   xc       yc    ',
     &(' mag     merr  ',i=1,nfilter)
c     0st file header##################
c     1st file header##################
      write(62,*) " F : Filrer"
      write(62,*) " T : UT"
      write(62,*) " X : Airmass"
c     1st file header##################
c     2nd file header##################
      write(63,*) ""
      write(63,*) " M : Master image"
      write(63,*) " F : Filrer"
      write(63,*) " T : UT"
      write(63,*) " X : Airmass"
      write(63,*) " *"
      write(63,*) ' ID(of all detected stars),   ID(each master filter)'
      write(63,*) ' xc(on the Xid master image), Mag(combined)'
      write(63,*) ' yc(on the Xid master image), Merr(combined error)'
      write(63,*) ' xcerror(combined),           nobs (each filter)'
      write(63,*) ' ycerror(combined)'
      write(63,*) ' Double (0=single,1=dounle)                     '
      write(63,*) ""
      write(63,*) ""
      write(63,*) ""
      write(63,*) ""
c     2nd file header##################
      write(62,'(t1,a1,a6,30a10)')'F','',(filter(master(i)),i=1,nfilter)
      write(62,'(t1,a1,30a10)')   'M',(input(master(i)),i=1,nfilter)
      write(62,'(t1,a1,30f10.6)')'T',(ut(master(i)),i=1,nfilter)
      write(62,'(t1,a1,30f10.7)')'X',(air(master(i)),i=1,nfilter)
      write(62,'(t1,a1)') "*"
      write(63,'(t1,a1,t11,a6,30a10)')'F',''
     *,(filter(master(i)),i=1,nfilter)
      write(63,'(t1,a1,t11,30a10)')   'M',(input(master(i)),i=1,nfilter)
      write(63,'(t1,a1,t11,30f10.6)')'T',(ut(master(i)),i=1,nfilter)
      write(63,'(t1,a1,t11,30f10.7)')'X',(air(master(i)),i=1,nfilter)
      write(63,'(t1,a1)') "*"   

      do i=1,nstars
       nid=0
       do j=1,nfilter
        if(wcmag(j,i).eq.0) then
         nid=nid+1
        endif
       enddo
       if(nid.eq.nfilter) then
        goto 500
       endif
       write(61,'(t1,i6,2f9.3,20f7.3)') i,wxc(i),wyc(i),
     & (wcmag(j,i),wcmerr(j,i),j=1,nfilter)
       write(62,'(30i10)')   (IDF(i,master(j)),j=1,nfilter)
       write(62,'(30I10)')   (IDW(j,i),j=1,nfilter)
       write(62,'(30f10.3)') (wcmag(j,i),j=1,nfilter)
       write(62,'(30f10.3)') (wcmerr(j,i),j=1,nfilter)
       write(62,*) ""
       write(63,'(t1,31i10)')   i,(IDF(i,master(j)),j=1,nfilter)
       write(63,'(t1,31f10.3)') wxc(i),(wcmag(j,i),j=1,nfilter)
       write(63,'(t1,31f10.3)') wyc(i),(wcmerr(j,i),j=1,nfilter)
       write(63,'(t1,f10.3,t11,30i10)') wxe(i),(nobs(j,i),j=1,nfilter)
       write(63,'(t1,f10.3)') wye(i)
       write(63,'(t1,i10)') nwdouble(i)
       write(63,*) ""
       write(63,*) ""
       write(63,*) ""
       write(63,*) ""
      enddo
500   close(63)
      close(62)
      return
      end

      subroutine selectmaster(nf,mas,master,air,nstars)
      parameter(ns=100000,ni=30)
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      real*8 air(ni),psf(ni)
      integer nimage,nstars,nf,master(ni),mas(ni)     
      do i=1,ni
       if(ims(i).eq.0) then
        goto 100
       endif
      enddo
100   nimage=i-1
      master(nf)=ims(1)
      do i=2,nimage
       if(psf(ims(i)).gt.psf(ims(i-1))) then
        master(nf)=ims(i)
        elseif(psf(ims(i)).lt.psf(ims(i-1))) then
        master(nf)=ims(i-1)
       endif
      enddo
      do i=1,nimage
       if(mas(ims(i)).eq.1) then
        master(nf)=ims(i)
       endif
      enddo
      return
      end

      subroutine readfilter(nimage,filter)
      parameter(ns=100000,ni=30)
      character H*1,filter(ni)*10,ch(10)*1,tf*1
      integer nimage,nfi,nco,nf
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      open(51,file='filterpar.wc_cid_3rd')
      nf=0
      do i=nimage,ni
       nfc(i,1)=0
       nfc(i,2)=0
      enddo
      do i=1,20
       read(51,'(t1,a1)',end=500) H
       if(H.eq.'*') then
        nf=nf+1
        backspace(51)
        read(51,'(t2,a1,t4,i5,t9,i16)',end=500) tf,nfi,nco
        read(filter(nimage),'(10a1)') (ch(j),j=1,10)
        do j=1,10
         if     (ch(j).eq.tf) then
          nfc(nimage,1)=nfi
          nfc(nimage,2)=nco
         endif
        enddo
       endif
      enddo
500   close(51)
      return
      end


      subroutine read_out_fs(inwc,nimage,input,filter,IDF,MAGF,
     &mas,n_apcor,air,psf,ut,minmax,nstars)
      parameter (ns=100000,ni=30)
      integer nimage,nstars,IDF(ns,ni),mas(ni),n_apcor(ni)
      character input(ni)*10,filter(ni)*10,inwc*32,H*1,line*400
      real*8 MAGF(ns,ni,8),air(ni),psf(ni),ut(ni)
      real minmax(ni,2)
      open(21,file=inwc)
      read(21,'(t1,30a10)') (input(j),j=1,ni)
      do i=1,ni
       psf(i)=0.d0
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
       elseif(H.eq.'A') then
        read(line,'(30I10)') (n_apcor(j),j=1,ni)
       endif
      enddo
c99    print *,'       nimage      master '
99    do i=1,ni
       minmax(i,1)=90.
       minmax(i,2)=0.
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
       if(idf(i,j).ne.0.and.magf(i,j,5).lt.minmax(j,1)) 
     & minmax(j,1)=magf(i,j,5)
      enddo
199   nstars=i-1
      do i=1,nimage
       minmax(i,2)=minmax(i,1)-8.
      enddo
      print *,nstars,' stars, ',nimage,' images read.'
      return
      end
      subroutine read_dub()
      parameter (ns=100000,ni=30)
      common/double/ndub(ni),imnum_d(ni),xc_d(ni,ns),yc_d(ni,ns),
     &dub,nim_dub
      character H*1,dub(ni)*32
      open(52,file='input.dub.wc')
      read(52,'(t1,a1)') H
      do i=1,ni
       read(52,*,end=99) imnum_d(i),dub(i)
       open(53,file=dub(i))
       do j=1,ns
        read(53,*,end=89) xc_d(i,j),yc_d(i,j)
       enddo
89     close(53) 
       ndub(i)=j-1
      enddo
99    nim_dub=i-1
      close(52)
      return
      end

