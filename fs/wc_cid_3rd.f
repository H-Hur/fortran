c  Made for fortran95.
c  Weight combine part only.
c  Last modified on 12 Sep. 2009 bu Hyeon-Oh Hur.
c  gfortran -o wc_cid_3rd ~/work/fortran/subroutine/fs/wc_cid_3rd.f
c  21,51,62,63,64,65,66
      parameter (ns=100000,ni=30)
      integer nimage,nstars,IDF(ns,ni),mas(ni)
      character input(ni)*10,filter(ni)*10,inwc*16,fs*16
      real*8 MAGF(ns,ni,8),air(ni),ut(ni),psf(ni)
      print *,'input output from fs file=  ? ( 1 = out2.cid_3rd )'
      read(*,'(a16)') inwc
      if(inwc.eq.'1') write(inwc,'(a16)') 'out2.cid_3rd          '
      call read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstars)
      do i=1,nimage
       call readfilter(i,filter)     
      enddo
      open(65,file='dm.wc')
      write(65,'(t1,a45)')
     &'#ref  tar d_mag   Ncal calmag1 calmag2 Filter                   '
      close(65)
      call wccoo (nimage,nstars,MAGF,IDF)
      call wcombpar(nimage,mas,IDF,MAGF,air,ut,filter,nstars)
      call writewcomb(input,filter,ut,air,IDF,nstars)

      stop
      end

      subroutine wcombpar(ninput,mas,IDF,MAGF,air,ut,filter,nstars)
      parameter (ns=100000,ni=30)
      character filter(ni)*10
      integer ninput,IDF(ns,ni),nfnum(20,ni),n,nf(ni),nc,nstars
      integer mas(ni)
      real*8 MAGF(ns,ni,8),ut(ni),air(ni)             
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      common/wwc/nfilter
      common/dmm/master(ni),dm(ni,ni)
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
       call wcomb(nf(i),mas,IDF,MAGF,air,filter,nstars)
100    nc=nc
      enddo   
      return
      end 
 
      subroutine cal_dm(nf,nfimage,ims,IDF,MAGF,filter,nstars)
      parameter (ns=100000,ni=30)
      common/dmm/master(ni),dm(ni,ni)
      real*8 wsum,dmsum,mmax1,mmax2,mmin1,mmin2
      real*8 MAGF(ns,ni,8),wt,ch1,ch2,wt1,wt2,err,err_st,err_psf
      real dm2(ns),mag_mas(ns),ww(ns),avg
      integer nfimage,IDF(ns,ni),nstars,ims(ni),nf,it1,it2,nit,nc,nt(ns)
      character H*1,ps*25,filter(ni)*10
      open(64,file='dm.wc')
      do i=1,ns
       read(64,'(t1,a1)',end=9) H
      enddo
9     backspace(64) 
      write(64,'(t1,a1)') ' '
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
     &   MAGF(j,ims(i),5).lt.mmax2-1.5d0.and. 
     &   MAGF(j,master(nf),5).gt.mmin1+0.5d0.and.
     &   MAGF(j,master(nf),5).lt.mmax1-1.5d0) then
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
         wt1=0.0001d0/MAGF(j,ims(i),6)/MAGF(j,ims(i),6)
         wt2=0.0001d0/MAGF(j,master(nf),6)/MAGF(j,master(nf),6)
         if(wt1.ge.1.) wt1=1.d0
         if(wt2.ge.1.) wt2=1.d0
         avg = (MAGF(j,ims(i),5)*wt1+MAGF(j,master(nf),5)*wt2)/(wt1+wt2)
c         err_ st=sqrt((wt1*(MAGF(j,ims(i),5)-avg)**2+
c     &   wt2*(MAGF(j,master(nf),5)-avg)**2)/(wt1+wt2))
c         err_psf=sqrt(MAGF(j,master(nf),6)**2+MAGF(j,ims(i),6)**2)
         err=sqrt(MAGF(j,master(nf),6)**2+MAGF(j,ims(i),6)**2)
c         err=max(err_st,err_psf)
         wt=0.0001d0/err/err
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
       write(64,'(2i4,f9.5,i5,2f8.3,a12)') master(nf),ims(i),dm(nf,i),
     & nit,mmin2+0.5d0,mmax2-2.0d0,filter(ims(i))
       if(master(nf).ne.ims(i)) then 
        write(ps,'(t1,a7,i1,a1,i1,a3)')
     &  'delmag_',master(nf),'-',ims(i),'.wc                           '
        if(master(nf).ge.10) write(ps,'(t1,a7,i2,a1,i1,a3)')
     &  'delmag_',master(nf),'-',ims(i),'.wc                           '
        if(ims(i).ge.10) write(ps,'(t1,a7,i1,a1,i2,a3)')
     &  'delmag_',master(nf),'-',ims(i),'.wc                           '
        if(master(nf).ge.10.and.ims(i).ge.10)
     &  write(ps,'(t1,a7,i2,a1,i2,a3)')
     &  'delmag_',master(nf),'-',ims(i),'.wc                           '
        open(66,file=ps)
        write(66,'(a43)') '#mag_ref  mag_tar  mag_trs   dmag_weight   '
        do j=1,nc
         write(66,'(4f9.3)') mag_mas(j),dm2(j),ww(j)
        enddo
        close(66)
       endif  
      enddo
      write(64,'(t1,a1)') ' '
      close(64)
      return
      end

      subroutine wcomb(ninput,mas,IDF,MAGF,air,filter,nstars)
      parameter (ns=100000,ni=30)
      integer ninput,IDF(ns,ni),nstars,nfilt
      integer npoint,mas(ni)
      character filter(ni)*10
      real*8 MAGF(ns,ni,8),air(ni)
      real wsum,wmsum,wesum,wt,wer_er,wer_st
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      common/dmm/master(ni),dm(ni,ni)
      do i=1,ni
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
      call cal_dm(ninput,nfilt,ims,IDF,MAGF,filter,nstars)
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
        npoint=npoint+1
        wsum=wsum+wt
        wmsum=wmsum+(MAGF(i,ims(j),5)+dm(ninput,j))*wt
        wesum=wesum+MAGF(i,ims(j),6)*wt
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
         wer_st=wer_st+
     &   real(npoint)*wt*
     &   (MAGF(i,ims(j),5)+dm(ninput,j)-wcmag(ninput,i))**2!###check
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
       if(npoint.ne.0.and.wcmerr(ninput,i).le.0.001) then
c        write(*,'(i3,15f9.3)') npoint,
c     &  (MAGF(i,ims(j),5)+dm(ninput,j),MAGF(i,ims(j),6),j=1,nfilt)
c     &  ,wcmag(ninput,i),wesum,wsum,wcmerr(ninput,i)
       endif
      enddo
      return
      end
     
      subroutine wccoo (nimage,nstars,MAGF,IDF)
      parameter (ns=100000,ni=30)
      integer nimage,nstars,npoint,IDF(ns,ni)
      real*8 wt,wtsum,MAGF(ns,ni,8)
      common/rcoo/wxc(ns),wyc(ns),wxe(ns),wye(ns)
      do i=1,nstars
       wxc(i)=0.d0
       wyc(i)=0.d0
       wxe(i)=0.d0
       wye(i)=0.d0
       wtsum=0.d0
       npoint=0
       do j=1,nimage
        if(IDF(i,j).ne.0) then
         npoint=npoint+1
         if(MAGF(i,j,6).le.0.01d0) then
          wt=1.d0
         else
          wt=(0.01d0/MAGF(i,j,6))**2.d0
         endif
         wtsum=wtsum+wt
         wxc(i)=wxc(i)+wt*MAGF(i,j,3)
         wyc(i)=wyc(i)+wt*MAGF(i,j,4)    
        endif    
       enddo
       wxc(i)=wxc(i)/wtsum
       wyc(i)=wyc(i)/wtsum
       if(npoint.eq.1) then
        wxe(i)=-1.d0
        wye(i)=-1.d0
       else
        do j=1,nimage
         if(IDF(i,j).ne.0) then
          if(MAGF(i,j,6).le.0.01d0) then
           wt=1.d0
          else
           wt=(0.01d0/MAGF(i,j,6))**2.d0
          endif
          wxe(i)=wxe(i)+(wxc(i)-MAGF(i,j,3))**2.d0
          wye(i)=wye(i)+(wyc(i)-MAGF(i,j,4))**2.d0
         endif      
        enddo 
        wxe(i)=sqrt(wxe(i)/real(npoint)/wtsum)
        wye(i)=sqrt(wye(i)/real(npoint)/wtsum)
       endif
      enddo

      return
      end

      subroutine writewcomb(input,filter,ut,air,IDF,nstars)
      parameter(ns=100000,ni=30)
      real*8 air(ni),ut(ni)
      character filter(ni)*10,input(ni)*10
      integer IDF(ns,ni),nid,nstars
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),nobs(ni,ns),
     &nfc(ni,2),IDW(ni,ns),ims(ni)
      common/rcoo/wxc(ns),wyc(ns),wxe(ns),wye(ns)
      common/wwc/nfilter
      common/dmm/master(ni),dm(ni,ni)
      open(62,file='wcout1.wc')
      open(63,file='wcout2.wc')
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
      write(63,*) ""
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
       write(63,*) ""
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


      subroutine read_out_fs(inwc,nimage,input,filter,IDF,MAGF,mas,air,
     &psf,ut,nstars)
      parameter (ns=100000,ni=30)
      integer nimage,nstars,IDF(ns,ni),mas(ni)
      character input(ni)*10,filter(ni)*10,inwc*16,H*1,line*400
      real*8 MAGF(ns,ni,8),air(ni),psf(ni),ut(ni)
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
       endif
      enddo
c99    print *,'       nimage      master '
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
       if(IDF(i,2).ne.0) 
     & write(*,'(3f9.3)') MAGF(i,2,1),MAGF(i,2,2),MAGF(i,2,5)
      enddo
199   nstars=i-1
      print *,nstars,' stars, ',nimage,' images read.'
      return
      end
