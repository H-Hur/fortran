c     gfortran -o aex_lint_mac ~/work/fortran/subroutine/aex/aex_lint/aex_lint_mac.f  
c     21,22,23,25,31,32,42
      character input*16,output1*16,output2*16,coeff*16
      call read_ZAMS()
      write(coeff,'(a16)') 'coeff.input                                '
      call read_coeff(coeff)
      write(input,'(a16)') 'wcout2.wc                                  '
      call readwchead(input)
      write(output1,'(a16)') 'out1.aex                                 '
      write(output2,'(a16)') 'out2.aex                                 '
      call write_output(input,coeff,output1,output2)
      stop
      end

      subroutine write_output(input,coeff,output1,output2)
      parameter(nf=30,ns=100000)
      common /zams/zMv(41),zVI(41),zBV(41),zUB(41)
      common /ncw/nc_U,nc_B,nc_V,nc_R,nc_I,nc_Ha,nc_K,nfilter,nc_H,nc_J
     &,nstars
      common/temp_col/temp_uub,temp_bbv,temp_vbv,temp_vvi,temp_ivi,
     &temp_iri,dt_uub,dt_bbv,dt_vbv,dt_vvi,dt_ivi,dt_iri,
     &nt_uub,nt_bbv,nt_vbv,nt_vvi,nt_ivi,nt_iri
      common /wcdata/ID(ns),xc(ns),yc(ns),xe(ns),ye(ns),nobs(nf,ns),
     &cmag(nf,ns),cmerr(nf,ns),id_f(nf,ns),ndouble(ns),
     &airmass(nf),UT(nf),filter
      common/coeff_const/num_filter(nf),num_color(nf),aeK1(nf),aeK2(nf)
     &,zero(nf),xcenter(nf),ycenter(nf),color_dp,ut_dp,x_dp,y_dp,r_dp,
     &theta_dp,add1_dp,add2_dp,add3_dp,nconst
      character color_dp(nf)*16,ut_dp(nf)*16,x_dp(nf)*16,y_dp(nf)*16
     &,r_dp(nf)*16,theta_dp(nf)*16,add1_dp(nf)*16,add2_dp(nf)*16
     &,add3_dp(nf)*16
      character filter(nf)*10,input*16,coeff*16,output2*16,double*1
      character output1*16
      real o_U,o_B,o_V,o_Vvi,o_Vbv,o_r,o_I,o_Iri,o_Ivi,o_Ha,mean_iter
      real x_U,x_B,x_Vvi,x_Vbv,x_r,x_Iri,x_Ivi,x_Ha,std_iter
      real UB,BV,VR,VI,RI,e_ub,e_bv,e_vi,e_ri,wt_bv,wt_vi,wt_ri
      real o_Icol,e_V,e_temp_high,e_temp_low
      integer nUB,nBV,nVR,nVI,nRI,nback,niter,nid,nit(ns),ncal_V
      integer ncal_U,ncal_B,ncal_R,ncal_vvi,ncal_vbv,ncal_iri,ncal_ivi
      integer ncal_I,nV
      data ncal_U,ncal_B,ncal_V,ncal_R,ncal_I,ncal_vvi,ncal_vbv,ncal_iri
     &,ncal_ivi,ncal_Ha,nid,nmax,mean_iter/0,0,0,0,0,0,0,0,0,0,0,0,0/
      open(41,file=output1)
      write(41,'(t1,a10,a16)') '# input = ',input
      write(41,'(t1,a10,a16)') '# coeff = ',coeff
      write(41,'(t1,a7,t11,30a10)')'#Filter',(filter(i),i=1,nfilter)
      write(41,'(t1,a5,t11,30f10.6)')'#Time',(UT(i),i=1,nfilter)
      write(41,'(t1,a8,t11,30f10.7)')'#Airmass',(airmass(i),i=1,nfilter)
      write(41,'(t1,a140)')  '#rawID    ID    Xc       Yc      U       B
     &       V       R      I        Ha     e_U    e_B    e_V    e_R    
     &e_I    e_Ha   nobs------------ D                                 '
      open(42,file=output2)
      write(42,'(t1,a10,a16)') '# input = ',input
      write(42,'(t1,a10,a16)') '# coeff = ',coeff
      write(42,'(t1,a7,t11,30a10)')'#Filter',(filter(i),i=1,nfilter)
      write(42,'(t1,a5,t11,30f10.6)')'#Time',(UT(i),i=1,nfilter)
      write(42,'(t1,a8,t11,30f10.7)')'#Airmass',(airmass(i),i=1,nfilter)
      write(42,'(t1,a153)')  '#rawID    ID    Xc       Yc       V       
     &I      U-B    B-V    V-I    R-I    Ha     eV     eI    e(U-B) e(B-
     &V) e(V-I) e(R-I) e(Ha) nobs--------------- D                     '
50    format(t1,2i6,2f9.3,2f8.3,12f7.3,7i3,1x,a1)
51    format(t1,2i6,2f9.3,6f8.3,6f7.3,6i3,1x,a1)
      nid=0
      nmax=0
      do i=1,nconst
       if(num_filter(i).eq.1.and.num_color(i).eq.1) ncal_U=1
       if(num_filter(i).eq.2.and.num_color(i).eq.2) ncal_B=1
       if(num_filter(i).eq.4.and.num_color(i).eq.5) ncal_R=1
       if(num_filter(i).eq.3.and.num_color(i).eq.2) ncal_vbv=1
       if(num_filter(i).eq.3.and.num_color(i).eq.4) ncal_vvi=1
       if(ncal_vbv.ne.0.or.ncal_vvi.ne.0) ncal_V=1
       if(num_filter(i).eq.5.and.num_color(i).eq.4) ncal_ivi=1
       if(num_filter(i).eq.5.and.num_color(i).eq.5) ncal_iri=1
       if(ncal_ivi.ne.0.or.ncal_iri.ne.0) ncal_I=1
       if(num_filter(i).eq.4.and.num_color(i).eq.5) ncal_rri=1
       if(num_filter(i).eq.6) ncal_Ha=1
      enddo
      do i=1,nstars
       e_ub=sqrt(cmerr(nc_U,i)**2+cmerr(nc_B,i)**2)
       e_bv=sqrt(cmerr(nc_B,i)**2+cmerr(nc_V,i)**2)
       e_vr=sqrt(cmerr(nc_V,i)**2+cmerr(nc_R,i)**2)
       e_vi=sqrt(cmerr(nc_V,i)**2+cmerr(nc_I,i)**2)
       e_ri=sqrt(cmerr(nc_R,i)**2+cmerr(nc_I,i)**2)
       nUB=min(nobs(nc_U,i),(nobs(nc_B,i)))
       nBV=min(nobs(nc_B,i),(nobs(nc_V,i)))
       nVR=min(nobs(nc_V,i),(nobs(nc_R,i)))
       nVI=min(nobs(nc_V,i),(nobs(nc_I,i)))
       nRI=min(nobs(nc_R,i),(nobs(nc_I,i)))
       nV=max(nBV,nVI)
       if(nUB.eq.0) e_ub=0.
       if(nBV.eq.0) e_bv=0.
       if(nVR.eq.0) e_ur=0.
       if(nVI.eq.0) e_vi=0.
       if(nRI.eq.0) e_ri=0.
       if(nV.eq.0) e_v=0.
       o_U=cmag(nc_U,i)
       o_B=cmag(nc_B,i)
       o_Vvi=cmag(nc_V,i)
       o_Vbv=cmag(nc_V,i)
       o_R=cmag(nc_R,i)
       o_Ivi=cmag(nc_I,i)
       o_Iri=cmag(nc_I,i)
       o_Ha=cmag(nc_Ha,i)
       UB=cmag(nc_U,i)-cmag(nc_B,i)
       BV=cmag(nc_B,i)-cmag(nc_V,i)
       VR=cmag(nc_V,i)-cmag(nc_R,i)
       VI=cmag(nc_V,i)-cmag(nc_I,i)
       RI=cmag(nc_R,i)-cmag(nc_I,i)
       niter=0
101    nback=0
       x_U=o_U
       x_B=o_B
       x_Vvi=o_Vvi
       x_Vbv=o_Vbv
       x_R=o_R
       x_Ivi=o_Ivi
       x_Iri=o_Iri
       x_Ha=o_Ha
       niter=niter+1
       do j=1,nconst
c       U : U-B
        if(nobs(nc_U,i).ne.0.and.nobs(nc_B,i).ne.0.and.
     &  num_filter(j).eq.1.and.num_filter(j).eq.1.and.ncal_U.ne.0) then
         call standard(1,1,cmag(nc_U,i),UB,id(i),x_Vbv,
     &   UB,BV,VR,VI,RI,airmass(nc_U),UT(nc_U),xc(i),yc(i),o_U)
         if(abs(o_U-x_U).ge.0.00001) nback=nback+1        
        endif
c       B : B-V
        if(nobs(nc_B,i).ne.0.and.nobs(nc_V,i).ne.0.and.ncal_B.ne.0.and.
     &  num_filter(j).eq.2.and.num_color(j).eq.2) then
         call standard(2,2,cmag(nc_B,i),BV,id(i),x_Vbv,
     &   UB,BV,VR,VI,RI,airmass(nc_B),UT(nc_B),xc(i),yc(i),o_B)
         if(abs(o_B-x_B).ge.0.00001) nback=nback+1
        endif
c       V      
        if(nobs(nc_V,i).ne.0.and.num_filter(j).eq.3) then
         if(nobs(nc_B,i).ne.0.and.num_color(j).eq.2.and.ncal_vbv.ne.0
     &   ) then ! V : B-V
          call standard(3,2,cmag(nc_V,i),BV,id(i),x_Vbv,
     &    UB,BV,VR,VI,RI,airmass(nc_V),UT(nc_V),xc(i),yc(i),o_Vbv)
          if(abs(o_Vbv-x_Vbv).ge.0.00001) nback=nback+1
         endif
         if(nobs(nc_I,i).ne.0.and.num_color(j).eq.4.and.ncal_vvi.ne.0)
     &   then ! V : V-I
          call standard(3,4,cmag(nc_V,i),VI,id(i),x_Vvi,
     &    UB,BV,VR,VI,RI,airmass(nc_V),UT(nc_V),xc(i),yc(i),o_Vvi)
          if(abs(o_Vvi-x_Vvi).ge.0.00001) nback=nback+1
         endif
        endif
c       R : R-I
        if(nobs(nc_R,i).ne.0.and.nobs(nc_I,i).ne.0.and.ncal_R.ne.0.and.
     &  num_filter(j).eq.4.and.num_color(j).eq.5) then
         call standard(4,5,cmag(nc_R,i),RI,id(i),x_Vvi,
     &   UB,BV,VR,VI,RI,airmass(nc_R),UT(nc_R),xc(i),yc(i),o_R)
         if(abs(o_R-x_R).ge.0.00001) nback=nback+1
        endif
c       I
        if(nobs(nc_I,i).ne.0.and.num_filter(j).eq.5) then
         if(nobs(nc_V,i).ne.0.and.num_color(j).eq.4.and.ncal_ivi.ne.0) 
     &   then ! I : V-I
          call standard(5,4,cmag(nc_I,i),VI,id(i),x_Vvi,
     &    UB,BV,VR,VI,RI,airmass(nc_I),UT(nc_I),xc(i),yc(i),o_Ivi)
          if(abs(o_Ivi-x_Ivi).ge.0.00001) nback=nback+1
         endif 
         if(nobs(nc_R,i).ne.0.and.num_color(j).eq.5.and.ncal_iri.ne.0) 
     &   then ! I : R-I
          call standard(5,5,cmag(nc_I,i),RI,id(i),x_Vvi,
     &    UB,BV,VR,VI,RI,airmass(nc_I),UT(nc_I),xc(i),yc(i),o_Iri)
          if(abs(o_Iri-x_Iri).ge.0.00001) nback=nback+1
         endif
        endif
c       Ha : V-I
        if(nobs(nc_Ha,i).ne.0.and.num_filter(j).eq.6.and.ncal_Ha.ne.0) 
     &  then
         call standard(6,0,cmag(nc_Ha,i),0.,id(i),x_Vvi,
     &   UB,BV,VR,VI,RI,airmass(nc_Ha),UT(nc_Ha),xc(i),yc(i),o_Ha)
         if(abs(o_Ha-x_Ha).ge.0.00001) nback=nback+1
        endif
       enddo
       UB=o_U-o_B
       BV=o_B-o_Vbv
       VR=o_Vvi-o_R
       VI=o_Vvi-o_Ivi
       if(ncal_vvi.eq.0.and.ncal_vbv.ne.0) VI=o_Vbv-o_Ivi
       if(ncal_Ivi.eq.0.and.ncal_Iri.ne.0) VI=o_Vvi-o_Iri
       RI=o_R-o_Iri
       if(ncal_Iri.eq.0.and.ncal_Ivi.ne.0) RI=o_R-o_Ivi
       if(niter.ge.20) goto 121
       if(nback.eq.0) goto 121
       goto 101
121    if(niter.ge.20) nmax=nmax+1
       wt_ri=0.
       wt_vi=0.
       if(nobs(nc_R,i).ne.0.and.nobs(nc_I,i).ne.0.and.ncal_iri.ne.0) 
     & then 
        wt_ri=0.01/e_ri/e_ri 
        if(e_ri.eq.0.) wt_ri=1.00
        if(wt_ri.ge.100.) wt_ri=1.00
       endif
       if(nobs(nc_V,i).ne.0.and.nobs(nc_I,i).ne.0.and.ncal_ivi.ne.0) 
     & then
        wt_vi=0.01/e_vi/e_vi
        if(e_vi.eq.0.) wt_vi=100.
        if(wt_vi.ge.100.) wt_vi=100.
       endif
       o_I=(wt_vi*o_Ivi+wt_ri*o_Iri)/(wt_vi+wt_ri)
       wt_bv=0.
       wt_vi=0.
       if(nobs(nc_B,i).ne.0.and.nobs(nc_V,i).ne.0.and.ncal_vbv.ne.0) 
     & then
        wt_bv=0.01/e_bv/e_bv
        if(e_bv.eq.0.) wt_bv=100.
        if(wt_bv.ge.100.) wt_bv=100.
       endif 
       if(nobs(nc_I,i).ne.0.and.nobs(nc_V,i).ne.0.and.ncal_vvi.ne.0) 
     & then
        wt_vi=0.01/e_vi/e_vi
        if(e_vi.eq.0.) wt_vi=100.
        if(wt_vi.ge.100.) wt_vi=100.
       endif
       o_V=(wt_vi*o_Vvi+wt_bv*o_Vbv)/(wt_vi+wt_bv)
       if(nobs(nc_V,i).eq.0.or.nV.eq.0) o_V=0.
       if(nobs(nc_B,i).ne.0.and.nobs(nc_V,i).ne.0.and.
     & ncal_B.ne.0.and.(ncal_vbv.ne.0.or.ncal_vvi.ne.0)) BV=o_B-o_V
       if(nobs(nc_V,i).ne.0.and.nobs(nc_I,i).ne.0.and.(ncal_vbv.ne.0
     & .or.ncal_vvi.ne.0).and.(ncal_iri.ne.0.or.ncal_ivi.ne.0)) 
     & VI=o_V-o_I
       if(nobs(nc_R,i).ne.0.and.nobs(nc_I,i).ne.0.and.ncal_R.ne.0.and.
     & (ncal_iri.ne.0.or.ncal_ivi.ne.0)) RI=o_R-o_Ivi!@o_Iri
       if(nUB.eq.0.or.ncal_U.eq.0.or.ncal_B.eq.0) UB=0.
       if(nUB.eq.0.or.ncal_U.eq.0.or.ncal_B.eq.0) e_ub=0.
       if(nBV.eq.0.or.ncal_B.eq.0.or.ncal_vbv.eq.0) BV=0.
       if(nBV.eq.0.or.ncal_B.eq.0.or.ncal_vbv.eq.0) e_bv=0.
       if(nVR.eq.0.or.ncal_R.eq.0.or.ncal_vri.eq.0) VR=0.
       if(nVR.eq.0.or.ncal_R.eq.0.or.ncal_vri.eq.0) e_vr=0.
       if(nVI.eq.0.or.ncal_vvi.eq.0.or.ncal_ivi.eq.0) VI=0.
       if(nVI.eq.0.or.ncal_vvi.eq.0.or.ncal_ivi.eq.0) e_vi=0.
       if(nRI.eq.0.or.ncal_iri.eq.0) RI=0.
       if(nRI.eq.0.or.ncal_iri.eq.0) e_ri=0
       if(ncal_iri.eq.0) nRI=0
       if(ncal_ivi.eq.0.or.(ncal_iri.eq.0.and.ncal_v.eq.0)) nVI=0
       if(ncal_B.eq.0.or.ncal_V.eq.0) nBV=0
       if(ncal_B.eq.0.or.ncal_U.eq.0) nUB=0
       if(ncal_U.eq.0) o_U=0.
       if(ncal_B.eq.0) o_B=0.
       if(ncal_V.eq.0) o_V=0.
       if(ncal_R.eq.0) o_R=0.
       if(ncal_I.eq.0) o_I=0.
       if(ncal_Ha.eq.0) o_Ha=0.
       write(double,'(a1)') ' '
       if(ndouble(i).ne.0) write(double,'(a1)') 'D'
       if(nobs(nc_V,i).ne.0.and.nV.ne.0) then
        nid=nid+1
        if(nobs(nc_I,i).eq.0) o_I=0.
        write(41,51)i,nid,xc(i),yc(i),o_U,o_B,o_V,o_R,o_I,o_Ha,
     &  cmerr(nc_U,i),cmerr(nc_B,i),cmerr(nc_V,i),cmerr(nc_R,i),
     &  cmerr(nc_I,i),cmerr(nc_Ha,i),nobs(nc_U,i),nobs(nc_B,i),
     &  nobs(nc_V,i),nobs(nc_R,i),nobs(nc_I,i),nobs(nc_Ha,i),double
        write(42,50)i,nid,xc(i),yc(i),o_V,o_I,UB,BV,VI,RI,o_Ha,
     &  cmerr(nc_V,i),cmerr(nc_I,i),e_UB,e_BV,e_VI,e_RI,cmerr(nc_Ha,i),
     &  nobs(nc_V,i),nobs(nc_I,i),nUB,nBV,nVI,nRI,nobs(nc_Ha,i),double
        nit(nid)=niter
        mean_iter=mean_iter+niter
       elseif(nobs(nc_I,i).ne.0.and.nRI.ne.0) then
        nid=nid+1
        o_V=0.
        write(41,51)i,nid,xc(i),yc(i),o_U,o_B,o_V,o_R,o_I,o_Ha,
     &  cmerr(nc_U,i),cmerr(nc_B,i),cmerr(nc_V,i),cmerr(nc_R,i),
     &  cmerr(nc_I,i),cmerr(nc_Ha,i),nobs(nc_U,i),nobs(nc_B,i),
     &  nobs(nc_V,i),nobs(nc_R,i),nobs(nc_I,i),nobs(nc_Ha,i),double
        write(42,50)i,nid,xc(i),yc(i),o_V,o_I,UB,BV,VI,RI,o_Ha,
     &  cmerr(nc_V,i),cmerr(nc_I,i),e_UB,e_BV,e_VI,e_RI,cmerr(nc_Ha,i),
     &  nobs(nc_V,i),nobs(nc_I,i),nUB,nBV,nVI,nRI,nobs(nc_Ha,i),double
        nit(nid)=niter
        mean_iter=mean_iter+niter
       elseif(nobs(nc_I,i).ne.0.and.nt_iri.ne.0) then
        nid=nid+1
        call standard(5,5,cmag(nc_I,i),temp_iri,id(i),o_V,
     &  UB,BV,VR,VI,temp_iri,airmass(nc_I),UT(nc_I),xc(i),yc(i),o_Icol)
        call standard(5,5,cmag(nc_I,i),temp_iri+dt_iri,id(i),o_V,
     &  UB,BV,VR,VI,temp_iri+dt_iri
     &  ,airmass(nc_I),UT(nc_I),xc(i),yc(i),e_temp_high)
        call standard(5,5,cmag(nc_I,i),temp_iri-dt_iri,id(i),o_V,
     &  UB,BV,VR,VI,temp_iri-dt_iri
     &  ,airmass(nc_I),UT(nc_I),xc(i),yc(i),e_temp_low)
        e_temp_low=abs(e_temp_low-o_Icol)
        e_temp_high=abs(e_temp_high-o_Icol)
        e_temp=sqrt(cmerr(nc_I,i)**2+max(e_temp_low,e_temp_high)**2)
        write(41,51)i,nid,xc(i),yc(i),o_U,o_B,o_V,o_R,o_Icol,o_Ha,
     &  cmerr(nc_U,i),cmerr(nc_B,i),cmerr(nc_V,i),cmerr(nc_R,i),
     &  e_temp,cmerr(nc_Ha,i),nobs(nc_U,i),nobs(nc_B,i),
     &  nobs(nc_V,i),nobs(nc_R,i),nobs(nc_I,i),nobs(nc_Ha,i),double
        write(42,50)i,nid,xc(i),yc(i),o_V,o_Icol,UB,BV,VI,RI,o_Ha,
     &  cmerr(nc_V,i),e_temp,e_UB,e_BV,e_VI,e_RI,cmerr(nc_Ha,i),
     &  nobs(nc_V,i),nobs(nc_I,i),nUB,nBV,nVI,nRI,nobs(nc_Ha,i),double
        nit(nid)=niter
        mean_iter=mean_iter+niter
       elseif(nobs(nc_I,i).ne.0.and.nt_ivi.ne.0) then
        call standard(5,4,cmag(nc_I,i),temp_ivi,id(i),o_V,
     &  UB,BV,VR,temp_ivi,RI,airmass(nc_I),UT(nc_I),xc(i),yc(i),o_Icol)
        call standard(5,4,cmag(nc_I,i),temp_ivi+dt_ivi,id(i),o_V,
     &  UB,BV,VR,temp_ivi+dt_ivi
     &  ,RI,airmass(nc_I),UT(nc_I),xc(i),yc(i),e_temp_high)
        call standard(5,4,cmag(nc_I,i),temp_ivi-dt_ivi,id(i),o_V,
     &  UB,BV,VR,temp_ivi-dt_ivi
     &  ,RI,airmass(nc_I),UT(nc_I),xc(i),yc(i),e_temp_low)
        e_temp_low=abs(e_temp_low-o_Icol)
        e_temp_high=abs(e_temp_high-o_Icol)
        e_temp=sqrt(cmerr(nc_I,i)**2+max(e_temp_low,e_temp_high)**2)
        nid=nid+1
        write(41,51)i,nid,xc(i),yc(i),o_U,o_B,o_V,o_R,o_Icol,o_Ha,
     &  cmerr(nc_U,i),cmerr(nc_B,i),cmerr(nc_V,i),cmerr(nc_R,i),
     &  e_temp,cmerr(nc_Ha,i),nobs(nc_U,i),nobs(nc_B,i),
     &  nobs(nc_V,i),nobs(nc_R,i),nobs(nc_I,i),nobs(nc_Ha,i),double
        write(42,50)i,nid,xc(i),yc(i),o_V,o_Icol,UB,BV,VI,RI,o_Ha,
     &  cmerr(nc_V,i),cmerr(nc_I,i),e_UB,e_BV,e_VI,e_RI,cmerr(nc_Ha,i),
     &  nobs(nc_V,i),nobs(nc_I,i),nUB,nBV,nVI,nRI,nobs(nc_Ha,i),double
        nit(nid)=niter
        mean_iter=mean_iter+niter
       endif
      enddo
      std_iter=0.
      mean_iter=real(mean_iter)/real(nid)
      do i=1,nid
       std_iter=std_iter+real((nit(i)-mean_iter)*(nit(i)-mean_iter))
      enddo
      std_iter=sqrt(std_iter/real(nid))
      print *,nmax,nstars,mean_iter,std_iter
      close(41)
      close(42)
      return
      end
c     nfil = filter number
c     nc = color number
c     Filter number 1 = U, 2 = B, 3 = V, 4 = R, 5 = I,6 = Ha
c     Color number  1 = U-B, 2 = B-V, 3 = V-R, 4 = V-I, 5 = R-I  
c     Other number 100 = Ut   101=xc,      102=yc,    103=rad,   104=theta
      subroutine standard(nfil,nc,cmag1,col,id,v,ub,bv,vr,vi,ri,X,UT,
     &xc,yc,omag)
      parameter(nf=30,ns=100000)
      common/coeff_const/num_filter(nf),num_color(nf),aeK1(nf),aeK2(nf)
     &,zero(nf),xcenter(nf),ycenter(nf),color_dp,ut_dp,x_dp,y_dp,r_dp,
     &theta_dp,add1_dp,add2_dp,add3_dp,nconst
      character color_dp(nf)*16,ut_dp(nf)*16,x_dp(nf)*16,y_dp(nf)*16
     &,r_dp(nf)*16,theta_dp(nf)*16,add1_dp(nf)*16,add2_dp(nf)*16
     &,add3_dp(nf)*16
      integer nfil,nc,nt
      real cmag1,omag,r,th,ub,bv,vr,vi,ri,X,UT,xc,yc,pi,col,v
      do i=1,nf
       if(nfil.eq.num_filter(i).and.nc.eq.num_color(i)) then
        nt=i
        goto 100
       endif
100   i100=i100
      enddo
      r=sqrt((xc-xcenter(nt))**2+(yc-ycenter(nt))**2)
      pi=3.141593
      if(yc.ge.CY) then
       th=asin(y/r)*180./pi
      else
       th=asin(-y/r)*180./pi
      endif
      call cal_trans(color_dp(nt),id,v,ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,
     &dZ_col)
      call cal_trans(ut_dp(nt),id,v,ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,
     &dZ_ut)
      call cal_trans(x_dp(nt),id,v,ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_xc)
      call cal_trans(y_dp(nt),id,v,ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_yc)
      call cal_trans(r_dp(nt),id,v,ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,
     &dZ_rad)
      call cal_trans(theta_dp(nt),id,v,ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,
     &dZ_th)
      call cal_trans(add1_dp(nt),id,v,ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,
     &dZ_a1)
      call cal_trans(add2_dp(nt),id,v,ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,
     &dZ_a2)
      call cal_trans(add3_dp(nt),id,v,ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,
     &dZ_a3)
      omag=cmag1-(aeK1(nt)-aeK2(nt)*col)*X+zero(nt)
     &+dZ_col+dZ_ut+dZ_xc+dZ_yc+dZ_rad+dZ_th+dZ_a1+dZ_a2+dZ_a3
c      print *,nt,color_dp(nt),ut_dp(nt),x_dp(nt),y_dp(nt),r_dp(nt),
c     &theta_dp(nt),add1_dp(nt),add2_dp(nt),add3_dp(nt)
      return
      end
 
c     Color number  1 = U-B, 2 = B-V, 3 = V-R, 4 = V-I, 5 = R-I  
c     Other number 100 = Ut   101=xc,      102=yc,    103=rad,   104=theta
      subroutine cal_trans(input,id,v,ub,bv,vr,vi,ri,Air,UT,xc,yc,r,th,
     &yy)
      parameter(ns=100000)
      common/red_inp/ebv_fg,ebv_cl,ebv_bg,rv_fg,rv_cl,rv_bg,id_wc(ns),
     &ebv_wc(ns),nredip
      character input*16,H*1
      real x(100),y(100),ub,bv,vr,vi,ri,Air,UT,xc,yc,r,th,xxx,yy,q,v
      integer n,npoints,n_depend,id
c      write(*,'(t1,a16)') input
      open(31,file=input)
c      write(*,'(t1,a16)') input
      rewind(31)
      read(31,*) npoints
      yy=0. 
      if(npoints.eq.0) goto 900   
      read(31,*) n_depend
      n=0
      do i=1,100
       read(31,'(t1,a1)',end=99) H
       if(H.ne.'#') then
        backspace(31)
        n=n+1
        read(31,*) x(n),y(n) 
       endif
      enddo
99    if(n.ne.npoints) print *,input,' has different points!!!'
      close(31)
      if(n_depend.eq.0) yy=0.
      if(n_depend.eq.1) xxx=ub
      if(n_depend.eq.2) xxx=bv
      if(n_depend.eq.3) xxx=vr
      if(n_depend.eq.4) xxx=vi
      if(n_depend.eq.5) xxx=ri
      if(n_depend.eq.100) xxx=UT
      if(n_depend.eq.101) xxx=xc
      if(n_depend.eq.102) xxx=yc
      if(n_depend.eq.103) xxx=r
      if(n_depend.eq.104) xxx=theta
      q = ub - 0.72*bv
      if(n_depend.eq.105) then
       do i=1,nredip
        if(id.eq.id_wc(i)) then
         xxx=bv-ebv_wc(i)
         goto 800
        endif
       enddo
       if(q.gt.-1.05.and.q.lt.-0.5.and.bv.lt.1.8.and.v.lt.17.5) then
        call cor_red(bv,ub,xxx,ub0)
        q = ub0 - 0.72*xxx
c         print *,xxx,ub0,q, '(B-V)0,(U-B)0'
        if(q.gt.-1.05.and.q.lt.-0.5.and.bv.lt.1.8.and.v.lt.17.5) then
         call cor_red(bv,ub,xxx,ub0)
c         print *,xxx,ub0, '(B-V)0,(U-B)0'
         goto 800
        else
         xxx=bv-ebv_fg
c         print *,xxx,ub0, '(B-V)0,(U-B)0'
         goto 800
        endif
       endif 
      endif
800   if(n_depend.ne.0) call lint(xxx,x,y,yy,n,n)
900   return
      end

      subroutine read_coeff(input)
      parameter(nf=30,ns=100000)
      character H*1,input*16,line*50
      common/coeff_const/num_filter(nf),num_color(nf),aeK1(nf),aeK2(nf)
     &,zero(nf),xcenter(nf),ycenter(nf),color_dp,ut_dp,x_dp,y_dp,r_dp,
     &theta_dp,add1_dp,add2_dp,add3_dp,nconst
      common/temp_col/temp_uub,temp_bbv,temp_vbv,temp_vvi,temp_ivi,
     &temp_iri,dt_uub,dt_bbv,dt_vbv,dt_vvi,dt_ivi,dt_iri,
     &nt_uub,nt_bbv,nt_vbv,nt_vvi,nt_ivi,nt_iri
      character color_dp(nf)*16,ut_dp(nf)*16,x_dp(nf)*16,y_dp(nf)*16
     &,r_dp(nf)*16,theta_dp(nf)*16,add1_dp(nf)*16,add2_dp(nf)*16
     &,add3_dp(nf)*16
      data temp_uub,temp_bbv,temp_vbv,temp_vvi,temp_ivi,temp_iri,
     &dt_uub,dt_bbv,dt_vbv,dt_vvi,dt_ivi,dt_iri
     &/0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0., 0./ 
      open(21,file=input)
      nconst=0 
      do i=1,ns
       read(21,'(t1,a1,a50)',end=99) H,line
       if(H.ne.'#'.and.H.ne.'C'.and.H.ne.'D'.and.H.ne.'N') then
        backspace(21)
        nconst=nconst+1
        read(21,*) num_filter(nconst),num_color(nconst),aeK1(nconst),
     &  aeK2(nconst),zero(nconst),xcenter(nconst),ycenter(nconst),
     &  color_dp(nconst),ut_dp(nconst),x_dp(nconst),y_dp(nconst),
     &  r_dp(nconst),theta_dp(nconst),add1_dp(nconst),add2_dp(nconst),
     &  add3_dp(nconst)
       elseif(H.eq.'C') then
        read(line,*) 
     &  temp_uub,temp_bbv,temp_vbv,temp_vvi,temp_ivi,temp_iri
       elseif(H.eq.'D') then
        read(line,*) dt_uub,dt_bbv,dt_vbv,dt_vvi,dt_ivi,dt_iri
       elseif(H.eq.'N') then
        read(line,*) nt_uub,nt_bbv,nt_vbv,nt_vvi,nt_ivi,nt_iri
       endif
c      print *,nconst,add3_dp(nconst)
      enddo
99    close(21)
      return
      end
! linear interpolation program for fortran.
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

      subroutine readwchead(input)
      parameter (nf=30,ns=100000)
      character fil(nf)*10,input*16,H*1,coll(10)*1,filter(nf)*10
      common /ncw/nc_U,nc_B,nc_V,nc_R,nc_I,nc_Ha,nc_K
     &,nfilter,nc_H,nc_J,nstars
      common /wcdata/ID(ns),xc(ns),yc(ns),xe(ns),ye(ns),nobs(nf,ns),
     &cmag(nf,ns),cmerr(nf,ns),id_f(nf,ns),ndouble(ns),
     &airmass(nf),UT(nf),filter
      data nc_U,nc_B,nc_V,nc_R,nc_I,nc_Ha,nc_K,nfilter,nc_H,nc_J
     &     /0  ,0   ,0   ,0   ,0   ,0    ,0   ,0      ,0   ,0  / 
      open(32,file=input)
      do i=1,ns
       read(32,'(a1)') H
       if(H.eq.'F') then
        backspace(32)
        read(32,'(t11,30a10)') (filter(m),m=1,nf)
        do j=1,nf
         read(filter(j),'(10a1)') (coll(m),m=1,10)
         do k=1,10
          if(coll(k).eq.'U'.or.coll(k).eq.'u') then
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
     &    ' ') then
           nfilter=nfilter+1
           nc_Ha=nfilter
          endif
         enddo
        enddo
       elseif(H.eq.'X') then
        backspace(32)
        read(32,'(t11,30a10)') (fil(m),m=1,nf)
        do j=1,nfilter
         read(fil(j),*) airmass(j)
        enddo
       elseif(H.eq.'T') then
        backspace(32)
        read(32,'(t11,30a10)') (fil(m),m=1,nf)
        do j=1,nfilter
         read(fil(j),*) UT(j)
        enddo
       elseif(H.eq.'*') then
        goto 100
       endif
      enddo
100   do i=1,ns
       read(32,*,end=300) ID(i),(id_f(j,i),j=1,nfilter)
       read(32,*,end=300) xc(i),(cmag(j,i),j=1,nfilter)
       read(32,*,end=300) yc(i),(cmerr(j,i),j=1,nfilter)
       read(32,*,end=300) xe(i),(nobs(j,i),j=1,nfilter)
       read(32,*,end=300) ye(i)
       read(32,*,end=300) ndouble(i)
       read(32,'(a1)',end=300) H
       read(32,'(a1)',end=300) H
       read(32,'(a1)',end=300) H
       read(32,'(a1)',end=300) H
      enddo
300   close(32)
      nstars=i-1
      return
      end

      subroutine read_ZAMS()
      character H*1
      common /zams/zMv(41),zVI(41),zBV(41),zUB(41)
      open(23,file='Sung2013.ZAMS.cc')
      read(23,'(a1)') H
      do i=1,41
       read(23,*,end=100) zMv(i),zVI(i),zBV(i),zUB(i)
      enddo
100   close(23)
      return
      end
      subroutine reddening_input()
      parameter (ns=100000)
      common/red_inp/ebv_fg,ebv_cl,ebv_bg,rv_fg,rv_cl,rv_bg,id_wc(ns),
     &ebv_wc(ns),nredip
      character H*1
      open(22,file='reddening.aex')
      read(22,*,end=99) ebv_fg,ebv_cl,ebv_bg
      read(22,*,end=99) rv_fg,rv_cl,rv_bg
      nredip=0
      do i=1,ns
       read(22,'(t1,a1)',end=99) H
       if(H.ne.'#') then
        backspace(22)
        nredip=nredip+1
        read(22,*) id_wc(nredip),ebv_wc(nredip)
       endif
      enddo
99    close(22)
      return
      end

      subroutine cor_red(bv,ub,bv0,ub0)
      common /zams/zMv(41),zVI(41),zBV(41),zUB(41)
      real bv,ub,bv0,ub0,s,ss,bv1,ub1,ebv_sum
      s=0.72
c      ss=0.025
      ebv_sum=0.
      bv1=bv
      ub1=ub
100   call lint(ub1,zUB,zBV,bv0,41,41)
      ebv_sum=ebv_sum+(bv1-bv0)
      ub0=ub-s*ebv_sum-ss*ebv_sum*ebv_sum
      bv0=bv-ebv_sum
      if((ub0.le.zub(1).or.bv0.le.zbv(1))) goto 500
      if(abs(ub0-ub1).gt.0.0001.or.bv1-bv0.gt.0.) then
       ub1=ub0
       bv1=bv0
       goto 100
      endif
500   return
      end
