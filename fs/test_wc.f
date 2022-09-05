c     f77 -o test_wc /home/gjgusdh/fortran/subroutine/fs/test_wc.f -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 
      character input*16

      write(input,'(a16)') 'wcout2.wc                 '
      call readwchead(input)



      stop
      end


      subroutine draw_combined()
      parameter (nf=30,ns=100000)
      character filter(nf)*10,ps*19,fname*3
      real*8 airmass(nf),UT(nf)
      common /ncw/nc_U,nc_B,nc_V,nc_R,nc_I,nc_Ha,nc_K
     &,nfilter,nc_H,nc_J,nstars
      common /wcdata/ID(ns),xc(ns),yc(ns),xe(ns),ye(ns),nobs(nf,ns),
     &cmag(nf,ns),cmerr(nf,ns),id_f(nf,ns),airmass,UT,filter
      real x1,x2

      do j=1,nfilter
       if(j.eq.nc_U) write(fname,'(a3)') 'dU_'
       if(j.eq.nc_B) write(fname,'(a3)') 'dB_'
       if(j.eq.nc_V) write(fname,'(a3)') 'dV_'
       if(j.eq.nc_R) write(fname,'(a3)') 'dR_'
       if(j.eq.nc_I) write(fname,'(a3)') 'dI_'
       if(j.eq.nc_Ha) write(fname,'(a3)') 'dHa'
       if(j.eq.nc_J) write(fname,'(a3)') 'dJ_'
       if(j.eq.nc_H) write(fname,'(a3)') 'dH_'
       if(j.eq.nc_K) write(fname,'(a3)') 'dK_'
       write(ps,'(t1,a13,2a3)') 'postlandfile ',fname,'.ps' 
       x1=0.
       x2=0.
       do i=1,nstars
        if(id_f(j,i).ne.0) then
         x1=min(x1,cmag(j,i))
         x2=min(x2,cmag(j,i))
        endif 
       enddo

       call sm_device (ps)
       call sm_graphics
       call sm_erase()
       call sm_limits(x1-0.5,x2+0.5,-0.3,0.3)
       call sm_box(1,2,0,0)
c       call sm_xlabel('Mag')
cd       call sm_                 
       call sm_gflush()
       call sm_hardcopy
       call sm_alpha

      enddo


      return
      end


      subroutine readwchead(input)
      parameter (nf=30,ns=100000)
      character fil(nf)*10,input*16,H*1,coll(10)*1
      character filter(nf)*10
      real*8 airmass(nf),UT(nf)
      common /ncw/nc_U,nc_B,nc_V,nc_R,nc_I,nc_Ha,nc_K
     &,nfilter,nc_H,nc_J,nstars
      common /wcdata/ID(ns),xc(ns),yc(ns),xe(ns),ye(ns),nobs(nf,ns),
     &cmag(nf,ns),cmerr(nf,ns),id_f(nf,ns),airmass,UT,filter
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
      enddo
300   close(31)
      return
      end

