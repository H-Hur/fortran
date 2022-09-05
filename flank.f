c for only fortran77 to draw color-color diagram.
c f77 -o 'flank.exe' 'flank.f' -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 

      parameter(ninput=1000000,ntem=3)
      real logwav1(ninput),logflux1(ninput)
      real logwav2(ninput),logflux2(ninput)
      real logwav3(ninput),logflux3(ninput)
      real*8 T,B,wmax,peak,wav,B0,Bsol
      real xstart(ntem),ystart(ntem),xlimit(ntem),ylimit(ntem),x1,x2,y1
      real y2
      integer i,j,tran,ndata(ntem),nstart,nend
      character output(ninput)*13
      
c      c=2.99792458*10.**10.
c      h=6.6260755*10.**27.
c      k=1.380658*10.**16.
      call sm_device ('postportfile flank.ps')
      call sm_graphics
      call sm_erase()
      call sm_location(3000,30000,3000,30000)

      wmax=0.d0
      peak=0.d0
      Bsol=0.
      T=6000.d0
      do i=1,ntem
       xstart(i)=0.
       xlimit(i)=0.
       ystart(i)=0.
       ylimit(i)=0.
      enddo
      tran=1
      open(21,file='600.flank')
      open(22,file='6000.flank')
      open(23,file='60000.flank')

      do j=1,ntem
       T=600.d0*10.d0**real(j-1)
       peak=0.d0
       wmax=0.d0

       do i=1,ninput,1
        wav=real(i)/10.d0**8.d0
        call flank(wav,T,B)
        if(tran.eq.1) then
         call flank(wav,6000.d0,B0)
         if(B0.gt.Bsol) then
          Bsol=B0
          tran=1
          elseif(B0.lt.Bsol) then
          tran=0
         endif      
        endif
        if(B.gt.peak) then
         if(peak.eq.0.d0) then
          xstart(j)=wav
          nstart=i
         endif
         peak=B
         wmax=wav
         
         elseif(B.lt.peak.and.B.eq.0.d0) then
         goto 100
        endif
        if(j.eq.1) then
         logwav1(i)=log10(wav)
         logflux1(i)=log10(B)
         elseif(j.eq.2) then
         logwav2(i)=log10(wav)
         logflux2(i)=log10(B)
         elseif(j.eq.3) then
         logwav3(i)=log10(wav)
         logflux3(i)=log10(B)
        endif
       if(B.eq.0.d0) then
        goto 50
       endif
       write(20+j,*) j,T,wav*10.d0**8.d0,B,log10(wav),log10(B)
50     B=B       
       enddo
100    if(wav.gt.xlimit(j)) then
        xlimit(j)=wav
       endif
       nend=i
       ndata(j)=nend-nstart
      print *,j,T,log10(peak),wmax*10.d0**8.d0,xstart(j),xlimit(j)
      ylimit(j)=log10(peak*1.01d0)
      enddo
      do i=2,ntem
       if(xstart(i).lt.xstart(i-1)) then
        x1=xstart(i)
        else
        x1=xstart(i-1)
       endif
       if(xlimit(i).gt.xlimit(i-1)) then
        x2=xlimit(i)
        else
        x2=xlimit(i-1)
       endif
       if(ylimit(i).gt.ylimit(i-1)) then
        y2=ylimit(i)
        else
        y2=ylimit(i-1)
       endif
      enddo
      print *,'end of compute'
      x1=log10(x1)
      x2=log10(x2)
      y1=-10.
      print *,x2,x1,y1,y2  

      call sm_limits(x1,x2,y1,y2)
      print *,'end of limits'

      call sm_ctype('black')
      print *,'end of ctype'

      call sm_box(1,2,3,3)
      print *,'drawing 1'
      call sm_conn(logwav1,logflux1,ndata(1))
      print *,'drawing 2'

      call sm_conn(logwav2,logflux2,ndata(2))
      print *,'drawing 3'
      call sm_conn(logwav3,logflux3,ndata(3))
      call sm_xlabel('wavelength(cm) : logscale')
      call sm_ylabel('B(lamda) : logscale')

      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      close(21)
      close(22)
      close(23)
      stop
      end
 
      subroutine flank(wav,temp,B)
      real*8 wav,temp,c,h,k,B
      c=2.99792458d0*10.d0**10.d0
      h=6.6260755d0*10.d0**(-27.d0)
      k=1.380658d0*10.d0**(-16.d0)
      B=(2.d0*h*c**2.d0/wav**5.d0)/(exp(h*c/(wav*k*temp))-1.d0)
      return
      end

