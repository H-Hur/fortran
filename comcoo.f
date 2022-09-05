      parameter(ns=50000)
      real*8 xc_in(ns),yc_in(ns),xc_compare(ns),yc_compare(ns),fitrad
      real*8 dx,dy,d!,xc_out(ns),yc_out(ns)
      integer nstar_in,nstar_compare,nstar_new,nfitrad
      character input*16,compare*16,output*7,H*1

      print *,'Input coordinate file=? (from all star file)'
      read(*,'(a16)') input
       
      print *,'Coordinate file to compare to the first input file=?(from
     & rem.coo)'
      read(*,'(a16)') compare

      print *,'Fitrad=? (generally, psfrad(=4FWHM))'
      read(*,*) fitrad
      Print *,'Output file name is [rem.coo].'
      write(output,'(a7)') 'rem.coo'
      
      open(31,file=input)
      open(32,file=compare)
      do i=1,ns
        read(31,*,end=100) xc_in(i),yc_in(i)
      enddo
100   nstar_in=i-1
      close(31)

      do i=1,ns
        read(32,*,end=200) xc_compare(i),yc_compare(i)
      enddo
200   nstar_compare=i-1
      close(32)
      open(33,file=output)
      nstar_new=0

      do j=1,nstar_compare
       nfitrad=0
       do i=1,nstar_in
        dx=xc_in(i)-xc_compare(j)
        dy=yc_in(i)-yc_compare(j)
        d=sqrt(dx**2.d0+dy**2.d0)
        if(d.le.fitrad) then
         nfitrad=nfitrad+1
        endif
       enddo
       if(nfitrad.eq.0) then
        nstar_new=nstar_new+1
        write(33,*) xc_compare(j),yc_compare(j)
       endif
      enddo
      close(33)
      close(32)
      close(31) 

      print *,nstar_new,'stars found.'
      stop
      end
