c     f95 -o make_eqcoo make_eqcoo.f
      integer inum
      character input*32

      print *,'type the input coordinate file name from the WEBDA 
     &(1 = ad2000.coo)'
      read(*,'(t1,a32)') input
      read(input,*) inum
      if(inum.eq.1) then
       open(21,file='ad2000.coo')
      else
       open(21,file=input) 
      endif
       
      call read_coo_webda(21)
      stop
      end


      subroutine read_coo_webda(input)
      parameter (ns=100000)
      integer input,nstar,ncoo,no(ns),ra1,ra2,dec1,dec2,numb
      integer no_last,num_ref,nobs(ns)
      real*8 ra(ns),dec(ns),ra3,dec3,pi,ra_r,dec_r
      character H*1
      pi=3.14159265358979d0
      nstar=0
      ncoo=1
      no_last=0
      do i=1,ns
       read(input,'(t1,a1)') H
       if(H.eq.'-') goto 10
      enddo
10    do i=1,ns     
       read(input,*,end=99) numb,num_ref,ra1,ra2,ra3,dec1,dec2,dec3
       ra_r=15.d0*(real(ra1)+real(ra2)/60.d0+ra3/3600.d0)
       if(dec1.ge.0) then
        dec_r=real(dec1)+real(dec2)/60.d0+dec3/3600.d0
       else
        dec_r=real(dec1)-real(dec2)/60.d0-dec3/3600.d0
       endif
       if(numb.ne.no_last) then
        no_last=numb
        nstar=nstar+1
        nobs(nstar)=1
        no(nstar)=numb
        ra(nstar)=ra_r
        dec(nstar)=dec_r
       else
        nobs(nstar)=nobs(nstar)+1
        ra(nstar)=ra(nstar)+ra_r
        dec(nstar)=dec(nstar)+dec_r
       endif
      enddo
99    close(input)
      open(22,file='eqcoo.dat')
      write(22,'(t1,a75)') '#   NO R.A.(degree)     D.E.C.(degree) R.A.(
     &hour)   D.E.C(degree) Nobs
     &                                                                 '
      do i=1,nstar
       ra(i)=ra(i)/real(nobs(i))
       dec(i)=dec(i)/real(nobs(i))
       ra_r=ra(i)/15.d0
       ra1=int(ra_r)
       ra2=int((ra_r-real(ra1))*60.d0)
       ra3=(ra_r-real(ra1)-real(ra2)/60.d0)*3600.d0
       if(dec(i).ge.0.d0) then
        dec1=int(dec(i))
        dec2=int((dec(i)-real(dec1))*60.d0)
        dec3=(dec(i)-real(dec1)-real(dec2)/60.d0)*3600.d0      
       else
        dec1=int(dec(i))
        dec2=-int((dec(i)-real(dec1))*60.d0)
        dec3=-(dec(i)-real(dec1)+real(dec2)/60.d0)*3600.d0
       endif      
       write(22,'(t1,i6,2f16.11,2i3,f6.2,i4,i3,f6.2,i4)') 
     & no(i),ra(i),dec(i),ra1,ra2,ra3,dec1,dec2,dec3,nobs(i)
      enddo
      close(22)
      return
      end




