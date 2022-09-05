c      gfortran -o cal_dra ~/work/fortran/subroutine/caltools/cal_dra.f
c      f95 -o cal_dra ~/fortran/subroutine/caltools/cal_dra.f
       integer ra1,ra2,dec1,dec2
       real*8 ra,dec,cra,cdec,ra3,dec3,dra,ddec,pi
       print *,'R.A. and D.E.C of the Center =? (h,m,s,d,m,s)'
       read(*,*) ra1,ra2,ra3,dec1,dec2,dec3
       cra=(real(ra1)+real(ra2)/60.d0+ra3/3600.d0)*15.d0
       if(dec1.ge.0.d0) then
        cdec=real(dec1)+real(dec2)/60.d0+dec3/3600.d0
       else
        cdec=real(dec1)-real(dec2)/60.d0-dec3/3600.d0
       endif
       print *,'R.A. and D.E.C of the target =? (h,m,s,d,m,s)'
       read(*,*) ra1,ra2,ra3,dec1,dec2,dec3
       ra=(real(ra1)+real(ra2)/60.d0+ra3/3600.d0)*15.d0
       if(dec1.ge.0.d0) then
        dec=real(dec1)+real(dec2)/60.d0+dec3/3600.d0
       else
        dec=real(dec1)-real(dec2)/60.d0-dec3/3600.d0
       endif
 
       pi=3.14159265358979d0
       dra=(ra-cra)*60.d0*cos(dec*pi/180.d0)
       ddec=(dec-cdec)*60.d0
       if(dra.gt.1.and.ddec.gt.1) then  
        print *,'delta R.A and delta D.E.C. in arcmin are, '
        write(*,'(2f10.3)') dra,ddec
       else
        print *,'delta R.A and delta D.E.C. in arcsec are, '
        write(*,'(2f10.3)') dra*60.d0,ddec*60.d0
       endif

       stop
       end
