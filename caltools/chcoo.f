c     input ra,dec,cra,cdec must be in degrees!!!
c     output dra,ddec is in arcmin
      subroutine eq_del(ra,dec,cra,cdec,dra,ddec)
      real*8 ra,dec,cra,cdec,dra,ddec,pi
      pi=3.14159265358979d0
      ddec=60.d0*(dec-cdec)
      dra=60.d0*(ra-cra)*dcos(dec*pi/180.d0)
      return
      end

c     dra,ddec,cra,cdec must be in degrees!!!
c     dra,ddec is in arcmin
      subroutine del_eq(dra,ddec,cra,cdec,ra,dec)
      real*8 ra,dec,cra,cdec,dra,ddec,pi
      pi=3.14159265358979d0
      dec=cdec+ddec/60.d0
      ra=cra+dra/60.d0/dcos(dec*pi/180.d0)
      return
      end
