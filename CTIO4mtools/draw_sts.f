c     f77 -o draw_sts draw_sts.f -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 
      parameter (ns=10000)
      real*8 ra(ns),dec(ns),cra,cdec,ram,ramax,decmax,decmin,pi,dr,dd
      real dra,ddec,ra1,ra2,dec1,dec2,v(ns),cv,pp,cut,rac,decc
      real radd,decdd,racs,deccs
      integer nstar,nmi,nbox,nfill
      character region*12,input*16,c(16)*1, H*1,sname(ns)*15,ps*28
      character outp*15,chname(15)*1,chnum(ns)*15,H*1,fillname(ns)*15
      pi=3.14159265358979d0
      open(31,file='draw_sts.fill')
      open(32,file='draw_sts.input')
      do i=1,ns
       read(31,'(t1,a15)',end=9) fillname(i)
      enddo
9     nfill=i-1
      close(31)
c      print *,'Input region name=?'
      write(region,'(a16)')'                             '
      write(input,'(a16)')'                             '
      read(32,'(a12)') region
      read(region,'(12a1)') (c(k),k=1,12)
      do k=1,16
       if(c(k).eq.' ') then
        goto 10
       endif
      enddo
10    write(c(k),'(a1)') '.'
      write(c(k+1),'(a1)') 'd'
      write(c(k+2),'(a1)') 'a'
      write(c(k+3),'(a1)') 't'
      write(input,'(16a1)') (c(i),i=1,k+3)
      write(c(k+1),'(a1)') 'p'
      write(c(k+2),'(a1)') 's'
      write(outp,'(16a1)') (c(i),i=1,k+2)

      open(21,file=input)
      read(21,'(t1,a1)') H
      nstar=0
      do i=1,ns
       nstar=nstar+1
       read(21,'(t1,f8.3,t105,f15.11,t121,f15.11,t137,a15)',end=99) 
     & cv,cra,cdec,sname(nstar)
 
       read(sname(nstar),'(15a1)') (chname(j),j=1,15)
       do j=1,15
        if(chname(j).eq.'-') then
         nmi=j 
         goto 20
        endif
       enddo
20     write(chnum(nstar),'(a15)') '                                  '
       write(chnum(nstar),'(15a1)') (chname(j),j=nmi+1,15)
       ra(nstar)=cra
       dec(nstar)=cdec
       v(nstar)=cv
       if(nstar.eq.1) then
        ramax=cra!ra(nstar)
        decmax=dec(nstar)
        ram=cra!ra(nstar)
        decmin=dec(nstar)

       elseif(ra(nstar).gt.ramax) then
        ramax=cra!ra(nstar)
       elseif(ra(nstar).lt.ram) then
        ram=cra!ra(nstar)
       endif
       if(dec(nstar).gt.decmax) then
        decmax=dec(nstar)
       elseif(dec(nstar).lt.decmin) then
        decmin=dec(nstar)
       endif

      enddo 
99    close(21)
      nstar=nstar-1

      cdec=(decmax+decmin)/2.d0
      cra=(ramax+ram)/2.d0
c      dr=(ramax-ram)/2.d0
c      dd=(decmax-decmin)/2.d0
c      dr=08.d0/60.d0
c      dd=08.d0/60.d0
c      print *,ra(nstar),ramax*60.*cos(cdec*pi/180.)+1.
c     &,ram*60.*cos(cdec*pi/180.)-1.
c      print *,'FOV, cut magnitude=?'
      read(32,*) dd,cut
c      print *,'ra shift,dec shift=? (arcmin)'
      read(32,*) radd,decdd   
      read(32,*) racs,deccs
c      dd=18.8
c      cut=18.
      dr=dd

      dec2=(dd )+decdd
      dec1=(-dd )+decdd
c      ra2=(dr)*cos((cdec+decdd)*pi/180.)+1.+radd
c      ra1=(-dr)*cos((cdec+decdd)*pi/180.)-1.+radd
      ra2=dr+radd
      ra1=-dr+radd

      print *,ra1,ra2

      write(ps,'(a13,a15)') 'postportfile ',outp

      call sm_device (ps)
      call sm_graphics
      call sm_erase()
      call sm_location(3000,30000,27000,30000)
      call sm_limits(1.,100.,1.,100.)
      call sm_relocate(38.,50.)
      call sm_expand(3.)
      call sm_label(region)
      call sm_expand(1.)

      call sm_limits(ra2,ra1,dec1,dec2)
      call sm_location(2500,32000,3000,27000)
      call sm_box(1,2,1,2)
      call sm_xlabel('\\gDR.A.')
      call sm_ylabel('\\gDD.E.C.')

c      print *,'box=? (1=CTIO mosaic II)'
c      read(*,*) nbox
      nbox=1
      call sm_ctype('red')
      call draw_box(nbox,racs,deccs)
      call sm_ctype('black')
      call sm_lweight(0.01)
      do i=1,nstar
       pp=(18.-v(i))/10. 
       if(pp.gt.0.99)  pp=0.99
       if(pp.lt.0.01)  pp=0.01
       ddec=(dec(i)-cdec)*60. 
       dra=(ra(i)-cra)*cos(dec(i)*pi/180.)*60.
       call sm_ptype(990.+pp,1)
       call sm_expand(2.)
       call sm_ctype('black')
       do j=1,nfill
        if(sname(i).eq.fillname(j)) then
         call sm_ptype(993.+pp/0.7,1)
         if(pp/0.7.gt.0.99) call sm_ptype(993.99,1)
         call sm_ctype('blue')
         print *,'The following stars are marked as a different color.'
         print *,sname(i)
        endif
       enddo
       call sm_points(dra,ddec,1)
       if(v(i).lt.cut) then
        call sm_expand(0.4)
        call sm_relocate(dra,ddec)
        call sm_label(chnum(i))
       endif
      enddo
      call sm_expand(1.)


      call sm_location(3000,30000,500,1500)
      call Csize(18.,8.)

      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      close(32)

      stop
      end

      subroutine draw_box(nbox,cra,cdec)
      integer nbox
      real cra,cdec,cx,cy,scl
      real x(8),y(4),dra(4),ddec(8)
      data x/1.,2052.,2119.,4170.,4239.,6292.,6359.,8410./
      data y/1.,4098.,4131.,8228./
      cx=4205.5
      cy=4114.5
      scl=0.268
      do i=1,8
       if(i.le.4) dra(i)=((cy-y(i))*scl)/60.+cra
       ddec(i)=((x(i)-cx)*scl)/60.+cdec
      enddo
      call sm_ltype(0)
      call sm_lweight(0.5)
      if(nbox.eq.1) then
       call box
     & (dra(1),ddec(1),dra(2),ddec(1),dra(2),ddec(2),dra(1),ddec(2))
       call box
     & (dra(1),ddec(3),dra(2),ddec(3),dra(2),ddec(4),dra(1),ddec(4))
       call box
     & (dra(3),ddec(1),dra(4),ddec(1),dra(4),ddec(2),dra(3),ddec(2))
       call box
     & (dra(3),ddec(3),dra(4),ddec(3),dra(4),ddec(4),dra(3),ddec(4))
       call box
     & (dra(1),ddec(5),dra(2),ddec(5),dra(2),ddec(6),dra(1),ddec(6))
       call box
     & (dra(1),ddec(7),dra(2),ddec(7),dra(2),ddec(8),dra(1),ddec(8))
       call box
     & (dra(3),ddec(5),dra(4),ddec(5),dra(4),ddec(6),dra(3),ddec(6))
       call box
     & (dra(3),ddec(7),dra(4),ddec(7),dra(4),ddec(8),dra(3),ddec(8))


      endif

      return
      end 

      subroutine box(x1,y1,x2,y2,x3,y3,x4,y4)
      real x1,y1,x2,y2,x3,y3,x4,y4
      call sm_relocate(x1,y1)
      call sm_draw(x2,y2)
      call sm_draw(x3,y3)
      call sm_draw(x4,y4)
      call sm_draw(x1,y1)
      return
      end

      subroutine Csize(maxmag,bright)
      real bright,f,pp(1),x(1),y(1),p,t,maxmag
      character mags*5
      common /id/idw
      t=0.1*(maxmag-bright)
      call sm_limits(0.,100.,0.,100.)
      p=400.

      call sm_expand(2.)
      f =-0.01/(maxmag-bright)+1.
      pp(1)=p+f
      x(1)=70.
      y(1)=60.
      call sm_ptype(pp,1)
      call sm_expand(2.)
      call sm_points(x(1),y(1),1)
      call sm_expand(1.)
      write(mags,'(f4.1)') bright+0.01
      call sm_relocate(x(1)-2.0,y(1)-55.)
      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(2.)

      f =-2.5*t/(maxmag-bright)+1.
      pp(1)=p+f
      x(1)=80.
      call sm_ptype(pp,1)
      call sm_expand(2.)
      call sm_points(x(1),y(1),1)
      call sm_expand(1.)
      write(mags,'(f4.1)') bright+2.5*t
      call sm_relocate(x(1)-2.0,y(1)-55.)
      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(2.)

      f =-5.*t/(maxmag-bright)+1.
      pp(1)=p+f
      x(1)=85.
      call sm_ptype(pp,1)
      call sm_expand(2.)
      call sm_points(x(1),y(1),1)
      call sm_expand(1.)
      write(mags,'(f4.1)') bright+5.*t
      call sm_relocate(x(1)-2.0,y(1)-55.)
      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(2.)

      f =-7.5*t/(maxmag-bright)+1.
      pp(1)=p+f
      x(1)=90.
      call sm_ptype(pp,1)
      call sm_expand(2.)
      call sm_points(x(1),y(1),1)
      call sm_expand(1.)
      write(mags,'(f4.1)') bright+7.5*t
      call sm_relocate(x(1)-2.0,y(1)-55.)
      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(2.)

      f =-9.90*t/(maxmag-bright)+1.
      pp(1)=p+f
      x(1)=95.
      call sm_ptype(pp,1)
      call sm_expand(2.)
      call sm_points(x(1),y(1),1)
      call sm_expand(1.)
      write(mags,'(f4.1)') bright+9.9*t
      call sm_relocate(x(1)-2.0,y(1)-55.)

      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(1.)

      return
      end

