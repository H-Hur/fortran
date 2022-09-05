c       f77 -o draw_webda draw_webda.f -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 
c     g77 -o draw_webda draw_webda.f -L/Applications/scisoft/i386/Packages/sm/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 -L/Users/apple/utils/sm_mac/scisoft/lib -laquaterm

      character coodat*32
      open(31,file='draw_webda.input')
c      print *,'input data file name='
c      print *, 'input coo file name='
      call read_combine()

      read(31,'(a32)') coodat
      call read_coo(coodat)
   
      call draw_map(31)

      stop
      end



      subroutine draw_map(input)
      parameter (ns=100000,nc=4)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      common/coo_dat/no_coo(ns),ra,dec,nobs_coo(ns),ncoo
      real*8 ra(ns),dec(ns),pi
      real vlimit,siz,ra_d,dec_d,cra,cdec
      integer n_color,input
      character region*16,lid*6
      pi=3.14159265358979d0
      read(input,'(a16)') region
      read(input,*) cra,cdec,dra,ddec
      dra=dra/2.
      ddec=ddec/2.
      read(input,*) vlimit,n_color

      call sm_device ('postportfile map_webda.ps')
      call sm_graphics
      call sm_erase()

      call sm_location(3000,32000,26000,30000)
      call sm_limits(0.,100.,0.,100.)
      call sm_relocate(0.,100.)
      call write_color(n_color)

      
      call sm_location(3000,32000,3000,25000)
      call sm_limits(dra,-dra,-ddec,ddec)
      call sm_xlabel('\\gD(R.A.)\\o|')
      call sm_ylabel('\\gD(D.E.C.)\\o|')
      call sm_box(1,2,0,0)
      do i=1,nstar
       do j=1,ncoo
        if(no(i).eq.no_coo(j)) then
         if(n_color.ge.0) call ctype(n_color,p(i))
         if(n_color.ge.1) call ctype(n_color,com(i,n_color)) ! V=1, B-V=2, U-B=3, V-I=4
         siz=(vlimit-5.)/(com(i,1)-5.)/3.
         if(siz.gt.0.99) siz=0.99
         if(siz.lt.0.01) siz=0.01
         call sm_ptype(400.+siz,1)
         ra_d=cos(dec(j)*pi/180.)*(ra(j)-cra*15.)*60.
         dec_d=(dec(j)-cdec)*60.
         call sm_expand(2.)
         if(com(i,1).gt.vlimit) goto 10
         call sm_points(ra_d,dec_d,1)
         call sm_relocate(ra_d,dec_d)
         write(lid,'(i6)') no(i) 
         call sm_expand(0.5)
         call sm_label(lid)
         goto 10
        endif
       enddo
10     i10=i10
      enddo
      call sm_expand(1.)
      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      return
      end

      subroutine ctype(n_color,color)
      integer n_color
      real color
      call sm_lweight(1.)
      if(n_color.eq.0) then
       if(color.ge.0.75) then
        call sm_ctype('blue')
        call sm_lweight(2.5)
       elseif(color.ge.0.5) then
        call sm_lweight(1.5)
        call sm_ctype('cyan')
       elseif(color.ge.0.25) then
        call sm_ctype('black')
        call sm_lweight(1.0)
       elseif(color.ge.0.) then
        call sm_ctype('green')
        call sm_lweight(0.7)
       else
        call sm_ctype('yellow')
        call sm_lweight(0.2)
       endif
      elseif(n_color.eq.2) then
       if(color.lt.-0.3) then
        call sm_ctype('blue')
        call sm_lweight(1.2)
       elseif(color.lt.-0.05) then
        call sm_ctype('cyan')
       elseif(color.lt.0.3) then
        call sm_ctype('black')
        call sm_lweight(0.7)
       elseif(color.lt.0.5) then
        call sm_ctype('green')
       elseif(color.lt.0.85) then
        call sm_ctype('yellow')
        call sm_lweight(5.0)
       elseif(color.lt.1.5) then
        call sm_ctype('magenta')
        call sm_lweight(0.7)
       else
        call sm_ctype('red')
        call sm_lweight(1.2)
       endif
      elseif(n_color.eq.3) then
       if(color.lt.-1.05) then
        call sm_ctype('blue')
        call sm_lweight(1.2)
       elseif(color.lt.-0.5) then
        call sm_ctype('cyan')
       elseif(color.lt.0.) then
        call sm_ctype('black')
        call sm_lweight(0.7)
       elseif(color.lt.0.5) then
        call sm_ctype('green')
       elseif(color.lt.1.) then
        call sm_ctype('yellow')
        call sm_lweight(5.0)
       elseif(color.lt.1.5) then
        call sm_ctype('magenta')
        call sm_lweight(0.7)
       else
        call sm_ctype('red')
        call sm_lweight(1.2)
       endif
      elseif(n_color.eq.3) then
       call sm_ctype('black')
      elseif(n_color.eq.4) then
       if(color.lt.-.2) then
        call sm_ctype('blue')
        call sm_lweight(1.2)
       elseif(color.lt.0.) then
        call sm_ctype('cyan')
       elseif(color.lt.0.5) then
        call sm_ctype('black')
        call sm_lweight(0.7)
       elseif(color.lt.1.) then
        call sm_ctype('green')
       elseif(color.lt.1.5) then
        call sm_ctype('yellow')
        call sm_lweight(5.0)
       elseif(color.lt.2.5) then
        call sm_ctype('magenta')
        call sm_lweight(0.7)
       else
        call sm_ctype('red')
        call sm_lweight(1.2)
       endif
      endif
      return
      end

      subroutine write_color(n_color)
      integer n_color
      real crange(6)
      character cname*8,label*23
      if(n_color.eq.0) then
       write(cname,'(a8)') 'P\\d\\gm'
       crange(1)= 1.01
       crange(2)= 0.75
       crange(3)= 0.5
       crange(4)= 0.25
       crange(5)= 0.
       crange(6)= -1.00
      elseif(n_color.eq.2) then
       crange(1)=-0.3
       crange(2)=-0.05
       crange(3)=0.3
       crange(4)=0.5
       crange(5)=0.85
       crange(6)=1.5
       write(cname,'(a8)') 'B-V'
      elseif(n_color.eq.3) then
       crange(1)=-1.05
       crange(2)=-0.5
       crange(3)=0.
       crange(4)=0.5
       crange(5)=1.
       crange(6)=1.5
       write(cname,'(a8)') 'U-B'
      elseif(n_color.eq.5) then
       crange(1)=0.
       crange(2)=0.
       crange(3)=0.
       crange(4)=0.
       crange(5)=0.
       crange(6)=0.
       write(cname,'(a8)') 'V-R'
      elseif(n_color.eq.4) then
       crange(1)=-0.2
       crange(2)=0.
       crange(3)=0.5
       crange(4)=1.
       crange(5)=1.5
       crange(6)=2.5
       write(cname,'(a8)') 'V-I'
      endif
      call sm_expand(0.6)
      write(label,'(a8,a8,a7)') 'color : ',cname,'              '
      call sm_relocate(0.,70.)
      call sm_label(label)   
      call sm_relocate(0.,60.)
      write(label,'(a12,f5.2,a3)') 'Blue     : < ',crange(1),'         '
      call sm_label(label)
      call sm_relocate(0.,50.)
      write(label,'(a12,f5.2,a3)') 'Cyan     : < ',crange(2),'         '
      call sm_label(label)
      call sm_relocate(0.,40.)
      write(label,'(a12,f5.2,a3)') 'Black    : < ',crange(3),'         '
      call sm_label(label)
      call sm_relocate(0.,30.)
      write(label,'(a12,f5.2,a3)') 'Green   : < ',crange(4),'         '
      call sm_label(label)
      call sm_relocate(0.,20.)
      write(label,'(a12,f5.2,a3)') 'Yellow   : < ',crange(5),'         '
      call sm_label(label)
      call sm_relocate(0.,10.)
      write(label,'(a12,f5.2,a3)') 'Magenta : < ',crange(6),'          '
      call sm_label(label)
      call sm_relocate(0.,0.)
      call sm_expand(1.)
      return
      end

      subroutine read_coo(input)
      parameter (ns=100000,nc=4)
      character H*1,input*32,line*99
      integer ih
      real rh
      common/coo_dat/no_coo(ns),ra,dec,nobs_coo(ns),ncoo
      real*8 ra(ns),dec(ns),ramin,ramax,cra,cdec,rmax,d,pi
      open(22,file=input)
      ncoo=0
      do i=1,ns
       read(22,'(t1,a1,t2,a99)',end=99) H,line
       if(H.ne.'#') then
        ncoo=ncoo+1    
        read(line,*) no_coo(ncoo),ra(ncoo),dec(ncoo),
     &  ih,ih,rh,ih,ih,rh,nobs_coo(ncoo)
        if(ncoo.eq.1) then
         ramin=ra(ncoo)
         ramax=ra(ncoo)
         decmin=dec(ncoo)
         decmax=dec(ncoo)
        else
         ramin=min(ra(ncoo),ramin)
         ramax=max(ra(ncoo),ramix)
         decmin=min(dec(ncoo),decmin)
         decmax=max(dec(ncoo),decmax)
        endif
       endif
      enddo
99    close(22)
      pi=3.14159265358979d0
      cra=(ramin+ramax)/2.d0
      cdec=(decmin+decmax)/2.d0
      rmax=0.d0 
      do i=1,ncoo
       d=sqrt((cos(dec(i)*pi/180.d0)*(cra-ra(i))/15.d0)**2
     & +(cdec-dec(i))**2)
       rmax=max(d,rmax)
      enddo
      print *,'R.A.min,max,  D.E.C. min,max,   Center,  Radius ='
      write(*,'(4f10.3)') ramin/15.d0,ramax/15.d0,decmin,decmax, 
     &cra/15.d0, cdec, rmax*60.d0
      return
      end
      subroutine read_combine()
      parameter (ns=100000)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      character H*1
      open(21,file='combine.dat')
      read(21,'(t1,a1)') H
      do i=1,ns
       read(21,*,end=99) no(i),(com(i,j),j=1,4),(nobs(i,j),j=1,4),p(i)
      enddo
99    nstar=i-1
      close(21)
      return
      end
