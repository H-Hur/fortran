c for only fortran77 to draw color-magnitude daigram.
c f77 -o 'cmd_t' 'cmd_t.f' -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11
      parameter (ns=100000,nf=10,nc=10)
      character input*16,output*16,input_cmd*16
      integer nmag,ncolor,nsize,nerror,nx,nir,nh,nneg,ncon,nzdraw
      real x1,x2,y1,y2,del_v,del_c
      write(input,'(a16)') 'Sung2001.ZAMS.cc'
      call read_ZAMS(input)


c      print *,'Input file name=?'
c      read(*,'(a16)') input_cmd            
      open(21,file='input.cmd')
      read(21,'(t1,a16)') input
      print *,input
      open(22,file=input)
      call passhead(22)
      call read_table(22)
      close(22) 
      
      do i=1,20
       read(21,'(t1,a16)') output
       read(21,*)nmag,ncolor,nsize,nerror,nx,nir,nh,nneg,nmp,x1,x2,y1,y2
     & ,del_v,del_c,nzdraw
       call draw_cmd(output,nmag,ncolor,nsize,nerror,nx,nir,nh,
     & nneg,nmp,x1,x2,y1,y2,del_v,del_c,nzdraw)
       read(21,*) ncon
       if(ncon.eq.0) then
        goto 100
       endif
      enddo
100   close(21) 
c      write(output,'(a16)') 'V.B-V.V16.ps    '
c      call draw_cmd(output,1,3,1,1,0,0,0,1,-0.2,2.,16.,7.)

c     write(output,'(a16)') 'V.V-I.ps        '
c     call draw_cmd(output,1,2,2,0,2,1,1,1,-0.2,3.5,20.,6.)

c     write(output,'(a16)') 'V.V-I.x.ps      '
c     call draw_cmd(output,1,2,2,2,1,0,0,1,-0.2,3.5,20.,6.)

c     write(output,'(a16)') 'V.B-V.ps        '
c     call draw_cmd(output,1,3,2,0,2,1,1,1,-0.2,2.5,20.,6.)

c     write(output,'(a16)') 'V.R-Ha.ps       '
c     call draw_cmd(output,1,5,2,0,2,1,1,1,-5.5,-2.5,20.,6.)


      stop
      end
!     nmag=1 -> mag : V
!     nmag=2 -> mag : I

!     ncolor : col(ns,ncolor)

!     nsise=1 -> A4
!     nsize=2 -> B5

!     nerror=0 -> draw all     
!     nerror=1 -> draw combine error> 0.1 only

!     nx=0  -> no mark x-ray
!     nx=1  -> mark x-ray XS only
!     nx=2  -> makr x-ray xc also

!     nir=1 -> mark ir-excess
      subroutine draw_cmd(output,nmag,ncolor,nsize,nerror,nx,nir,nh,
     &nneg,nmp,x1,x2,y1,y2,del_v,del_c,nzdraw)
      parameter (ns=100000,nf=10,nc=10)
      common/dat/id(ns),ra(ns),dec(ns),aV(ns),aI(ns),
     & col(ns,5),eV(ns),eI(ns),cole(ns,5),nobs(ns,7),xex,iex,double,Haex
     & ,mp(ns),pname,aJ(ns),aJH(ns),aHK(ns),ej(ns),ejh(ns),ehk(ns),qir,
     & name2mass,nameHD,ids(ns),nstar
      common /zams/zMv(34),zVI(34),zBV(34),zUB(34)
      character name2mass(ns)*16,qir(ns)*3,xex(ns)*1,iex(ns)*1
      character double(ns)*1,Haex(ns)*1,pname(ns)*8,nameHD(ns)*8
      character output*16,psname*29
      character xlabel*10,ylabel*10
      real amag(ns),color(ns),pp(1),x1,x2,y1,y2,redV(34),redc(34)
      real del_v,del_c
      real*8 e1,e2,e(ns)
      integer mag,ncolor,nerror,nx,nir,nsize,nh,nneg,nmp,nzdraw
      do i=1,34
       redV(i)=zMv(i)+del_v
       if(ncolor.eq.3) then
        redc(i)=zBV(i)+del_c
       elseif(ncolor.eq.2) then
        redc(i)=zVI(i)+del_c
       elseif(ncolor.eq.4) then
        redc(i)=zUB(i)+del_c
       endif
      enddo       
      write(psname,'(a13,a16)') 'postportfile ',output
      call sm_device (psname)
      call sm_graphics
      call sm_erase()
      if(nsize.eq.1) then
       call sm_location(3000,30000,3000,30000)
      elseif(nsize.eq.2) then
       call sm_location(3000,25000,4000,26000)
      endif
      call sm_limits(x1,x2,y1,y2)
      call sm_ctype('black')
      call sm_expand(1.)
      call sm_box(1,2,0,0)
      if(nmag.eq.1) then
       write(ylabel,'(a10)') 'V         '
       do i=1,nstar
        amag(i)=aV(i)
        color(i)=col(i,ncolor)
        e1=eV(i)
        e2=cole(i,ncolor)
        call cal_error(e1,e2,e(i))
       enddo
      elseif(nmag.eq.2) then
       write(ylabel,'(a10)') 'I        '
       do i=1,nstar
        amag(i)=aI(i)
        color(i)=col(i,ncolor)
        e1=eI(i)
        e2=cole(i,ncolor)
        call cal_error(e1,e2,e(i))
       enddo
      endif
      call w_label(ncolor,xlabel)
      call sm_xlabel(xlabel)
      call sm_ylabel(ylabel)
      if(nmp.eq.1) then
       do i=1,nstar
        if(mp(i).gt.80) then
         pp(1)=993.+0.9
         call sm_expand(2.)
         call sm_points(color(i),amag(i),1)
         call sm_expand(1.)
        endif
       enddo
       goto 900
      endif


      do i=1,nstar
       if(id(i).lt.0.and.nneg.eq.1.and.amag(i).ne.0.
     & .and.color(i).ne.0.) then
c        if((xex(i).ne.'X').and.(Iex(i).ne.'I').and.(Haex(i).ne.'H').and.
c     &  (nx.ne.0).and.(nir.ne.0).and.(nh.ne.0)) then
         pp(1)=52.+0.9
         call sm_expand(2.)
         call sm_ptype(pp,1)
         call sm_points(color(i),amag(i),1)
         call sm_expand(1.)
c        endif
       endif

       if(nobs(i,nmag).ne.0.and.nobs(i,ncolor+2).ne.0) then
        if((nerror.eq.1).and.(e(i).ge.0.d0).and.(e(i).lt.0.1d0)
     &  .and.(id(i).gt.0).and.(ids(i).ge.0)) then

         if(nmp.ge.2.and.mp(i).gt.70) then
          pp(1) = 993. + 0.6
          call sm_ptype(pp,1)
          if(nmp.eq.2) then
           call sm_ctype('blue')
          endif
          call sm_points(color(i),amag(i),1)
          call sm_ctype('black')
         elseif(xex(i).eq.' '.and.Iex(i).eq.' '.and.Haex(i).eq.' ')
     &   then
          pp(1) = 991. + 0.1
          call sm_ptype(pp,1)
          call sm_points(color(i),amag(i),1)
         endif
        endif
        if((nx.ne.0).and.(xex(i).eq.'X').and.(e(i).ge.0.d0)
     &  .and.(e(i).lt.0.1d0).and.(ids(i).ge.0)) then
         if(id(i).gt.0) then
          if(Haex(i).eq.'H') then
           pp(1) = 61. + 0.7
           call sm_ctype('red')
          elseif(Iex(i).eq.'I') then
           pp(1) = 43.+ 0.7
           call sm_ctype('red')
          else
           pp(1) = 41. + 0.5
           call sm_ctype('blue')
          endif
           call sm_expand(2.)
           call sm_ptype(pp,1)
           call sm_lweight(3.)
           call sm_points(color(i),amag(i),1)
           call sm_lweight(1.)
           call sm_ctype('black')
           call sm_expand(1.)
         endif
        endif
        if((nir.eq.1).and.(iex(i).eq.'I').and.(e(i).ge.0.d0)
     &  .and.(e(i).lt.0.1d0).and.(ids(i).ge.0).and.(xex(i).eq.' ')
     &  .and.mp(i).lt.70)then
         if(id(i).gt.0) then
          pp(1) = 990. + 0.7
          call sm_ctype('red')
          call sm_ptype(pp,1)
          call sm_lweight(2.)
          call sm_points(color(i),amag(i),1)
          call sm_ctype('black')
          call sm_expand(1.)
          call sm_lweight(1.)
         endif
        endif

        if((nh.eq.1).and.(Haex(i).eq.'H').and.(e(i).ge.0.d0)
     &  .and.(e(i).lt.0.1d0).and.(ids(i).ge.0).and.(xex(i).eq.' ')
     &  .and.mp(i).lt.70) then
         if(id(i).gt.0) then
          pp(1) = 40. + 0.8
          call sm_expand(2.5)
          call sm_ctype('magenta')
          call sm_ptype(pp,1)
          call sm_expand(1.)
          call sm_lweight(2.)
          call sm_points(color(i),amag(i),1)
          call sm_lweight(1.)
          call sm_expand(1.)
          call sm_ctype('black')
         endif
        endif
       endif
500   enddo

      if(nzdraw.eq.1) then     
       call sm_lweight( 4.0 )
c       call sm_ctype('blue')
       call sm_conn(redc,redV,34)
c       call sm_ctype('black')
       call sm_lweight( 0. )
      endif

900   call sm_gflush()
      call sm_hardcopy
      call sm_alpha


      return
      end 

      subroutine w_label(ncolor,label)
      character label*10
      integer ncolor
      if(ncolor.eq.1) then
       write(label,'(a10)') 'R-I       '
      elseif(ncolor.eq.2) then
       write(label,'(a10)') 'V-I       '
      elseif(ncolor.eq.3) then
       write(label,'(a10)') 'B-V       '
      elseif(ncolor.eq.4) then
       write(label,'(a10)') 'U-B       '
      elseif(ncolor.eq.5) then
       write(label,'(a10)') 'R-Ha      '
      endif
      return
      end

      subroutine read_table(input)
      integer input
      parameter (ns=100000)
      common/dat/id(ns),ra(ns),dec(ns),aV(ns),aI(ns),
     & col(ns,5),eV(ns),eI(ns),cole(ns,5),nobs(ns,7),xex,iex,double,Haex
     & ,mp(ns),pname,aJ(ns),aJH(ns),aHK(ns),ej(ns),ejh(ns),ehk(ns),qir,
     & name2mass,nameHD,ids(ns),nstar
      character name2mass(ns)*16,qir(ns)*3,xex(ns)*1,iex(ns)*1
      character double(ns)*1,Haex(ns)*1,pname(ns)*8,nameHD(ns)*8
1     format(t1,i6,2f8.3,7f7.3,7f6.3,7i2,1x,4a1,1x,i2,1x,a8,6f7.3,
     &1x,a3,1x,a16,1x,a8,13x,i6)
      do i=1,ns
       ids(i)=0
       read(input,1,end=100) id(i),ra(i),dec(i),aV(i),aI(i)
     & ,(col(i,j),j=1,5),eV(i),eI(i),(cole(i,j),j=1,5),(nobs(i,j),j=1,7)
     & ,double(i),xex(i),iex(i),Haex(i)
     & ,mp(i),pname(i)
     & ,aJ(i),ajh(i),ahk(i),ej(i),ejh(i),ehk(i),qir(i),name2mass(i)
     & ,nameHD(i)!,ids(i)
      enddo
100   nstar=i-1
      return
      end

      subroutine passhead(ninst)
      parameter (ns=100000)
      integer ninst
      character H*2,H1
10    format(t1,a2,t1,a1)
      do i=1,ns
       read(ninst,10) H,H1
       if(H.eq.'* ') then 
        goto 100
       elseif(H.eq.'#*') then
        goto 100
       elseif(H.eq.'--') then
        goto 100
       elseif(H1.eq.'*') then
        goto 100
       endif
      enddo
100   return
      end

!     combine two errors.     
      subroutine cal_error(e1,e2,e3)
      real*8 e1,e2,e3
      if(e1.lt.0.d0.or.e2.lt.0.d0) then
       e3=-1.d0
      else
       e3=sqrt(e1**2.d0+e2**2.d0)
      endif
      return
      end

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
12    y=yt1(i)
      goto 14
10    continue
11    y=yt1(k)+(yt1(k+1)-yt1(k))*((x-xt1(k))/(xt1(k+1)-xt1(k)))
14    return
      end

      subroutine read_ZAMS(input)
      character H*1,input*16
      common /zams/zMv(34),zVI(34),zBV(34),zUB(34)
      open(31,file=input)
      read(31,'(a1)') H
      do i=1,34
       read(31,*,end=100) zMv(i),zVI(i),zBV(i),zUB(i)
      enddo
100   nzams=i-1
      close(31)
      return
      end
