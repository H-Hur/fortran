c for only fortran77 to draw color-magnitude daigram.
c f77 -o 'ccd_t' 'ccd_t.f' -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11
      parameter (ns=100000,nf=10,nc=10)
      character input*16,output*16,input_ccd*16
      integer ncolorx,ncolory,nsize,nerror,nx,nir,nh,nz,nneg,ncon,nmp
      real x1,x2,y1,y2 
c      print *,'Input file name=?'
c      read(*,'(a16)') input_ccd            
      write(input,'(a16)') 'tr16_table_n.dat'
      open(21,file='input.ccd')

      read(21,'(a16)') input
      open(22,file=input)
      call passhead(22)
      call read_table(22)
      close(22) 

c      write(output,'(a16)') 'U-B.B-V.ps      '
c      call draw_ccd(output,3,4,2,0,2,1,1,-0.2,2.5,2.,-1.,1)
c      print *,'Input file name=?'
c      read(*,'(a16)') input     
c      open(22,file=input)
      do i=1,20
       print *,'0'
       read(21,'(t1,a16)') output
       print *,output
       read(21,*) ncolorx,ncolory,nsize,nerror,nx,nir,nh,nneg,
     & x1,x2,y1,y2,nz,nmp
       print *,ncolorx,ncolory,nsize,nerror,nx,nir,nh,nneg,
     & x1,x2,y1,y2,nz
       call draw_ccd(output,ncolorx,ncolory,nsize,nerror,nx,nir,
     & nh,nneg,x1,x2,y1,y2,nz,nmp)
       goto 100
       read(21,*) ncon
       print *,'ncon',ncon
       if(ncon.eq.0) then
        print *,'2'
        goto 100
       endif
       print *,'3'      
      enddo
100   close(21) 

c      write(output,'(a16)') 'R-Ha.V-I.ps      '
c      call draw_ccd(output,2,5,2,1,0,0,1,-0.,4.5,-5.5,-3.5,2)


c      write(output,'(a16)') 'V.B-V           '
c      call draw_cmd(output,1,3,2,0,2,1,-0.2,2.5,20.,7.)


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
!     nx=1  -> mark x-ray X only

!     nir=1 -> mark ir-excess

!     nh=1 -> mark Ha excess stars

!     nz=0 -> Do ont draw ZAMS
!     nz=1 -> Draw zams - U-B:B-V
!     nz=2 -> Draw zams - R-Ha:V-I

      subroutine draw_ccd(output,ncolorx,ncolory,nsize,nerror,nx,nir,
     &nh,nneg,x1,x2,y1,y2,nz,nmp)
      parameter (ns=100000,nf=10,nc=10)
      character output*16,psname*29
      character xlabel*10,ylabel*10,zamsfile*16,datafile*16
      common/dat/id(ns),ra(ns),dec(ns),aV(ns),aI(ns),
     & col(ns,5),eV(ns),eI(ns),cole(ns,5),nobs(ns,7),xex,iex,double,Haex
     & ,mp(ns),pname,aJ(ns),aJH(ns),aHK(ns),ej(ns),ejh(ns),ehk(ns),qir,
     & name2mass,nameHD,ids(ns),nstar
      character name2mass(ns)*16,qir(ns)*3,xex(ns)*1,iex(ns)*1
      character double(ns)*1,Haex(ns)*1,pname(ns)*8,nameHD(ns)*8
      real xcolor(ns),ycolor(ns),pp(1),x1,x2,y1,y2,ebv(5),dx,dy,qqq,bv0
      real*8 e1,e2,e(ns)
      integer ncolorx,ncolory,nerror,nx,nir,nsize,nz,nreden,nh,nneg,nmp
      write(psname,'(a13,a16)') 'postlandfile ',output
c      call write_dataname(output,datafile)
c      open(21,file=datafile)
c      call write_datahead(21,ncolorx,ncolory)

      call sm_device (psname)
      call sm_graphics
      call sm_erase()
      if(nsize.eq.1) then
       call sm_location(3000,30000,3000,30000)
      elseif(nsize.eq.2) then
       call sm_location(4000,26000,9000,30000)
      endif
      call sm_limits(x1,x2,y1,y2)
      call sm_ctype('black')
      call sm_expand(1.)
      call sm_box(1,2,0,0)
      do i=1,nstar
       xcolor(i)=col(i,ncolorx)
       ycolor(i)=col(i,ncolory)
       e1=cole(i,ncolorx)
       e2=cole(i,ncolory)
       call cal_error(e1,e2,e(i))
      enddo
      call w_label(ncolorx,xlabel)
      call w_label(ncolory,ylabel)
      call sm_xlabel(xlabel)
      call sm_ylabel(ylabel)
1     format(t1,i6,3f7.3,1x,4a1)
      do i=1,nstar

       if(xcolor(i).ne.0..and.ycolor(i).ne.0) then

        if(id(i).lt.0.and.nneg.eq.1.and.xcolor(i).ne.0.
     & .and.ycolor(i).ne.0) then
         call sm_expand(2.)
         pp(1) = 52.+0.9
         call sm_ptype(pp,1)
         call sm_points(xcolor(i),ycolor(i),1)
         call sm_expand(1.)
        endif


        if((nerror.eq.1).and.(e(i).ge.0.d0).and.(e(i).lt.0.1d0)
     &  .and.(id(i).gt.0).and.(ids(i).ge.0)) then
         qqq=col(i,4)-0.72*col(i,3)
         if((nmp.eq.2.and.mp(i).gt.70).and.
     &   ((qqq.lt.-0.4.and.aV(i).lt.15.5).or.(col(i,4).lt.0.1.and.
     &   aV(i).lt.14.5))) then
          pp(1) = 993. + 0.6
          call sm_ptype(pp,1)
          call sm_ctype('blue')
          call sm_points(xcolor(i),ycolor(i),1)
          call sm_ctype('black')
         elseif(xex(i).eq.' '.and.Iex(i).eq.' '.and.Haex(i).eq.' ')
     &   then
          pp(1) = 991. + 0.1
          call sm_ptype(pp,1)
          call sm_points(xcolor(i),ycolor(i),1)
         endif
        elseif((nerror.eq.0).and.(e(i).ge.0.d0).and.(e(i).lt.0.1d0)
     &  .and.(id(i).gt.0).and.(ids(i).ge.0)) then
         pp(1) = 991. + 0.1
         call sm_ptype(pp,1)
         call sm_points(xcolor(i),ycolor(i),1)

        elseif((nerror.eq.0).and.(e(i).ge.0.d0)
     &  .and.(id(i).gt.0).and.(ids(i).ge.0)) then
         pp(1) = 991.+ 0.1
         call sm_ptype(pp,1)
         call sm_points(xcolor(i),ycolor(i),1)
       
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
           call sm_points(xcolor(i),ycolor(i),1)
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
          call sm_points(xcolor(i),ycolor(i),1)
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
          call sm_points(xcolor(i),ycolor(i),1)
          call sm_lweight(1.)
          call sm_expand(1.)
          call sm_ctype('black')
         endif
        endif
       endif
      enddo

      if(nz.eq.1) then
       write(zamsfile,'(a16)') 'Sung2001.ZAMS.cc'
       call readhead(zamsfile)
       call read_ZAMS(zamsfile)
       call draw_ZAMScc(ncolorx,ncolory,0.,0.,0,0,1)
       print *,'Drawing ',output
       print *,'Number of redden line of ZAMS?'
       read(*,*) nreden
       if(nreden.eq.0) then
        goto 200
       endif
       do i=1,nreden
        print *,'dx,dy of ZAMS line=? (end -> type negative number) '
        read(*,*) dx,dy
        if(dx.le.0..and.dy.le.0.) then
         goto 200
        endif
        call draw_zamscc(ncolorx,ncolory,dx,dy,1,0,1)
       enddo

      elseif(nz.eq.2) then
       write(zamsfile,'(a16)') 'Sung1997.ZAMS.cc'
       call read_Ha(zamsfile)
       call draw_ZAMScc(ncolorx,ncolory,0.,0.,0,0,2)
       print *,'Drawing ',output
       print *,'Number of redden line of ZAMS?'
       read(*,*) nreden
       if(nreden.eq.0) then
        goto 200
       endif
       do i=1,nreden
        print *,'dx,dy of ZAMS line=? (end -> type negative number) '
        read(*,*) dx,dy
        if(dx.le.0..and.dy.le.0.) then
         goto 200
        endif
        call draw_zamscc(ncolorx,ncolory,dx,dy,1,2,2)
       enddo
      endif

     

200   call sm_gflush()
      call sm_hardcopy
      call sm_alpha
c      close(21)
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
     & ,nameHD(i),ids(i)
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

      subroutine draw_ZAMScc(ncolorx,ncolory,dx,dy,ltype,nctype,nz)
      parameter (ns=1000)
      integer ltype,ncolorx,ncolory,nz,nctype
      real x(ns),y(ns),ebv,dx,dy
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),zVR(ns),zRI(ns),nzams
      common/zHa/zRH(ns),zrhVI(ns),nzHa
      character ctype*8
 
      if(ncolorx.eq.4) then
       do i=1,nzams
        x(i)=zUB(i)+dx
       enddo
       goto 100

      elseif(ncolorx.eq.3) then
       do i=1,nzams
        x(i)=zBV(i)+dx
       enddo
       goto 100

      elseif(ncolorx.eq.2) then
       do i=1,nzams
        x(i)=zVI(i)+dx
       enddo
       if(ncolory.eq.5) then
        do i=1,nzHa
         x(i)=zrhVI(i)+dx
        enddo 
       endif
       goto 100

      elseif(ncolorx.eq.5) then
       do i=1,nzams
        x(i)=zRH(i)+dx
       enddo

       goto 100
      endif

100   if(ncolory.eq.4) then
       do i=1,nzams
        y(i)=zUB(i)+dy
       enddo
       goto 200

      elseif(ncolory.eq.3) then
       do i=1,nzams
        y(i)=zBV(i)+dy
       enddo
       goto 200

      elseif(ncolory.eq.2) then
       do i=1,nzams
        y(i)=zVI(i)+dy
       enddo
       goto 200

      elseif(ncolory.eq.5) then
       do i=1,nzHa
        y(i)=zRH(i)+dy
       enddo
       goto 200

      endif

200   call select_ctype(nctype,ctype)
      call sm_ctype(ctype)
      if(nz.eq.1) then
       call sm_lweight( 3.0 )
       call sm_ltype(ltype)
       call sm_conn(x,y,nzams)
      elseif(nz.eq.2) then
       call sm_lweight( 3.0 )
       call sm_ltype(ltype)
       call sm_conn(x,y,nzHa)
       call sm_ctype(ctype)
      endif
      call sm_lweight( 0. )
      call sm_ctype('black')
900   return
      end
    
      subroutine select_ctype(nctype,ctype)
      integer nctype
      character ctype*8
      if(nctype.eq.0) then
       write(ctype,'(a8)') 'black   '
      elseif(nctype.eq.1) then
       write(ctype,'(a8)') 'write   '
      elseif(nctype.eq.2) then
       write(ctype,'(a8)') 'red     '
      elseif(nctype.eq.3) then
       write(ctype,'(a8)') 'green   '
      elseif(nctype.eq.4) then
       write(ctype,'(a8)') 'blue    '
      elseif(nctype.eq.5) then
       write(ctype,'(a8)') 'cyan    '
      elseif(nctype.eq.6) then
       write(ctype,'(a8)') 'magenta '
      elseif(nctype.eq.7) then
       write(ctype,'(a8)') 'yellow  '
      endif
      return
      end

      subroutine read_ZAMS(input)
      parameter (nf=10,ns=1000,nmcol=20)
      character line*200,H(200)*1,input*16,coll(nmcol)*32
      character Hv(32)*1,lastcoll1*32
      integer numcol,nc1,nc2,ncol,nc10,nc20
      common /nc/nc_filter,nc_star,nc_mag,nc_merr,nc_UT,nc_air,nc_weight
     &,ncollum,nc_V,nc_BV,nc_VR,nc_RI,nc_VI,nc_UB,nc_Mv
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),zVR(ns),zRI(ns),nzams
      do i=1,32
       write(Hv(i),'(a1)') ' '
      enddo
      open(31,file=input)
      write(lastcoll1,'(32a1)') (Hv(k),k=1,32)
      nzams=0
      do i=1,ns
       read(31,'(a200)',end=500) line
       nzams=nzams+1
       read(line,'(200a1)') (H(j),j=1,200)
       if(H(1).eq.'#') then
        nzams=nzams-1
        goto 100
       endif
       do j=1,16
        if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &  '^')then
         nc2=j
         goto 10
        endif
       enddo

10     do numcol=1,ncollum
        if(numcol.eq.1) then
         nc1=1
         do j=1,nc2+1+32
          if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &    '^')  then
           nc2=j
           goto 20
          endif
         enddo
        else
         do j=nc2+1,nc2+1+32
          if(H(j).eq.' '.and.H(j+1).ne.' ') then
           nc1=j+1
          endif
          if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &    '^')  then
           nc2=j
           goto 20
          endif
         enddo
        endif
20      write(coll(numcol),'(32a1)') (Hv(k),k=1,32)
        write(coll(numcol),'(32a1)') (H(k),k=nc1,nc2)
       enddo

       read(coll(nc_Mv),*) zMv(nzams)
       read(coll(nc_BV),*) zBV(nzams)
c       read(coll(nc_VR),*) zVR(nzams)
c       read(coll(nc_RI),*) zRI(nzams)
       read(coll(nc_VI),*) zVI(nzams)
       read(coll(nc_UB),*) zUB(nzams)

100   enddo
500   close(31)
      return
      end
      subroutine readhead(input)
      parameter (nf=10,ns=1000)
      character header*200,H(200)*1,input*16,coll*3
      common /nc/nc_filter,nc_star,nc_mag,nc_merr,nc_UT,nc_air,nc_weight
     &,ncollum,nc_V,nc_BV,nc_VR,nc_RI,nc_VI,nc_UB,nc_Mv
      open(31,file=input)
      read(31,'(a200)') header
      ncollum=0
      read(header,'(200a1)') (H(i),i=1,200)

      do i=1,200
       if(H(i).eq.' '.or.H(i).eq.'#') then
        if(H(i+1).ne.' ') then
         write(coll,'(3a1)') H(i+1),H(i+2),H(i+3)
         if(coll.eq.'Sta'.or.coll.eq.'sta') then
          ncollum=ncollum+1
          nc_star=ncollum
          elseif(coll.eq.'Mag'.or.coll.eq.'mag') then
          ncollum=ncollum+1
          nc_mag=ncollum
          elseif(coll.eq.'Mer'.or.coll.eq.'mer') then
          ncollum=ncollum+1
          nc_merr=ncollum
          elseif(coll.eq.'V  '.or.coll.eq.'v  ') then
          ncollum=ncollum+1
          nc_V=ncollum
          elseif(coll.eq.'B-V'.or.coll.eq.'b-v') then
          ncollum=ncollum+1
          nc_BV=ncollum
          elseif(coll.eq.'U-B'.or.coll.eq.'u-b') then
          ncollum=ncollum+1
          nc_UB=ncollum
          elseif(coll.eq.'V-R'.or.coll.eq.'v-r') then
          ncollum=ncollum+1
          nc_VR=ncollum
          elseif(coll.eq.'R-I'.or.coll.eq.'r-i') then
          ncollum=ncollum+1
          nc_RI=ncollum
          elseif(coll.eq.'V-I'.or.coll.eq.'v-i') then
          ncollum=ncollum+1
          nc_VI=ncollum
          elseif(coll.eq.'UT '.or.coll.eq.'ut '.or.coll.eq.'UT('.or.
     &    coll.eq.'ut(') then
          ncollum=ncollum+1
          nc_ut=ncollum
          elseif(coll.eq.'Air'.or.coll.eq.'air'.or.coll.eq.'X  '.or.
     &    coll.eq.'x  ') then
          ncollum=ncollum+1
          nc_air=ncollum
          elseif(coll.eq.'wei'.or.coll.eq.'Wei') then
          ncollum=ncollum+1
          nc_weight=ncollum
          elseif(coll.eq.'fil'.or.coll.eq.'Fil') then
          ncollum=ncollum+1
          nc_filter=ncollum
          elseif(coll.eq.'Mv') then
          ncollum=ncollum+1
          nc_Mv=ncollum
         endif
        endif
       endif
      enddo
      close(31)
      return
      end
      
      subroutine read_Ha(input)
      parameter(ns=1000)
      character input*16,H*1
      common/zHa/zRH(ns),zrhVI(ns),nzHa
      open(21,file=input)
      nzHa=0
      do i=1,ns
       read(21,'(t1,a1)',end=100) H
       if(H.ne.'#') then
        backspace(21)
        nzHa=nzHa+1
        read(21,*) zrhVI(nzHa),zRH(nzHa)
       endif
      enddo
100   close(21)
      return
      end

      subroutine write_dataname(dname,dataname)
      character dname*10, dataname*16, ch(16)*1
      read(dname,'(10a1)') (ch(i),i=1,10)
      do i=1,10
       if(ch(i).eq.'.') then
        goto 100
       endif
      enddo
100   write(ch(i+1),'(a1)') 'c'
      write(ch(i+2),'(a1)') 'c'
      write(ch(i+3),'(a1)') 'd'
      write(ch(i+4),'(a1)') '_.'
      write(ch(i+5),'(a1)') 't'
      write(dataname,'(16a1)') (ch(j),j=1,i+5)
      return
      end

      subroutine write_datahead(nopen,ncolorx,ncolory)
      character xc*4,yc*4
      integer ncolorx,ncolory,nopen

      if(ncolorx.eq.1) then
       write(xc,'(a4)') 'R-I '
      elseif(ncolorx.eq.2) then
       write(xc,'(a4)') 'V-I '
      elseif(ncolorx.eq.3) then
       write(xc,'(a4)') 'B-V '
      elseif(ncolorx.eq.4) then
       write(xc,'(a4)') 'U-B '
      elseif(ncolorx.eq.5) then
       write(xc,'(a4)') 'R-Ha'
      endif

      if(ncolory.eq.1) then
       write(yc,'(a4)') 'R-I '
      elseif(ncolory.eq.2) then
       write(yc,'(a4)') 'V-I '
      elseif(ncolory.eq.3) then
       write(yc,'(a4)') 'B-V '
      elseif(ncolory.eq.4) then
       write(yc,'(a4)') 'U-B '
      elseif(ncolory.eq.5) then
       write(yc,'(a4)') 'R-Ha'
      endif
 
      write(nopen,'(t1,a13,2x,a4,2x,a4)')'#   ID    V  ',xc,yc 

      return
      end
