c for only fortran77 to draw color-color diagram.
c f77 -o ccd_o.exe /home/gjgusdh/fortran/subroutine/aex/ccd_o.f -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 
c     Last updated 06,30,2008

      parameter(ns=100000,ni=30)
      integer ID(ns,ni),nim(ni),nimage,yaxis1,xaxis1,xaxis2,nstar,nc
      integer yaxis2,IDW(ns),inf_out(ns,ni),answer,nreden
      integer num_f(ni),num_c(ni),ntransform
      real*8 mag(ns,ni),merr(ns,ni),inmag(ns,ni),xc(ns,ni),yc(ns,ni)
      real*8 meany1,meany2,meanx1,meanx2
      real pp(1),lx,mx,ly,my,x(ns),y(ns),xw(ns),yw(ns),ebv(ns)
      real zeropoint(ni)
      character input*16,output*23,fil(ni)*10,ch(10),filter(ni)*1
      character filter1*1,filter2*1,filtery1*1,filtery2*1,xlabel*5
      character ylabel*5,output2*14,Ha*2,H*1
      write(Ha,'(a2)') 'Ha'
      print *,'Input file name = ? ( 1 = out1.aex )'
      read(*,'(a16)')input
      if(input.eq.'1') then
       print *,'Input file name = out1.aex'
       open(21,file='out1.aex')
      else
       open(21,file=input)
      endif
      do nhead=1,50
       read(21,51) H,(fil(i),i=1,30)
51     format(t1,a1,t11,30a10)
       read(fil(1),'(10a1)') (ch(j),j=1,10)
       if(H.eq.'*') then
        goto 50
       elseif(H.eq.'F') then
        do i=1,30
         read(fil(i),'(10a1)') (ch(j),j=1,10)
         nim(i)=0
         k=j-1
         do j=1,k
          if(ch(j).eq.'I') then
           nim(i)=1
           write(filter(i),'(a1)') 'I'
          elseif(ch(j).eq.'V') then
           nim(i)=1
           write(filter(i),'(a1)') 'V'
          elseif(ch(j).eq.'B') then
           nim(i)=1
           write(filter(i),'(a1)') 'B'
          elseif(ch(j).eq.'U') then
           nim(i)=1
           write(filter(i),'(a1)') 'U'
          elseif(ch(j).eq.'R') then
           nim(i)=1
           write(filter(i),'(a1)') 'R'
          elseif(ch(j).eq.'H') then
           nim(i)=1
           write(filter(i),'(a1)') 'H'
          elseif(ch(j).eq.'G') then
           nim(i)=1
           write(filter(i),'(a1)') 'G'
          elseif(ch(j).eq.'J') then
           nim(i)=1
           write(filter(i),'(a1)') 'J'
          elseif(ch(j).eq.'K') then
           nim(i)=1
           write(filter(i),'(a1)') 'K'
          elseif(ch(j).eq.'L') then
           nim(i)=1
           write(filter(i),'(a1)') 'L'
          endif
         enddo
        enddo
       endif
      enddo

50    do i=1,30
       if(nim(i).eq.0) then
        goto 60
       endif
      enddo
60    nimage=i-1
      print *,'Input filters for x-axis color =? (Ha=H)'
      read(*,'(a1)') filter1
      read(*,'(a1)') filter2
      print *,'Input filters fot Y-axis color  =? (Ha=H)'
      read(*,'(a1)') filtery1
      read(*,'(a1)') filtery2

      do i=1,nimage
       if(filter1.eq.filter(i)) then
        xaxis1=i
       elseif(filter2.eq.filter(i)) then 
        xaxis2=i
       endif
       if(filtery1.eq.filter(i)) then
        yaxis1=i
       elseif(filtery2.eq.filter(i)) then
        yaxis2=i
       endif
      enddo
      meany1=0.d0
      meany2=0.d0
      meanx1=0.d0
      meanx2=0.d0
c     print *,xaxis1,xaxis2,yaxis,nimage
c      read(21,'(a13)') output
70    do i=1,ns
       read(21,*,end=100) i1,(ID(i,j),j=1,nimage)
       read(21,*,end=100) r1,(inf_out(i,j),j=1,nimage)
       read(21,*,end=100) r1,(mag(i,j),j=1,nimage)
       read(21,*,end=100) r1,(inmag(i,j),j=1,nimage)
       read(21,*,end=100) r1,(merr(i,j),j=1,nimage)
c       read(21,*,end=100) (xc(i,j),j=1,nimage)
c       read(21,*,end=100) (yc(i,j),j=1,nimage)
       read(21,'(a19)',end=100) output
       read(21,'(a19)',end=100) output
       read(21,'(a19)',end=100) output
       read(21,'(a19)',end=100) output
       read(21,'(a19)',end=100) output
      enddo

100   write(xlabel,'(5a1)') filter1,' ','-',' ',filter2
      write(ylabel,'(5a1)') filtery1,' ','-',' ',filtery2
      if(filter1.eq.'R'.and.filter2.eq.'H') then
       write(xlabel,'(4a1,a2)') filter1,' ','-',' ',Ha
      elseif(filtery1.eq.'R'.and.filtery2.eq.'H') then
       write(ylabel,'(4a1,a2)') filtery1,' ','-',' ',Ha
      endif
      nstar=i-1
      nc=0
      do i=1,nstar
       if(mag(i,xaxis1).gt.0.d0.and.mag(i,xaxis2).gt.0.d0.and.
     * mag(i,yaxis1).gt.0.d0.and.mag(i,yaxis2).gt.0.d0) then
        nc=nc+1
        meanx1=meanx1+mag(i,xaxis1)
        meanx2=meanx2+mag(i,xaxis2)
        meany1=meany1+mag(i,yaxis1)
        meany2=meany2+mag(i,yaxis2)
       endif
      enddo

      meanx1=meanx1/real(nc)
      meanx2=meanx2/real(nc)
      meany1=meany1/real(nc)
      meany2=meany2/real(nc)
      lx=meanx1-meanx2
      mx=lx
      ly=meany1-meany2
      my=ly    

      if(filter1.eq.'R'.and.filter2.eq.'H') then
       print *,'Input coeff file=?'
       read(*,'(a16)') input
       call readhead(input)
       call readcoeff(input,num_f,num_c,ntransform,zeropoint)    
      elseif(filtery1.eq.'R'.and.filtery2.eq.'H') then
       print *,'Input coeff file=?'
       read(*,'(a16)') input
       call readhead(input)
       call readcoeff(input,num_f,num_c,ntransform,zeropoint)
      endif


      do i=1,nstar
       if(mag(i,xaxis1).gt.0.d0.and.mag(i,xaxis2).gt.0.d0.and.
     * mag(i,yaxis1).gt.0.d0.and.mag(i,yaxis2).gt.0.d0) then
        x(i)=mag(i,xaxis1)-mag(i,xaxis2)
        y(i)=mag(i,yaxis1)-mag(i,yaxis2)
        if(filter1.eq.'R'.and.filter2.eq.'H') then
         do j=1,ntransform
          if(num_f(j).eq.4) then
           x(i)=mag(i,xaxis1)-zeropoint(j)-mag(i,xaxis2)
          endif
         enddo
        elseif(filtery1.eq.'R'.and.filtery2.eq.'H') then
         do j=1,ntransform
          if(num_f(j).eq.4) then
           y(i)=mag(i,yaxis1)-zeropoint(j)-mag(i,yaxis2)
          endif
         enddo
        endif
       endif
      enddo
c      print *,lx,mx,ly,my
      write(output,101) 'postlandfile ',filtery1,'-',filtery2,'.',
     *filter1,'-',filter2,'.ps'
      if(filter1.eq.'R'.and.filter2.eq.'H') then
       write(output,'(a13,6a1,a2,a3)') 'postportfile',filtery1,'-',
     & filtery1,'.',filter1,'-',Ha,'.ps'
      endif
101   format(a13,7a1,a3)
      write(output2,102) filter1,'_',filtery1,'-',filtery2,'.',
     *filter1,'-',filter2,'_.ccd'
      if (filtery1.eq.'R'.and.filtery2.eq.'H') then
       write(output2,'(4a1,a2,4a1,a5)') filter1,'_',filtery1,'-',Ha,'.',
     * filter1,'-',filter2,'_.ccd'
      endif
102   format(9a1,a5)
      open(31,file=output2)
      write(31,*) 'Catalog of color-color diagram ',output,'.'
      write(31,*) ' '
      write(31,*) 'ID(',filter1,')   ',filter1,'-',filter2,'     ',
     *filtery1,'-',filtery2
      call sm_device (output)
      call sm_graphics
      call sm_erase()
      call sm_location(3000,30000,3000,30000)
      print *,'X-axis range=?'
      read(*,*) lx,mx
      print *,'Y-axis range=?'
      read(*,*) ly,my
      call sm_limits(lx,mx,ly,my)
      call sm_ctype('black')
      call sm_box(1,2,1,2)
      call sm_expand(0.8)
      call sm_xlabel(xlabel)
      call sm_expand(1.)
      call sm_ylabel(ylabel)
      call sm_expand(1.)
      pp(1)=403.+0.15
      call sm_ptype(pp,1)
      nprint=0
      do i=1,nstar
       if(mag(i,xaxis1).gt.0.d0.and.mag(i,xaxis2).gt.0.d0.and.
     * mag(i,yaxis1).gt.0.d0.and.mag(i,yaxis2).gt.0.d0) then
        call sm_points(x(i),y(i),1)
        
c        write(31,110) ID(i,xaxis1),x(i),y(i)
        nprint=nprint+1
        xw(nprint)=x(i)
        yw(nprint)=y(i)
        IDW(nprint)=ID(i,xaxis1)
       endif
      enddo
      call sort(xw,yw,IDw,nprint)
      do i=1,nprint 
        write(31,110) IDW(i),xw(i),yw(i)
      enddo
110   format(I6,2f8.3)

      print *,'Do yo want to draw the ZAMS line ? ( 1=yes,0=no)'
      read(*,*) answer
      if(answer.eq.1) then
       write(input,'(a16)') 'Sung2001.ZAMS.cc'
       call readhead(input)
       call read_ZAMS(input)
       call draw_zams(filter1,filter2,filtery1,filtery2,0.,0)
       print *,'Number of redden line of ZAMS?'
       read(*,*) nreden
       if(nreden.eq.0) then
        goto 200
       endif
       do i=1,nreden
        print *,'e(B-V)=? (end -> type negative number) '
        read(*,*) ebv(i)
        if(ebv(i).lt.0.) then
         goto 200
        endif 
        call draw_zams(filter1,filter2,filtery1,filtery2,ebv(i),1)
       enddo
      endif
      goto 200
200   call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      write(31,*) nprint,' stars draw on the ',output,'.'
      write(31,*) ' '
      close(31)
      close(21)
      print *,nprint,' stars draw on the ',output,'.'
      print *,'drawn catalog wrote on the ',output2,'.'
      stop
      end

     

      subroutine readcoeff(coeff,num_f,num_c,ntrans,zeropoint)
      parameter (nf=30,ns=100000)
      character coeff*16,H*1
      real zeropoint(nf),coeff1(nf),coeff2(nf)
      integer ntrans,num_f(nf),num_c(nf),numf
      common /co/num_fil(nf),num_col(nf),transcoeff(nf),
     &delta_ut(nf),colorran(nf,2),timeran(nf,2)
     &,timevar(nf),transzero(nf),timezero(nf),Lfilter(nf)
      common/dat1/dmag(nf,ns),dmerr(nf,ns),ut(nf,ns),air(nf,ns),
     &weight(nf,ns),nall
      common/dat2/V(nf,ns),BV(nf,ns),VR(nf,ns),VI(nf,ns),UB(nf,ns),
     &nstar(nf),nfilter,mfilter(nf)
      open(31,file=coeff)
      ntrans=0
      numf=1

      do i=1,50
       read(31,'(a1)',end=500) H
       if(H.eq.'@'.or.H.eq.'#') then
        goto 100
        else
        backspace(31)
        ntrans=ntrans+1
        read(31,*) num_fil(ntrans),num_col(ntrans),colorran(ntrans,1),
     &  colorran(ntrans,2),transcoeff(ntrans),transzero(ntrans),
     &  coeff1(ntrans),coeff2(ntrans),timevar(ntrans),timeran(ntrans,1),
     &  timeran(ntrans,2),timezero(ntrans),delta_ut(ntrans),
     &  zeropoint(ntrans)
        Lfilter(mfilter(ntrans))=ntrans
        num_f(ntrans)=num_fil(ntrans)
        num_c(ntrans)=num_col(ntrans)
        num_f(ntrans)=num_fil(ntrans)
        num_c(ntrans)=num_col(ntrans)
       endif
100   enddo
500   close(31)

      return
      end



      subroutine draw_ZAMS(c1,c2,y1,y2,reden,ltype)
      parameter (ns=1000)
      character c1*1,c2*1,y1*1,y2*1
      integer ltype
      real x(ns),y(ns),slope,reden,dx,dy
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),zVR(ns),zRI(ns),nzams
      dx=0.
      slope=0.
      dy=0.
      if (c1.eq.'B'.or.c1.eq.'b'.and.c2.eq.'V'.or.c2.eq.'v'.and.
     &y1.eq.'U'.or.y1.eq.'u'.and.y2.eq.'B'.or.y2.eq.'b') then
       slope=0.72
       dx=reden
       dy=reden*slope
      endif 


      if(c1.eq.'U'.or.c1.eq.'u'.and.c2.eq.'B'.or.c2.eq.'b') then
       do i=1,nzams
        x(i)=zUB(i)+dx
       enddo
       goto 100
      elseif(c1.eq.'B'.or.c1.eq.'b'.and.c2.eq.'V'.or.c2.eq.'v') then
       do i=1,nzams
        x(i)=zBV(i)+dx
       enddo
       goto 100
      elseif(c1.eq.'V'.or.c1.eq.'v'.and.c2.eq.'I'.or.c2.eq.'i') then
       do i=1,nzams
        x(i)=zVI(i)+dx
       enddo
       goto 100
      elseif(c1.eq.'V'.or.c1.eq.'v'.and.c2.eq.'R'.or.c2.eq.'r') then
       print *,'This color cannot draw ZAMS line.'
       goto 900
      elseif(c1.eq.'R'.or.c1.eq.'R'.and.c2.eq.'I'.or.c2.eq.'i') then
       print *,'This color cannot draw ZAMS line.'
       goto 900
      elseif(c1.eq.'R'.or.c1.eq.'R'.and.c2.eq.'h') then
       print *,'This color cannot draw ZAMS line.'
       goto 900
      endif

100   if(y1.eq.'U'.or.y1.eq.'u'.and.y2.eq.'B'.or.y2.eq.'b') then
       do i=1,nzams
        y(i)=zUB(i)+dy
       enddo
       goto 200
      elseif(y1.eq.'B'.or.y1.eq.'b'.and.y2.eq.'V'.or.y2.eq.'v') then
       do i=1,nzams
        y(i)=zBV(i)+dy
       enddo
       goto 200
      elseif(y1.eq.'V'.or.y1.eq.'v'.and.y2.eq.'I'.or.y2.eq.'i') then
       do i=1,nzams
        y(i)=zVI(i)+dy
       enddo
       goto 200
      elseif(y1.eq.'V'.or.y1.eq.'v'.and.y2.eq.'R'.or.y2.eq.'r') then
       print *,'This color cannot draw ZAMS line.'
       goto 900
      elseif(y1.eq.'R'.or.y1.eq.'R'.and.y2.eq.'I'.or.y2.eq.'i') then
       print *,'This color cannot draw ZAMS line.'
       goto 900
      elseif(y1.eq.'R'.or.y1.eq.'R'.and.y2.eq.'h') then
       print *,'This color cannot draw ZAMS line.'
       goto 900
      endif

200   call sm_lweight( 3.0 )
      call sm_ltype(ltype)
      call sm_conn(x,y,nzams)


900   return
      end




      subroutine read_ZAMS(input)
      parameter (nf=10,ns=1000,nmcol=20)
      character line*200,H(200)*1,input*16,coll(nmcol)*32
      character Hv(32)*1,lastcoll1*32
      integer numcol,nc1,nc2                
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





      subroutine sort(x,y,ID,nstar)
      parameter (ns=100000,ni=30)
      integer ID(ns),IDs
      real x(ns),y(ns),xs,ys
      do i=1,nstar
       do j=1,nstar
        if (x(i).gt.x(j)) then
          xs=x(i)
          ys=y(i)
          IDs=ID(i)
          x(i)=x(j)    
          y(i)=y(j)
          ID(i)=ID(j)
          x(j)=xs
          y(j)=ys
          ID(j)=IDs
        endif
       enddo
      enddo
      return
      end 
      
