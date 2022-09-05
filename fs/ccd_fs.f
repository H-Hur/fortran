c for only fortran77 to draw color-color diagram.
c f77 -o ccd /home/gjgusdh/fortran/subroutine/fs/ccd_fs.f -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 
      parameter(ns=100000,ni=30)
      integer ID(ns,ni),nim(ni),nimage,yaxis1,xaxis1,xaxis2,nstar,nc
      integer ID2(ns,ni),yaxis2,IDW(ns)
      real*8 mag(ns,ni),merr(ns,ni)
      real*8 meany1,meany2,meanx1,meanx2
      real pp(1),lx,mx,ly,my,x(ns),y(ns),xw(ns),yw(ns)
      character input*16,output*23,fil(ni)*10,ch(10),filter(ni)*1
      character filter1*1,filter2*1,filtery1*1,filtery2*1,xlabel*5
      character ylabel*5,output2*14
      print *,'Input file name = ? ( 1 = wcout1.fs )'
      read(*,'(a16)')input
      if(input.eq.'1') then
       print *,'Input file name = wcout1.fs'
       open(21,file='wcout1.fs')
      else
       open(21,file=input)
      endif
      do nhead=1,10
       read(21,51) (fil(i),i=1,30)
51     format(30a10)
       read(fil(1),'(10a1)') (ch(j),j=1,10)
c       print *,nhead,ch(1),fil(1)
       if(ch(1).eq.'*') then 
        goto 50
       elseif(ch(1).eq.'F') then 
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
          elseif(ch(j).eq.'E') then
           nim(i)=1
           write(filter(i),'(a1)') 'E'
          elseif(ch(j).eq.'F') then
           nim(i)=1
           write(filter(i),'(a1)') 'F'
          elseif(ch(j).eq.'O') then
           nim(i)=1
           write(filter(i),'(a1)') 'O'
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
      do i=1,ns
       read(21,*,end=100) (ID(i,j),j=1,nimage)
       read(21,*,end=100) (ID2(i,j),j=1,nimage)
       read(21,*,end=100) (mag(i,j),j=1,nimage)
       read(21,*,end=100) (merr(i,j),j=1,nimage)
       read(21,'(a19)',end=100) output
c       print *,ID(i,1),mag(i,1),merr(i,1)
      enddo
100   write(xlabel,'(5a1)') filter1,' ','-',' ',filter2
      write(ylabel,'(5a1)') filtery1,' ','-',' ',filtery2
      nstar=i-1
      print *,nstar,nimage,(filter(i),i=1,nimage)
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
     
      do i=1,nstar
       if(mag(i,xaxis1).gt.0.d0.and.mag(i,xaxis2).gt.0.d0.and.
     * mag(i,yaxis1).gt.0.d0.and.mag(i,yaxis2).gt.0.d0) then
        x(i)=mag(i,xaxis1)-mag(i,xaxis2)
        y(i)=mag(i,yaxis1)-mag(i,yaxis2)
         lx=min(lx,x(i))
         mx=max(mx,x(i))
         ly=min(ly,y(i))
         my=max(my,y(i))
       endif
      enddo
      print *,lx,mx,ly,my
      write(output,101) 'postlandfile ',filtery1,'-',filtery2,'.',
     *filter1,'-',filter2,'.ps'
101   format(a13,7a1,a3)
      write(output2,102) filter1,'_',filtery1,'-',filtery2,'.',
     *filter1,'-',filter2,'_.ccd'
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
      call sm_limits(lx-0.1,mx+0.1,my+0.1,ly-0.1)
c      call sm_limits(-0.5,4.,6.,0.)
      call sm_ctype('black')
      call sm_box(1,2,0,0)
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
      call sm_gflush()
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
      