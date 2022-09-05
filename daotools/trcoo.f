! for fortran95
! final Update : 01,04,2008
c     gfortran -o tr ~/fortran/subroutine/daotools/trcoo.f
      parameter(mnstar=999999)
      real*8 coo(2,mnstar),xcoo(100),ycoo(100),x,y
      real*8 dx,dy,xm,ym,xt,yt,xmin,ymin,xmax,ymax
      character input(100)*15,output(100)*15,a*1,master*15,targets*15
      integer ninput,nstar,m
100   format(a15,t17,f7.2,t25,f7.2) 
101   format(a15)
      xmin=0.d0
      ynim=0.d0
      xmax=4096.d0
      ymax=4096.d0
      print *,''
      print *,'Coordinate transformming program for IRAF photometry.(204
     &8pixels).'
      print *,''
      print *,'V1.12, updated at Sep,19, 2007  by Hyeon-Oh Her'
      print *,''
      print *,'input file name : dcoo.tr'
      print *,'input format : filename,xcoo,ycoo(a15,t17,f7.2,t25,f7.2)'
      print *,'O.K? (Write the master file name, or press (0) to change 
     *single image.)'
      read(*,'(a15)') master
      open(31,file='dcoo.tr')
      do i=1,100
       read(31,100,end=110) input(i),xcoo(i),ycoo(i)
       write(output(i),101) input(i)
      enddo
110   close(31)
      ninput=i-1
      print *,ninput,' input files read.'

      read(master,'(a1)') a
      if(a.eq.'0') then
       goto 300
      elseif(a.eq.' ') then
       write(master,'(a15)') input(1)
      endif      


210   format (2f9.3)
211   format(t8,a15,t25,f7.2,t34,f7.2,t42,f6.2,t49,f6.2)
212   format(t1,a7,t8,a15,t25,f7.2,t34,f7.2)
213   format(t8,a15,t26,a4,t35,a4,t44,a2,t51,a2)


      write(*,213) 'Input file name','Xcoo','Ycoo','dx','dy'
      do i=1,ninput
       if(input(i).eq.master) then
        m=i
       endif 
      enddo
      do i=1,ninput
       dx=xcoo(m)-xcoo(i)
       dy=ycoo(m)-ycoo(i)
       if(input(i).eq.master) then
        write(*,212) 'Master-',input(i),xcoo(i),ycoo(i)
        goto 230
       endif 
       write(*,211) input(i),xcoo(i),ycoo(i),dx,dy
       open(200,file=master)
       open(201,file=input(i))
       do j=1,mnstar
        read(200,*,end=220) x,y
        x=x-dx
        y=y-dy
       if(x.gt.xmin.and.y.gt.ymin.and.x.lt.xmax.and.y.lt.ymax) then
        write(201,210) x,y
       endif
       enddo
220    nstar=j-1
       close(200)
       close(201)
230    nstar=nstar  
      enddo
      goto 800
 

300   xm=99999.d0
      xt=0.d0
      print *,'Master coordinate file name=?'
      read(*,'(a15)') master
      do i=1,ninput
       if(input(i).eq.master) then
        xm=xcoo(i)
        ym=ycoo(i)
        write(*,302) 'Master : ', input(i),xcoo(i),ycoo(i)
        goto 301
       endif
      enddo
      print *,'File not found.'

301   print *,''
      print *,'Target coordinate filen name=?'
      read(*,'(a15)') targets
      do i=1,ninput
       if(input(i).eq.targets) then
        xt=xcoo(i)
        yt=ycoo(i)
        dx=xm-xt
        dy=ym-yt
        write(*,302) 'Target : ',input(i),xcoo(i),ycoo(i)
        write(*,303) 'dx = ',dx,', dy = ',dy 
        goto 350
       endif
      enddo
      print *,'File not found.'

302   format(t8,a9,t17,a15,t34,f7.2,t43,f7.2)
303   format(t8,a5,f6.2,a7,f6.2)    

350   open(200,file=master)
      open(201,file=targets)
      do i=1,mnstar
       read(200,*,end=400) x,y
       x=x-dx
       y=y-dy
       if(x.gt.xmin.and.y.gt.ymin.and.x.lt.xmax.and.y.lt.ymax) then
        write(201,210) x,y
       endif
      enddo
400   nstar=i-1
      close(200)
      close(201)



800   print *,''
      if(nstar.eq.0.or.dx.eq.99999.d0) then
       print *,'File not found. Check input coordinate file name.'
      else
       print *,nstar,' stars transformed coordinate.'
      endif
      stop
      end
