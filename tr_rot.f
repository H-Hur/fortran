      integer nstar
      real*8 xc,yc,x,y,dcx,dcy,dx,dy,r,l,the1,the2,del,pi,x1,y1
      pi=3.14159265358979
      del=-0.23*pi/180.d0
      xc=1355.353d0
      yc=1386.246d0
   
      open(31,file='1.coo')
      open(32,file='4.coo')
      do i=1,100000
       read(31,*,end=500) x,y
       dcx=xc-x
       dcy=yc-y
       r=sqrt(dcx**2.d0+dcy**2.d0)
       the1=acos(dcx/r)
        the2=the1+del
       x1=xc-r*cos(the2)
       y1=yc-r*sin(the2)
       if(y.gt.yc) then
        the2=the1-del
        x1=xc-r*cos(the2)
        y1=yc+r*sin(the2)

       endif

       write(32,'(2f10.3)') x1,y1
      enddo
500   close(31)    
      close(32)
 
      stop
      end  
