c     f95 -o cross /home/hhur/fortran/subroutine/daotools/cross.f
      parameter (ns=100000)
      character input*16,output*16
      integer nsat,ncrt
      real*8 xcen,ycen,xsat(ns),ysat(ns),cx1,cx2,cx3,cx4,cy1,cy2,cy3,cy4
      
      write(input,'(a16)') 'sat.coo                                    '
      write(output,'(a16)') 'cross.c                                   '

      xcen=2032.5d0
      ycen=2032.5d0
      open(21,file=input)
      do i=1,ns
       read(21,*,end=99) xsat(i),ysat(i)
      enddo
99    close(21)
      nsat=i-1
 
      open(22,file=output)
      do i=1,nsat
       if(xsat(i).lt.xcen.and.ysat(i).lt.ycen) then !3
        cx1=2.d0*xcen -xsat(i)
        cy1=2.d0*ycen -ysat(i)
        cx2=2.d0*xcen -xsat(i)
        cy2=ysat(i)
        cx3=xsat(i)
        cy3=ysat(i)
        cx4=xsat(i)
        cy4=2.d0*ycen -ysat(i)
        write(22,'(t1,2f9.3)') cx1,cy1
        write(22,'(t1,2f9.3)') cx2,cy2
        write(22,'(t1,2f9.3)') cx4,cy4
       elseif(xsat(i).ge.xcen.and.ysat(i).lt.ycen) then !4
        cx1=xsat(i)
        cy1=2.d0*ycen -ysat(i)
        cx2=2.d0*xcen -xsat(i)
        cy2=2.d0*ycen -ysat(i)
        cx3=2.d0*xcen -xsat(i)
        cy3=ysat(i)
        cx4=xsat(i)
        cy4=ysat(i)
        write(22,'(t1,2f9.3)') cx1,cy1
        write(22,'(t1,2f9.3)') cx2,cy2
        write(22,'(t1,2f9.3)') cx3,cy3
       elseif(xsat(i).lt.xcen.and.ysat(i).ge.ycen) then !2
        cx1=2.d0*xcen -xsat(i)
        cy1=ysat(i)
        cx2=xsat(i)
        cy2=ysat(i)
        cx3=xsat(i)
        cy3=2.d0*ycen -ysat(i)
        cx4=2.d0*xcen -xsat(i)
        cy4=2.d0*ycen -ysat(i)
        write(22,'(t1,2f9.3)') cx1,cy1
        write(22,'(t1,2f9.3)') cx3,cy3
        write(22,'(t1,2f9.3)') cx4,cy4
       elseif(xsat(i).ge.xcen.and.ysat(i).ge.ycen) then !1
        cx1=xsat(i)
        cy1=ysat(i)
        cx2=2.d0*xcen -xsat(i)
        cy2=ysat(i)
        cx3=2.d0*xcen -xsat(i)
        cy3=2.d0*ycen -ysat(i)
        cx4=xsat(i)
        cy4=2.d0*ycen -ysat(i)
        write(22,'(t1,2f9.3)') cx2,cy2
        write(22,'(t1,2f9.3)') cx3,cy3
        write(22,'(t1,2f9.3)') cx4,cy4
       endif
      enddo
      close(22)
      stop
      end
