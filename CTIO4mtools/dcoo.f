c     f95 -o dcoo dcoo.f
      real xc,yc,dx,dy,x1,x2,y1,y2
      character sname*12,line*21,input*40,H*1
      data x1,x2,y1,y2/25.,2023.,25.,4071./
      print *,'input file name=?'
      read(*,'(t1,a40)') input
      print *,'dx,dy=?'
      read(*,*) dx,dy
      open(21,file=input)
      read(21,'(t1,a1)') H
      do i=1,100000
       read(21,'(t1,a21,a12)',end=99) line,sname
       read(line,*) xc,yc
       xc=xc+dx
       yc=yc+dy
       if(xc.ge.x1.and.xc.le.x2.and.yc.ge.y1.and.yc.le.y2) 
     & write(*,'(t1,a12,2f9.3)')sname,xc,yc
      enddo
99    close(21)
      stop
      end
