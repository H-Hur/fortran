      real angle,rad

      call cal_angle(1.,1.,1.7,2.2,9000.,5000.,210.,297.,angle,rad)
      print *,angle,rad
      stop
      end

c     dx,dy : slope of the line
c     range : data range in the figure
c     pos : location position in the sheet
c     A4 postport : xleng=210,yleng=297  landscape : xleng=297,yleng=210
      subroutine cal_angle(dx,dy,xrange,yrange,xpos,ypos,xleng,yleng,
     &angle,rad)
      real dx,dy,xrange,yrange,xpos,ypos,xleng,yleng,angle,rad
      real*8 pi
      pi=3.14159265358979d0
      rad=atan((xrange/yrange)*(ypos/xpos)*(yleng/xleng)*(dy/dx))
      angle=rad*180./pi
      return
      end


