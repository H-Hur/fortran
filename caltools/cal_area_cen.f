      subroutine cal_cent_weight_sqare(x1,y1,x2,y2,x3,y3,x4,y4,x,y)
      real x1,y1,x2,y2,x3,y3,x4,y4,xc1,yc1,xc2,yc2,a1,a2,d2,d3,d4
      real xxx(4),yyy(4)
      d2=sqrt((x2-x1)**2+(y2-y1)**2)
      d3=sqrt((x3-x1)**2+(y3-y1)**2)
      d4=sqrt((x4-x1)**2+(y4-y1)**2)
      xxx(1)=x1
      yyy(1)=y1
      if(d2.le.d3.and.d2.le.d4.and.d3.le.d4) then
       xxx(2)=x2
       yyy(2)=y2
       xxx(3)=x3
       yyy(3)=y3
       xxx(4)=x4
       yyy(4)=y4
      elseif(d2.le.d3.and.d2.le.d4.and.d4.le.d3) then
       xxx(2)=x2
       yyy(2)=y2
       xxx(3)=x4
       yyy(3)=y4
       xxx(4)=x3
       yyy(4)=y3
      elseif(d3.le.d2.and.d3.le.d4.and.d2.le.d4) then
       xxx(2)=x3
       yyy(2)=y3
       xxx(3)=x2
       yyy(3)=y2
       xxx(4)=x4
       yyy(4)=y4
      elseif(d3.le.d2.and.d3.le.d4.and.d4.le.d2) then
       xxx(2)=x3
       yyy(2)=y3
       xxx(3)=x4
       yyy(3)=y4
       xxx(4)=x2
       yyy(4)=y2
      elseif(d4.le.d2.and.d4.le.d2.and.d2.le.d3) then
       xxx(2)=x4
       yyy(2)=y4
       xxx(3)=x2
       yyy(3)=y2
       xxx(4)=x3
       yyy(4)=y3
      elseif(d4.le.d2.and.d4.le.d3.and.d3.le.d2) then
       xxx(2)=x4
       yyy(2)=y4
       xxx(3)=x3
       yyy(3)=y3
       xxx(4)=x2
       yyy(4)=y2
      endif
      call cal_cent_weight_triangle(xxx(1),yyy(1),xxx(2),yyy(2),
     &xxx(3),yyy(3),xc1,yc1)
      call cal_cent_weight_triangle(xxx(4),yyy(4),xxx(2),yyy(2),
     &xxx(3),yyy(3),xc2,yc2)
      call cal_area_triangle(xxx(1),yyy(1),xxx(2),yyy(2),
     &xxx(3),yyy(3),a1)
      call cal_area_triangle(xxx(4),yyy(4),xxx(2),yyy(2),
     &xxx(3),yyy(3),a2)
      x=(xc1*a1+xc2*a2)/(a1+a2)
      y=(yc1*a1+yc2*a2)/(a1+a2)
      return
      end

      subroutine cal_cent_weight_triangle(x1,y1,x2,y2,x3,y3,x,y)
      real x1,y1,x2,y2,x3,y3,x,y
      x=(x1+x2+x3)/3.
      y=(y1+y2+y3)/3.
      return
      end

      subroutine cal_area_triangle(x1,y1,x2,y2,x3,y3,a)
      real x1,y1,x2,y2,x3,y3,a
      a=abs((x1-x2)*y3 + (x2-x3)*y1 + (x3-x1)*y2 )/2.
      return
      end

