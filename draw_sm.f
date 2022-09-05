c     f77 -o name name.f ~/work/fortran/subroutine/draw_sm.f -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 
      subroutine draw_box(xc,yc,siz)
      real xc,yc,siz
      call sm_relocate(xc-siz/2.,yc-siz/2.)
      call sm_draw(xc-siz/2.,yc+siz/2.)
      call sm_draw(xc+siz/2.,yc+siz/2.)
      call sm_draw(xc+siz/2.,yc-siz/2.)
      call sm_draw(xc-siz/2.,yc-siz/2.)
      return
      end

      subroutine draw_circle(xc,yc,rad)
      parameter (n=720)
      real xc,yc,rad,x(n),y(n),pi
      pi=3.14159265358979
      do i=1,n
       x(i)=xc+rad*cos((real(i)/real(n))*(2.*pi))
       y(i)=yc+rad*sin((real(i)/real(n))*(2.*pi))
      enddo
      call sm_conn(x,y,n)
      return
      end

      subroutine draw_plot(x,y,s,n,c1,c2)
      integer n
      real x(n),y(n),s(n),siz
      do i=1,n
       siz=(max(c1,c2)-s(i))/abs(c1-c2)
       if(siz.gt.0.9) siz=0.9
       if(siz.lt.0.1) siz=0.1
       call sm_ptype(400.+siz,1)
       call sm_points(x(i),y(i),1)
      enddo 
      return
      end

