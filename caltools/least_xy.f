      subroutine least_xy(xl,yl,n,slope,yp)
      real xl(n),yl(n),slope,yp
      integer n
      real xx,yy,xy,x,y,a,b,c,d,xc,yc,xyc,dslope,sq,rt
      xx=0.
      yy=0.
      xy=0.
      x=0.
      y=0.
      do i=1,n
       xx=xx+x(i)*x(i)
       yy=yy+y(i)*y(i)
       xy=xy+x(i)*y(i)
       x=x+x(i)
       y=y+y(i)
      enddo
      xc=real(n)*xx - x*x
      yc=real(n)*yy - y*y
      xyc = real(n)*xy - x*y  
      a = 2.*xc
      b = 2.*xc*yc - xc
      c = x - 2.*xc*yc
      d = 2.*xc*yc
  
      rt = 2.*b*b*b - 9.*a*b*c + 27.*a*a*d  
      sq = sqrt(rt*rt - 4.*(b*b-3.*a*c)**3) 
      dslope = (-1./(3.*a))*
     &(b + ((rt+sq)/2.)**(1./3.) + ((rt-sq)/2.)**(1./3.))

      slope = real(dslope)
      yp = - (dslope*x + y)/real(n)
      return
      end
