c      gfortran -o cal_lin ~/work/fortran/subroutine/caltools/cal_lin.f
       real x1,y1,x2,y2,s,y
       print *,'x1,y1,x2,y2=?'
       read(*,*) x1,y1,x2,y2
       call cal_lin(x1,y1,x2,y2,s,y)
       print *,s,y
       stop
       end
       subroutine cal_lin(x1,y1,x2,y2,s,y)
       real x1,y1,x2,y2,s,y
       s=(y2-y1)/(x2-x1)
       y=y1-s*x1
       return
       end
