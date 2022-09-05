      subroutine cal_det(n,mat,det)
      parameter (nx=100) 
      integer n,n1,n2
      real*8 mat(nx,nx),cd1(n),cd2(n),det,d1,d2            
      d1=0.d0
      d2=0.d0
      do i=1,n
       cd1(i)=1.d0
       cd2(i)=1.d0
       do j=1,n
        n1=j
        n2=i+j
        if(n2.ge.n+1) then
         n2=n2-n
        endif
        cd1(i)=cd1(i)*mat(n1,n2)
        n2=i-j
        if(n2.le.0) then
         n2=n2+n
        endif
        cd2(i)=cd2(i)*mat(n1,n2)
       enddo
       d1=d1+cd1(i)
       d2=d2+cd2(i)
      enddo
      det=d1-d2 
      return
      end
