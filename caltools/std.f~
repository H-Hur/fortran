      subroutine STD(n,dat,average,stddev)
      real dat(n),average,stddev
      average=0.
      stddev=0.
      if(n.eq.0) then
       print*,'Subroutine STD received no data, ndata=0'
       goto 100
      elseif(n.eq.1) then
       average=dat(1)
       stddev=0.
c       print*,'Subroutine STD received single data, ndata=1 '
      else
       do i=1,n
        average=average+dat(i)
       enddo
       average=average/real(n)
       do i=1,n
        stddev=stddev+(average-dat(i))*(average-dat(i))
       enddo
       stddev=sqrt(stddev/real(n-1))
      endif
100   return
      end

      subroutine wSTD_dp(n,dat,w,average,std)
      real*8 dat(n),average,std,w(n),wsum
      integer n,np
      average=0.d0
      std=0.d0
      wsum=0.d0
      np=0
      if(n.eq.0) then
       print*,'Subroutine wSTD_dp received no data'
       goto 100
      elseif(n.eq.1) then
       average=dat(1)
       std=0.d0
      else
       do i=1,n
        average=average+dat(i)*w(i)
        wsum=wsum+w(i)
        if(w(i).ne.0.d0) np=np+1
       enddo
       average=average/wsum
       do i=1,n
        std=std+w(i)*(average-dat(i))*(average-dat(i))
       enddo
       std=sqrt(std/(real(n-1)*wsum/real(np)))
      endif
100   return
      end

      subroutine wSTD(n,dat,w,average,std)
      real dat(n),average,std,w(n),wsum
      integer n,np
      average=0.
      std=0.
      wsum=0.
      np=0
      if(n.eq.0) then
       print*,'Subroutine wSTD received no data'
       goto 100
      elseif(n.eq.1) then
       average=dat(1)
       std=0.
      else
       do i=1,n
        average=average+dat(i)*w(i)
        wsum=wsum+w(i)
        if(w(i).ne.0.) np=np+1
       enddo
       if(np.ne.0) average=average/wsum
       if(np.eq.0) average=0.
       do i=1,n
        std=std+w(i)*(average-dat(i))*(average-dat(i))
       enddo
       if(np.ge.2) std=sqrt(std/(real(n-1)*wsum/real(np)))
       if(np.le.1) std=0.
      endif
100   return
      end
