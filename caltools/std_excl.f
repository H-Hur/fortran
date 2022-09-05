
c     x,n : input data
c     th : threshold sigma
c     xx,nn : calculated data (within th)
c     av,st : output
c     this subroutine require the subroutin STD(n,dat,average,stddev,ncall)
      subroutine std_ex(x,n,th,nn,av,st)
      integer n,nn,n_ex,n_ex_last,niter
      real x(n),xx(100000),av,st,th
      do i=1,n
       xx(i)=x(i)
      enddo
      nn=n
      n_ex_last=0
      n_ex=0 
      niter=0
10    call std(nn,xx,av,st)
      n_ex_last=n_ex
      nn=0
      n_ex=0
      niter=niter+1
      do i=1,n
       if(abs(x(i)-av).le.st*th) then
        nn=nn+1
        xx(nn)=x(i)
       else
        n_ex=n_ex+1
       endif
      enddo
      if(n_ex.eq.n_ex_last) goto 100 
      if(nn.lt.10.and.niter.ge.2) goto 100
      goto 10
100   return
      end

      subroutine std_ex_w(x,w,n,th,nmin,nmax,nn,av,st)
      integer n,nn,n_ex,n_ex_last,nmin,nmax,ncal
      real x(100000),w(100000),xx(100000),ww(100000),av,st,th
      do i=1,n
       xx(i)=x(i)
       ww(i)=w(i)
      enddo
      nn=n
      n_ex_last=0
      n_ex=0
      ncal=0
10    call wSTD(nn,xx,ww,av,st) 
      n_ex_last=n_ex
      ncal=ncal+1
      nn=0
      n_ex=0
      do i=1,n
       if(abs(x(i)-av).le.st*th) then
        nn=nn+1
        xx(nn)=x(i)
        ww(nn)=w(i)
       else
        n_ex=n_ex+1
       endif
      enddo
      if(n_ex.gt.n_ex_last.and.nn.gt.nmin.and.ncal.lt.nmax) goto 10
100   return
      end

      subroutine STD(n,dat,average,stddev)
      real dat(n),average,stddev
      integer n 
      average=0.
      stddev=0.
      if(n.eq.0) then
       goto 100
      elseif(n.eq.1) then
       average=dat(1)
       stddev=0.
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

      subroutine wSTD(n,dat,w,average,std)
      real dat(n),average,std,w(n),wsum
      integer n,np
      np=0
      do i=1,n
       if(w(i).ne.0.) np=np+1
      enddo
      average=0.
      std=0.
      wsum=0.
      if(n.eq.0) then
       goto 100
      elseif(n.eq.1) then
       average=dat(1)
       std=0.
      else
       do i=1,n
        average=average+dat(i)*w(i)
        wsum=wsum+w(i)
       enddo
       average=average/wsum
       do i=1,n
        std=std+w(i)*(average-dat(i))*(average-dat(i))
       enddo
       std=sqrt(std/(real(np-1)*wsum/real(np)))
      endif
100   return
      end
      subroutine wSTD_dp(n,dat,w,average,std)
      real*8 dat(n),average,std,w(n),wsum
      integer n,np
      np=0
      do i=1,n
       if(w(i).ne.0.d0) np=np+1
      enddo
      average=0.d0
      std=0.d0
      wsum=0.d0
      if(n.eq.0) then
       goto 100
      elseif(n.eq.1) then
       average=dat(1)
       std=0.d0
      else
       do i=1,n
        average=average+dat(i)*w(i)
        wsum=wsum+w(i)
       enddo
       average=average/wsum
       do i=1,n
        std=std+w(i)*(average-dat(i))*(average-dat(i))
       enddo
       std=sqrt(std/(dble(np-1)*wsum/dble(np)))
      endif
100   return
      end
