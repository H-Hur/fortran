c      f95 -o delcoo /home/gjgusdh/fortran/subroutine/daotools/del_coo.f
       parameter (ns=100000)
       integer ncoo,nstar,nclose
       real xcoo(ns),ycoo(ns),xstar(ns),ystar(ns),dist,fitrad
       character coo*16,star*16
       print *,'coo file =?'
       read(*,'(t1,a16)') coo
       print *,'star=?'
       read(*,'(t1,a16)') star
       print *,'Erasing rad=?'
       read(*,*) fitrad
       open(21,file=coo)
       do i=1,ns
        read(21,*,end=99) xcoo(i),ycoo(i)
       enddo
99     ncoo=i-1
       print *,ncoo
       close(21)
       open(22,file=star)
       do i=1,ns
        read(22,*,end=199) xstar(i),ystar(i)
       enddo
199    nstar=i-1
       print *,nstar
       close(22)
       
       open(31,file='coo.out')
       do i=1,ncoo
        nclose=0
        do j=1,nstar        
         dist=sqrt((xcoo(i)-xstar(j))**2+(ycoo(i)-ystar(j))**2)
         if(dist.le.fitrad) nclose=nclose+1
        enddo
        if(nclose.eq.0) write(31,'(t1,2f8.3)')  xcoo(i),ycoo(i)
       enddo
       close(31)
       stop
       end

