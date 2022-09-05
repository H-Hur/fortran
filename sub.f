      parameter (ns=100000)
      integer nstar1,nstar2,IH,npart(ns)
      real*8 xc1(ns),yc1(ns),xc2(ns),yc2(ns),mag1(ns),mag2(ns)
      real*8 dist,dmag,RH
      character main*16,targeta*16,H*1
      do i=1,ns
       npart(i)=0
      enddo

      print *,'Main allstar file=?'
      read(*,'(a16)') main
     
      print *,'Target allstar file=?'
      read(*,'(a16)') targeta
 
      
      open(31,file=main)
      open(32,file=targeta)
      open(33,file='sub.coo')
      open(34,file='part.coo') 

      nstar1=0
      nstar2=0
      do i=1,ns
       read(31,'(a1)',end=500) H
       if(H.ne.'#') then
        nstar1=nstar1+1
        backspace(31)
        read(31,*) IH,xc1(nstar1),yc1(nstar1),mag1(nstar1)
        read(31,*) RH
       endif
       read(32,'(a1)',end=500) H
       if(H.ne.'#') then
        nstar2=nstar2+1
        backspace(32)
        read(32,*) IH,xc2(nstar2),yc2(nstar2),mag2(nstar2)
        read(32,*) RH
       endif
      enddo 
  
500   close(31)
      close(32)

      do i=1,nstar1
       do j=1,nstar2
        dist=sqrt((xc1(i)-xc2(j))**2.d0+(yc1(i)-yc2(j))*2.d0)
        dmag=mag1(i)-mag2(j)      
        if(dist.le.1.d0.and.dmag.le.1.d0) then
         npart(j)=npart(j)+1
        endif
       enddo
      enddo
 
      do i=1,nstar2
       if(npart(i).eq.0) then
        write(33,'(10f11.3)') xc2(i),yc2(i)
       else
        write(34,'(10f11.3)') xc2(i),yc2(i)
       endif
      enddo
      close(33)
      close(34)

      stop
      end        
        
