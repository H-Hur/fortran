c     f95 -o list_stsstd list_stsstd.f 
      parameter (ns=100000)
      integer night(5),nobs(5), nstar, npos
      real*8 mag(5),err(5),vary, v(ns), ra(ns),dec(ns), sra,sdec
      real*8 bve,vie,rie
      character region*16, output*16, c(16)*1,H*1,sname(ns)*15
      character line(ns)*135, input*16, name2(ns)*15
      
      print *,'Input region name=?'
      write(region,'(a16)')'                             '
      write(output,'(a16)')'                             '
      write(input,'(a16)')'                             '
      read(*,'(a12)') region
      read(region,'(12a1)') (c(k),k=1,12)   
      do k=1,16
       if(c(k).eq.' ') then
        goto 10
       endif
      enddo
10    write(c(k),'(a1)') '.'
      write(c(k+1),'(a1)') 'd'
      write(c(k+2),'(a1)') 'a'
      write(c(k+3),'(a1)') 't'
      write(output,'(16a1)') (c(j),j=1,k+3)
      write(c(k+1),'(a1)') 'p'
      write(c(k+2),'(a1)') 'h'
      write(c(k+3),'(a1)') 'o'
      write(input,'(16a1)') (c(j),j=1,k+3)


      open(21,file=input)
      read(21,'(t2,a1)') H
      do i=1,ns
       read(21,'(t1,a15,a135)',end=99) sname(i),line(i) 
c       print *,line(i)
       backspace(21)
       read(21,'(t66,f6.3)') v(i)
      enddo
99    nstar=i-1
      close(21)  

      call datasort(v,line,sname,nstar)

      write(c(k+1),'(a1)') 'p'
      write(c(k+2),'(a1)') 'o'
      write(c(k+3),'(a1)') 's'
      write(input,'(16a1)') (c(j),j=1,k+3)
      open(21,file=input)
      read(21,'(t2,a1)') H
      do i=1,ns
       read(21,'(t1,f15.11,t18,f15.11,t100,a15)',end=199) 
     & ra(i),dec(i),name2(i)
      enddo
199   close(21)  
      npos=i-1

201   format(t1,8f8.3,8i4,f7.3,2f16.11,a15)

      open(22,file=output)
      write(22,'(t1,a112)') '# V        B-V     V-I     R-I    Nnight---
     &-----  nobs----------  Vary? R.A.             D.E.C.         Name 
     &                                                               '
      do i=1,nstar
       read(line(i),*) (mag(j),err(j),night(j),nobs(j),j=1,5),vary
       call cal_error(err(3),err(2),bve)
       call cal_error(err(3),err(5),vie)
       call cal_error(err(4),err(5),rie)
       if(nobs(3).eq.0) then
        goto 220
       endif
       sra=0.d0
       sdec=0.d0
       do j=1,npos
        if(name2(j).eq.sname(i)) then
         sra=ra(j)
         sdec=dec(j)
         goto 210
        endif
       enddo
       goto 220
210    write(22,201) mag(3),mag(2)-mag(3),mag(3)-mag(5), mag(4)-mag(5),
     & err(3),bve,vie,rie,
     & night(3),min(night(2),night(3)),min(night(3),night(5)),
     & min(night(4),night(5)), nobs(3),min(nobs(2),nobs(3)),
     & min(nobs(3),nobs(5)),min(nobs(4),nobs(5)),vary,sra,sdec,sname(i)
220    l=l
      enddo
      close(22)       
      print *,nstar,npos
      stop
      end
 

      subroutine datasort(val,ch,ch2,nstar)
      parameter (ns=100000)
      integer nstar
      real val(ns),ref
      character ch(ns)*135,chref*135,ch2(ns)*15,chref2*15
      do i=1,nstar
       do j=1,nstar-1
        if(val(j+1).lt.val(j)) then
         ref=val(j)
         write(chref,'(a135)') ch(j)
         write(chref2,'(a15)') ch2(j)

         val(j)=val(j+1)
         write(ch(j),'(a135)') ch(j+1)
         write(ch2(j),'(a15)') ch2(j+1)

         val(j+1)=ref
         write(ch(j+1),'(a135)') chref
         write(ch2(j+1),'(a15)') chref2

        endif
       enddo
      enddo
      return
      end

!     combine two errors.     
      subroutine cal_error(e1,e2,e3)
      real*8 e1,e2,e3
      if(e1.lt.0.d0.or.e2.lt.0.d0) then
       e3=-1.d0
      else
       e3=sqrt(e1**2.d0+e2**2.d0)
      endif
      return
      end

