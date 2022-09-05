c     f95 -o dstar /home/gjgusdh/fortran/subroutine/fs/dstar_id.f
      parameter (ns=100000)
      integer nstar,nhead,id(ns),did(ns),ndel
      character line1(ns)*80,line2(ns)*80,head(ns)*80,inp*16,outp*17
      character ch(17)*1,H*1
      print *,'input file name=?'
      open(23,file='input.dstar')
      do i=1,ns
       read(23,*,end=5) did(i)
      enddo
5     ndel=i-1
      print *,'input als file=?'
      read(*,'(a16)') inp
      open(21,file=inp)
      open(22,file='out.dstar')
      nhead=0
      nstar=0
      do i=1,ns
       read(21,'(t1,a1)',end=100) H
       if(H.eq.'#') then
        nhead=nhead+1
        backspace(21)
        read(21,'(t1,a80)') head(nhead)
       else
        nstar=nstar+1
        backspace(21)
        read(21,'(t1,a80)') line1(nstar)
        read(21,'(t1,a80)') line2(nstar)
        read(line1(nstar),'(t1,i6)') id(nstar)
       endif
      enddo
100   print *,nhead,'lines head,'
      print *,nstar,'stars read.'
      do i=1,nhead
       write(22,'(t1,a80)') head(i)
      enddo
      do i=1,nstar
       do j=1,ndel
        if(id(i).eq.did(j)) then
         goto 250
        endif
       enddo
       write(22,'(t1,a80)') line1(i)
       write(22,'(t1,a80)') line2(i)
250    ndel=ndel
      enddo
      print *,ndel,'stars deleted.'
      close(23)
      close(21)
      close(22)
      stop
      end
