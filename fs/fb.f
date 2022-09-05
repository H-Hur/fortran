c     for fortran95 only!
c This program pickes up only coordinates to output files.
c    f95 -o fb /home/gjgusdh/fortran/subroutine/fs/fb.f
      integer line,refID,targetID,nread,nt
      real*8 refcoo(2),tarcoo(2)                
      character H*1,H2*1,Hh*1,input*16      
      print *,'Input file name = ? ( 1 = binary.fs )'
      read(*,'(a16)')input
      print *,'Start line on the input file = ?'
      read(*,*) line          

      if(input.eq.'1') then
       print *,'Input file name = binary.fs'
       open(21,file='binary.fs')
       write(input,'(a16)') 'binary.fs              '
      else
       open(21,file=input)
      endif
 
      open(22,file='ref.coo')
      open(23,file='tar.coo')
      nread=0
      nt=0
      do i=1,line-1
       read(21,'(a1)') H
      enddo
10    format (I6,T8,F8.3,T17,F8.3)
11    format (T27,I6,T36,F8.3)      
      do i=1,1000
       read(21,'(t2,a1,t6,a1,t32,a1)',end=50) Hh,H,H2   
c       print *,i,Hh,H,H2
       if(Hh.eq.'T'.or.Hh.eq.'R') goto 50
       if(H2.eq.' '.and.H.ne.' ') then
        backspace(21)
        nread=nread+1
        read(21,*) refID,refcoo(1),refcoo(2)
        write(22,'(2f9.3)') refcoo(1),refcoo(2)
       elseif(H.eq.' '.and.H2.ne.' ') then
        backspace(21)
        nt=nt+1
        read(21,*) targetID,tarcoo(1),tarcoo(2)
        write(23,'(2f9.3)') tarcoo(1),tarcoo(2)
       endif
      enddo
50    print *,'Wrote to [ref.coo], [tar.coo].'
 
      close(21)
      close(22)
      close(23)
      stop
      end     
