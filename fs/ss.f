C   For fortran95.
C   Last update : June,26th, 2008  by Hyeon-oh Hur
c   f95 -o ss /home/gjgusdh/fortran/subroutine/fs/ss.f
      parameter (ns=100000)
      real xc(ns),yc(ns),mag(ns),merr(ns),chi(ns)
      real sharp(ns),sx(3),sy(3),cx,cy,dat(7),answer,apcrad,dfit
      integer ID(ns),ninput,nstar,sid(3),cids(5),np,cid
      character input*16,filename*16,apcname*16,h(16)*1
      open(26,file='input.ss')
      open(25,file='out.ss')
c50    print *,'Input allstar file name.'
50    read (26,'(a16)') filename
c      print *,'aperture correction radius=?'
      read(26,*) apcrad
c      print *,'fitrad to run cross-identification program=?'
      read(26,*) dfit
      open(21,file=filename)
      call reading(nstar,canid,dat,mag,xc,yc,ID,filename)
c      print *,'searching radius=',dat(7)
      do i=1,5
       cids(i)=0
      enddo
 
      do i=2,4
c60     print *,'input x,y coordinate on the image=?'
60     read(26,*) cx,cy
        if (cx.lt.0..or.cy.lt.0.) then
         goto 100
        endif   
        call fcan(cx,cy,nstar,dat,xc,yc,ID,mag,np,cid)       
        if(np.eq.0) then
         goto 60
         elseif(np.eq.1.and.cid.eq.cids(i-1)) then
c         print *,''
         print *,filename,'Choiced the same star. Select other star.'
c         print *,''
         goto 60
         elseif(np.eq.1) then
         cids(i)=cid
c         print *,'O.K'
        endif
       enddo
      read(filename,'(16a1)') (h(i),i=1,16)
      do i=1,13
       if(h(i).eq.'.'.and.h(i+1).eq.'a'.and.h(i+2).eq.'l'.and.
     &  h(i+3).eq.'s') then
        write(h(i+2),'(a1)') 'p'
        write(h(i+3),'(a1)') 'c'
        write(apcname,'(16a1)') (h(j),j=1,i+3)
        goto 80
       endif
      enddo
80    write(25,501) filename,apcname,cids(2),cids(3),cids(4),apcrad,dfit

c      print *,''
c100   print *,'Do you want to find same stars on other images?'
c      print *,'yes=1,no=0'
100   read(26,*) answer
      if (answer.eq.1) then
       close(21)
       goto 50
      endif
      close(21)
      close(25)
      close(26)
501   format(t1,a16,t17,a16,t33,i7,t40,i7,t47,i7,t54,f7.3,t62,f4.1)

      stop
      end




c This subroutine made to read and write the result of task allstar in IRAF 
      subroutine reading(n,canid0,dat,mag,xc,yc,ID,filename)
      parameter(ns=100000)
      real xc(ns),yc(ns),mag(2,ns),merr(ns),chi(ns),dat(7)
      character*2 H
      character filename*16,filter*16,input*16,INimage,Pimage
      integer ns0,ni0,n,id(ns),i,cid
  

      open(22,file='result',err=22)
      goto 30
22    print *,'Failed to open ',filename
      

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
23    format(A16)                                   !chatacter format
24    format(T17,A16)                               !character format 
25    format (T17,F7.1)                         ! data min,max format
26    format (T17,F5.2)                            !gain,noise format
27    format(I4,2F9.3,F7.2,F7.3,F10.5,I3,F8.3,F7.3,I3)!writing format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


30    print *,filename       
              do i=1,6!건너뛰기             1~6
              read(21,'(a2)',err=31) H
              enddo

      read(21,24,err=31) INimage!인풋이미지 7
      read(21,'(a2)',err=31) H!건너뛰기         8
      read(21,24,err=31) Pimage!PSF이미지   9
              do i=10,13!건너뛰기           10~13
              read(21,'(a2)',err=31) H
              enddo
       do i=14,15!datamin,datamax           14~15
       read(21,25,err=31) dat(i-13)
       enddo
       do i=16,17!gain,noise                16~17
       read(21,25,err=31) dat(i-13)
       enddo
              do i=18,19!건너뛰기           18~19
              read(21,'(a2)',err=31) H
              enddo

      read(21,24,err=31) filter!필터   20
              do i=21,23
              read(21,'(a2)',err=31) H!건너뛰기 21~23
              enddo

       do i=24,26!PSFMAG,PSFRAD,FWHM        24~26
       read(21,25,err=31) dat(i-19)
       enddo
              do i=27,44!건너뛰기           27~44
              read(21,'(a2)',err=31) H
              enddo
      goto 32

31    print *,'Failed to start reading from ',filename


32    do i=1,ns   !                     45~
       read(21,*,err=33,end=34) ID(i),xc(i),yc(i),mag(1,i),mag(2,i)
       read(21,*,err=33,end=34) H
      enddo
33    print *,'Failed reading from ',filename
34    print *,'Succesfully read from data. number of stars=',i-1
      n=i-1

      do i=1,n
       Write(22,27,err=35) ID(i),xc(i),yc(i),mag(1,i),mag(2,i)
      enddo
      goto 36

35    print *,'Failed to write the result.'
36    print *,' '

      close(21)
      close(22)
      return
      end



      subroutine fcan(cx,cy,n,dat,xc,yc,ID,mag,np,cid)
      parameter (ns=100000)
      integer ID(ns),np,n,cid
      real cx,cy,dat(7),xc(ns),yc(ns),mag(2,ns)
      np=0
      do i=1,n
       if(sqrt((cx-xc(i))**2.+(cy-yc(i))**2.).le.dat(7)) then
        if(cx-xc(i).lt.dat(7)*0.75d0.and.cy-yc(i).lt.dat(7)*0.75d0) then
         np=np+1
c        print *,'ID,Xcenter,Ycenter,Magnitude,Mag error'
c        print *,ID(i),xc(i),yc(i),mag(1,i),mag(2,i)
         cid=ID(i)
        endif
       endif
      enddo
      if(np.eq.0) then
       write(*,'(a19,2f8.3)')'no star found near ',cx,cy
      endif
      return
      end






       subroutine sort(Nm11,EID11)
       integer Nm11,CID
       real Crad,EID11(999999,2)
       CID=0
       Crad=0.

       do i1=1,Nm11
        do j1=1,Nm11-1
         if (EID11(j1+1,2).lt.EID11(j1,2)) then
          Crad=EID11(j1,2)
          EID11(j1,2)=EID11(j1+1,2)
          EID11(j1+1,2)=Crad
          CID=EID11(j1,1)
          EID11(j1,1)=EID11(j1+1,1)
          EID11(j1+1,1)=CID
         endif
        enddo
       enddo
      return
      end




