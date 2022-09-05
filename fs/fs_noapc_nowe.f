c  Made for fortran95.
c  No aperture correction.
c  No weight combine term with difference exposed images.
      parameter (ns=100000,ni=30)
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),merr(ni,ns),chi(ni,ns)
      real*8 sharp(ni,ns),dx(ni,3),dy(ni,3),dat(ni,7),fitten(ni,ns)
      integer ID(ni,ns),nstar(ni),sid(ni,3),nimage,nstars
      character input(ni)*16,filter(ni)*16
      integer mpart(ns,ni,10),tpart(ns,ni,10),IDF(ns,ni)
      real*8 MAGF(ns,ni,8)

      open(21,file='inputfs.dat')
22    format(a16,t17,i7,t24,i7,t31,i7) 
      do n=1,ni
      read(21,22,end=50) input(n),sid(n,1),sid(n,2),sid(n,3)
      enddo
50    nimage=n-1
      close(21)

      do n=1,nimage
       call reading(n,xc,yc,mag,merr,dat,filter,input,ID,sharp,chi,nstar
     *)
      enddo 
      nstars=nstar(1)

      do i=1,nstar(1)
       IDF(i,1)=ID(1,i)
       MAGF(i,1,1)=xc(1,i)
       MAGF(i,1,2)=yc(1,i)
       MAGF(i,1,3)=0.d0
       MAGF(i,1,4)=0.d0
       MAGF(i,1,5)=mag(1,i)
       MAGF(i,1,6)=merr(1,i)
       MAGF(i,1,7)=sharp(1,i)
       MAGF(i,1,8)=chi(1,i)
      enddo
      do n=2,nimage
       call fs(n,input,filter,xc,yc,dx,dy,mag,merr,mpart,tpart,ID,sid,MA
     *GF,IDF,sharp,chi,nstars,nstar)
      enddo
       
      print *,'Writing binary candidates...' 
      call binary(nimage,input,ID,xc,yc,dx,dy,mpart,tpart,mag,merr,nstar
     *s,nstar)
      print *,'Writing result files [out1.fs], [out2.fs].'
500   open(32,file='out1.fs')
      open(33,file='out2.fs')
      write(32,'(30a10)') (input(i),i=1,nimage)
      write(32,*) ""
      write(33,'(30a10)') (input(i),i=1,nimage)
      write(33,*) ""
      write(32,'(a7,30a10)') '  ',(filter(i),i=1,nimage)
      write(32,*) ""
      write(33,'(a7,30a10)') '  ',(filter(i),i=1,nimage)
      write(33,*) ""

       do i=1,nstars
       write(32,'(30I10)') (IDF(i,j),j=1,nimage)
       write(32,'(30f10.3)') (MAGF(i,j,5),j=1,nimage)
       write(32,'(30f10.3)') (MAGF(i,j,6),j=1,nimage)
       write(32,*) ""
       write(33,'(30I10)') (IDF(i,j),j=1,nimage)
       do k=1,4
        write(33,'(30f10.3)') (MAGF(i,j,k),j=1,nimage)
       enddo
       k=5
       write(33,'(30f10.3)') (MAGF(i,j,k),j=1,nimage)       
       do k=6,8
        write(33,'(30f10.3)') (MAGF(i,j,k),j=1,nimage)
       enddo
       write(33,*) ""
      enddo
      close(32)
      close(33)

      stop
      end

      subroutine fs(ninput,input,filter,xc,yc,dx,dy,mag,merr,mpart,tpart
     *,ID,sid,MAGF,IDF,sharp,chi,nstars,nstar)
      parameter(ns=100000,ni=30)
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),dat(ni,7),fitten(ni,ns)
      integer ID(ni,ns),nstar(ni),ninput,sid(ni,3),i,j,k
      real*8 fitrad,dx(ni,3),dy(ni,3),tarx,tary,refx,refy,fit,tar(4)
      real*8 merr(ni,ns),MAGF(ns,ni,8),sharp(ni,ns),chi(ni,ns),dfit
      real*8 xmat(3,3),ymat(3,3),xtar(3,1),ytar(3,1),xref(3,1),yref(3,1)
      integer npart(ni,3),npartner(ns),used(ns),pu(ns),npu,IDF(ns,ni)
      integer tpart(ns,ni,10),nfit,tarID,mpart(ns,ni,10),nstars
      character input(ni)*16,filter(ni)*16,H*16,output*20,ch(20)*1

      do i=1,3
       xmat(i,1)=1.d0
       ymat(i,1)=1.d0
       npart(ninput,i)=0
      enddo

      do i=1,nstars
       do j=1,10
        mpart(i,ninput,j)=0
       enddo
       do j=1,3
        if(ID(1,i).eq.sid(1,j)) then
         xref(j,1)=xc(1,i)
         yref(j,1)=yc(1,i)
        endif
       enddo
      enddo

      do i=1,nstar(ninput)
       do j=1,10
       tpart(i,ninput,j)=0
       enddo
       do j=1,3
        if(ID(ninput,i).eq.sid(ninput,j)) then
         xtar(j,1)=xref(j,1)
         ytar(j,1)=yref(j,1)
         xmat(j,2)=xc(ninput,i)
         xmat(j,3)=yc(ninput,i)
         ymat(j,2)=yc(ninput,i)
         ymat(j,3)=xc(ninput,i)
        endif
       enddo
      enddo

c      call printing(xmat,xtar)
c      call printing(ymat,ytar)
      call cramer(ninput,xmat,xtar,dx)
      call cramer(ninput,ymat,ytar,dy)  
!1=dx,2=a,3=b,ref(x)=dx+a*tar(x)+b*tar(y)
      fitrad=2.5d0
      read(input(ninput),'(16a1)',end=22) (ch(i),i=1,16)
22    j=i-1

      do i=1,j
       if(ch(i).eq.'.') then
      goto 23
      endif
      enddo

23    write(ch(i+1),'(a1)') 't'
      write(ch(i+2),'(a1)') 'r'
      write(ch(i+3),'(a1)') 's'
      write(output,'(20a1)') (ch(j),j=1,i+3)

      open(31,file=output)
      write(31,*) 'input (target) als file : ',input(ninput)
      write(31,*) 'reference all star file : ',input(1)
      write(31,*) 'filter : ', filter(ninput)
      write(31,*) 'refx = ',dx(ninput,1),' + ',dx(ninput,2),' X tar(x)
     * + ',dx(ninput,3),' Xtar(y)'
      write(31,*) 'refy = ',dy(ninput,1),' + ',dy(ninput,2),' X tar(y)
     * + ',dy(ninput,3),' Xtar(x)'
      write(31,*) ''
90    format('fitrad=',f4.1,'...fitting...')
      print *,'Start finding counterparts : ',input(ninput)
      open(32,file='ID.fs')
      open(33,file='mag.fs') 

      do i=1,nstar(ninput)
       fitten(ninput,i)=0.d0
       npartner(i)=0
       pu(i)=0
      enddo

      do i=1,nstars
       used(i)=0
      enddo
      npu=0
      dfit=0.2d0
      fit=dfit
      do k=1,50,3
       nfit=0
       fit=fit+dfit
       if(fit.gt.fitrad) then
        fit=fitrad
       endif
       write(*,90) fit
       do i=1,nstar(ninput)
        refx=dx(ninput,1)+dx(ninput,2)*xc(ninput,i)+
     *  dx(ninput,3)*yc(ninput,i)
        refy=dy(ninput,1)+dy(ninput,2)*yc(ninput,i)+
     *  dy(ninput,3)*xc(ninput,i)
        do j=1,nstars
         if(sqrt((refx-xc(1,j))**2.d0+(refy-yc(1,j))**2.d0).le.fit.and.
     *    fitten(ninput,i).lt.0.1d0) then
          npartner(i)=npartner(i)+1
          used(j)=used(j)+1
          mpart(j,ninput,used(j))=ID(ninput,i)
          tpart(i,ninput,npartner(i))=ID(ninput,j) 
          pu(i)=used(j)
          fitten(ninput,i)=fit
          nfit=nfit+1
         endif  
        enddo
       enddo
       print *,'            ',nfit,'stars identified.'
       if(nfit.eq.nstar(ninput)) then
        goto 120
       elseif(fit.ge.fitrad) then
        goto 120
       endif 
      enddo

105   format(I6, T9,F8.3, T19,F8.3, T29,I3, T37,F4.1,T44,I1)

120   do i=1,nstars
       if(used(i).ge.2) then
        npu=npu+1
       endif
       IDF(i,ninput)=mpart(i,ninput,1)
       do j=1,nstar(ninput)
        if(ID(ninput,j).eq.IDF(i,ninput)) then
         goto 121
        endif
       enddo
121    MAGF(i,ninput,1)=xc(ninput,j)
       MAGF(i,ninput,2)=yc(ninput,j)
       MAGF(i,ninput,3)=xc(1,i)-xc(ninput,j)
       MAGF(i,ninput,4)=yc(1,i)-yc(ninput,j)
       MAGF(i,ninput,5)=mag(ninput,j)
       MAGF(i,ninput,6)=merr(ninput,j)
       MAGF(i,ninput,7)=sharp(ninput,j)
       MAGF(i,ninput,8)=chi(ninput,j)
      enddo

      do i=1,nstar(ninput)
       if(npartner(i).ge.2) then
        npart(ninput,2)=npart(ninput,2)+1
       elseif(npartner(i).eq.1) then
        npart(ninput,1)=npart(ninput,1)+1
       elseif(npartner(i).eq.0) then
        npart(ninput,3)=npart(ninput,3)+1
        nstars=nstars+1
        ID(1,nstars)=nstars
        mag(1,nstars)=9999.99d0
        merr(1,nstars)=9999.999d0
        xc(1,nstars)=dx(ninput,1)+dx(ninput,2)*xc(ninput,i)+dx(ninput,3)
     *  *yc(ninput,i)
        yc(1,nstars)=dy(ninput,1)+dy(ninput,2)*yc(ninput,i)+dy(ninput,3)
     *  *xc(ninput,i)
        sharp(1,nstars)=9999.999d0
        chi(1,nstars)=9999.999d0
        do j=1,ninput-1
         IDF(nstars,ninput-j)=0
         do k=1,8
          MAGF(nstars,ninput-j,k)=0.d0
         enddo
        enddo
        IDF(nstars,ninput)=ID(ninput,i) 
        do k=1,8
         MAGF(nstars,ninput,1)=xc(ninput,i)
         MAGF(nstars,ninput,2)=yc(ninput,i)
         MAGF(nstars,ninput,3)=9999.999d0
         MAGF(nstars,ninput,4)=9999.999d0
         MAGF(nstars,ninput,5)=mag(ninput,i)
         MAGF(nstars,ninput,6)=merr(ninput,i)
         MAGF(nstars,ninput,7)=sharp(ninput,i)
         MAGF(nstars,ninput,8)=chi(ninput,i)
        enddo
       endif
      enddo

      print *,nstars,'stars countedfrom all allstar files.'
      write(31,*) nstar(ninput),' stars read from allstar file.'
      write(31,*) npart(ninput,1),' stars found single partner.'
      write(31,*) npart(ninput,2),' stars found more than 2 partners.'
      write(31,*) npu,' stars have bad partner which has onother partner
     * on the same target image.' 
      write(31,*) npart(ninput,3),' stars found no partner.'
      write(31,*) ''
      write(31,*) '   ID     xc        yc   npartner fitrad refnpartner'

      do i=1,nstar(ninput)
       write(31,105)ID(ninput,i),xc(ninput,i),yc(ninput,i),
     * npartner(i),fitten(ninput,i),pu(i)
      enddo
      print *,'Wrote to [',output,'].' 
      close(31)       
      return
      end 
      
      
      subroutine binary(nimage,input,ID,xc,yc,dx,dy,mpart,tpart,mag,merr
     *,nstars,nstar)
      parameter (ns=100000,ni=30)
      integer nstar(ni),ID(ni,ns),mpart(ns,ni,10),i1,nstars
      integer tpart(ns,ni,10)
      character input(ni)*16
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),merr(ni,ns),tar(4),dx(ni,3)
      real*8 dy(ni,3),tx,ty
      open(31,file='binary.fs')
      write(31,*) 'This part shows master image stars who have many part
     *ners on target images.'
      write(31,*) 'Master allstar file : ',input(1)
      do k=2,nimage
       write(31,*) 'Target allstar file : ',input(k)
       write(31,*) ''
       write(31,*) 'RefID    Refxc    Refyc   TarID    Xcenter   Ycenter
     *   Mag     Merr    '
       do i=1,nstars
        if(mpart(i,k,2).eq.0) then
        else
        tx=(xc(1,i)-dx(k,1)-dx(k,3)*yc(1,i))/dx(k,2)
        ty=(yc(1,i)-dy(k,1)-dy(k,3)*xc(1,i))/dy(k,2)
         write(31,101) ID(1,i),tx,ty
         do j=1,10
          if(mpart(i,k,j).eq.0) then
           goto 50
          endif
           call idfs(k,mpart(i,k,j),ID,xc,yc,mag,merr,tar,nstar)
           write(31,100) mpart(i,k,j),(tar(l),l=1,4)
         enddo 
50       write(31,*) ''         
        endif
       enddo
      enddo

100   format(T27,I6,T36,F8.3,T46,F8.3,T57,F6.3,T65,F6.3)
101   format(I6,2f9.3)
      write(31,*) ''
      write(31,*) ''
      write(31,*) 'This part shows target image stars who have many part
     *ners on the master image.'
      write(31,*) 'Master allstar file : ',input(1)

      do k=2,nimage
       write(31,*) 'Target allstar file : ',input(k)
       write(31,*) ''
       write(31,*) 'RefID    Refxc    Refyc  TarID     Xcenter   Ycenter
     *  Mag   Merr   '
       do i=1,nstar(k)
        if(tpart(i,k,2).eq.0) then
        else
         write(31,101) ID(k,i),xc(k,i),yc(k,i)
         do j=1,10
          if(mpart(i,k,j).eq.0) then
           goto 150
          endif
           call idfs(k,tpart(i,k,j),ID,xc,yc,mag,merr,tar,nstar)
           write(31,100) tpart(i,k,j),(tar(l),l=1,4)
         enddo
150      write(31,*) '' 
        endif
       enddo
      enddo

      close(31)
      return
      end


      subroutine idfs(ninput,inID,ID,xc,yc,mag,merr,tar,nstar)
      parameter (ns=100000,ni=30)
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),merr(ni,ns),tar(4)
      integer ninput,nstar(ni),ID(ni,ns),inID,tarID
      do i=1,nstar(ninput)
       do j=1,nstar(1)
        if(ID(ninput,i).eq.inID) then
         tarID=ID(ninput,i)
         tar(1)=xc(ninput,i)
         tar(2)=yc(ninput,i)
         tar(3)=mag(ninput,i)
         tar(4)=merr(ninput,i)
        endif
       enddo
      enddo
      return
      end

           
      subroutine printing(m2,s1)
      real*8 m2(3,3),s1(3,1)
      print *,'|',m2(1,1),m2(1,2),m2(1,3),'|  ',s1(1,1)
      print *,'|',m2(2,1),m2(2,2),m2(2,3),'| =',s1(2,1)
      print *,'|',m2(3,1),m2(3,2),m2(3,3),'|  ',s1(3,1)
      print *,''
      return
      end


c This subroutine made to read and write the result of task allstar in IRAF 
      subroutine reading(ninput,xc,yc,mag,merr,dat,filter,input,ID,sharp
     *,chi,nstar)
      parameter(ns=100000,ni=30)
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),merr(ni,ns),dat(ni,7),a1
      real*8 sharp(ni,ns),chi(ni,ns)
      character filename*16,filter(ni)*16,input(ni)*16,INimage,Pimage,H
      integer ID(ni,ns),nstar(ni),ninput,i1
      write(filename,'(a16)') input(ninput)
      open(21,file=filename,ERR=22)
      print *,' '
      open(22,file='result')
      goto 30
22    print *,'Failed to open ',filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
23    format(A16)                                   !chatacter format
24    format(T17,A16)                               !character format 
25    format (T17,F7.1)                         ! data min,max format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
30    print *,'Running subroutine reading...'
              do i=1,6!건너뛰기             1~6
              read(21,23,err=31) H
              enddo
      read(21,24,err=31) INimage!인풋이미지 7
      read(21,23,err=31) H!건너뛰기         8
      read(21,24,err=31) Pimage!PSF이미지   9
              do i=10,13!건너뛰기           10~13
              read(21,23,err=31) H
              enddo
       do i=14,15!datamin,datamax           14~15
       read(21,25,err=31) dat(ninput,i-13)
       enddo
       do i=16,17!gain,noise                16~17
       read(21,25,err=31) dat(ninput,i-13)
       enddo
              do i=18,19!건너뛰기           18~19
              read(21,23,err=31) H
              enddo
      read(21,24,err=31) filter(ninput)!필터   20
              do i=21,23
              read(21,23,err=31) H!건너뛰기 21~23
              enddo
       do i=24,26!PSFMAG,PSFRAD,FWHM        24~26
       read(21,25,err=31) dat(ninput,i-19)
       enddo
              do i=27,44!건너뛰기           27~44
              read(21,23,err=31) H
              enddo
      goto 32
31    print *,'Failed to start reading from ',filename

32    do i=1,ns   !                     45~
      read(21,*,err=33,end=34) ID(ninput,i),xc(ninput,i),yc(ninput,i),ma
     *g(ninput,i),merr(ninput,i)
      read(21,*,err=33,end=34) sharp(ninput,i),chi(ninput,i)
      enddo

33    print *,'Failed reading from ',filename
34    print *,'Succesfully read from data. number of stars=',i-1,':',
     *filename
      nstar(ninput)=i-1
      close(21)
      close(22)
      return
      end


      subroutine sort(Nm11,EID11)
      integer Nm11,CID
      real*8 Crad,EID11(999999,2)
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

       subroutine cramer(n,matrix,res1,eq)
       parameter(ni=30)
       real*8 matrix(3,3),m(3,3),res1(3,1),res2(3,1),eq(ni,3),slot(3,1)
       real*8 mastereq
       integer i,l,n
       do i=1,3
        res2(i,1)=res1(i,1)
        do j=1,3
        m(i,j)=matrix(i,j)
        enddo
       enddo
       mastereq=m(1,1)*m(2,2)*m(3,3)+m(1,2)*m(2,3)*m(3,1)+m(1,3)*m(2,1)*
     * m(3,2)-m(1,2)*m(2,1)*m(3,3)-m(1,1)*m(2,3)*m(3,2)-m(1,3)*m(2,2)*
     * m(3,1)

       do l=1,3
        do i=1,3
         slot(i,1)=m(i,l)
         m(i,l)=res2(i,1)
         res2(i,1)=slot(i,1)
        enddo
       eq(n,l)=m(1,1)*m(2,2)*m(3,3)+m(1,2)*m(2,3)*m(3,1)+m(1,3)*m(2,1)*
     * m(3,2)-m(1,2)*m(2,1)*m(3,3)-m(1,1)*m(2,3)*m(3,2)-m(1,3)*m(2,2)*
     * m(3,1)
       eq(n,l)=eq(n,l)/mastereq
        do i=1,3
         res2(i,1)=res1(i,1)
         do j=1,3
         m(i,j)=matrix(i,j)
         enddo
        enddo
       enddo
       return
       end

