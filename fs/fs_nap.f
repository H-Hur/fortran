c  Made for fortran95.
c  Aperture correction part realized on 09. Jan. 2008 by Hyeon-oh Hur
c  Weight combine realized on 01. Feb. 2008 by Hyeon-oh Hur
c  Last modified on 09. Feb. 2008 bu Hyeon-Oh Hur.

      parameter (ns=100000,ni=30)
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),merr(ni,ns),chi(ni,ns)
      real*8 sharp(ni,ns),dx(ni,3),dy(ni,3),fitten(ni,ns),ut(ni)
      integer ID(ni,ns),nstar(ni),sid(ni,3),nimage,nstars
      character input(ni)*16,filter(ni)*16,infs*16,inapc(ni)*16,h(16)*1
      integer mpart(ns,ni,10),tpart(ns,ni,10),IDF(ns,ni)
      real*8 MAGF(ns,ni,8),cmag(ni,ns),dfit(ni),apcrad(ni),dat(ni,8)

C Reading input.
      print *,'Input file=?, (1:default=input.fs)'
      read(*,'(a16)') infs
      do i=1,16
       read(infs,'(a1)') h(i)
       if(h(i).eq.'1') then
        open(21,file='input.fs')
        print *,'Input file name = input.fs'
        goto 30
       endif
      enddo
      open(21,file=infs)
      print *,'Input file name = ',infs
22    format(t1,a16,t17,a16,t33,i7,t40,i7,t47,i7,t54,f7.3,t62,f4.1)
30    do n=1,ni
       read(21,22,end=50) input(n),inapc(n),sid(n,1),sid(n,2),sid(n,3)
     & ,apcrad(n),dfit(n)
      enddo
50    nimage=n-1
      close(21)

C reading als and apc files.
      do n=1,nimage
       call reading(n,xc,yc,mag,merr,filter,input,ID,sharp,chi,ut,dat,
     * nstar) 
       call readapc(inapc(n),n,apcrad)
       call apcxid(n,nstar,ID,xc,yc,mag,merr)
       call apcorrecting(n,nstar,xc,yc,mag,merr,cmag)
       call readfilter(n,filter)
      enddo 

      nstars=nstar(1)
 
      do i=1,nstar(1)
       IDF(i,1)=ID(1,i)
       MAGF(i,1,1)=xc(1,i)
       MAGF(i,1,2)=yc(1,i)
       MAGF(i,1,3)=0.d0
       MAGF(i,1,4)=0.d0
       MAGF(i,1,5)=mag(1,i) !!!!!!!!!!!chamge to mag->cmag
       MAGF(i,1,6)=merr(1,i)
       MAGF(i,1,7)=sharp(1,i)
       MAGF(i,1,8)=chi(1,i)
      enddo
C Finding same stars on each images.
      do n=2,nimage
       call fs(n,input,filter,xc,yc,dx,dy,mag,merr,mpart,tpart,ID,
     *sid,MAGF,IDF,sharp,chi,nstars,nstar,dfit)
      enddo
      call wcombpar(nimage,IDF,MAGF,dat,ut,filter,nstar)
      call writewcomb(input,filter,ut,dat,IDF,nstar)

       
      print *,'Writing binary candidates...' 

C Finding binary candidates.
      call binary(nimage,input,ID,xc,yc,dx,dy,mpart,tpart,mag,merr
     *,nstars,nstar)
      print *,'Wrote result files [out1.fs], [out2.fs].'
      print *,'Writing fake star candidates...'
C Finding fake candidates.
      call writefake(nimage,filter,ID,xc,yc,dx,dy,tpart,mag,merr,
     *nstars,nstar,input)
      print *,'Wrote fake star candidates to [each input filename.fake]'
     
C Writing result.
500   open(32,file='out1.fs')
      open(33,file='out2.fs')
      write(32,'(t2,30a10)') (input(i),i=1,nimage)
      write(33,'(t2,30a10)') (input(i),i=1,nimage)
c     1st file header##################
      write(32,*) " F : Filrer"
      write(32,*) " T : UT"
      write(32,*) " X : Airmass"
c     1st file header##################
c     2nd file header##################
      write(33,*) ""
      write(33,*) " F : Filrer"
      write(33,*) ' P : PSFMAG'
      write(33,*) " T : UT"
      write(33,*) " X : Airmass"
      write(33,*) " *"
      write(33,*) ' ID'
      write(33,*) ' X center'
      write(33,*) ' Y cente'
      write(33,*) ' dX to the Master image'
      write(33,*) ' dY to the Master image'
      write(33,*) ' Mag(aperture corrected)'
      write(33,*) ' Merr(combined PSF photometry & aperture correctoin)'
      write(33,*) ' Sharpness'
      write(33,*) ' Chi'
      write(33,*) ""
c     2nd file header##################
      write(32,'(t1,a1,a6,30a10)') 'F','  ',(filter(i),i=1,nimage)
      write(32,'(t1,a1,30f10.6)') 'P',(dat(i,6),i=1,nimage)
      write(32,'(t1,a1,30f10.6)') 'T',(ut(i),i=1,nimage)
      write(32,'(t1,a1,30f10.7)') 'X',(dat(i,5),i=1,nimage)
      write(32,'(t1,a1)') "*"
      write(33,'(t1,a1,a6,30a10)') 'F','  ',(filter(i),i=1,nimage)
      write(33,'(t1,a1,30f10.6)') 'P',(dat(i,6),i=1,nimage)
      write(33,'(t1,a1,30f10.6)') 'T',(ut(i),i=1,nimage)
      write(33,'(t1,a1,30f10.7)') 'X',(dat(i,5),i=1,nimage)
      write(33,'(t1,a1)') "*"   

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

c To find same stars on other images.
      subroutine fs(ninput,input,filter,xc,yc,dx,dy,mag,merr,mpart,tpart
     *,ID,sid,MAGF,IDF,sharp,chi,nstars,nstar,dfit)
      parameter(ns=100000,ni=30)
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),fitten(ni,ns)
      integer ID(ni,ns),nstar(ni),ninput,sid(ni,3),i,j,k
      real*8 fitrad,dx(ni,3),dy(ni,3),tarx,tary,refx,refy,fit,tar(4)
      real*8 merr(ni,ns),MAGF(ns,ni,8),sharp(ni,ns),chi(ni,ns),dfit(ni)
      real*8 xmat(3,3),ymat(3,3),xtar(3,1),ytar(3,1),xref(3,1),yref(3,1)
      integer npart(ni,3),npartner(ns),used(ns),pu(ns),IDF(ns,ni)
      integer tpart(ns,ni,10),nfit,tarID,mpart(ns,ni,10),nstars
      character input(ni)*16,filter(ni)*16,H*16,output*20,ch(20)*1
      real*8 xcen,ycen,disdance,tx(ni,ns),ty(ni,ns)
      common/rcoo/rx(ni,ns),ry(ni,ns) 
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
c      print *,'Start finding counterparts : ',input(ninput)

      do i=1,nstar(ninput)
       fitten(ninput,i)=0.d0
       npartner(i)=0
       pu(i)=0
      enddo
      do i=1,nstars
       used(i)=0
      enddo
      fit=dfit(ninput)
      distance=1000.d0
c      write(*,90) fit
      nfit=0
      do i=1,nstar(ninput)
       refx=dx(ninput,1)+dx(ninput,2)*xc(ninput,i)+
     * dx(ninput,3)*yc(ninput,i)
       refy=dy(ninput,1)+dy(ninput,2)*yc(ninput,i)+
     * dy(ninput,3)*xc(ninput,i)
       tx(ninput,i)=refx
       ty(ninput,i)=refy
       do j=1,nstars
        if(sqrt((refx-xc(1,j))**2.d0+(refy-yc(1,j))**2.d0).le.fit.and.
     *   fitten(ninput,i).lt.0.1d0) then
         if(sqrt((refx-xc(1,j))**2.d0+(refy-yc(1,j))**2.d0).le.
     *    distance) then
          xcen=xc(1,j)
          ycen=yc(1,j)        
          distance=(sqrt((refx-xc(1,j))**2.d0+(refy-yc(1,j))**2.d0))
     *    **0.5d0
         endif
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
       print *,'            '
       print *,'          ',nfit,'of',nstar(ninput),'stars identified on
     * ',input(ninput)
       if(nfit.eq.nstar(ninput)) then
        goto 120
        elseif(fit.ge.fitrad) then
        goto 120
       endif 

105   format(I6, T9,F8.3, T19,F8.3, T29,I3, T37,F4.1,T44,I1)
120   do i=1,nstars
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
        rx(nstars,ninput)=tx(ninput,i)
        ry(nstars,ninput)=ty(ninput,i)
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

c      print *,nstars,'stars counted from all of allstar files.'
      write(31,*) nstar(ninput),' stars read from allstar file.'
      write(31,*) npart(ninput,1),' stars found single partner.'
      write(31,*) npart(ninput,2),' stars found more than 2 partners.'
      write(31,*) npart(ninput,3),' stars found no partner.'
      write(31,*) ''
      write(31,*) '   ID     xc        yc   npartner fitrad refnpartner'

      do i=1,nstar(ninput)
       write(31,105)ID(ninput,i),xc(ninput,i),yc(ninput,i),
     * npartner(i),fitten(ninput,i),pu(i)
      enddo
c      print *,'Wrote to [',output,'].' 
      close(31)       
      return
      end 

c To write fake candidates which have no partner on other images.
      subroutine writefake(nimage,filter,ID,xc,yc,dx,dy,tpart,mag,merr,n
     *stars,nstar,input)
      parameter (ns=100000,ni=30)
      integer nstar(ni),ID(ni,ns),nstars
      integer tpart(ns,ni,10),nimage,nout
      character input(ni)*16,filter(ni)*16,output*17,ch(17)*1
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),merr(ni,ns),tar(4),dx(ni,3)
      real*8 dy(ni,3),tx,ty
      open(31,file='fakestar.fs')
      write(31,*) 'This part shows target image stars who have no partne
     *rs on the master images.'
      write(31,*) 'Master filter : ',filter(1)
      write(ch(14),'(a1)') 'f'
      write(ch(15),'(a1)') 'a'
      write(ch(16),'(a1)') 'k'
      write(ch(17),'(a1)') 'e'
      do k=2,nimage
       write(31,*) ''
       write(31,*) 'Target filter : ',filter(k)
       read(input(k),'(13a1)') (ch(j),j=1,13)
       do i=1,13
        if(ch(i).eq.'.') then
         goto 50
        endif
       enddo
50     write(ch(i+1),'(a1)') 'f'
       write(ch(i+2),'(a1)') 'a'
       write(ch(i+3),'(a1)') 'k'
       write(ch(i+4),'(a1)') 'e'
       write(output,'(17a1)') (ch(j),j=1,i+4)
       open(32,file=output)
       write(31,*) 'Target allstar file : ',input(k)
       write(31,*) ''
       write(31,*) 'TerID    Xcenter   Ycenter  Refxc   Refyc    Mag  
     *Merr    '
       do i=1,nstar(k)
        if(tpart(i,k,1).eq.0) then
         tx=dx(k,1)+dx(k,2)*xc(k,i)+dx(k,3)*yc(k,i)
         ty=dy(k,1)+dy(k,2)*yc(k,i)+dy(k,3)*xc(k,i) 
         write(31,101) ID(k,i),xc(k,i),yc(k,i),tx,ty,mag(k,i),merr(k,i)
         write(32,102) tx,ty
        endif 
       enddo
       close(32)
      enddo
101   format(I6,T9,F8.3,T18,F8.3,T27,F8.3,T36,F8.3,T46,F6.3,T53,F6.3)
102   format(2f9.3)
      close(31)
      return
      end
 
c To find and write the binary candidates.       
      subroutine binary(nimage,input,ID,xc,yc,dx,dy,mpart,tpart,mag,merr
     *,nstars,nstar)
      parameter (ns=100000,ni=30)
      integer nstar(ni),ID(ni,ns),mpart(ns,ni,10),i1,nstars
      integer tpart(ns,ni,10),nimage
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
          if(tpart(i,k,j).eq.0) then
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

c To find stars using ID.
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

c To print matrix.           
      subroutine printing(m2,s1)
      real*8 m2(3,3),s1(3,1)
      print *,'|',m2(1,1),m2(1,2),m2(1,3),'|  ',s1(1,1)
      print *,'|',m2(2,1),m2(2,2),m2(2,3),'| =',s1(2,1)
      print *,'|',m2(3,1),m2(3,2),m2(3,3),'|  ',s1(3,1)
      print *,''
      return
      end


c This subroutine made to read and write the result of task allstar in IRAF 
      subroutine reading(ninput,xc,yc,mag,merr,filter,input,ID,sharp
     *,chi,ut,dat,nstar)
      parameter(ns=100000,ni=30)
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),merr(ni,ns),a1
      real*8 sharp(ni,ns),chi(ni,ns),sec
      character filename*16,filter(ni)*16,input(ni)*16,H*1
      character rhead*10,INimage*16,Pimage*16,cut*16
      integer ID(ni,ns),nstar(ni),ninput,i1,hour,mini
      real*8 dat(ni,8),ut(ni)
      write(filename,'(a16)') input(ninput)
      open(21,file=filename,ERR=22)
      print *,' '
      open(22,file='result')
      goto 30
22    print *,'Failed to open ',filename
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
24    format(T17,A16)                               !character format 
25    format (T17,F7.1)                         ! data min,max format
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
c30    print *,'Running subroutine reading...'
30    do i=1,100
       read(21,'(t1,a1)',err=31) H
       if(H.eq.'#') then
        backspace(21)
        read(21,'(t4,a10)',err=31) rhead
        if(rhead.eq.'IMAGE     ') then
         backspace(21)
         read(21,24,err=31) INimage
         elseif(rhead.eq.'DATAMIN   ') then
         backspace(21)
         read(21,25,err=31) dat(ninput,1)

         elseif(rhead.eq.'DATAMAX   ') then
         backspace(21)
         read(21,25,err=31) dat(ninput,2)

         elseif(rhead.eq.'OTIME     ') then
         backspace(21)
         read(21,24,err=31) cut
         read(cut,'(I2,a1,I2,a1,f6.3)',err=31) hour,H,mini,H,sec
         ut(ninput)=real(hour)+real(mini)/60.d0+real(sec)/3600.d0  
         elseif(rhead.eq.'GAIN      ') then
         backspace(21)
         read(21,25,err=31) dat(ninput,3)

         elseif(rhead.eq.'READNOISE ') then
         backspace(21)
         read(21,25,err=31) dat(ninput,4)

         elseif(rhead.eq.'XAIRMASS  ') then
         backspace(21)
         read(21,25,err=31) dat(ninput,5)

         elseif(rhead.eq.'IFILTER   ') then
         backspace(21)
         read(21,24,err=31) filter(ninput)

         elseif(rhead.eq.'PSFMAG    ') then
         backspace(21)
         read(21,25,err=31) dat(ninput,6)

         elseif(rhead.eq.'PSFRAD    ') then
         backspace(21)
         read(21,25,err=31) dat(ninput,7)

         elseif(rhead.eq.'FITRAD    ') then
         backspace(21)
         read(21,25,err=31) dat(ninput,8)
        endif

        else 
        backspace(21)
        goto 32 
       endif
      enddo


31    print *,'Failed to start reading from ',filename
32    do i=1,ns 
      read(21,*,end=34) ID(ninput,i),xc(ninput,i),yc(ninput,i),ma
     *g(ninput,i),merr(ninput,i)
      read(21,*,end=34) sharp(ninput,i),chi(ninput,i)
      enddo
      
c34    print *,'Succesfully read from data. number of stars=',i-1,':',
c     *filename
34    nstar(ninput)=i-1
      close(21)
      close(22)
      return
      end

c To sort the data
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

c To perform the cramer matrix.
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

c To reading the aperture correction photometry file.
      subroutine readapc(input,nimage,apcrad)
      parameter (ns=100000,ni=30)
      character input*16,temp(20)*3,title*7,cmer*5
      integer nimage,nhead,naperture,numap,nbad,nindef
      real*8 fitrad,ann,dann,ap(16),corap,readap,skip,apcrad(ni)
      common/ap/amerr(ns),amag(ns),axc(ns),ayc(ns),nap,pmag(ns)
      corap=apcrad(nimage)
      open(31,file=input)
      nhead=0
      do i=1,100
       read(31,'(t1,12a1)') (temp(j),j=1,12)
       if(temp(1).eq.'#') then
        nhead=nhead+1
        if(temp(2).eq.'K') then
         write(title,'(7a1)') (temp(j),j=4,10)
         if(title.eq.'ANNULUS') then
          backspace(31)
          read(31,'(t17,f7.3)') ann
          elseif(title.eq.'DANNULU') then
          backspace(31)
          read(31,'(t17,f7.3)') dann
         endif
        endif
        else
        backspace(31)
        goto 100
       endif
      enddo
100   fitrad=ann+dann

      read(31,'(t1,a3)') temp(1)
      do i=2,20
       read(31,'(t1,a3)') temp(i)
       if(temp(i).eq.temp(1)) then
        goto 200
       endif
      enddo
200   naperture=i-5

      do i=1,naperture+1
       backspace(31)
      enddo
      do i=1,naperture
       read(31,*) ap(i)
       if(ap(i).eq.corap) then
        numap=i
        goto 250
       endif
      enddo
250   rewind(31)
      nbad=0
      nindef=0
      do i=1,nhead
       read(31,'(a7)') title
      enddo

      do 300 i=1,ns
       read(31,'(a7)',end=500) title
       read(31,*) axc(i),ayc(i)
       read(31,'(a7)') title
       read(31,'(a7)') title
       read(31,*,end=500) readap,skip,skip,skip,pmag(i)
       do j=2,naperture
        if(j.eq.numap) then
         read(31,'(t52,a5)',end=500) cmer
         if (cmer.eq.'INDEF') then
          nindef=nindef+1
          goto 300 
          else
          backspace(31)
         endif
         read(31,*,end=500) readap,skip,skip,skip,amag(i),amerr(i)
         if(amerr(i).gt.0.01d0) then
          nbad=nbad+1
         endif
         else
         read(31,'(a7)',end=500) title
        endif
       enddo 
300   continue
500   close(31)
      nap=i-1
      if(nbad.gt.0) then
       print *,'Warring. ',input,'has ',nbad,' of high merr stars higher
     * than 0.01.'
      endif
      return
      end

c     Find aperture correction stars on allstar file.
      subroutine apcxid(nimage,nstar,ID,xc,yc,mag,merr)
      parameter (ns=100000,ni=30)
      integer nimage,ID(ni,ns),np,nstar(ni),nbp,ndp,ncor
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns)
      real*8 merr(ni,ns),dist,dmag
      common/apc/apmerr(ni,ns),apmag(ni,ns),apxc(ni,ns),apyc(ni,ns),napc
     *(ni),psfmag(ni,ns),cormean(ni),apcor(ni,ns),apcerr(ni,ns)
      common/ap/amerr(ns),amag(ns),axc(ns),ayc(ns),nap,pmag(ns)
      write(41,*)nap,'apc stars found from ',nimage,'th apc file.'
      nbp=0
      ndp=0
      ncor=0
      cormean(nimage)=0.
      open(41,file='apcxid.fs')
      do i=1,nap
       np=0
       do j=1,nstar(nimage)
        dist=sqrt((xc(nimage,j)-axc(i))**2.d0+
     *  (yc(nimage,j)-ayc(i))**2.d0)
        dmag=abs(mag(nimage,j)-pmag(i))
        if(dist.lt.0.5d0.and.dmag.lt.0.5d0) then
         np=np+1
         ncor=ncor+1
         apxc(nimage,ncor)=axc(i)
         apyc(nimage,ncor)=ayc(i)
         psfmag(nimage,ncor)=pmag(i)
         apmag(nimage,ncor)=amag(i)
         apmerr(nimage,ncor)=amerr(i)
         apcor(nimage,ncor)=apmag(nimage,ncor)-mag(nimage,j)
         apcerr(nimage,ncor)=sqrt(apmerr(nimage,ncor)**2.
     *   +merr(nimage,j)*2.)
         cormean(nimage)=cormean(nimage)+apcor(nimage,ncor)
        endif
       enddo
       if(np.eq.1) then
        elseif(np.eq.0) then
        nbp=nbp+1
        elseif(np.ge.2) then
        ndp=ndp+1
       endif
      enddo

      napc(nimage)=ncor
      cormean(nimage)=cormean(nimage)/real(ncor)
      i=i-1
      if(nbp.ge.1) then
       write(41,*) nbp,'of apc correction star found no partner!!!'
       elseif(nbp.eq.0.and.ndp.eq.0.and.nap.gt.16) then
      endif
      if(ndp.ge.1) then
       write(41,*)ndp,'of apc correction star found bad partners!!!'
      endif
      close(41)
      return
      end

c performing aperture correction.
      subroutine apcorrecting(nimage,nstar,xc,yc,mag,merr,cmag)
      parameter (ns=100000,ni=30)
      integer nimage,nstar(ni)
      real*8 xc(ni,ns),yc(ni,ns),mag(ni,ns),cmag(ni,ns)
      real*8 dist,wd,we,wsum,apsum
      common/apc/apmerr(ni,ns),apmag(ni,ns),apxc(ni,ns),apyc(ni,ns)
     *,napc(ni),psfmag(ni,ns),cormean(ni),apcor(ni,ns),apcerr(ni,ns)
      do i=1,nstar(nimage)
       wsum=0.d0
       apsum=0.d0
       do j=1,napc(nimage)
        dist=sqrt((xc(nimage,i)-apxc(nimage,j))**2.d0+
     *  (yc(nimage,i)-apyc(nimage,j))**2.d0)
        we=apcerr(nimage,j)
        wd=exp(-dist/150.d0)
        apsum=apsum+apcor(nimage,j)*we*wd
        wsum=wsum+we*wd
       enddo
       cmag(nimage,i)=mag(nimage,i)+apsum/wsum
      enddo
      return
      end


      subroutine wcombpar(ninput,IDF,MAGF,dat,ut,filter,nstar)
      parameter (ns=100000,ni=30)
      character filter(ni)*16
      integer ninput,IDF(ns,ni),nfnum(20,ni),n,nstar(ni),nwc,nf(ni),nc
      real*8 MAGF(ns,ni,8)             
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),wcxc(ni,ns),
     *wcyc(ni,ns),nfc(ni,2),IDW(ni,ns),ims(ni)
      common/wwc/nfilter
      common/dmm/master(ni),dm(ni,ni)
 
      do i=1,ni
       ims(i)=0
       nf(i)=0
       do j=1,20
        nfnum(j,i)=0
       enddo
      enddo    
         
      nfilter=0
      do j=1,20 
       n=0
       do i=1,ninput
        if(nfc(i,1).eq.j) then
         n=n+1
         nfnum(j,n)=i
        endif
       enddo
       if(nfnum(j,1).ne.0) then
        nfilter=nfilter+1
       endif
      enddo
    
      nc=0
      do i=1,10
       n=0
       do j=1,ninput
        if(nfnum(i,j).ne.0) then
         n=n+1
         ims(n)=nfnum(i,j)
         if(n.eq.1) then
          nc=nc+1
         endif
        endif
       enddo
       nf(i)=nc
       if(n.eq.0) then
        nf(i)=0
        goto 100
       endif
       call wcomb(nf(i),ims,IDF,MAGF,dat,nstar)
100   nc=nc
      enddo   
      
      return
      end 

 
      subroutine cal_dm(nf,nfimage,ims,IDF,MAGF,dat,nstar)
      parameter (ns=100000,ni=30)
      common/dmm/master(ni),dm(ni,ni)
      real*8 wsum,dmsum,mmax1,mmax2,mmin1,mmin2
      real*8 MAGF(ns,ni,8),wt,ch1,ch2,dat(ni,8)
      integer nfimage,IDF(ni),nstar(ni),ims(ni),nf,it1,it2,nit

      do i=1,nstar(master(nf))
       if(MAGF(i,master(nf),5).ne.0.d0) then
        mmax1=MAGF(i,master(nf),5)
        mmin1=MAGF(i,master(nf),5)
        goto 10
       endif
      enddo
10    do i=1,nstar(master(nf))
       if(MAGF(i,master(nf),5).gt.mmax1) then
        mmax1=MAGF(i,master(nf),5)
       endif
       if(MAGF(i,master(nf),5).lt.mmin1.and.
     & MAGF(i,master(nf),5).ne.0.d0)then
        mmin1=MAGF(i,master(nf),5)
c        print *,MAGF(i,master(nf),1),MAGF(i,master(nf),2)
       endif
      enddo


      do i=1,nfimage 
       nit=0
        wt=0.
        wsum=0.
        dmsum=0.
        do j=1,nstar(ims(i))
         if(MAGF(j,ims(i),5).ne.0.d0) then
          mmax2=MAGF(j,ims(i),5)
          mmin2=MAGF(j,ims(i),5)
          goto 50
         endif
        enddo
50      do j=1,nstar(ims(i))
         if(MAGF(j,ims(i),5).gt.mmax2) then
          mmax2=MAGF(j,ims(i),5)
         endif
         if(MAGF(j,ims(i),5).lt.mmin2.and.MAGF(j,ims(i),5).ne.0.d0) then
          mmin2=MAGF(j,ims(i),5)
c          print *,MAGF(j,ims(i),1),MAGF(j,ims(i),2)
         endif
        enddo
        write(*,51) nf,i,mmax1,mmin1,mmax2,mmin2
51    format(I2,I2,4f8.4)
        do j=1,nstar(ims(i))
         it1=0
         it2=0
c         if(MAGF(j,ims(i),5).ne.0.d0.and.
c     &    MAGF(j,master(nf),5).ne.0.d0.and.
c     &    MAGF(j,ims(i),6).lt.0.05d0.and.
c     &    MAGF(j,master(nf),6).lt.0.05d0) then
         if(MAGF(j,ims(i),5).ne.0.d0.and.
     &    MAGF(j,master(nf),5).ne.0.d0.and.
     &    MAGF(j,ims(i),5).gt.mmin2+0.35d0.and.
     &    MAGF(j,ims(i),5).lt.mmax2-1.0d0.and. 
     &    MAGF(j,master(nf),5).gt.mmin1+0.35d0.and.
     &    MAGF(j,master(nf),5).lt.mmax1-1.0d0) then
          nit=nit+1
          if(MAGF(j,ims(i),6).lt.0.005d0) then
           it1=1
           ch1=MAGF(j,ims(i),6) 
           MAGF(j,ims(i),6)=0.005d0
          endif
          if(MAGF(j,master(nf),6).lt.0.005d0) then
           it2=1
           ch2=MAGF(j,master(nf),6)
           MAGF(j,master(nf),6)=0.005d0
          endif
          wt=sqrt(0.005d0**2.d0+0.005d0*2.d0)/
     &    sqrt(MAGF(j,ims(i),6)**2.d0+MAGF(j,master(nf),6)**2.d0)
          wsum=wsum+wt
          dmsum=dmsum+wt*(MAGF(j,master(nf),6)-MAGF(j,ims(i),6))
          if(it1.eq.1) then
           MAGF(j,ims(i),6)=ch1
           it1=0
          endif
          if(it2.eq.2) then
           MAGF(j,master(nf),6)=ch2
           it2=0
          endif
         endif
        enddo
        dm(nf,i)=dmsum/wsum
c       write(*,'(i3,i3,f9.5,i5)'),nf,i,dm(nf,i),nit
      enddo
      return
      end


      subroutine wcomb(ninput,ims,IDF,MAGF,dat,nstar)
      parameter (ns=100000,ni=30)
      integer ninput,IDF(ns,ni),wgt(ni),nstar(ni),nfilt
      integer ims(ni)
      real*8 MAGF(ns,ni,8),dat(ni,8)
      real wsum,wmsum,wesum,wxcsum,wycsum,wt
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),wcxc(ni,ns),
     *wcyc(ni,ns),nfc(ni,2),IDW(ni,ns)
      common/rcoo/rx(ni,ns),ry(ni,ns)
      common/dmm/master(ni),dm(ni,ni)
      common/wwc/nfilter
      do i=1,ni
       if(ims(i).eq.0) then
        goto 100
       endif
      enddo
100   nfilt=i-1
      do i=1,ni
       do j=1,ni
        dm(i,j)=0.
       enddo
      enddo
      call selectmaster(ninput,master,ims,dat,ut,nstar)
      call cal_dm(ninput,nfilt,ims,IDF,MAGF,dat,nstar)

      do i=1,ns
       wsum=0.
       wmsum=0.
       wesum=0.
       wxcsum=0.
       wycsum=0.
       do j=1,nfilt  
        if(MAGF(i,ims(j),6).le.0.005d0) then  
         wt=1.
         elseif(MAGF(i,ims(j),6).gt.90.d0.or.
     &   MAGF(i,ims(j),5).gt.90.d0) then
         wt=0.
         elseif(MAGF(i,ims(j),6).gt.0.005d0) then
         wt=0.005/MAGF(i,ims(j),6) 
        endif
        if(IDF(i,ims(j)).eq.0) then
         wt=0.
        endif
        wsum=wsum+wt
        wmsum=wmsum+(MAGF(i,ims(j),5)+dm(j,ims(j)))*wt
        wesum=wesum+MAGF(i,ims(j),6)*wt
        wxcsum=wxcsum+MAGF(i,ims(j),1)*wt
        wycsum=wycsum+MAGF(i,ims(j),2)*wt
       enddo
       if(wsum.gt.0.) then
        IDW(ninput,i)=IDF(i,master(ninput))
        wcmag(ninput,i)=wmsum/wsum
        wcmerr(ninput,i)=wesum/wsum
        wcxc(ninput,i)=wxcsum/wsum
        wcyc(ninput,i)=wycsum/wsum
        elseif (wsum.eq.0.) then
        IDW(ninput,i)=0.
        wcmag(ninput,i)=0.
        wcmerr(ninput,i)=0.
        wcxc(ninput,i)=0.
        wcyc(ninput,i)=0.
       endif
      enddo

500   return
      end


      subroutine writewcomb(input,filter,ut,dat,IDF,nstar)
      parameter(ns=100000,ni=30)
      real*8 dat(ni,8),ut(ni)
      character filter(ni)*16,input(ni)*16
      integer IDF(ns,ni),nid,nfnum(20,ni),ims(ni),nstar(ni)
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),wcxc(ni,ns),
     *wcyc(ni,ns),nfc(ni,2),IDW(ni,ns)
      common/rcoo/rx(ni,ns),ry(ni,ns)
      common/wwc/nfilter
      common/dmm/master(ni),dm(ni,ni)
      open(62,file='wcout1.fs')
      open(63,file='wcout2.fs')
c     1st file header##################
      write(62,*) " F : Filrer"
      write(62,*) " T : UT"
      write(62,*) " X : Airmass"
c     1st file header##################
c     2nd file header##################
      write(63,*) ""
      write(63,*) " M : Master image"
      write(63,*) " F : Filrer"
      write(63,*) " T : UT"
      write(63,*) " X : Airmass"
      write(63,*) " *"
      write(63,*) ' ID(of all detected stars)'
      write(63,*) ' ID(on the master image of each filters)'
      write(63,*) ' X center(on the Xid master image)'
      write(63,*) ' Y center(on the Xid master image)'
      write(63,*) ' Mag(aperture corrected)'
      write(63,*) ' Merr(combined PSF photometry & aperture correctoin)'
      write(63,*) ""
c     2nd file header##################
      write(62,'(t1,a1,a6,30a10)')'F','',(filter(master(i)),i=1,nfilter)
      write(62,'(t1,a1,30a10)')   'M',(input(master(i)),i=1,nfilter)
      write(62,'(t1,a1,30f10.6)')'T',(ut(master(i)),i=1,nfilter)
      write(62,'(t1,a1,30f10.7)')'X',(dat(master(i),5),i=1,nfilter)
      write(62,'(t1,a1)') "*"
      write(63,'(t1,a1,a6,30a10)')'F','',(filter(master(i)),i=1,nfilter)
      write(63,'(t1,a1,30a10)')   'M',(input(master(i)),i=1,nfilter)
      write(63,'(t1,a1,30f10.6)')'T',(ut(master(i)),i=1,nfilter)
      write(63,'(t1,a1,30f10.7)')'X',(dat(master(i),5),i=1,nfilter)
      write(63,'(t1,a1)') "*"   

      do i=1,ns
       nid=0
       do j=1,nfilter
        if(wcmag(j,i).eq.0) then
         nid=nid+1
        endif
       enddo
       if(nid.eq.nfilter) then
        goto 500
       endif
       write(62,'(30i10)')   (IDF(i,master(j)),j=1,nfilter)
       write(62,'(30I10)')   (IDW(j,i),j=1,nfilter)
       write(62,'(30f10.3)') (wcmag(j,i),j=1,nfilter)
       write(62,'(30f10.3)') (wcmerr(j,i),j=1,nfilter)
       write(62,*) ""
       write(63,'(30i10)')   (IDF(i,master(j)),j=1,nfilter)
       write(63,'(30I10)')   (IDW(j,i),j=1,nfilter)
       write(63,'(30f10.3)') (wcxc(j,i),j=1,nfilter)
       write(63,'(30f10.3)') (wcyc(j,i),j=1,nfilter)
       write(63,'(30f10.3)') (wcmag(j,i),j=1,nfilter)
       write(63,'(30f10.3)') (wcmerr(j,i),j=1,nfilter)
       write(63,*) ""
       write(63,*) ""
       write(63,*) ""
       write(63,*) ""
      enddo
500   close(63)
      close(62)
      return
      end

      subroutine selectmaster(nf,master,ims,dat,ut,nstar)
      parameter(ns=100000,ni=30)
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),wcxc(ni,ns),
     *wcyc(ni,ns),nfc(ni,2),IDW(ni,ns)
      real*8 dat(ni,8),ut(ni)
      integer ims(ni),nimage,nstar(ni),nf,master(ni)
     
      do i=1,ni
       if(ims(i).eq.0) then
        goto 100
       endif
      enddo
100   nimage=i-1
      master(nf)=ims(1)
      do i=2,nimage
       if(nstar(ims(i)).gt.nstar(master(nf))) then
        master(nf)=ims(i)
       endif
      enddo
      do i=2,nimage
       if(real(nstar(ims(i))).gt.real(nstar(ims(i-1)))*0.9.or.
     & real(nstar(ims(i))).gt.real(nstar(ims(i-1)))*1.1) then
        if(dat(ims(i),8).lt.dat(ims(i-1),8)) then
         master(nf)=ims(i)
         elseif(dat(ims(i),8).gt.dat(ims(i-1),8)) then
         master(nf)=ims(i-1)
        endif
       endif
      enddo
      return
      end


      subroutine readfilter(nimage,filter)
      parameter(ns=100000,ni=30)
      character H*1,filter(ni)*16,ch(16)*1,tf*1
      integer nimage,nfi,nco,nf
      common/wc/wcmag(ni,ns),wcmerr(ni,ns),wcxc(ni,ns),
     *wcyc(ni,ns),nfc(ni,2),IDW(ni,ns)
      open(51,file='filterpar.fs')
      nf=0
      do i=nimage,ni
       nfc(i,1)=0
       nfc(i,2)=0
      enddo
      do i=1,20
       read(51,'(t1,a1)',end=500) H
       if(H.eq.'*') then
        nf=nf+1
        backspace(51)
        read(51,'(t2,a1,t4,i5,t9,i16)',end=500) tf,nfi,nco
        read(filter(nimage),'(16a1)') (ch(j),j=1,16)
        do j=1,16
         if     (ch(j).eq.tf) then
          nfc(nimage,1)=nfi
          nfc(nimage,2)=nco
         endif
        enddo
       endif
      enddo
500   close(51)
      return
      end
