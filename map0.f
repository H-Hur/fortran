      real coo1(2,999999),mag1(2,999999),SKY1(999999),Nmag1(2,36)
      real P1(3,999999),dat1(7),chi1(36),sharp1(36),Emag1(36),m0,mm
      integer n1,ITER1(999999),ID1(999999),start1,end1,idw,maxm1
      character*16 filename1,INimage1,Pimage1,filter1           
      call sm_device ('postportfile smt.ps')
      call sm_graphics
      call sm_erase()

      call reading(coo1,mag1,SKY1,P1,dat1,n1,ITER1,ID1,filename1,INimage
     *1,Pimage1,filter1)

      call dist(mag1,n1,dat1,P1,start1,end1,maxm1,Nmag1,Emag1,chi1,sharp
     *1)

      call sm_location(11000,17000,29800,30000)
      call sm_limits(0.,100.,0.,100.)
      call sm_relocate(0.,0.)  
      call sm_expand(2.5)
      call sm_label(INimage1)   
 
      call sm_location(3000,30000,25800,28500)
      call header(filename1,INimage1,Pimage1,filter1,dat1,n1)

      call sm_location(3000,30000,6000,25000)
      call maindraw(coo1,mag1,n1,start1,end1,maxm1,ID1)

      call sm_location(3000,11000,1500,4000)
      call Nmagdraw(Nmag1,start1,end1)

      call sm_location(12500,20500,1500,4000)
      call Emagdraw(mag1,n1,m0,mm)

c      call sm_location(22000,30000,1500,5000)
c      call chidraw(mag1,n1,P1,m0,mm)

      call sm_location(22000,30000,1500,4000)
      call sharpdraw(mag1,n1,P1,m0,mm)



      call sm_gflush()
      call sm_hardcopy
      call sm_alpha



      stop
      end



c This subroutine made to read and write the result of task allstar in IRAF 
      subroutine reading(coo,mag,SKY,P,dat,n,ITER,ID,filename,INimage,Pi
     *mage,filter)

      real coo(2,999999),mag(2,999999),SKY(999999)
      real P(3,999999),dat(7)
      character*2 H
      character*16 filename,INimage,Pimage,filter
      integer n,ITER(999999),ID(999999),i  


      print *,'Input allstar file name.'
      read (*,*) filename
      open(21,file=filename,ERR=22)
      print *,' '
      open(22,file='result',err=22)
      goto 30
22    print *,'Failed to open ',filename


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
23    format(A16)                                   !chatacter format
24    format(T17,A16)                               !character format 
25    format (T17,F7.1)                         ! data min,max format
26    format (T17,F5.2)                            !gain,noise format
27    format(I5,2F9.3,F7.2,F7.3,F10.5,I3,F8.3,F7.3,I3)!writing format
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
       read(21,25,err=31) dat(i-13)
       enddo
       do i=16,17!gain,noise                16~17
       read(21,25,err=31) dat(i-13)
       enddo
              do i=18,19!건너뛰기           18~19
              read(21,23,err=31) H
              enddo

      read(21,24,err=31) filter!필터        20
              do i=21,23
              read(21,23,err=31) H!건너뛰기 21~23
              enddo

       do i=24,26!PSFMAG,PSFRAD,FWHM        24~26
       read(21,25,err=31) dat(i-19)
       enddo
              do i=27,44!건너뛰기           27~44
              read(21,23,err=31) H
              enddo
      goto 32

31    print *,'Failed to start reading from ',filename


32    do i=1,999999   !                     45~
      read(21,*,err=33,end=34) ID(i),coo(1,i),coo(2,i),mag(1,i),mag(2,i)
     *,SKY(i),ITER(i)
      read(21,*,err=33,end=34) P(1,i),P(2,i),P(3,i)
      enddo
33    print *,'Failed reading from ',filename
34    print *,'Succesfully read from data. number of stars=',i-1
      n=i-1

      do i=1,n
      Write(22,27,err=35) ID(i),coo(1,i),coo(2,i),mag(1,i),mag(2,i),SKY(
     *i),ITER(i),P(1,i),P(2,I),P(3,i)
      enddo
      goto 36

35    print *,'Failed to write the result.'     
36    print *,'Writed to the [result.dat],succesfully.'    
      print *, ' ' 

      close(21)
      close(22)
      return
      end


c This subroutine made to find max and min of magnitudes, and to count of magnitudes and magnitude errors.
       subroutine dist(mag,n,dat,P,startmag,endmag,maxmag,Nmag,Emag,chi,
     *sharp)
       integer i,j,startmag,endmag,PSF,n,maxmag
       real dat(7),mag(2,999999),Emag(36),chi(36),P(3,999999),sharp(36)
       real Nmag(2,36)

       print *,'Runnibg subroutine dist...'
       open(40,file='mag.dat')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
41     format(F4.1,F8.1,3F10.5)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

       do i=1,36
        Nmag(1,i)=real(i)
        Nmag(2,i)=0.
        Emag(i)=0.
        sharp(i)=0.
        chi(i)=0.
       enddo

       do i=1,n
        do j=5,30
         if (mag(1,i).ge.j.and.mag(1,i).lt.j+1) then
          Nmag(2,j)=Nmag(2,j)+1.
          Emag(j)=Emag(j)+mag(2,i)
          sharp(j)=sharp(j)+P(1,i)
          chi(j) = chi(j)+  P(2,i)
         endif
        enddo
       enddo

       do i=5,30
        Emag(i)=Emag(i)/Nmag(2,i)
        sharp(i)=sharp(i)/Nmag(2,i)
        chi(i)=chi(i)/Nmag(2,i)
       enddo

       PSF=dat(5)

       do i=PSF,30
        if(Nmag(2,i).eq.0.) then
         goto 50
        endif
       enddo
50     endmag=i-1

       do i=PSF,5,-1
        if(Nmag(2,i).eq.0.) then
         goto 51
        endif
       enddo
51     startmag=i+1
       print *,'Counted magnitudes of data..., mag=',startmag,'~',endmag
       maxmag=0
       do i=startmag,endmag
       write (40,41) Nmag(1,i)+0.5,Nmag(2,i),Emag(i),sharp(i),chi(i)
        if(Nmag(2,i).ge.maxmag) then
        maxmag=Nmag(1,i)
        endif
       enddo

       print *,'Writed Mag,Nmag,Emag,charpness,chi to [mag.dat].'
       print *,' '
       close(40)
       return
       end

!This subroutine will draw starrymap in the field of view.
      subroutine maindraw(coo,mag,n,startmag,endmag,maxmag,ID)
      integer n,i,j,startmag,endmag,ID(999999),np,idw,maxmag,tf,tb
      real coo(2,999999),mag(2,999999),xm,ym,xcc
      real pp(1),f,siz(1),bright,faint,bt
      character*6 wid
      print *,'input brightper limit of magniture, highest=',startmag
      read(*,*) bright
       if(bright.lt.startmag) then
       bt=startmag-0.99
       bright=startmag
       else
       bt=bright
       endif
      print *,'input fainter limit of magnitude, lowest=',endmag
      read(*,*) faint
       if (faint.gt.endmag) then
       maxmag=endmag+0.99
       elseif(faint.lt.endmag) then
       maxmag=faint
       endif
      print *,'Do you wand to put ID on the map? 1=yes, 0=no'
      read(*,*) idw
      print *,'drawing main map...'
      print *,' ' 
      xm=coo(1,1)
      ym=coo(2,1)
      
      do i=1,n
       if (coo(1,i).gt.xm) then
        xm=coo(1,i)
       endif
       if (coo(2,i).gt.ym) then
        ym=coo(2,i)
       endif
      enddo
      ym=ym+5.
      xm=xm+5.
            
      call sm_limits(0.,xm,0.,ym)
      call sm_ctype('black')
      call sm_expand(0.6)
      call sm_xlabel('Xcenter(pixel)')
      call sm_ylabel('Ycenter(pixel)')
      call sm_box(1,2,1,2)
      call sm_expand(1.)
      np=0
      tf=0
      tb=0
      do i=1,n
        if(mag(1,i).le.faint.and.mag(1,i).ge.bt) then
          np=np+1
          f =(-mag(1,i)+real(bright) )/real(maxmag-bright)+1.
          if (f.ge.0.999)then
           f=0.99
           tb=tb+1
          elseif (f.le.0.)then
           f=0.01
           tf=tf+1
          endif 
          write(wid,'(I6)') ID(i)
          call sm_expand(2.)
           if(idw.eq.1) then
           pp(1) = 403. + f
           call sm_ptype(pp,1)
           else
           pp(1)=400. + f
           call sm_ptype(pp,1)
           endif
          call sm_points(coo(1,i),coo(2,i),1)
          siz(1)=1.-(mag(1,i)-real(bright))/real(faint-bright)
          call sm_expand(siz)
          xcc=coo(1,i)-80.*(xm-5.)/2010.15295
          call sm_relocate(xcc,coo(2,i))
           if (idw.eq.1) then
           call sm_label(wid)
           endif 
        goto 100
        endif
100   np=np
      enddo
      print *,'nprint=',np,', nbrightest=',tb,', nfaintest=',tf
      print *,'Created the main map to file [smt.ps].'
      print *,' '
        
      return
      end    

!This subroutine will draw a chart which shows that how many start observated per a magnitude.
      subroutine Nmagdraw(Nmag,startmag,endmag)
      real Nmag(2,36),ylimits,pp(1)
      integer startmag,endmag
      ylimits=0.

      do i=startmag,endmag
      if (Nmag(2,i).gt.ylimits) then
      ylimits=Nmag(2,i)
      endif
      enddo
 
      print *,'Drawing Number of stars per magnitude histogram...' 

      call sm_ctype('black')
      call sm_limits(real(startmag),real(endmag),0.,ylimits)
      call sm_expand(0.4)
      call sm_box(1,2,3,3)
      call sm_xlabel('mag')
      call sm_ylabel('Nstar')
      call sm_expand(22.)
      pp(1)=23.+0.1
      call sm_ptype(pp,1)

      do i=startmag,endmag
      call sm_points(real(i),Nmag(2,i),1)
      enddo
      call sm_expand(1.)

      print *,'Created Number of stars per magnitude histogram to [smt.p
     *s].'
      print *,' '

      return
      end
  
!This subroutine will show that there are how much errors in magnitudes.
      subroutine Emagdraw(mag,n,m0,mm)
      real ylimits,pp(1),mag(2,999999),m0,mm
      integer n
      ylimits=0.

      do i=1,n
      if (mag(2,i).gt.ylimits) then
      ylimits=mag(2,i)
      endif
      enddo
      m0=mag(1,1)
      mm=m0
      do i=1,n
      if (mag(1,i).gt.mm) then
      mm=mag(1,i)
      elseif (mag(1,i).lt.m0) then
      m0=mag(1,i)
      endif
      enddo

      print *,'Drawing magerror per magnitude histogram...'

      call sm_ctype('black')
      call sm_limits(m0,mm,0.,ylimits)
      call sm_expand(0.4)
      call sm_xlabel('mag')
      call sm_ylabel('magerror')
      call sm_box(1,2,3,3)
      call sm_expand(1.)
      pp(1)=403.+0.1
      call sm_ptype(pp,1)

      do i=1,n
      call sm_points(mag(1,i),mag(2,i),1)
      enddo
      call sm_expand(1.)

      print *,'Created magerror per magnitude histogram to [smt.ps].'
      print *,' '

      return
      end

!This subroutine will show that there is how many chi in magnitudes.
      subroutine chidraw(mag,n,P,m0,mm)
      real chi(36),ylimits,pp(1),mag(2,999999),P(3,999999),m0,mm,ylim
      integer n
      ylimits=0.
      ylim=0.

      do i=1,n
      if (P(2,i).gt.ylimits) then
      ylimits=P(2,i)
      elseif (P(2,i).lt.ylim) then
      ylim=P(2,i)
      endif
      enddo

      print *,'Drawing chi per magnitude histogram...'

      call sm_ctype('black')
      call sm_limits(m0,mm,ylim,ylimits)
      call sm_expand(0.4)
      call sm_xlabel('mag')
      call sm_ylabel('chi')
      call sm_box(1,2,3,3)
      call sm_expand(1.)
      pp(1)=403.+0.1
      call sm_ptype(pp,1)

      do i=1,n
      call sm_points(mag(1,i),P(2,i),1)
      enddo
      call sm_expand(1.)

      print *,'Created chi per magnitude histogram to [smt.ps].'
      print *,' '

      return
      end

!This subroutine will show that sharpness of stars are much or not.
      subroutine sharpdraw(mag,n,P,m0,mm)
      real sharp(36),ylimits,pp(1),mag(2,999999),P(3,999999),m0,mm,ylim
      integer n
      ylimits=0.
      ylim=0.
      do i=1,n
      if (P(1,i).gt.ylimits) then
      ylimits=P(1,i)
      elseif (P(1,i).lt.ylim) then
      ylim=P(1,i)
      endif
      enddo

      print *,'Drawing sharpness per magnitude histogram...'

      call sm_ctype('black')
      call sm_limits(m0,mm,ylim,ylimits)
      call sm_expand(0.4)
      call sm_xlabel('mag')
      call sm_ylabel('sharpness')
      call sm_box(1,2,3,3)
      call sm_expand(1.)
      pp(1)=403.+0.1
      call sm_ptype(pp,1)

      do i=1,n
      call sm_points(mag(1,i),P(1,i),1)
      enddo
      call sm_expand(1.)

      print *,'Created sharpness per magnitude histogram to [smt.ps].'
      print *,' '

      return
      end


      subroutine header(filename,INimage,Pimage,filter,dat,n)
      real dat(7)!,dmin,dmax,gain,nois,pmag,prad,fwhm
      integer n
      character*16 filename,INimage,Pimage,filter,nums
      character*16 dmin,dmax,gain,nois,pmag,prad,fwhm
    
      write(dmin,'(f6.1)') dat(1)
      write(dmax,'(f8.1)') dat(2)
      write(gain,'(f5.2)') dat(3)
      write(nois,'(f5.1)') dat(4)
      write(pmag,'(f5.2)') dat(5)
      write(prad,'(f5.2)') dat(6)
      write(fwhm,'(f5.2)') dat(7)
      write(nums,'(I6)') n
      print *,'Writing labels...'

      call sm_limits(0.,100.,0.,100.)
      call sm_box(3,3,3,3)
      call sm_expand(0.6)

      call sm_relocate(1.,85.)
      call sm_label('Filename : ')
      call sm_relocate(22.,85.)
      call sm_label(filename)

      call sm_relocate(1.,70.)
      call sm_label('Input image : ')
      call sm_relocate(22.,70.)
      call sm_label(INimage)

      call sm_relocate(1.,55.)
      call sm_label('datamin : ')
      call sm_relocate(22.,55.)
      call sm_label(dmin)

      call sm_relocate(1.,40.)
      call sm_label('datamax : ')
      call sm_relocate(22.,40.)
      call sm_label(dmax)

      call sm_relocate(1.,25.)
      call sm_label('Gain :')
      call sm_relocate(22.,25.)
      call sm_label(gain)

      call sm_relocate(1.,10.)
      call sm_label('Readout Noise :')
      call sm_relocate(22.,10.)
      call sm_label(nois)

      call sm_relocate(50.,85.)
      call sm_label('Filter : ')
      call sm_relocate(72.,85.)
      call sm_label(filter)

      call sm_relocate(50.,70.)
      call sm_label('Number of stars : ')
      call sm_relocate(70.,70.)
      call sm_label(nums)

      call sm_relocate(50.,55.)
      call sm_label('PSF image : ')
      call sm_relocate(72.,55.)
      call sm_label(Pimage)

      call sm_relocate(50.,40.)
      call sm_label('PSF mag : ')
      call sm_relocate(72.,40.)
      call sm_label(pmag)

      call sm_relocate(50.,25.)
      call sm_label('PSF rad : ')
      call sm_relocate(72.,25.)
      call sm_label(prad)

      call sm_relocate(50.,10.)
      call sm_label('FWHM : ')
      call sm_relocate(72.,10.)
      call sm_label(fwhm)

      return
      end
 
