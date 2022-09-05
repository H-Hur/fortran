c for only fortran77 to draw color-color diagram.
c f77 -o 'map.exe' 'map.f' -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 


      real coo1(2,999999),mag1(2,999999),SKY1(999999),Nmag1(2,36)
      real P1(3,999999),dat1(7),chi1(36),sharp1(36),Emag1(36),m0,mm
      real fai1,bri1
      integer n1,ITER1(999999),ID1(999999),start1,end1,idw,maxm1
      character*16 filename1,INimage1,Pimage1,filter1
      character output*15,ps*28

      print *,'Ver. 1.3 ,   Modified at Sep,27, 2007, by Hyeon-Oh Hur'
      print *,''
      call reading(coo1,mag1,SKY1,P1,dat1,n1,ITER1,ID1,filename1,INimage
     *1,Pimage1,filter1,output)
      write(ps,'(a13,a15)') 'postportfile ',output
      print *,ps
      call sm_device (ps)
      call sm_graphics
      call sm_erase()

      call dist(mag1,n1,dat1,P1,start1,end1,maxm1,Nmag1,Emag1,chi1,sharp
     *1)
 
      call sm_location(3000,30000,27300,30600)
      call header(filename1,INimage1,Pimage1,filter1,dat1,n1,output)

      call sm_location(3000,30000,7550,26550)
      call maindraw(coo1,mag1,n1,start1,end1,maxm1,bri1,fai1,ID1,output)
 
      call sm_location(3000,30000,6500,7500)
      call Csize(maxm1,bri1) 
 
      call sm_location(3000,15000,1500,3600)
      call Nmagdraw(Nmag1,start1,end1,output)

      call sm_location(18000,30000,1500,3600)
      call Emagdraw(mag1,n1,m0,mm,output)

      call sm_location(3000,15000,4110,6250)
      call chidraw(mag1,n1,P1,m0,mm,output)

      call sm_location(18000,30000,4110,6250)
      call sharpdraw(mag1,n1,P1,m0,mm,output)


      call sm_gflush()
      call sm_hardcopy
      call sm_alpha


      stop
      end


c This subroutine made to read and write the result of task allstar in IRAF 
      subroutine reading(coo,mag,SKY,P,dat,n,ITER,ID,filename,INimage,Pi
     *mage,filter,output)

      real coo(2,999999),mag(2,999999),SKY(999999)
      real P(3,999999),dat(7)
      character H*2,output*15,ch(15)*1
      character*16 filename,INimage,Pimage,filter
      integer n,ITER(999999),ID(999999),i  


      print *,'Input allstar file name.'
      read (*,*) filename
      open(21,file=filename,ERR=22)
      read(filename,'(15a1)') (ch(i),i=1,15)
      do i=1,15
       if(ch(i).eq.'.') then
        goto 10
       endif
      enddo
10    write(ch(i+1),'(a1)') 'p'
      write(ch(i+2),'(a1)') 's'
      write(output,'(15a1)') (ch(j),j=1,i+2)
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
      subroutine maindraw(coo,mag,n,startmag,endmag,maxmag,bright,faint,
     *ID,output)
      integer n,i,j,startmag,endmag,ID(999999),np,idw,maxmag,tf,tb
      real coo(2,999999),mag(2,999999),xm,ym,xcc
      real pp(1),f,siz(1),bright,faint,bt,xl,yl,x1(1),y1(1)
      character*6 wid
      character output*15
      common /id/idw
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
c          print *,mag(1,i),f
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
          siz(1)=1.-(mag(1,i)-bright)/(faint-bright)
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
      print *,'Created the main map to file [',
     *output,'].'
      print *,' '
      return
      end    

      subroutine Csize(maxmag,bright)
      integer maxmag
      real bright,f,pp(1),x(1),y(1),p,t
      character mags*5
      common /id/idw
      t=0.1*(real(maxmag)-bright)
      call sm_limits(0.,100.,0.,100.)
      if(idw.eq.1) then
       p=403.
      else
       p=400.
      endif

      call sm_expand(2.)
      f =-0.01/(real(maxmag)-bright)+1.
      pp(1)=p+f
      x(1)=70.
      y(1)=60.
      call sm_ptype(pp,1)
      call sm_points(x(1),y(1),1)    
      write(mags,'(f4.1)') bright+0.01
      call sm_relocate(x(1)-2.0,y(1)-55.)
      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(2.)        

      f =-2.5*t/(real(maxmag)-bright)+1.
      pp(1)=p+f
      x(1)=80.
      y(1)=60.
      call sm_ptype(pp,1)
      call sm_points(x(1),y(1),1)
      write(mags,'(f4.1)') bright+2.5*t
      call sm_relocate(x(1)-2.0,y(1)-55.)
      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(2.)

      f =-5.*t/(real(maxmag)-bright)+1.
      pp(1)=p+f
      x(1)=85.
      y(1)=60.
      call sm_ptype(pp,1)
      call sm_points(x(1),y(1),1)
      write(mags,'(f4.1)') bright+5.*t
      call sm_relocate(x(1)-2.0,y(1)-55.)
      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(2.)

      f =-7.5*t/(real(maxmag)-bright)+1.
      pp(1)=p+f
      x(1)=90.
      y(1)=60.
      call sm_ptype(pp,1)
      call sm_points(x(1),y(1),1)
      write(mags,'(f4.1)') bright+7.5*t
      call sm_relocate(x(1)-2.0,y(1)-55.)
      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(2.)

      f =-9.90*t/(real(maxmag)-bright)+1.
      pp(1)=p+f
      x(1)=95.
      y(1)=60.
      call sm_ptype(pp,1)
      call sm_points(x(1),y(1),1)
      write(mags,'(f4.1)') bright+9.9*t
      call sm_relocate(x(1)-2.0,y(1)-55.)
      call sm_expand(0.5)
      call sm_label(mags)
      call sm_expand(1.)

      return
      end

!This subroutine will draw a chart which shows that how many start observated per a magnitude.
      subroutine Nmagdraw(Nmag,startmag,endmag,output)
      real Nmag(2,36),ylimits,pp(1)
      integer startmag,endmag
      character output*15
      ylimits=0.

      do i=startmag,endmag
      if (Nmag(2,i).gt.ylimits) then
      ylimits=Nmag(2,i)
      endif
      enddo
 
      print *,'Drawing Number of stars per magnitude histogram...' 

      call sm_ctype('black')
      call sm_limits(real(startmag),real(endmag),0.,ylimits+50.)
      call sm_expand(0.44)
      call sm_box(1,2,3,3)
      call sm_expand(0.5)
      call sm_xlabel('mag')
      call sm_ylabel('Nstar')
      call sm_expand(22.)
      pp(1)=23.+0.1
      call sm_ptype(pp,1)

      do i=startmag,endmag
      call sm_points(real(i),Nmag(2,i),1)
      enddo
      call sm_expand(1.)
      print *,'Created Number of stars per magnitude histogram to [',
     *output,'].'
      print *,' '

      return
      end
  
!This subroutine will show that there are how much errors in magnitudes.
      subroutine Emagdraw(mag,n,m0,mm,output)
      real ylimits,pp(1),mag(2,999999),m0,mm
      integer n
      character output*15
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
      call sm_expand(0.5)
      call sm_xlabel('mag')
      call sm_ylabel('merr')
      call sm_expand(0.45)
      call sm_box(1,2,3,3)
      call sm_expand(1.)
      pp(1)=403.+0.1
      call sm_ptype(pp,1)
      do i=1,n
      call sm_points(mag(1,i),mag(2,i),1)
      enddo
      call sm_expand(1.)
      print *,'Created magerror per magnitude histogram to [',
     *output,'].'
      print *,' '

      return
      end

!This subroutine will show that there is how many chi in magnitudes.
      subroutine chidraw(mag,n,P,m0,mm,output)
      real chi(36),ylimits,pp(1),mag(2,999999),P(3,999999),m0,mm,ylim
      integer n
      character output*15
      ylimits=0.
      ylim=0.
      print *,'1'
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
      call sm_expand(0.5)
c      call sm_xlabel('mag')
      call sm_ylabel('chi')
      call sm_box(1,2,3,3)
      call sm_expand(1.)
      pp(1)=403.+0.1
      call sm_ptype(pp,1)

      do i=1,n
      call sm_points(mag(1,i),P(2,i),1)
      enddo
      call sm_expand(1.)

      print *,'Created chi per magnitude histogram to [',
     *output,'].'
      print *,' '

      return
      end

!This subroutine will show that sharpness of stars are much or not.
      subroutine sharpdraw(mag,n,P,m0,mm,output)
      real sharp(36),ylimits,pp(1),mag(2,999999),P(3,999999),m0,mm,ylim
      integer n
      character output*15
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
      call sm_expand(0.5)
c      call sm_xlabel('mag')
      call sm_ylabel('sharpness')
      call sm_box(1,2,3,3)
      call sm_expand(1.)
      pp(1)=403.+0.1
      call sm_ptype(pp,1)

      do i=1,n
      call sm_points(mag(1,i),P(1,i),1)
      enddo
      call sm_expand(1.)

      print *,'Created sharpness per magnitude histogram to [',
     *output,'].'
      print *,' '

      return
      end


      subroutine header(filename,INimage,Pimage,filter,dat,n,output)
      real dat(7)!,dmin,dmax,gain,nois,pmag,prad,fwhm
      integer n
      character*16 filename,INimage,Pimage,filter,nums
      character*16 dmin,dmax,gain,nois,pmag,prad,fwhm
      character output*15    
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

      call sm_relocate(1.,89.)
      call sm_label('Filename : ')
      call sm_relocate(22.,89.)
      call sm_label(filename)

      call sm_relocate(1.,76.)
      call sm_label('Input image : ')
      call sm_relocate(22.,76.)
      call sm_label(INimage)

      call sm_relocate(1.,63.)
      call sm_label('datamin : ')
      call sm_relocate(22.,63.)
      call sm_label(dmin)

      call sm_relocate(1.,50.)
      call sm_label('datamax : ')
      call sm_relocate(22.,50.)
      call sm_label(dmax)

      call sm_relocate(1.,37.)
      call sm_label('Gain :')
      call sm_relocate(22.,37.)
      call sm_label(gain)

      call sm_relocate(1.,24.)
      call sm_label('Readout Noise :')
      call sm_relocate(22.,24.)
      call sm_label(nois)

c      call sm_relocate(50.,11.)
c      call sm_label('')
c      call sm_relocate(72.,11.)
c      call sm_label()


      call sm_relocate(50.,89.)
      call sm_label('Filter : ')
      call sm_relocate(72.,89.)
      call sm_label(filter)

      call sm_relocate(50.,76.)
      call sm_label('Number of stars : ')
      call sm_relocate(70.,76.)
      call sm_label(nums)

      call sm_relocate(50.,63.)
      call sm_label('PSF image : ')
      call sm_relocate(72.,63.)
      call sm_label(Pimage)

      call sm_relocate(50.,50.)
      call sm_label('PSF mag : ')
      call sm_relocate(72.,50.)
      call sm_label(pmag)

      call sm_relocate(50.,37.)
      call sm_label('PSF rad : ')
      call sm_relocate(72.,37.)
      call sm_label(prad)

      call sm_relocate(50.,24.)
      call sm_label('FWHM : ')
      call sm_relocate(72.,24.)
      call sm_label(fwhm)

c      call sm_relocate(50.,11.)
c      call sm_label('')
c      call sm_relocate(72.,11.)
c      call sm_label()



      return
      end
 
