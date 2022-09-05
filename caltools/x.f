      real*8 shrate,s(30),lat,long,height,X(12,90),ra3,dec3,hour,mini
      real*8 startt,endt,dt,pt,cz,sz,rad,dec(90),ra(90),st,time(12)
      integer nsts,ld,ra1,ra2,dec1,dec2,nday,nprint,nprinted,npage,nt
      character H*1,starname(90)*12,output(9)*11
      open(21,file='input.air')
      read(21,'(a1)') H
      read(21,'(a1)') H
10    format(t1,a7,t14,i3,t18,i3,t22,f7.3,t30,i3,t34,f7.3)
      do i=1,200
       read(21,10,end=50) starname(i),ra1,ra2,ra3,dec1,dec3
       ra(i)=real(ra1)+real(ra2)/60.d0+ra3/3600.d0
       dec(i)=real(dec1)+real(dec3)/60.d0 
      enddo
50    nsts=i-1
      print *,'data=?'
      read(*,*) id
               
55    close(21)
      rad=3.14159265358979d0/180.0d0 ! transform from degree to radian
      lat=38.67332d0 ! latitude,longitude of observatory
      long=66.89641d0
      height=2593.d0

      shrate=24.d0/(23.d0+56.d0/60.d0+4.09d0/3600.d0)! rate of siderial day to helios day
      s(1)=23.05d0/60.d0!siderial time of 0h heilos day
      s(2)=27.81d0/60.d0
      s(3)=31.78d0/60.d0
      s(4)=35.71d0/60.d0
      s(5)=59.65d0/60.d0 
      nday=5

      dt=1.d0/1.d0 ! time resolution to calculate airmass
      print *,'time to start observing?'
      read(*,*) hour,mini
      startt=real(hour)+real(mini)/60.d0
      print *,'time to end observing?'
      read(*,*) hour,mini
      endt=real(hour)+real(mini)/60.d0
      do i=1,12
       time(i)=startt+real(i-1)*dt
        if(time(i).gt.endt) then
         goto 80
        endif
      enddo
80    nt=i-1

      write(*,'(a12,12f7.3)') 'starname',(time(i),i=1,nt)

      do j=1,nsts
       do i=1,nt
        pt=startt+real(i-1)*dt
        if(pt.gt.12.d0) then
         st=s(id)-(24.d0-pt)*shrate
         else
         st=s(ld)+(pt)*shrate
        endif
        if(st.lt.0.d0) then
         st=st+24.d0
         elseif(st.gt.24.d0) then
         st=st-24.d0
        endif
        if(pt.gt.endt) then
         goto 100
        endif
        cz=sin(rad*lat)*sin(rad*dec(j)) + cos(rad*lat)*cos(rad*dec(j))
     *  *cos(rad*15.d0*(ra(j)-st))
        sz=1.d0/cz
        X(i,j)=sz-0.0018167d0*(sz-1.d0)-0.002885d0*(sz-1.d0)**2.d0
     *  -0.0008083d0*(sz-1.d0)**3.d0
        if(X(i,j).gt.2.5d0.or.X(i,j).lt.0.d0) then
         X(i,j)=9.999999d0
        endif
c        write(*,'(a12,12f8.4)') starname(j),(X(i,j),i=1,nt)
       enddo
100    write(*,'(a12,12f7.3)') starname(j),(X(i,j),i=1,nt)
      enddo
 


      stop
      end

