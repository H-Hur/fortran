      parameter (ns=100000)
      integer n2mass,nstar,idx(ns),id2m(ns),nobx(ns,7),nobs(ns,7)
      integer nnp,mtotr(ns),trto2m(ns),nsame,nig,ncount,idtem
      real*8 drax(ns),ddecx(ns),vx(ns),ix(ns),cx(ns,5),cex(ns,5)
      real*8 vex(ns),iex(ns),jx(ns),xc(ns),yc(ns),vt(ns),it(ns)
      real*8 ct(ns,5),vet(ns),iet(ns),cet(ns,5),dra,ddex,dist
      real*8 pi,era,edec,ra2m(ns),dec2m(ns),fitrad(3)
      character namex(ns)*16,chx(ns)*61,name2m(ns)*16,H*1
      character nigname(ns)*16
      pi=3.14159265358979d0
      era=(10.d0 + 45.d0/60.d0 + 3.59d0 / 3600d0) * 15.d0
      edec= -59.d0 -41.d0/60.d0 -4.3d0/3600.d0
      open(21,file='2mass_25mbox.dat')
      open(22,file='tr16_302.2mx')
      open(23,file='tr16_302.tr')
      do i=1,ns
       read(21,'(t1,a1)') H
       if(H.eq.'*') then
        goto 10
       endif
      enddo
10    do i=1,ns
       read(22,'(t1,a1)') H
       if(H.eq.'*') then
        goto 20
       endif
      enddo
20    do i=1,ns
       read(23,'(t1,a1)') H
       if(H.eq.'*') then
        goto 30
       endif
      enddo
 
30    do i=1,ns
       read(21,'(t58,a1)',end=39) H
       backspace(21)
       if(H.ne.' ') then
        read(21,31) id2m(i),ra2m(i),dec2m(i),name2m(i)
       else   
        read(21,32) id2m(i),ra2m(i),dec2m(i),name2m(i)
       endif
      enddo
31    format(t1,i7,t34,f10.6,t46,f10.6,t58,a16)
32    format(t1,i7,t35,f10.6,t47,f10.6,t59,a16)
39    n2mass=i-1

      nig=0
      do i=1,ns
       jx(i)=0.d0
       write(namex(i),'(a16)') '                '
       write(chx(i),'(a61)') '                                         '
       read(22,'(t130,f6.3)',end=49) jx(i)
       backspace(22)
       read(22,*) idx(i),drax(i),ddecx(i),vx(i),ix(i),(cx(i,j),j=1,5)
     & ,vex(i),iex(i),(cex(i,j),j=1,5),(nobx(i,j),j=1,7)
       if(jx(i).ne.0.d0) then
        nig=nig+1
        backspace(22)
        read(22,'(t173,a16)') namex(i)
        write(nigname(i),'(a16)') namex(i)
        backspace(22)
        read(22,'(t128,a61)') chx(i)
       endif
      enddo
49    nstar=i-1

      do i=1,nstar
       read(23,*,end=59) xc(i),yc(i),vt(i),it(i),(ct(i,j),j=1,5),vet(i)
     & ,iet(i),(cet(i,j),j=1,5),(nobs(i,j),j=1,7)
       do j=1,nstar
        nsame=0
        do k=1,5
         if((ct(i,k).eq.cx(j,k)).and.(cet(i,k).eq.cex(j,k))) then
          nsame=nsame+1
         endif
        enddo
        do k=1,7
         if(nobs(i,k).eq.nobx(j,k)) then
          nsame=nsame+1
         endif
        enddo
        if(nsame.eq.12) then
         mtotr(i)=j!2mx to tr
         trto2m(j)=i!tr to 2m
         goto 55
        endif
       enddo
55    enddo
59    close(23)
      close(22)
      close(21)


      fitrad(1)=1.5d0
      fitrad(2)=2.d0
      fitrad(3)=2.5d0

      open(21,file='coo.1') 
      open(22,file='coo.2')
      open(23,file='nopart17.coo')
      open(24,file='nopart18.coo')
      open(25,file='nopart17id.dat')
      open(26,file='nopart18id.dat')
      open(27,file='doubleid.dat')
      write(23,'(t1,a23)') '#   xc      yc      idx'
      write(24,'(t1,a23)') '#   xc      yc      idx'
      write(21,'(t1,a86)') '#    xc2m      yc2m       xc       yc       
     & dx          dy      dist    idx   id2mass'
      write(22,'(t1,a86)') '#    xc        yc         xc2m     yc2m     
     & dx          dy      dist    id2mass idx'
      write(25,'(t1,a91)') '#   xc2m     yc2m  rad  dist dist(s)  xc  
     &    yc      dx       dy       idx id2mass            '    
      write(26,'(t1,a91)') '#   xc2m     yc2m  rad  dist dist(s)  xc  
     &    yc      dx       dy       idx id2mass            ' 
      write(27,'(t1,a60)')'# xc        yc      2mxid1  2mxid2   dist  2m
     &ass-name                                                         '
      nnp=0
      do i=1,nstar
       if(namex(i).eq.' '.and.vx(i).gt.1.d0.and.vx(i).le.17.d0) then
        write(23,102) xc(trto2m(i)),yc(trto2m(i)),idx(i)        
       elseif(namex(i).eq.' '.and.vx(i).gt.1.d0.and.vx(i).le.18.d0)then
        write(24,102) xc(trto2m(i)),yc(trto2m(i)),idx(i)    
       idtem=0
       endif
       do j=1,n2mass 
         ddec=(ddecx(i)/60.d0 + edec - dec2m(j))*3600.d0/0.6015d0
         dra=(drax(i)/60.d0/cos((ddecx(i)/60.d0+edec)*pi/180.d0)+era
     &  -ra2m(j))*cos((ddecx(i)/60.d0+edec)*pi/180.d0)*3600.d0/0.6015d0
         dist=sqrt(ddec**2.d0+dra**2.d0)
        if(namex(i).eq.' '.and.vx(i).gt.1.d0.and.vx(i).le.17.d0) then
         do k=1,3
          if(dist.le.fitrad(k)) then
           ncount=0
           do l=1,nig
            if(name2m(j).ne.nigname(l)) then
             ncount=ncount+1
            endif
           enddo
           if((ncount.eq.nig).and.(idtem.ne.id2m(j))) then
            idtem=id2m(j)
            write(25,70) xc(trto2m(i))-ddec,yc(trto2m(i))-dra,fitrad(k),
     &      dist,dist/0.6015d0,xc(trto2m(i)),yc(trto2m(i)),ddec,dra
     &      ,idx(i),id2m(j)
           endif
          endif
         enddo 
        elseif(namex(i).eq.' '.and.vx(i).gt.1.d0.and.vx(i).le.18.d0)then
         do k=1,3
          if(dist.le.fitrad(k)) then
           ncount=0
           do l=1,nig
            if(name2m(j).ne.nigname(l)) then
             ncount=ncount+1
            endif
           enddo
           if(ncount.eq.nig.and.(idtem.ne.id2m(j))) then
            idtem=id2m(j)
            write(26,70) xc(trto2m(i))-ddec,yc(trto2m(i))-dra,fitrad(k),
     &      dist,dist/0.6015d0,xc(trto2m(i)),yc(trto2m(i)),ddec,dra
     &      ,idx(i),id2m(j)
           endif
          endif
         enddo
        elseif(namex(i).eq.name2m(j)) then
         nnp=nnp+1
         write(21,101) xc(trto2m(i))-ddec,yc(trto2m(i))-dra,
     &   xc(trto2m(i)),yc(trto2m(i)),ddec,dra,dist,idx(i),id2m(j)
         write(22,101) xc(trto2m(i)),yc(trto2m(i)),xc(trto2m(i))-ddec,
     &   yc(trto2m(i))-dra,ddec,dra,dist,id2m(j),idx(i)
         goto 69
        endif
69     enddo
       do j=1,nstar
        if((namex(i).ne.' ').and.(namex(i).eq.namex(j)).and.(i.ne.j))
     &  then
         dist=sqrt((xc(trto2m(i))-xc(trto2m(j)))**2.d0+(yc(trto2m(i))
     &   -yc(trto2m(j)))**2.d0)
         write(27,71) xc(trto2m(i)),yc(trto2m(i)),idx(i),idx(j),dist,
     &   namex(i)
        endif
       enddo
      enddo
70    format(2f9.3,f4.1,2f6.2,4f9.3,2i7)
71    format(t1,2f9.3,2i8,f7.3,2x,a16)
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)
      close(26)
      close(27)
101   format(t1,7f10.3,2i7)
102   format(t1,2f10.3,i7)
      print *,nstar,n2mass,nnp




        
      stop
      end
