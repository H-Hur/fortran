c for only fortran77.
c f77 -o cor_red.UBV cor_red.UBV.f -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11
      parameter (nf=10,ns=1000)
      character filter(nf)*2,star(nf,ns)*16,input*24,coeff*16,output*16
      real*8 transy(nf,ns),airy(nf,ns),timey(nf,ns),colx(nf,ns)
      real*8 timex(nf,ns)
      integer num_f(nf),num_c(nf),ntransform,lastnum_f,ne,nne
      write(input,'(a24)') 'Sung2001.ZAMS.cc        '
      write(output,'(a16)') 'out.cor_red.UBV '
      
      call readhead(input)
      call read_ZAMS(input)
 
c      write(input,'(a24)') 'point.cc                '
c      call readhead(input)
      call cor_red()


      stop
      end
   

      subroutine cor_red()
      parameter (nf=10,ns=1000,maxiter=1000)
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),zVR(ns),zRI(ns),nzams
      integer niter(ns),nunred
      real fx,fy,dy,unredx,unredy,sBV,sUB,slope,eBV,eUB,eVI,Rv,Av,Ab,Au
      real c_bv(21),c_u(21),f
      data c_bv/-.35d0,-.33d0,-0.3d0,-0.2d0,-0.1d0, 0.0d0, 0.1d0, 0.2d0,
     &0.3d0, 0.4d0,0.5d0,0.6d0,0.7d0,0.8d0,0.9d0,1.0d0,1.2d0,1.4d0,1.6d0
     &,1.8d0,2.0d0/
      data c_u/-.132d0,-.132d0,-.126d0,-.074d0,-.030d0,0.002d0,0.024d0,
     &0.026d0,0.010d0,-.017d0,-.06d0,-.104d0,-.132d0,-.132d0,-.132d0,
     &-.132d0,-.132d0,-.132d0,-.132d0,-.132d0,-.132d0/

      print *,'B-V,U-B=?'
      read(*,*) sBV,sUB
c      print *,'Slope=?'
c      read(*,*) slope
c      sBV=-0.271
c      sUB=-1.081
      slope=0.72
      Rv=3.1
      nunred=0
      do i=1,1  
       call lint(sUB,zUB,zBV,fx,ns,ns)
       dx=fx-sBV
       call lint(sBV,zBV,zUB,fy,ns,ns)
       dy=sUB-fy
c       print *,dx,dy
c       if(dx.gt.0.05.and.dy.gt.0.05) then
        nunred=nunred+1
        niter(i)=0  
        unredx=sBV
        unredy=sUB
        do j=1,1000
         niter(i)=niter(i)+1
         call lint(unredy,zUB,zBV,fx,ns,ns) 
         dx=fx-unredx
         dy=dx*slope
         unredx=unredx+dx
         unredy=unredy+dy
         if(abs(dx).lt.0.00001.and.abs(dy).lt.0.00001) then
          goto 100
         endif
        enddo
100     eBV=sBV-unredx
        eUB=sUB-unredy
        Av=Rv*eBV
        Ab=eBV+Av
        call lint(unredx,c_bv,c_u,f,21,21)

c       endif
c       write(*,*) i,sBV,sUB,unredx,unredy,eBV,eUB,Av,Ab
        write(*,'(a5,f7.4)') 'U-B= ',sUB
        write(*,'(a5,f7.4)') 'B-V= ',sBV
        write(*,'(a7,f7.4)') 'e(U-B)=  ',eUB
        write(*,'(a7,f7.4)') 'e(B-V)=  ',eBV
        write(*,'(a13,f7.4)') 'unreden U-B= ',unredy
        write(*,'(a13,f7.4)') 'unrened B-V= ',unredx
        write(*,'(a7,f7.4)') 'slope= ',slope
        write(*,'(a6,i4)') 'niter=  ',niter(i)
        write(*,'(a11,f7.4)')'f[(B-V)0]= ',f
        write(*,*) 'Last dx,dy=  ',dx,dy
       
       enddo

      return
      end    


   


      subroutine read_ZAMS(input)
      parameter (nf=10,ns=1000,nmcol=20)
      character line*200,H(200)*1,input*24,coll(nmcol)*32
      character Hv(32)*1,lastcoll1*32
      integer numcol,nc1,nc2,ncol,nc10,nc20
      common /nc/nc_filter,nc_star,nc_mag,nc_merr,nc_UT,nc_air,nc_weight
     &,ncollum,nc_V,nc_BV,nc_VR,nc_RI,nc_VI,nc_UB,nc_Mv
      common /zams/zMv(ns),zVI(ns),zBV(ns),zUB(ns),zVR(ns),zRI(ns),nzams
      do i=1,32
       write(Hv(i),'(a1)') ' '
      enddo
      open(31,file=input)
      write(lastcoll1,'(32a1)') (Hv(k),k=1,32)
      nzams=0
      do i=1,ns
       nzams=nzams+1
       read(31,'(a200)',end=500) line
       read(line,'(200a1)') (H(j),j=1,200)
       if(H(1).eq.'#') then
        nzams=nzams-1
        goto 100
       endif
       do j=1,16
        if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &  '^')then
         nc2=j
         goto 10
        endif
       enddo

10     do numcol=1,ncollum
        if(numcol.eq.1) then
         nc1=1
         do j=1,nc2+1+32
          if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &    '^')  then
           nc2=j
           goto 20
          endif
         enddo
        else
         do j=nc2+1,nc2+1+32
          if(H(j).eq.' '.and.H(j+1).ne.' ') then
           nc1=j+1
          endif
          if(H(j).ne.' '.and.H(j+1).eq.' '.or.H(j).ne.' '.and.H(j+1).eq.
     &    '^')  then
           nc2=j
           goto 20
          endif
         enddo
        endif
20      write(coll(numcol),'(32a1)') (Hv(k),k=1,32)
        write(coll(numcol),'(32a1)') (H(k),k=nc1,nc2)
       enddo

       read(coll(nc_Mv),*) zMv(nzams)
       read(coll(nc_BV),*) zBV(nzams)
c       read(coll(nc_VR),*) zVR(nzams)
c       read(coll(nc_RI),*) zRI(nzams)
       read(coll(nc_VI),*) zVI(nzams)
       read(coll(nc_UB),*) zUB(nzams)

c       print *,zMv(nzams),zVI(nzams),zBV(nzams),zUB(nzams)
100   enddo
500   close(31)



      return
      end



      subroutine readhead(input)
      parameter (nf=10,ns=1000)
      character header*200,H(200)*1,input*24,coll*3
      common /nc/nc_filter,nc_star,nc_mag,nc_merr,nc_UT,nc_air,nc_weight
     &,ncollum,nc_V,nc_BV,nc_VR,nc_RI,nc_VI,nc_UB,nc_Mv
      open(31,file=input)
      read(31,'(a200)') header
      ncollum=0
      read(header,'(200a1)') (H(i),i=1,200)

      do i=1,200
       if(H(i).eq.' '.or.H(i).eq.'#') then
        if(H(i+1).ne.' ') then
         write(coll,'(3a1)') H(i+1),H(i+2),H(i+3)
         if(coll.eq.'Sta'.or.coll.eq.'sta') then
          ncollum=ncollum+1
          nc_star=ncollum
          elseif(coll.eq.'Mag'.or.coll.eq.'mag') then
          ncollum=ncollum+1
          nc_mag=ncollum
          elseif(coll.eq.'Mer'.or.coll.eq.'mer') then
          ncollum=ncollum+1
          nc_merr=ncollum
          elseif(coll.eq.'V  '.or.coll.eq.'v  ') then
          ncollum=ncollum+1
          nc_V=ncollum
          elseif(coll.eq.'B-V'.or.coll.eq.'b-v') then
          ncollum=ncollum+1
          nc_BV=ncollum
          elseif(coll.eq.'U-B'.or.coll.eq.'u-b') then
          ncollum=ncollum+1
          nc_UB=ncollum
          elseif(coll.eq.'V-R'.or.coll.eq.'v-r') then
          ncollum=ncollum+1
          nc_VR=ncollum
          elseif(coll.eq.'R-I'.or.coll.eq.'r-i') then
          ncollum=ncollum+1
          nc_RI=ncollum
          elseif(coll.eq.'Mv') then
          ncollum=ncollum+1
          nc_Mv=ncollum
          elseif(coll.eq.'V-I'.or.coll.eq.'v-i') then
          ncollum=ncollum+1
          nc_VI=ncollum
          elseif(coll.eq.'UT '.or.coll.eq.'ut '.or.coll.eq.'UT('.or.
     &    coll.eq.'ut(') then
          ncollum=ncollum+1
          nc_ut=ncollum
          elseif(coll.eq.'Air'.or.coll.eq.'air'.or.coll.eq.'X  '.or.
     &    coll.eq.'x  ') then
          ncollum=ncollum+1
          nc_air=ncollum
          elseif(coll.eq.'wei'.or.coll.eq.'Wei') then
          ncollum=ncollum+1
          nc_weight=ncollum
          elseif(coll.eq.'fil'.or.coll.eq.'Fil') then
          ncollum=ncollum+1
          nc_filter=ncollum
         endif
        endif
       endif
      enddo
      close(31)
      return
      end


      subroutine lint(x,xt,yt,y,n,m)
      parameter(nmax=10000)
      real xt(n),yt(n),xt1(nmax),yt1(nmax)
      real x,y
      integer i,k,m
      if(xt(m).lt.xt(1)) then
       do 100 i=1,m
        xt1(i)=xt(m+1-i)
        yt1(i)=yt(m+1-i)
100    continue
      else
       do 200 i=1,m
        xt1(i)=xt(i)
        yt1(i)=yt(i)
200    continue
      endif
      do 10 i=2,m
       k=i-1
       if(x-xt1(i).lt.0.d0) then
        goto 11
        elseif(x-xt1(i).eq.0.d0) then
        goto 12
        elseif (x-xt1(i).gt.0.d0) then
        goto 10
       endif
12     y=yt1(i)
       goto 14
10    continue
11    y=yt1(k)+(yt1(k+1)-yt1(k))*((x-xt1(k))/(xt1(k+1)-xt1(k)))
14    return
      end


