c     f95 -o ext_coo ext_coo.f ~/fortran/subroutine/caltools/mat.f
      parameter (nx=1000)
      real*8 raco(nx),deco(nx),cra,cded,x1,x2,y1,y2,ap
      integer nnn,norder
      character input*24,pho*24,pos*24,app*24,coo*24
      x1=0.d0
      x2=2048.d0
      y1=0.d0
      y2=4096.d0
      ap=18.6d0
      open(21,file='ext.input')
      read(21,*) norder
      read(21,*) cra
      read(21,*) cdec
      read(21,'(t1,a24)') input!standard coo, stad list with t1,a12,t13,*,*
      read(21,'(t1,a24)') app  !standard ap phot file *.ap
      read(21,'(t1,a24)') pho  !standard pho file
      read(21,'(t1,a24)') pos  !standard pos file
      read(21,'(t1,a24)') coo  !output coo file
      close(21)
      call read_cid(input)
      call read_ap(app)
      call read_std(pho,pos)
      call cid_star()
      call cal_1st(cra,cdec,nnn)
      call cal_2nd(cra,cdec,nnn)
      call cal_3rd(cra,cdec,nnn)
      if(norder.eq.1) call cal_1st(cra,cdec,nnn)
      if(norder.eq.2) call cal_2nd(cra,cdec,nnn)
      if(norder.eq.3) call cal_3rd(cra,cdec,nnn)
      call cal_coo(nnn,cra,cdec,x1,x2,y1,y2,ap,coo)


      stop
      end

      subroutine cal_coo(nnn,cra,cdec,x1,x2,y1,y2,ap,coo)
      parameter (ns=100000,nx=1000)
      character input(2)*24,name_std(ns)*12,star_cid(ns)*12,line*300
      character coo*24
      common/std/name_std,sv(ns),sub(ns),sbv(ns),svi(ns),sri(ns),
     &esv(ns),esub(ns),esbv(ns),esvi(ns),esri(ns),
     &nsv(ns),nsub(ns),nsbv(ns),nsvi(ns),nsri(ns), ! n_obs
     &msv(ns),msub(ns),msbv(ns),msvi(ns),msri(ns), ! nights
     &svar(ns),sra,sdec,nstd,input
      common/cid/star_cid,x_in(ns),y_in(ns),xc(ns),yc(ns),am(ns),ae(ns),
     &v(ns),bv(ns),vi(ns),ri(ns),ncid
      real*8 cra,cdec,xco(nx),yco(nx),sra(ns),sdec(ns),r,d,xx,yy
      real*8 pi,sc,x(ns),y(ns),ra(ns),dec(ns),xcal(ns),ycal(ns)
      real*8 x1,x2,y1,y2,ap
      integer nnn
      pi=3.14159265358979d0
      sc=1.d0
c      nnn=10
      do i=1,ncid
       x(i)=xc(i)
       y(i)=yc(i)
       do j=1,nstd
        if(name_std(j).eq.star_cid(i)) then
         ra(i)=(sra(j)-cra)*cos(sdec(j)*pi/180.d0)*60.d0
         dec(i)=(sdec(j)-cdec)*60.d0
        endif
       enddo
      enddo
      call make_3xy(nnn,ra,dec,x,ncid,xco,sc)
      call make_3xy(nnn,ra,dec,y,ncid,yco,sc)
      k=nnn+1
      if(k.lt.10) then
       do i=k,10
        xco(i)=0.d0
        yco(i)=0.d0
       enddo
      endif
      open(34,file='3rd_transform.coo')
      write(line,*) (xco(n),n=1,nnn)
      write(34,'(t1,a4,a300)') 'RA  ',line
      write(line,*) (yco(n),n=1,nnn)
      write(34,'(t1,a4,a300)') 'Dec ',line
      write(34,'(t1,a30)') '#t1,4f10.3,a13,2f18.12                     '
      write(34,'(t1,a135)') '#Pixel                                     
     &                Name           degree                           
     &                                                                 '
      write(34,'(t1,a135)') '# Xcal       Ycal       Xc        Yc       
     &  dx        dy    Star            R.A      Dec.                   
     &                                                                 '
      print *,'Running 3rd transform to extract pixel coordinates'
      do i=1,ncid
       ra(i)=ra(i)/sc
       dec(i)=dec(i)/sc
       xcal(i)=xco(1) + xco(2) * ra(i) + xco(3) * dec(i) +
     & xco(4) * ra(i)**2 + xco(5) * ra(i)*dec(i) + xco(6) * dec(i)**2
     & + xco(7) *( ra(i)**3 )+ xco(8) * (ra(i)**2)*dec(i)
     & + xco(9) * ra(i)*(dec(i)**2) + xco(10)*(dec(i)**3)
       ycal(i)=yco(1) + yco(2) * ra(i) + yco(3) * dec(i) +
     & yco(4) * ra(i)**2 + yco(5) * ra(i)*dec(i) + yco(6) * (dec(i)**2)
     & + yco(7) * (ra(i)**3) + yco(8) * (ra(i)**2)*dec(i)
     & + yco(9) * ra(i)*(dec(i)**2) + yco(10)*(dec(i)**3)
       write(34,'(t1,6f10.3,a13,2f18.12)')
     & xcal(i),ycal(i),xc(i),yc(i),xcal(i)-xc(i),ycal(i)-yc(i),
     & star_cid(i),ra(i),dec(i)
      enddo
      close(34)
      open(35,file=coo)
      write(35,'(a96)')'#   Xc        Yc     Star              V       U
     &-B     B-V     V-I    R-I    nV  nUB nBV nVI nRI                 '
      do i=1,nstd
       r=(sra(i)-cra)*cos(sdec(i)*pi/180.d0)*60.d0
       d=(sdec(i)-cdec)*60.d0
       xx=xco(1) + xco(2) * r + xco(3) * d +
     & xco(4) * r**2 + xco(5) * r*d + xco(6) * d**2
     & + xco(7) *( r**3 )+ xco(8) * (r**2)*d
     & + xco(9) * r*(d**2) + xco(10)*(d**3)
       yy=yco(1) + yco(2) * r + yco(3) * d +
     & yco(4) * r**2 + yco(5) * r*d + yco(6) * (d**2)
     & + yco(7) * (r**3) + yco(8) * (r**2)*d
     & + yco(9) * r*(d**2) + yco(10)*(d**3)
       if(xx.gt.x1-ap.and.xx.lt.x2+ap.and.
     & yy.gt.y1-ap.and.yy.lt.y2+ap) 
     & write(35,'(2f10.3,1x,a12,3x,5f8.3,5i4)') xx,yy,name_std(i),
     & sv(i),sub(i),sbv(i),svi(i),sri(i),
     & nsv(i),nsub(i),nsbv(i),nsvi(i),nsri(i) 
      enddo
      close(35)
      return
      end

      subroutine cal_3rd(cra,cdec,nnn)
      parameter (ns=100000,nx=1000)
      character input(2)*24,name_std(ns)*12,star_cid(ns)*12,line*300
      common/std/name_std,sv(ns),sub(ns),sbv(ns),svi(ns),sri(ns),
     &esv(ns),esub(ns),esbv(ns),esvi(ns),esri(ns),
     &nsv(ns),nsub(ns),nsbv(ns),nsvi(ns),nsri(ns), ! n_obs
     &msv(ns),msub(ns),msbv(ns),msvi(ns),msri(ns), ! nights
     &svar(ns),sra,sdec,nstd,input
      common/ap/id_ap(ns),x_ap(ns),y_ap(ns),amag(ns),amerr(ns),sky(ns)
     &,air,apert,exptime,nap
      common/cid/star_cid,x_in(ns),y_in(ns),xc(ns),yc(ns),am(ns),ae(ns),
     &v(ns),bv(ns),vi(ns),ri(ns),ncid
      real*8 rac(ns),decc(ns),cra,cdec,pi,dra(ns),ddec(ns),x(ns),y(ns)
      real*8 sra(ns),sdec(ns),ra(ns),dec(ns),sc,raco(nx),deco(nx)
      integer nnn
      nnn=10
      sc=1.d0
      pi=3.14159265358979d0
      do i=1,ncid
       x(i)=xc(i)
       y(i)=yc(i)
       do j=1,nstd
        if(name_std(j).eq.star_cid(i)) then
         dra(i)=(sra(j)-cra)*cos(sdec(j)*pi/180.d0)*60.d0
         ddec(i)=(sdec(j)-cdec)*60.d0
         ra(i)=sra(j)
         dec(i)=sdec(j)
        endif
       enddo
      enddo
      call make_3xy(nnn,x,y,ra,ncid,raco,sc)
      call make_3xy(nnn,x,y,dec,ncid,deco,sc)
      open(30,file='3rd_transform.ext')
      write(30,'(t1,a4,a24)') 'pho ',input(1)
      write(30,'(t1,a4,a24)') 'pos ',input(2)
      write(line,*) (raco(n),n=1,nnn)
      write(30,'(t1,a4,a300)') 'RA  ',line
      write(line,*) (deco(n),n=1,nnn)
      write(30,'(t1,a4,a300)') 'Dec ',line
      write(30,'(t1,a30)') '#t1,a12,2f9.3,4f18.12,4f9.3                '
      write(30,'(t1,a135)') '#Name        pixel               degree    
     &                                                            arcmin
     &        arcsec                                                   '
      write(30,'(t1,a135)') '#Star          Xc       Yc       R.A       
     &        Dec.            R.A (calculated)  Dec.(Calculated) del_ra 
     &del_dec ra_res dec_res                                           '
      print *,'Running 3rd transform...'
      do i=1,ncid
       rac(i)=raco(1) + raco(2) * xc(i) + raco(3) * yc(i) +
     & raco(4) * xc(i)**2 + raco(5) * xc(i)*yc(i) + raco(6) * yc(i)**2
     & + raco(7) *( xc(i)**3 )+ raco(8) * (xc(i)**2)*yc(i)
     & + raco(9) * xc(i)*(yc(i)**2) + raco(10)*(yc(i)**3)
       decc(i)=deco(1) + deco(2) * xc(i) + deco(3) * yc(i) +
     & deco(4) * xc(i)**2 + deco(5) * xc(i)*yc(i) + deco(6) * (yc(i)**2)
     & + deco(7) * (xc(i)**3) + deco(8) * (xc(i)**2)*yc(i)
     & + deco(9) * xc(i)*(yc(i)**2) + deco(10)*(yc(i)**3)
       write(30,'(t1,a12,2f9.3,4f18.12,4f11.3)') star_cid(i),xc(i),yc(i)
     & ,ra(i),dec(i),rac(i),decc(i),dra(i),ddec(i)
     & ,(ra(i)-rac(i))*cos(decc(i)*pi/180.d0)*60.d0*60.d0
     & ,(dec(i)-decc(i))*60.d0*60.d0
      enddo
      close(30)
      return
      end
 
      subroutine cal_2nd(cra,cdec,nnn)
      parameter (ns=100000,nx=1000)
      character input(2)*24,name_std(ns)*12,star_cid(ns)*12,line*300
      common/std/name_std,sv(ns),sub(ns),sbv(ns),svi(ns),sri(ns),
     &esv(ns),esub(ns),esbv(ns),esvi(ns),esri(ns),
     &nsv(ns),nsub(ns),nsbv(ns),nsvi(ns),nsri(ns), ! n_obs
     &msv(ns),msub(ns),msbv(ns),msvi(ns),msri(ns), ! nights
     &svar(ns),sra,sdec,nstd,input
      common/ap/id_ap(ns),x_ap(ns),y_ap(ns),amag(ns),amerr(ns),sky(ns)
     &,air,apert,exptime,nap
      common/cid/star_cid,x_in(ns),y_in(ns),xc(ns),yc(ns),am(ns),ae(ns),
     &v(ns),bv(ns),vi(ns),ri(ns),ncid
      real*8 rac(ns),decc(ns),cra,cdec,pi,dra(ns),ddec(ns),x(ns),y(ns)
      real*8 sra(ns),sdec(ns),ra(ns),dec(ns),sc,raco(nx),deco(nx)
      integer nnn
      nnn=6
      sc=1.d0
      pi=3.14159265358979d0
      do i=1,ncid
       x(i)=xc(i)
       y(i)=yc(i)
       do j=1,nstd
        if(name_std(j).eq.star_cid(i)) then
         dra(i)=(sra(j)-cra)*cos(sdec(j)*pi/180.d0)*60.d0
         ddec(i)=(sdec(j)-cdec)*60.d0
         ra(i)=sra(j)
         dec(i)=sdec(j)
        endif
       enddo
      enddo
      call make_3xy(nnn,x,y,ra,ncid,raco,sc)
      call make_3xy(nnn,x,y,dec,ncid,deco,sc)
      open(31,file='2nd_transform.ext')
      write(31,'(t1,a4,a24)') 'pho ',input(1)
      write(31,'(t1,a4,a24)') 'pos ',input(2)
      write(line,*) (raco(n),n=1,nnn)
      write(31,'(t1,a4,a300)') 'RA  ',line
      write(line,*) (deco(n),n=1,nnn)
      write(31,'(t1,a4,a300)') 'Dec ',line
      write(31,'(t1,a30)') '#t1,a12,2f9.3,4f18.12,4f9.3                '
      write(31,'(t1,a135)') '#Name        pixel               degree    
     &                                                            arcmin
     &        arcsec                                                   '
      write(31,'(t1,a135)') '#Star          Xc       Yc       R.A       
     &        Dec.            R.A (calculated)  Dec.(Calculated) del_ra 
     &del_dec ra_res dec_res                                           '
      print *,'Running 2nd transform...'
      do i=1,ncid
       rac(i)=raco(1) + raco(2) * xc(i) + raco(3) * yc(i) +
     & raco(4) * xc(i)**2 + raco(5) * xc(i)*yc(i) + raco(6) * yc(i)**2
       decc(i)=deco(1) + deco(2) * xc(i) + deco(3) * yc(i) +
     & deco(4) * xc(i)**2 + deco(5) * xc(i)*yc(i) + deco(6) * (yc(i)**2)
       write(31,'(t1,a12,2f9.3,4f18.12,4f11.3)') star_cid(i),xc(i),yc(i)
     & ,ra(i),dec(i),rac(i),decc(i),dra(i),ddec(i)
     & ,(ra(i)-rac(i))*cos(decc(i)*pi/180.d0)*60.d0*60.d0
     & ,(dec(i)-decc(i))*60.d0*60.d0
      enddo
      close(31)
      return
      end



      subroutine cal_1st(cra,cdec,nnn)
      parameter (ns=100000,nx=1000)
      character input(2)*24,name_std(ns)*12,star_cid(ns)*12,line*300
      common/std/name_std,sv(ns),sub(ns),sbv(ns),svi(ns),sri(ns),
     &esv(ns),esub(ns),esbv(ns),esvi(ns),esri(ns),
     &nsv(ns),nsub(ns),nsbv(ns),nsvi(ns),nsri(ns), ! n_obs
     &msv(ns),msub(ns),msbv(ns),msvi(ns),msri(ns), ! nights
     &svar(ns),sra,sdec,nstd,input
      common/ap/id_ap(ns),x_ap(ns),y_ap(ns),amag(ns),amerr(ns),sky(ns)
     &,air,apert,exptime,nap
      common/cid/star_cid,x_in(ns),y_in(ns),xc(ns),yc(ns),am(ns),ae(ns),
     &v(ns),bv(ns),vi(ns),ri(ns),ncid
      real*8 rac(ns),decc(ns),cra,cdec,pi,dra(ns),ddec(ns),x(ns),y(ns)
      real*8 sra(ns),sdec(ns),ra(ns),dec(ns),sc,raco(nx),deco(nx)
      integer nnn
      nnn=3
      print *,'Running 1st transform...'
      sc=1.d0
      pi=3.14159265358979d0
      do i=1,ncid
       x(i)=xc(i)
       y(i)=yc(i)
       do j=1,nstd
        if(name_std(j).eq.star_cid(i)) then
         dra(i)=(sra(j)-cra)*cos(sdec(j)*pi/180.d0)*60.d0
         ddec(i)=(sdec(j)-cdec)*60.d0
         ra(i)=sra(j)
         dec(i)=sdec(j)
        endif
       enddo
      enddo
      call make_3xy(nnn,x,y,ra,ncid,raco,sc)
      call make_3xy(nnn,x,y,dec,ncid,deco,sc)
      open(32,file='1st_transform.ext')
      write(32,'(t1,a4,a24)') 'pho ',input(1)
      write(32,'(t1,a4,a24)') 'pos ',input(2)
      write(line,*) (raco(n),n=1,nnn)
      write(32,'(t1,a4,a300)') 'RA  ',line
      write(line,*) (deco(n),n=1,nnn)
      write(32,'(t1,a4,a300)') 'Dec ',line
      write(32,'(t1,a30)') '#t1,a12,2f9.3,4f18.12,4f9.3                '
      write(32,'(t1,a135)') '#Name        pixel               degree    
     &                                                            arcmin
     &        arcsec                                                   '
      write(32,'(t1,a135)') '#Star          Xc       Yc       R.A       
     &        Dec.            R.A (calculated)  Dec.(Calculated) del_ra 
     &del_dec ra_res dec_res                                           '
      do i=1,ncid
       rac(i)=raco(1) + raco(2) * xc(i) + raco(3) * yc(i) 
       decc(i)=deco(1) + deco(2) * xc(i) + deco(3) * yc(i)  
       write(32,'(t1,a12,2f9.3,4f18.12,4f11.3)') star_cid(i),xc(i),yc(i)
     & ,ra(i),dec(i),rac(i),decc(i),dra(i),ddec(i)
     & ,(ra(i)-rac(i))*cos(decc(i)*pi/180.d0)*60.d0*60.d0
     & ,(dec(i)-decc(i))*60.d0*60.d0
      enddo
      close(32)
      return
      end

      subroutine cid_star()
      parameter (ns=100000)
      character input(2)*24,name_std(ns)*12,star_cid(ns)*12
      common/std/name_std,sv(ns),sub(ns),sbv(ns),svi(ns),sri(ns),
     &esv(ns),esub(ns),esbv(ns),esvi(ns),esri(ns),
     &nsv(ns),nsub(ns),nsbv(ns),nsvi(ns),nsri(ns), ! n_obs
     &msv(ns),msub(ns),msbv(ns),msvi(ns),msri(ns), ! nights
     &svar(ns),sra,sdec,nstd,input
      common/ap/id_ap(ns),x_ap(ns),y_ap(ns),amag(ns),amerr(ns),sky(ns)
     &,air,apert,exptime,nap
      common/cid/star_cid,x_in(ns),y_in(ns),xc(ns),yc(ns),am(ns),ae(ns),
     &v(ns),bv(ns),vi(ns),ri(ns),ncid
      real*8 sra(ns),sdec(ns),dist,dist_least,fitrad
      integer n_least
      fitrad=5.d0
      do i=1,ncid
       dist_least=fitrad
       n_least=0
       xc(i)=-99.999
       yc(i)=-99.999
       am(i)=-99.999
       v(i)=-99.999
       bv(i)=-99.999
       vi(i)=-99.999
       ri(i)=-99.999
       do j=1,nap
        dist=sqrt((x_ap(j)-x_in(i))**2+(y_ap(j)-y_in(i))**2)
        if(dist.lt.dist_least) then
         dist_least=dist
         n_least=j
        endif
        if(n_least.ne.0) then
         xc(i)=x_ap(n_least)
         yc(i)=y_ap(n_least)
         am(i)=-amag(n_least)
        endif
       enddo
200   i200=i200
      enddo
      return
      end
      subroutine read_std(pho,pos)
      parameter (ns=100000)
      character pho*24,pos*24,name_std(ns)*12,line*150,name2(ns)*12
      character input(2)*24
      common/std/name_std,sv(ns),sub(ns),sbv(ns),svi(ns),sri(ns),
     &esv(ns),esub(ns),esbv(ns),esvi(ns),esri(ns),
     &nsv(ns),nsub(ns),nsbv(ns),nsvi(ns),nsri(ns), ! n_obs
     &msv(ns),msub(ns),msbv(ns),msvi(ns),msri(ns), ! nights
     &svar(ns),sra,sdec,nstd,input
      real*8 sra(ns),sdec(ns),ra(ns),dec(ns)
      real b,eb,r,er,ai,ei,u,eu
      integer nb,nr,ni,mb,mr,mi,npos,nu,mu
      open(23,file=pho)
      open(24,file=pos)
      write(input(1),'(a24)') pho
      write(input(2),'(a24)') pos
      do i=1,ns
       read(24,'(t1,2f16.11,t101,a12)',err=49,end=50) 
     & ra(i),dec(i),name2(i) 
      enddo
49    print *,'Error in reading std position file.'
      goto 900
50    close(24)
      npos=i-1
      read(23,'(t1,a1)') H
      do i=1,ns
       read(23,'(t2,a12,t15,a150)',err=99,end=100) name_std(i),line
       read(line,*) u,eu,nu,mu,b,eb,nb,mb,sv(i),esv(i),nsv(i),msv(i),
     & r,er,nr,mr,ai,ei,ni,mi,svar(i)
       sub(i)=u-b
       nsub(i)=min(nu,nb)
       msub(i)=min(mu,nb)
       call cal_error(eu,eb,esub(i))
       sbv(i)=b-sv(i)
       nsbv(i)=min(nsv(i),nb)
       msbv(i)=min(msv(i),nb)
       call cal_error(esv(i),eb,esbv(i))
       svi(i)=sv(i)-ai
       nsvi(i)=min(nsv(i),ni)
       msvi(i)=min(msv(i),ni)
       call cal_error(esv(i),ei,esvi(i))
       sri(i)=r-ai
       nsri(i)=min(nr,ni)
       msri(i)=min(mr,ni)
       call cal_error(er,ei,esri(i))
       do j=1,npos
        if(name_std(i).eq.name2(j)) then 
         sra(i)=ra(j)
         sdec(i)=dec(j)
         goto 90
        endif
       enddo
       print *,i,'st star has no EQ coordinate informaiton!'
90    i90=i90
      enddo
99    print *,'Error in reading std photometry file.'
      goto 900
100   nstd=i-1
      close(23)
900   return
      end

      subroutine read_ap(input)
      parameter (ns=100000)
      character input*24,H*1
      common/ap/id_ap(ns),x_ap(ns),y_ap(ns),amag(ns),amerr(ns),sky(ns)
     &,air,apert,exptime,nap
      open(22,file=input)
      do i=1,ns
       read(22,'(t1,a1)',err=9) H
       if(H.ne.'#') goto 10
      enddo
      print *,'No star read from ',input
      goto 900
9     print *,'Errror!! No ap phot stars are found!!'
      goto 900
10    backspace(22)
      do i=1,ns
       read(22,'(t43,i7)',err=99,end=100) id_ap(i)
       read(22,*,err=99,end=100) x_ap(i),y_ap(i)
       read(22,*,err=99,end=100) sky(i)
       read(22,*,err=99,end=100) exptime,air
       read(22,*,err=99,end=100) apert
      enddo
99    print *,'Errror in reading ap file!!'
      goto 900
100   nap=i-1
      close(22)
      print *,nap,' ap stars read.'
900   return
      end
      subroutine read_cid(input)
      parameter (ns=100000)
      character input*24,H*1,star_cid(ns)*12,line*30
      common/cid/star_cid,x_in(ns),y_in(ns),xc(ns),yc(ns),am(ns),ae(ns),
     &v(ns),bv(ns),vi(ns),ri(ns),ncid
      open(21,file=input)
      do i=1,ns
       read(21,'(t1,a1)',err=9) H
       if(H.ne.'#') goto 10
      enddo
      print *,'No star read from ',input
      goto 900
9     print *,'Errror!! No input initial cross-id stars are found!!'
      goto 900
10    backspace(21)
      ncid=0
      do i=1,ns
       read(21,'(t1,a1)',end=100) H
       if(H.ne.'#') then
        backspace(21)
        ncid=ncid+1
        read(21,'(T1,A12,T13,a30)',err=99,end=100) 
     &  star_cid(ncid),line
        read(line,*) x_in(ncid),y_in(ncid)
       endif
      enddo
99    print *,'Errror in reading initial cross-id stars!!'
      goto 900
100   close(21)
      print *,ncid,' initial cross-id stars read.'
900   return
      end
      subroutine cal_error(e1,e2,e3)
      real e1,e2,e3
      if(e1.lt.0.d0.or.e2.lt.0.) then
       e3=-1.
      else
       e3=sqrt(e1**2.+e2**2.)
      endif
      return
      end
