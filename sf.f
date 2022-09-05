      real*8 ufi(27),bfi(27),vfi(27),rfi(27),ifi(27),Qwav(10)
      real*8 uw(27),bw(27),vw(27),rw(27),iw(27),bxfi(27),cc(6)
      real*8 wav(100),flux(100),pflux,test,w,R(10),Q(10),Rmir,Tfil,Qccd
      real*8 Stot(6),lamda0(6),moment(6),lameff(6),SS(6),SES(6),col(5)
      integer lamda,iskip,nwav
      character input*8,H*1

      open(21,file='u.sf')
      open(22,file='b.sf')
      open(23,file='v.sf')
      open(24,file='r.sf')
      open(25,file='i.sf')
      do i=1,27
       read(21,*) uw(i),ufi(i)
       ufi(i)=ufi(i)/1000.d0
       read(22,*) bw(i),bfi(i),bxfi(i)
       bfi(i)=bfi(i)/1000.d0
       bxfi(i)=bxfi(i)/1000.d0
       read(23,*) vw(i),vfi(i)
       vfi(i)=vfi(i)/1000.d0
       read(24,*) rw(i),rfi(i)
       rfi(i)=rfi(i)/1000.d0
       read(25,*) iw(i),ifi(i)
       ifi(i)=ifi(i)/1000.d0
      enddo
      close(21)
      close(22)
      close(23)
      close(24)
      close(25)

      open(21,file='q.sf') 
      do i=1,10
       read(21,*) Qwav(i),Q(i),R(i)
       Qwav(i)=Qwav(i)*1000.d0
       R(i)=R(i)/100.d0
       Q(i)=Q(i)/100.d0
      enddo 
      close(21)
 
      do i=1,100
       wav(i)=0.d0
       flux(i)=0.d0
      enddo
      
      do k=1,100
      print *,'input file name=? (0=end)'
      read(*,'(a8)') input
      write(H,'(a1)') input
      if(H.eq.'0') then
       goto 900
      endif

      open(21,file=input) ! 별 고르기
      do i=1,100
       read(21,*,end=10) iskip,iskip,wav(i),flux(i)
      enddo 
10    nwav=i-1

      do i=1,6
       lamda0(i)=0.d0
       moment(i)=0.d0
       lameff(i)=0.d0
      enddo

      do j=1,6  ! J do loof - 필터고르기
       SS(j)=0.d0
       SES(j)=0.d0
       cc(j)=0.d0
       do i=300,1200,1
        w=real(i)
        call lint(w,Qwav,R,Rmir,10,10) !미러의반응율
        call lint(w,Qwav,Q,Qccd,10,10) !양자효율
        if(j.eq.1) then                !필터별 투과함수
         call lint(w,uw,ufi,Tfil,27,27)
         elseif(j.eq.2) then
         call lint(w,bw,bfi,Tfil,27,27)
         elseif(j.eq.3) then
         call lint(w,vw,vfi,Tfil,27,27)
         elseif(j.eq.4) then
         call lint(w,rw,rfi,Tfil,27,27)
         elseif(j.eq.5) then
         call lint(w,iw,ifi,Tfil,27,27)
         elseif(j.eq.6) then
         call lint(w,bw,bxfi,Tfil,27,27)
        endif
        call lint(w,wav,flux,pflux,nwav,nwav)
        Stot(j)=Tfil*Qccd*Rmir**2.d0
        SS(j)=SS(j)+Stot(j)
        SES(j)=SES(j)+pflux*Stot(j)
        lamda0(j)=lamda0(j)+w*Stot(j)       
        lameff(j)=lameff(j)+w*pflux*Stot(j)     
        cc(j)=cc(j)+lameff(j)*Stot(j)
      enddo
      enddo
 
      do j=1,6
       lamda0(j)=lamda0(j)/SS(j)
       lameff(j)=lameff(j)/SES(j)
      enddo
 
      do j=1,6 ! J do loof - 필터고르기
       SS(j)=0.d0
       SES(j)=0.d0
       do i=300,1200,1
        w=real(i)
        call lint(w,Qwav,R,Rmir,10,10) !미러의반응율
        call lint(w,Qwav,Q,Qccd,10,10) !양자효율
        if(j.eq.1) then                !필터별 투과함수
         call lint(w,uw,ufi,Tfil,27,27)
         elseif(j.eq.2) then
         call lint(w,bw,bfi,Tfil,27,27)
         elseif(j.eq.3) then
         call lint(w,vw,vfi,Tfil,27,27)
         elseif(j.eq.4) then
         call lint(w,rw,rfi,Tfil,27,27)
         elseif(j.eq.5) then
         call lint(w,iw,ifi,Tfil,27,27)
         elseif(j.eq.6) then
         call lint(w,bw,bxfi,Tfil,27,27)
        endif
        call lint(w,wav,flux,pflux,nwav,nwav)
        Stot(j)=Tfil*Qccd*Rmir**2.d0
        SS(j)=SS(j)+Stot(j)
        SES(j)=SES(j)+pflux*Stot(j)
        moment(j)=moment(j)+(w-lamda0(j))**2.d0*Stot(j)
       enddo   
      enddo    
      do j=1,6
       moment(j)=sqrt(moment(j)/SS(j))
      enddo
      col(1)=-2.5d0*log10(cc(1)/cc(6))+0.79d0+0.
     *104d0
      col(2)=-2.5d0*log10(cc(2)/cc(3))-0.102d0-0
     *.008d0
      col(3)=-2.5d0*log10(cc(3)/cc(4))+0.008d0-0
     *.193d0
      col(4)=-2.5d0*log10(cc(3)/cc(5))+0.008d0-0
     *.443d0
      col(5)=-2.5d0*log10(cc(4)/cc(5))+0.193d0-0
     *.443d0
 
100   format(I1, T5, a7, T13, F7.3, T23, a7, T31, F7.3)
101   format(7f9.3)
102   format(T5,5f9.3)
      if(k.eq.1) then
       open(30,file='giant.sf')
       open(31,file='gub.sf')
       open(32,file='gbv.sf')
       open(33,file='gvr.sf')
       open(34,file='gvi.sf')
       open(35,file='gri.sf')
       write(30,*)'Filter 1=Ux, 2=B, 3=V, 4=R, 5=I, 6=Bx'
       do j=1,6
        write(30,100) j,'lamda0=',lamda0(j),'moment=',moment(j)
       enddo
       write(30,*) ''
      endif 
      write(30,*) ''
       do l=1,5
        write(30+l,101) col(l),(lameff(j),j=1,6)  
       enddo
      enddo

900   close(31)
      close(32)    
      close(33)
      close(34)
      close(35)
      stop
      end

      subroutine lint(x,xt,yt,y,n,m)
      parameter(nmax=10000)
      real*8 xt(n),yt(n),xt1(nmax),yt1(nmax)
      real*8 x,y
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
12    y=yt1(i)
      goto 14
10    continue
11    y=yt1(k)+(yt1(k+1)-yt1(k))*((x-xt1(k))/(xt1(k+1)-xt1(k)))
14    return
      end
    
