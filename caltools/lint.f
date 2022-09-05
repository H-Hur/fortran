! linear interpolation program for double precision.
      subroutine lint_dp(x,xt,yt,y,n,m)
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

! linear interpolation program for fortran.
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
12    y=yt1(i)
      goto 14
10    continue
11    y=yt1(k)+(yt1(k+1)-yt1(k))*((x-xt1(k))/(xt1(k+1)-xt1(k)))
14    return
      end
