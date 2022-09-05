      SUBROUTINE ks2d1s(x1,y1,n1,quadvl,d1,prob)
      INTEGER n1
      REAL d1,prob,x1(n1),y1(n1)
      EXTERNAL quadvl
CU    USES pearsn,probks,quadct,quadvl
      INTEGER j
      REAL dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,rr,sqen,probks
      d1=0.0
      do 11 j=1,n1
        call quadct(x1(j),y1(j),x1,y1,n1,fa,fb,fc,fd)
        call quadvl(x1(j),y1(j),ga,gb,gc,gd)
        d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
11    continue
      call pearsn(x1,y1,n1,r1,dum,dumm)
      sqen=sqrt(float(n1))
      rr=sqrt(1.0-r1**2)
      prob=probks(d1*sqen/(1.0+rr*(0.25-0.75/sqen)))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0"0)P.
      SUBROUTINE ks2d2s(x1,y1,n1,x2,y2,n2,d,prob)
      INTEGER n1,n2
      REAL d,prob,x1(n1),x2(n2),y1(n1),y2(n2)
CU    USES pearsn,probks,quadct
      INTEGER j
      REAL d1,d2,dum,dumm,fa,fb,fc,fd,ga,gb,gc,gd,r1,r2,rr,sqen,probks
      d1=0.0
      do 11 j=1,n1
        call quadct(x1(j),y1(j),x1,y1,n1,fa,fb,fc,fd)
        call quadct(x1(j),y1(j),x2,y2,n2,ga,gb,gc,gd)
        d1=max(d1,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
11    continue
      d2=0.0
      do 12 j=1,n2
        call quadct(x2(j),y2(j),x1,y1,n1,fa,fb,fc,fd)
        call quadct(x2(j),y2(j),x2,y2,n2,ga,gb,gc,gd)
        d2=max(d2,abs(fa-ga),abs(fb-gb),abs(fc-gc),abs(fd-gd))
12    continue
      d=0.5*(d1+d2)
      sqen=sqrt(float(n1)*float(n2)/float(n1+n2))
      call pearsn(x1,y1,n1,r1,dum,dumm)
      call pearsn(x2,y2,n2,r2,dum,dumm)
      rr=sqrt(1.0-0.5*(r1**2+r2**2))
      prob=probks(d*sqen/(1.0+rr*(0.25-0.75/sqen)))
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0"0)P.
      SUBROUTINE ksone(data,n,func,d,prob)
      INTEGER n
      REAL d,data(n),func,prob
      EXTERNAL func
CU    USES probks,sort
      INTEGER j
      REAL dt,en,ff,fn,fo,probks
      call sort(n,data)
      en=n
      d=0.
      fo=0.
      do 11 j=1,n
        fn=j/en
        ff=func(data(j))
        dt=max(abs(fo-ff),abs(fn-ff))
        if(dt.gt.d)d=dt
        fo=fn
11    continue
      en=sqrt(en)
      prob=probks((en+0.12+0.11/en)*d)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0"0)P.
      SUBROUTINE kstwo(data1,n1,data2,n2,d,prob)
      INTEGER n1,n2
      REAL d,prob,data1(n1),data2(n2)
CU    USES probks,sort
      INTEGER j1,j2
      REAL d1,d2,dt,en1,en2,en,fn1,fn2,probks
      call sort(n1,data1)
      call sort(n2,data2)
      en1=n1
      en2=n2
      j1=1
      j2=1
      fn1=0.
      fn2=0.
      d=0.
1     if(j1.le.n1.and.j2.le.n2)then
        d1=data1(j1)
        d2=data2(j2)
        if(d1.le.d2)then
          fn1=j1/en1
          j1=j1+1
        endif
        if(d2.le.d1)then
          fn2=j2/en2
          j2=j2+1
        endif
        dt=abs(fn2-fn1)
        if(dt.gt.d)d=dt
      goto 1
      endif
      en=sqrt(en1*en2/(en1+en2))
      prob=probks((en+0.12+0.11/en)*d)
      return
      END
C  (C) Copr. 1986-92 Numerical Recipes Software 0"0)P.
       subroutine sort(ns,EID11)
       integer CID
       real Crad,EID11(ns)
       CID=0
       Crad=0.
       do i=1,ns
        do j=1,ns-1
         if (EID11(j+1).lt.EID11(j)) then
         Crad=EID11(j)
         EID11(j)=EID11(j+1)
         EID11(j+1)=Crad
         endif
        enddo
       enddo
      return
      end


