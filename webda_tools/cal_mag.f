c f95 -o cal_mag cal_mag.f
c gfortran -o cal_mag ~/work/fortran/subroutine/webda_tools/cal_mag.f ~/work/fortran/subroutine/caltools/std.f 
      integer nvvi, nvi,nv,nbv,nub1,nub2,nadd
      character input*19,H*1,master*6
      open(25,file='draw_compare.input')
      nadd=0
      read(25,'(t1,a1)',end=99) H
      read(25,'(t1,a1)',end=99) H
      read(25,'(t1,a19,t20,6i3,t40)') input,nvvi,nvi,nv,nbv,nub1,
     &  nub2
      call read_master(input)
      do i=1,100
       read(25,'(t1,a19,t20,6i3)',end=99) input,nvvi,nvi,nv,nbv,
     & nub1,nub2
       call read_target(input,nadd,nvvi,nvi,nv,nbv,nub1,nub2)
      enddo
99    close(25)
      write(input,'(a19)') 'combine.dat                              '
      call read_proper()



      call write_table(input)
      stop
      end

      subroutine write_table(input)
      parameter (ns=100000)
      character input*19
      common/dat/no(ns),v(ns),bv(ns),ub(ns),vi(ns),
     &v_std(ns),bv_std(ns),ub_std(ns),vi_std(ns),
     &n_v(ns),n_bv(ns),n_ub(ns),n_vi(ns),nstarm,nstar
      common/proper/no_p(ns),nref_p(ns),pro(ns),npro
      real pm
      open(26,file=input)
      write(26,'(t1,a88)') '#   no   V      B-V      U-B     V-I   nobs 
     &           mem_p  V_std BV_std UB_std VI_std                   '
1     format(i6,4f8.3,4i4,f6.2,4f7.3)
      do i=1,nstar
       do j=1,npro
        if(no(i).eq.no_p(j)) then
         pm=pro(j)
         goto 10
        endif
       enddo
       pm=-1.
10     if(n_v(i).eq.0) then
        v(i)=0.
       endif
       if(n_bv(i).eq.0) then
        bv(i)=0.
       endif
       if(n_ub(i).eq.0) then
        ub(i)=0.
       endif
       if(n_vi(i).eq.0) then
        vi(i)=0.
       endif
       write(26,1) no(i),v(i),bv(i),ub(i),vi(i),n_v(i),n_bv(i),n_ub(i)
     & ,n_vi(i),pm,v_std(i),bv_std(i),ub_std(i),vi_std(i)
      enddo
      close(26)
      return
      end 


      subroutine read_target(input,nadd,nvvi,nvi,nv,nbv,nub1,nub2)
      parameter (ns=100000,nx=100)
      character input*19,H*1,num_ref*4
      common/dat/no(ns),v(ns),bv(ns),ub(ns),vi(ns),
     &v_std(ns),bv_std(ns),ub_std(ns),vi_std(ns),
     &n_v(ns),n_bv(ns),n_ub(ns),n_vi(ns),nstarm,nstar
      common/adat/av(ns),abv(ns),aub(ns),avi(ns),
     &rv(ns,nx),rbv(ns,nx),rub(ns,nx),rvi(ns,nx)
      real vr,tvi,tv,tbv,tub
      real dv,dbv,dub,dvi,cv(ns),cbv(ns),cub(ns),cvi(ns)
      integer no_t,nvvi,nvi,nv,nbv,nub1,nub2,nadd,no_a(ns)
      integer nvo,nbvo,nubo,nvio
      read(input,'(t1,a4)') num_ref
      call read_delta(num_ref,nvvi,nvi,nv,nbv,nub1,nub2)
      read(num_ref,*) nr
      open(24,file=input)
      read(24,'(t1,a1)') H
      do i=1,ns
       read(24,*,end=99) no_t,tv,tbv,tub,vr,tvi
       do j=1,nstarm
        if(no_t.eq.no(j)) then
         call cal_dv(num_ref,no_t,tvi,tv,tbv,tub1,tub2,dvi,dv,dbv,dub,
     &   nvvi,nvi,nv,nbv,nub1,nub2)
         call c_nobs(tv,tbv,tub,tvi,nvo,nbvo,nubo,nvio)
         av(j)=av(j)+(tv+dv)
         abv(j)=abv(j)+(tbv+dbv)
         aub(j)=aub(j)+(tub+dub)
         avi(j)=avi(j)+(tvi+dvi)
         n_v(j)=n_v(j)+nvo 
         n_bv(j)=n_bv(j)+nbvo
         n_ub(j)=n_ub(j)+nubo
         n_vi(j)=n_vi(j)+nvio
         if(nvo.eq.1) rv(j,n_v(j))=(tv+dv)
         if(nbvo.eq.1) rbv(j,n_bv(j))=(tbv+dbv)
         if(nubo.eq.1) rub(j,n_ub(j))=(tub+dub)
         if(nvio.eq.1) rvi(j,n_vi(j))=(tvi+dvi)
         goto 100
        endif
       enddo

       if(nadd.ne.0) then
        do j=nstarm,nstar
         if(no_t.eq.no(j)) then  
          call cal_dv(num_ref,no_t,tvi,tv,tbv,tub1,tub2,dvi,dv,dbv,dub,
     &    nvvi,nvi,nv,nbv,nub1,nub2)
          call c_nobs(tv,tbv,tub,tvi,nvo,nbvo,nubo,nvio)
          av(j)=av(j)+(tv+dv)
          abv(j)=abv(j)+(tbv+dbv)
          aub(j)=aub(j)+(tub+dub)
          avi(j)=avi(j)+(tvi+dvi)
          n_v(j)=n_v(j)+nvo
          n_bv(j)=n_bv(j)+nbvo
          n_ub(j)=n_ub(j)+nubo
          n_vi(j)=n_vi(j)+nvio
          if(nvo.eq.1) rv(j,n_v(j))=(tv+dv)
          if(nbvo.eq.1) rbv(j,n_bv(j))=(tbv+dbv)
          if(nubo.eq.1) rub(j,n_ub(j))=(tub+dub)
          if(nvio.eq.1) rvi(j,n_vi(j))=(tvi+dvi)
          goto 100
         endif
        enddo
       endif
       call cal_dv(num_ref,no_t,tvi,tv,tbv,tub1,tub2,dvi,dv,dbv,dub,
     & nvvi,nvi,nv,nbv,nub1,nub2)
       nadd=nadd+1
       nstar=nstar+1
       call c_nobs(tv,tbv,tub,tvi,nvo,nbvo,nubo,nvio)
       av(nstar)=(tv+dv)
       abv(nstar)=(tbv+dbv)
       aub(nstar)=(tub+dub)
       avi(nstar)=(tvi+dvi)
       n_v(nstar)=nvo
       n_bv(nstar)=nbvo
       n_ub(nstar)=nubo
       n_vi(nstar)=nvio
       if(nvo.eq.1) rv(j,n_v(nstar))=(tv+dv)
       if(nbvo.eq.1) rbv(j,n_bv(nstar))=(tbv+dbv)
       if(nubo.eq.1) rub(j,n_ub(nstar))=(tub+dub)
       if(nvio.eq.1) rvi(j,n_vi(nstar))=(tvi+dvi)
       no(nstar)=no_t
100    l=l
      enddo
99    close(24)
      do i=1,nstar
       v(i)=av(i)/real(n_v(i))
       bv(i)=abv(i)/real(n_bv(i))
       ub(i)=aub(i)/real(n_ub(i))
       vi(i)=avi(i)/real(n_vi(i))
       v_std(i)= 0.
       bv_std(i)=0.
       ub_std(i)=0.
       vi_std(i)=0.
       if(n_v(i).ge.1) then
        do j=1,n_v(i)
         cv(j)=rv(i,j)
        enddo
        call std(n_v(i),cv,v(i),v_std(i))
       endif
       if(n_bv(i).ge.1) then
        do j=1,n_bv(i)
         cbv(j)=rbv(i,j)
        enddo
        call std(n_bv(i),cbv,bv(i),bv_std(i))
       endif
       if(n_ub(i).ge.1) then
        do j=1,n_ub(i)
         cub(j)=rub(i,j)
        enddo
        call std(n_ub(i),cub,ub(i),ub_std(i))
       endif
       if(n_vi(i).ge.1) then
        do j=1,n_vi(i)
         cvi(j)=rvi(i,j)
        enddo
        call std(n_vi(i),cvi,vi(i),vi_std(i))
       endif
      enddo
      n=12
      return 
      end
      
      subroutine c_nobs(v,bv,ub,vi,nv,nbv,nub,nvi)
      integer nv,nbv,nub,nvi
      real v,bv,ub,vi
      if(v.lt.30.) then
       nv=1
      else
       nv=0
       v=0.
      endif
      if(bv.lt.30.) then
       nbv=1
      else
       nbv=0
       bv=0.
      endif
      if(ub.lt.30.) then
       nub=1
      else
       nub=0
       ub=0.
      endif
      if(vi.lt.30.) then
       nvi=1
      else
       nvi=0
       vi=0.
      endif
      return
      end 

      subroutine cal_dv(nref,no_t,tvi,tv,tbv,tub1,tub2,dvi,dv,dbv,dub,
     &nvvi,nvi,nv,nbv,nub1,nub2)
      parameter (nd=100)
      integer nvvi,nvi,nv,nbv,nub1,nub2,no_t
      common/delta/xvvi(nd),yvvi(nd),xvi(nd),yvi(nd),xv(nd),yv(nd)
     &,xbv(nd),ybv(nd),xub1(nd),yub1(nd),xub2(nd),yub2(nd),jvvi,jvi,jv,
     &jbv,jub1,jub2
      real tvi,tv,tbv,tub,dv,dbv,dub1,dub2,dvi,rbv,rub
      character nref*4
      if(nbv.eq.0) then
       dbv=0.
      else
       call lint(tbv,xbv,ybv,dbv,jbv,jbv)
       dbv=-dbv
      endif
      rbv=tbv+dbv 
      if(nvi.eq.0) then
       dvi=0.
      else
       call lint(tvi,xvi,yvi,dvi,jvi,jvi)
       dvi=-dvi
      endif
      if(nv.eq.0) then
       if(nvvi.eq.0) then
        dv=0.
       else
        call lint(tvi,xvvi,yvvi,dv,jvvi,jvvi)
       endif
      else
       call lint(rbv,xv,yv,dv,jv,jv)
       dv=-dv
      endif
      if(nub1.eq.0) then
       dub1=0.
      else
       call lint(rbv,xub1,yub1,dub1,jub1,jub1)
       dub1=-dub1
      endif
      if(nub2.eq.0) then
       dub2=0.
      else
       call lint(tub,xub2,yub2,dub2,jub2,jub2)
       dub2=-dub2
      endif
      if(tbv.ge.-0.4) then
       dub=dub1
      else
       dub=dub2
      endif 

      return
      end

      subroutine read_master(input)
      parameter (ns=100000,nx=100)
      character input*19,H*1
      common/dat/no(ns),v(ns),bv(ns),ub(ns),vi(ns),
     &v_std(ns),bv_std(ns),ub_std(ns),vi_std(ns),
     &n_v(ns),n_bv(ns),n_ub(ns),n_vi(ns),nstarm,nstar
      common/adat/av(ns),abv(ns),aub(ns),avi(ns),
     &rv(ns,nx),rbv(ns,nx),rub(ns,nx),rvi(ns,nx)
      real vr
      open(22,file=input)
      read(22,'(t1,a1)') H
      do i=1,ns
       read(22,*,end=99) no(i),av(i),abv(i),aub(i),vr,avi(i)
       n_v(i)=1
       n_bv(i)=1
       n_ub(i)=1
       n_vi(i)=1
       rv(i,n_v(i))=av(i)
       rbv(i,n_bv(i))=abv(i)
       rub(i,n_ub(i))=aub(i)
       rvi(i,n_vi(i))=avi(i)
c       if(no(i).eq.30) print *,av(i),n_v(i),rv(i,n_v(i)) !########
       if(av(i).gt.99.) then
        n_v(i)=0
        av(i)=0.
       endif
       if(abv(i).gt.99.) then
        n_bv(i)=0
        abv(i)=0.
       endif
       if(aub(i).gt.99.) then
        n_ub(i)=0
        aub(i)=0.
       endif
       if(avi(i).gt.99.) then
        n_vi(i)=0
        avi(i)=0.
       endif
c       if(no(i).eq.30) print *,av(i),n_v(i),rv(i,n_v(i)) !########
      enddo
99    nstarm=i-1
      nstar=nstarm
      close(22)
      return 
      end
      
      subroutine read_delta(num_ref,nvvi,nvi,nv,nbv,nub1,nub2)
      parameter (nd=100)
      integer nvvi,nvi,nv,nbv,nub1,nub2
      common/delta/xvvi(nd),yvvi(nd),xvi(nd),yvi(nd),xv(nd),yv(nd)
     &,xbv(nd),ybv(nd),xub1(nd),yub1(nd),xub2(nd),yub2(nd),jvvi,jvi,jv,
     &jbv,jub1,jub2
      character input*8,H*1,num_ref*4
      if(nvvi.ne.0) then
       write(input,'(2a4)') num_ref,'.vvi'
       open(23,file=input)
       read(23,'(t1,a1)') H
       do i=1,nd
        read(23,*,end=79) xvvi(i),yvvi(i)
       enddo
79     jvvi=i-1
       close(23)
      endif
      if(nvi.ne.0) then
       write(input,'(2a4)') num_ref,'.vi '
       open(23,file=input)
       read(23,'(t1,a1)') H
       do i=1,nd
        read(23,*,end=99) xvi(i),yvi(i)
       enddo
99     jvi=i-1
       close(23)
      endif
      if(nv.ne.0) then
       write(input,'(2a4)') num_ref,'.v  '
       open(23,file=input)
       read(23,'(t1,a1)') H
       do i=1,nd
        read(23,*,end=199) xv(i),yv(i)
       enddo
199    jv=i-1
       close(23)
      else
       jv=0
      endif
      if(nbv.ne.0) then
       write(input,'(2a4)') num_ref,'.bv '
       open(23,file=input)
       read(23,'(t1,a1)') H
       do i=1,nd
        read(23,*,end=299) xbv(i),ybv(i)
       enddo
299    jbv=i-1
       close(23)
      else
       jbv=0
      endif
      if(nub1.ne.0) then
       write(input,'(2a4)') num_ref,'.ub '
       open(23,file=input)
       read(23,'(t1,a1)') H
       do i=1,nd
        read(23,*,end=399) xub1(i),yub1(i)
       enddo
399    jub1=i-1
       close(23)
      else
       jub1=0
      endif
      if(nub2.ne.0) then
       write(input,'(2a4)') num_ref,'.ub2'
       open(23,file=input)
       read(23,'(t1,a1)') H
       do i=1,nd
        read(23,*,end=499) xub2(i),yub2(i)
       enddo
499    jub2=i-1
       close(23)
      else
       jub2=0
      endif
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
12    y=yt1(i)
      goto 14
10    continue
11    y=yt1(k)+(yt1(k+1)-yt1(k))*((x-xt1(k))/(xt1(k+1)-xt1(k)))
14    return
      end

      subroutine read_proper()
      parameter (ns=100000)
      common/proper/no_p(ns),nref_p(ns),pro(ns),npro
      character sub*40,tar*20,line*60,H*1,input*20,c(60)*1
      integer nprob_ref,no,nref
      real p
c      print *,'input file name=?'
c      read(*,'(a20)') input
      open(21,file='make_table.input')
      do i=1,100
       read(21,'(t1,a1)',end=99) H
       if(H.eq.'#') then
        backspace(21)
        read(21,'(t3,a60)') line
        read(line,'(60a1)') (c(j),j=1,60)
        do j=1,60
         if(c(j).eq.':') then 
          goto 50
         endif
        enddo
50      write(sub,'(40a1)') (c(k),k=1,j-2)
        if(sub.eq.'proper motion data') then
         write(input,'(20a1)') (c(k),k=j+2,j+21)
        elseif(sub.eq.'prob') then
         write(tar,'(20a1)') (c(k),k=j+2,j+21)
         if(tar.ne.'none') read(tar,*) nprob_ref
        endif
       endif
      enddo
99    open(22,file=input)
      read(22,'(t1,a1)',end=100) H
      read(22,'(t1,a1)') H
      npro=0
      do j=1,ns
       read(22,*,end=100) no,nref,p
       if(nref.eq.nprob_ref) then
        npro=npro+1
        no_p(npro)=no
        nref_p(npro)=nref
        pro(npro)=p
       endif
      enddo
100   close(22)                
      close(21)
      return
      end
