c     f95 -o make_std make_std.f
c     21,22,23,24
      character std*16,stdd*16,master*16,H*1
      open(21,file='input.make_std')
      read(21,'(t1,a1)') H
      do i=1,1000
       read(21,'(t1,3a16)',end=99) master,std,stdd
       call read_master(master)     
       call read_std(std)
       call write_std(stdd)
      enddo
99    close(21)
      stop
      end

      subroutine write_std(output)
      parameter(ns=10000)
      character output*16
      common/std/xc_s(ns),yc_s(ns),cmag(ns),cerr(ns),airmass(ns),ut,nstd
      character ut(ns)*12,star_m(ns)*13,line(ns)*135
      common/master/xc_m(ns),yc_m(ns),star_m,v_m(ns),ub_m(ns),bv_m(ns),
     &vi_m(ns),ri_m(ns),nv_m(ns),nub_m(ns),nbv_m(ns),nvi_m(ns),nri_m(ns)
     &,nmaster
      real dx,dy
      if(nstd.ne.nmaster) then
       print *,output,', has different numbers of stars'
       goto 900
      endif
      open(23,file=output)
      write(23,'(t1,a135)')'# Star             Xc        Yc      Mag    
     &  Merr  Airmass UT             V       U-B     B-V      V-I    R-I
     &  nV nUB nBV nVI nRI                                             '
      do i=1,nstd
       dx=abs(xc_s(i)-xc_m(i))
       dy=abs(yc_s(i)-yc_m(i))
       if(dx.le.30.and.dy.le.30) then
        write(line(i),'(t1,a15,2f10.3,2f8.3,f8.5,1x,a12,5f8.3,5i4)') 
     &  star_m(i),xc_s(i),yc_s(i),cmag(i),cerr(i),
     &  airmass(i),ut(i),v_m(i),ub_m(i),bv_m(i),vi_m(i),
     &  ri_m(i),nv_m(i),nub_m(i),nbv_m(i),nvi_m(i),nri_m(i)
       else
        print *,output,' ,',i,'st star has different coordinates!!!!'
        goto 900
       endif
      enddo
      call datasort(yc_s,line,nstd)
      do i=1,nstd
       write(23,'(t1,a135)') line(i)
      enddo
900   close(23)
      return
      end

      subroutine read_std(input)
      parameter(ns=10000)
      character input*16
      common/std/xc_s(ns),yc_s(ns),cmag(ns),cerr(ns),airmass(ns),ut,nstd
      character ut(ns)*12,c1*6,c2*5
      open(23,file=input)
      do i=1,ns
       cmag(i)=0.
       cerr(i)=0.
       read(23,*,end=99) xc_s(i),yc_s(i),c1,c2,airmass(i)
     & ,ut(i)
       if(c1.ne.'INDEF') read(c1,*) cmag(i)
       if(c2.ne.'INDEF') read(c2,*) cerr(i)
      enddo
99    nstd=i-1
      close(23)
      return
      end

      subroutine read_master(input)
      parameter(ns=10000)
      character input*16,H*1
      common/master/xc_m(ns),yc_m(ns),star_m,v_m(ns),ub_m(ns),bv_m(ns),
     &vi_m(ns),ri_m(ns),nv_m(ns),nub_m(ns),nbv_m(ns),nvi_m(ns),nri_m(ns)
     &,nmaster
      character star_m(ns)*13
      open(22,file=input)
      read(22,'(t1,a1)') H
      do i=1,ns
       read(22,*,end=99) xc_m(i),yc_m(i),star_m(i),v_m(i),ub_m(i),
     & bv_m(i),vi_m(i),ri_m(i),nv_m(i),nub_m(i),nbv_m(i),nvi_m(i),
     & nri_m(i)
      enddo
99    nmaster=i-1
      close(22)
      return
      end
      subroutine datasort(val,ch,nstar)
      integer nstar
      real val(nstar),ref
      character ch(nstar)*135,chref*135
      do i=1,nstar
       do j=1,nstar-1
        if(val(j+1).lt.val(j)) then
         ref=val(j)
         write(chref,'(a135)') ch(j)
         val(j)=val(j+1)
         write(ch(j),'(a135)') ch(j+1)
         val(j+1)=ref
         write(ch(j+1),'(a135)') chref
        endif
       enddo
      enddo
      return
      end

