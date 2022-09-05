c     f95 -o pcoo_fs /home/gjgusdh/fortran/subroutine/fs/pcoo_fs.f
      character input*30
      write(input,'(a30)') 'out2.fs                                    '
      call read_fs(input)

      write(input,'(a30)') 'input.pcoo                                 '
      call read_input(input)

      call write_coo()

      stop
      end



      subroutine write_coo()
      parameter (ns=100000,ni=30)
      common/fs/filter,pmag(ni),time(ni),airmass(ni),ID(ns,ni),
     &xc(ns,ni),yc(ns,ni),xc_ref(ns,ni),yc_ref(ns,ni),
     &smag(ns,ni),smerr(ns,ni),nobs(ns,ni),nstar,nimage
      character filter(ni)*10
      common/coo/output,noutput
      character output(ni)*30
      do k=2,nimage
       open(24,file=output(k))
       do i=1,nstar
        if(nobs(i,k).ne.0) then
         write(24,'(t1,2f10.3,i8)') xc_ref(i,k),yc_ref(i,k),id(i,k)
        endif
       enddo
       close(24)
      enddo

      return
      end

      subroutine read_input(input)
      parameter (ni=30)
      common/coo/output,noutput
      character output(ni)*30
      character input*30
      open(22,file=input)
      do i=1,ni
       read(22,'(t1,a30)',end=99) output(i)
       print *,output(i)
      enddo
99    noutput=i-1
      close(22)
      return
      end

      subroutine read_fs(input)
      parameter (ns=100000,ni=30)
      common/fs/filter,pmag(ni),time(ni),airmass(ni),ID(ns,ni),
     &xc(ns,ni),yc(ns,ni),xc_ref(ns,ni),yc_ref(ns,ni),
     &smag(ns,ni),smerr(ns,ni),nobs(ns,ni),nstar,nimage
      character filter(ni)*10
      character input*30,H*1,line*300
      real rh
      open(21,file=input)
      do j=1,ni
       airmass(j)=0
       time(j)=0
       pmag(j)=0
      enddo
      do 10  i=1,ns
       read(21,'(t1,a1,t2,a300)') H,line
       if(H.eq.'*') goto 98
       if(H.eq.'X') then
        read(line,*,end=10) (airmass(j),j=1,ni)
       elseif(H.eq.'T') then
        read(line,*,end=10) (time(j),j=1,ni)
       elseif(H.eq.'P') then
        read(line,*,end=10) (pmag(j),j=1,ni)
       endif
10    continue

98    do j=1,ni
       if(pmag(j).eq.0.and.airmass(i).eq.0.and.time(j).eq.0) then
        nimage=j-1
        goto 99
       endif
      enddo

99    do i=1,ns
       read(21,*,end=199) (id(i,j),j=1,nimage)
       read(21,*,end=199) (xc(i,j),j=1,nimage)
       read(21,*,end=199) (yc(i,j),j=1,nimage)
       read(21,*,end=199) (xc_ref(i,j),j=1,nimage)
       read(21,*,end=199) (yc_ref(i,j),j=1,nimage)
       read(21,*,end=199) (smag(i,j),j=1,nimage)
       read(21,*,end=199) (smerr(i,j),j=1,nimage)
       read(21,*,end=199) rh
       read(21,*,end=199) rh
       do j=1,nimage
        if(id(i,j).ne.0) nobs(i,j)=1
       enddo
      enddo
199   nstar=i-1
      close(21)
      return
      end
