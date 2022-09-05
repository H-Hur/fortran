c     f77 -o test_trs /home/gjgusdh/fortran/subroutine/fs/test_trs.f -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 
      parameter (ns=100000)
      integer nstar,id(ns),npartner(ns),nrefpartner(ns),nsingle     
      integer ndouble,nfake,ncall
      real xc(ns),yc(ns),dist_c(ns),fitrad(ns),x(3),y(3),std_d,av_d
      real xmax,ymax,dmax,rmax
      character cx(3)*22,cy(3)*22,input*30,H*1,als*30,ref*30,filter*30
      character cl(30)*1,output*30,ps*43
      print *,'input trs file=?'
      read(*,'(a30)') input
      open(21,file=input)
      read(21,'(t28,a30)') als
      read(21,'(t28,a30)') ref
      read(21,'(t11,a30)') filter
      read(21,'(t8,a22,t38,a22,t79,a22)') (cx(j),j=1,3)
      read(21,'(t8,a22,t38,a22,t79,a22)') (cy(j),j=1,3)
      do i=1,3
       read(cx(i),*) x(i)
       read(cy(i),*) y(i)
      enddo
      read(21,*) nstar
      read(21,*) nsingle
      read(21,*) ndouble
      read(21,*) nfake
      read(21,'(t5,a1)') H
      read(21,'(t5,a1)') H
      read(input,'(30a1)') (cl(j),j=1,30)
      do i=1,30
       if(cl(i).eq.'.') goto 100
      enddo
100   write(cl(i+1),'(a1)') 'p'
      write(cl(i+2),'(a1)') 's'
      write(output,'(30a1)') (cl(j),j=1,i+2)
      write(ps,'(a13,a30)') 'postportfile ',output
      print *,ps

      xmax=0.
      ymax=0.
      dmax=0.
      rmax=0.

      do i=1,nstar
       read(21,*) id(i),xc(i),yc(i),npartner(i),fitrad(i),dist_c(i),
     & nrefpartner(i)
       xmax=max(xmax,xc(i))
       ymax=max(ymax,yc(i))
       rmax=max(rmax,dist_c(i))
       dmax=max(dmax,fitrad(i))
      enddo

      call cal_std(nstar,fitrad,av_d,std_d,ncall)
      print *,'mean fitrad = ',av_d,' , stddev = ',std_d

      call sm_device (ps)
      call sm_graphics
      call sm_erase()
      
      call sm_location(3000,30000,3000,10000)
      call sm_limits(0.,xmax+10.,0.,dmax+0.1)
      call sm_box(1,2,0,0)
      call sm_relocate(0.,av_d)
      call sm_draw(xmax+10.,av_d)
      call sm_relocate(0.,av_d+2.5*std_d)
      call sm_draw(xmax+10.,av_d+2.5*std_d)
      call sm_xlabel('Xcenter')
      call sm_ylabel('Fitrad (pixel)')
      do i=1,nstar
       if(npartner(i).ge.1.and.nrefpartner(i).ge.1) then
        call sm_ptype(43.2,1)
        if(npartner(i).ge.2.or.nrefpartner(i).ge.2) then
         call sm_ptype(43.9,1)
        elseif(npartner(i).eq.0.or.nrefpartner(i).eq.0) then
         call sm_ptype(40.9,1)
        endif
        call sm_points(xc(i),fitrad(i),1)
       endif
      enddo

      call sm_location(3000,30000,13000,20000)
      call sm_limits(0.,ymax+10.,0.,dmax+0.1)
      call sm_box(1,2,0,0)
      call sm_xlabel('Ycenter')
      call sm_ylabel('Fitrad (pixel)')
      call sm_relocate(0.,av_d)
      call sm_draw(ymax+10.,av_d)
      call sm_relocate(0.,av_d+2.5*std_d)
      call sm_draw(ymax+10.,av_d+2.5*std_d)
      do i=1,nstar
       if(npartner(i).ge.1.and.nrefpartner(i).ge.1) then
        call sm_ptype(43.2,1)
        if(npartner(i).ge.2.or.nrefpartner(i).ge.2) then
         call sm_ptype(43.9,1)
        elseif(npartner(i).eq.0.or.nrefpartner(i).eq.0) then
         call sm_ptype(40.9,1)
        endif
        call sm_points(yc(i),fitrad(i),1)
       endif
      enddo

      call sm_location(3000,30000,23000,30000)
      call sm_limits(0.,rmax+10.,0.,dmax+0.1)
      call sm_box(1,2,0,0)
      call sm_xlabel('Distance from Center')
      call sm_ylabel('Fitrad (pixel)')
      call sm_relocate(0.,av_d)
      call sm_draw(rmax+10.,av_d)
      call sm_relocate(0.,av_d+2.5*std_d)
      call sm_draw(rmax+10.,av_d+2.5*std_d)
      do i=1,nstar
       if(npartner(i).ge.1.and.nrefpartner(i).ge.1) then
        call sm_ptype(43.2,1)
        if(npartner(i).ge.2.or.nrefpartner(i).ge.2) then
         call sm_ptype(43.9,1)
        elseif(npartner(i).eq.0.or.nrefpartner(i).eq.0) then
         call sm_ptype(40.9,1)
        endif
        call sm_points(dist_c(i),fitrad(i),1)
       endif
      enddo

      
      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      stop
      end
     
      subroutine cal_STD(n,dat,average,std,ncall)
      real dat(n),average,std
      integer ncall
      ncall=ncall+1
      average=0.
      std=0.
      if(n.eq.0) then
       print*,'Subroutine STD received no data, ndata=0 when ncall='
     & ,ncall, 'mean=0, stddev=0'
       goto 100
      elseif(n.eq.1) then
       average=dat(1)
       std=0.
       print*,'Subroutine STD received single data, ndata=1 when ncall='
     & ,ncall, 'mean=',average,' stddev=0'
      else
       do i=1,n
        average=average+dat(i)
       enddo
       average=average/real(n)
       do i=1,n
        std=std+(average-dat(i))*(average-dat(i))
       enddo
       std=sqrt(std/real(n-1))
      endif
100   return
      end
