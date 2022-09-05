c for only fortran77 to draw color-magnitude daigram.
c f77 -o 'draw_relation' 'draw_relation.f' -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11
      character datafile*16,input*16,xl*10,yl*10
      print *,'input file name=?'
      read(*,'(a16)') input
c      write(input,'(a16)') '348.rel               '
c      print *,'Input dat file=?'
c      read(*,'(a16)') datafile
      call read_input(input,datafile)
      open(21,file=datafile)
      call read_dmag(21)
      close(21)
      call sm_device ('postportfile relations.ps')
      call sm_graphics
      call sm_erase()
      call sm_location(4000,25000,27500,28500)
      call sm_limits(0.,100.,0.,100.)
      call sm_relocate(10.,50.)
      call sm_label(datafile)
      write(xl,'(a10)') 'R-I                        '
      write(yl,'(a10)') '\\gD(R-I)                        '
      call draw_dmag(2,2,-0.2,2.,-0.3,0.3,4000,25000,25000,27500
     &,xl,yl,1,2)
      write(xl,'(a10)') 'V-I                        '
      write(yl,'(a10)') '\\gD(V-I)                       '
      call draw_dmag(3,3,-0.2,3.,-0.3,0.3,4000,25000,19500,22000
     &,xl,yl,1,2)
      write(yl,'(a10)') '\\gDV                       '
      call draw_dmag(4,1,-0.2,2.5,-0.3,0.3,4000,25000,14000,16500
     &,xl,yl,0,2)
      write(yl,'(a10)') '\\gD(B-V)                        '
      call draw_dmag(4,4,-0.2,2.5,-0.3,0.3,4000,25000,11500,14000
     &,xl,yl,0,2)
      write(xl,'(a10)') 'B-V                        '
      write(yl,'(a10)') '\\gD(U-B)                        '
      call draw_dmag(4,5,-0.2,2.5,-0.5,0.5,4000,25000,9000,11500
     &,xl,yl,1,2)
      write(xl,'(a10)') 'U-B                      '
      write(yl,'(a10)') '\\gD(U-B)                        '
      call draw_dmag(5,5,-1.2,2.,-0.5,0.5,4000,25000,3500,6000
     &,xl,yl,1,2)
      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      stop
      end

      subroutine eject_highstd(nx,ny,slope,dy,xp,yp,np,dm,nin,
     &avin,stdin,avout,stdout,nstdout,nexout_all)
      parameter (ns=10000)
      common/dmag/ncl(ns),no(ns),id(ns),dcm(ns,5,3),nobs(ns,5),
     &ncal(5),nex(5),ndmag
      real*8 avin,stdin,avout,stdout
      integer nx,ny,np,nin(ns),nej(ns),nstdout,nexout,nroof,nexout_all
      real slope,dy,xp(ns),yp(ns),dm(ns)
      do i=1,ndmag
       nej(i)=0
      enddo
      nroof=0
      nexout_all=0
100   avout=0.d0
      stdout=0.d0
      nstdout=0
      nexout=0
      do i=1,ndmag 
c       if(nroof.eq.0.and.nx.eq.2.and.nin(i).eq.1) then
c        print *,dm(i),nin(i),abs(dm(i))-2.5d0*stdin
c       endif
       if((nin(i).ne.0).and.(abs(dm(i)).gt.2.5d0*stdin)) then
        nin(i)=0
        nej(i)=1
        nexout=nexout+1
       elseif(nin(i).ne.0) then
        nin(i)=1
        nej(i)=0
        avout=avout+dm(i)
        nstdout=nstdout+1
       endif
      enddo
      avout=avout/real(nstdout)
      do i=1,ndmag
       if(nin(i).ne.0) then
        stdout=stdout+(avout-dm(i))**2.d0
       endif
      enddo
      stdout=sqrt(stdout/real(nstdout-1))
      if(nexout.ne.0) then
c       print *,'100'
       nexout_all=nexout_all+nexout
       if(nx.eq.4.and.ny.eq.1) then
        print *,nx,ny,avin,stdin,avout,stdout,nstdout,nexout_all,nroof
       endif
       stdin=stdout
       nroof=nroof+1
       goto 100
      endif
      return
      end


      subroutine cal_mean_std(ncol,nrel,nx,ny,slope,dy,xp,yp,np,dm,nin,
     &average,std,nstd)
      parameter (ns=10000)
      common/dmag/ncl(ns),no(ns),id(ns),dcm(ns,5,3),nobs(ns,5),
     &ncal(5),nex(5),ndmag
      real*8 average,std
      real y_mean,slope,dy,xp(ns),yp(ns),dm(ns)
      integer nx,ny,nstd,np,ncol,nin(ns)
      nstd=0
      std=0.d0
      average=0.d0
      do i=1,ndmag
       nin(i)=0
       if((nx.eq.1).and.(nobs(i,ny).ne.0)) then
        nin(i)=1
        nstd=nstd+1
        if(nrel.eq.1) then
         y_mean=slope*dcm(i,nx,1)+dy
        elseif(nrel.eq.2) then
         call lint(dcm(i,nx,1),xp,yp,y_mean,np,np)
        endif
        dm(i)=dcm(i,ny,3)-y_mean
        average=average+dm(i)
       elseif((ny.eq.1).and.(nobs(i,nx)).ne.0) then
        nin(i)=1
        nstd=nstd+1
        if(nrel.eq.1) then
         y_mean=slope*dcm(i,nx,1)+dy
        elseif(nrel.eq.2) then
         call lint(dcm(i,nx,1),xp,yp,y_mean,np,np)
        endif
        dm(i)=dcm(i,ny,3)-y_mean
        average=average+dm(i)
       elseif((nobs(i,ny).ne.0).and.(nobs(i,nx).ne.0)) then
        nin(i)=1
        nstd=nstd+1
        if(nrel.eq.1) then
         y_mean=slope*dcm(i,nx,1)+dy
        elseif(nrel.eq.2) then
         call lint(dcm(i,nx,1),xp,yp,y_mean,np,np)
        endif
        dm(i)=dcm(i,ny,3)-y_mean
        average=average+dm(i)
       endif
      enddo
      average=average/real(nstd)

      do i=1,ndmag
       if((nx.eq.1).and.(nobs(i,ny).ne.0)) then
        std=std+(average-dm(i))**2.d0
       elseif((ny.eq.1).and.(nobs(i,nx)).ne.0) then
        std=std+(average-dm(i))**2.d0
       elseif((nobs(i,ny).ne.0).and.(nobs(i,nx).ne.0)) then
        std=std+(average-dm(i))**2.d0
       endif
      enddo
      std=sqrt(std/real(nstd-1))
      return
      end

      subroutine draw_line(nx,ny,x1,x2,y1,y2)
      parameter (ns=10000)
      common/rel/num_x(6),num_y(6),num_rel(6),dy(6),slope(6)
     &,relations
      character relations(6)*16,label*32
      real*8 average,std,avout,stdout
      real xp(ns),yp(ns),x1,x2,y1,y2,dm(ns)
      integer np,nin(ns),nstddev,neject
      do i=1,6
       if((nx.eq.num_x(i)).and.(ny.eq.num_y(i)).and.(num_rel(i).eq.1))
     & then
        call sm_relocate(-5.,-5*slope(i)+dy(i))
        call sm_ctype('red')
        call sm_draw(5.,5.*slope(i)+dy(i))
        call sm_ctype('black')
        call cal_mean_std(i,num_rel(i),nx,ny,slope(i),dy(i)
     &  ,xp,yp,np,dm,nin,average,std,nstd)
        write(*,'(2f8.4,3i2)')average,std,i,nx,ny
        call eject_highstd(nx,ny,slope(i),dy(i),xp,yp,np,dm,nin,
     &  average,std,avout,stdout,nstddev,neject) 
        write(label,'(f7.4,a5,f6.4,a5,i4,a1,i3,a1)') 
     &  avout,' \\g+',stdout,'mag, ',nstddev,'(',neject,')'  
        call sm_relocate(x1+0.59*(x2-x1),y1+0.12*(y2-y1))
        call sm_expand(0.5)
        call sm_label(label)
        call sm_expand(1.)
       elseif((nx.eq.num_x(i)).and.(ny.eq.num_y(i)).and.
     & (num_rel(i).eq.2)) then
        open(21,file=relations(i))
        call passhead(21)
        do j=1,ns
         read(21,*,end=99) xp(j),yp(j)
        enddo
99      np=j-1
        close(21)
        call sm_ctype('red')
        call sm_conn(xp,yp,np)
        call sm_ctype('black')
        call cal_mean_std(i,num_rel(i),nx,ny,slope(i),dy(i)
     &  ,xp,yp,np,dm,nin,average,std,nstd)
        write(*,'(2f8.4,3i2)')average,std,i,nx,ny
        call eject_highstd(nx,ny,slope(i),dy(i),xp,yp,np,dm,nin,
     &  average,std,avout,stdout,nstddev,neject)
        write(label,'(f7.4,a5,f6.4,a5,i4,a1,i3,a1)')
     &  avout,' \\g+',stdout,'mag, ',nstddev,'(',neject,')'      
        call sm_relocate(x1+0.59*(x2-x1),y1+0.12*(y2-y1))
        call sm_expand(0.5)
        call sm_label(label)
        call sm_expand(1.)
       endif
      enddo 
      return
      end

!     num_x,num_y : 1=dV 2=d(R-I) 3=d(V-I) 4=d(B-V) 5=d(U-B)
!     x1,x2,y1,y2 : xlange,ylange to draw
!     l1,l2,l3,l4 : location of figure
      subroutine draw_dmag(num_x,num_y,x1,x2,y1,y2,l1,l2,l3,l4,
     &xl,yl,mx,my)
      parameter (ns=10000)
      common/dmag/ncl(ns),no(ns),id(ns),dcm(ns,5,3),nobs(ns,5),
     &ncal(5),nex(5),ndmag
      integer num_x,num_y,l1,l2,l3,l4,mx,my
      real x1,x2,y1,y2,pp(1),xlocate,ylocate
      character xl*10,yl*10
      call sm_location(l1,l2,l3,l4)        !l1,l2,l3,l4
      call sm_limits(x1,x2,y1,y2)          !x1,x2,y1,y2
      call sm_ctype('black')
      call sm_expand(0.7)
      call sm_box(mx,my,0,0)               ! mx,my
      call sm_expand(1.)
      pp(1)=993.+0.3
      call sm_ptype(pp,1)
      do i=1,ndmag
       if((num_x.eq.1).and.(nobs(i,num_y).ne.0)) then
        call sm_points(dcm(i,num_x,1),dcm(i,num_y,3),1)
       elseif((num_y.eq.1).and.(nobs(i,num_x).ne.0)) then
        call sm_points(dcm(i,num_x,1),dcm(i,num_y,3),1)
       elseif((nobs(i,num_x).ne.0).and.(nobs(i,num_y).ne.0)) then
        call sm_points(dcm(i,num_x,1),dcm(i,num_y,3),1)
       endif
      enddo
      call sm_relocate(x1,0.)
      call sm_draw(x2,0.)
      call draw_line(num_x,num_y,x1,x2,y1,y2)
      call sm_expand(0.7)
      call sm_ylabel(yl)                   !yl
      call sm_xlabel(xl)
      call sm_expand(1.)
      write(xl,'(a10)') '                            '
      write(yl,'(a10)') '                            '
      return
      end

      subroutine read_input(input,datafile)
      common/rel/num_x(6),num_y(6),num_rel(6),dy(6),slope(6)
     &,relations
      character relations(6)*16
      character input*16,datafile*16
      open(21,file=input)
      call passhead(21)
      read(21,'(a16)') datafile
      do i=1,6
       read(21,'(t1,3i2,a16)') num_x(i),num_y(i),num_rel(i),relations(i)
       if(num_rel(i).eq.1) then 
        read(relations(i),*) dy(i),slope(i)
       endif 
      enddo
      close(21)
      return
      end

      subroutine read_dmag(input)
      parameter (ns=10000)
      integer input
      common/dmag/ncl(ns),no(ns),id(ns),dcm(ns,5,3),nobs(ns,5),
     &ncal(5),nex(5),ndmag
      call passhead(input)
      do i=1,5
       ncal(i)=0
       nex(i)=0
      enddo
      do i=1,ns
       read(input,*,end=99) ncl(i),no(i),id(i)
     & ,((dcm(i,k,j),j=1,3),k=1,5)
     & ,(nobs(i,j),j=1,5)
c       write(*,'(15f9.4)')((dcm(i,k,j),j=1,3),k=1,5)
       do j1=1,5
        if(nobs(i,j1).ne.0) then
         ncal(j1)=ncal(j1)+1
        endif
       enddo
      enddo
99    ndmag=i-1
      return
      end

      subroutine passhead(ninst)
      parameter (ns=100000)
      integer ninst
      character H*2,H1
10    format(t1,a2,t1,a1)
      do i=1,ns
       read(ninst,10) H,H1
       if((H.eq.'* ').or.(H.eq.'#*').or.(H.eq.'--')) then
        goto 100
       elseif(H1.eq.'*') then
        goto 100
       endif
      enddo
100   return
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

