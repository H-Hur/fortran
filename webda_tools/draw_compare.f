c f77 -o 'draw_compare' 'draw_compare.f' -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11
c g77 -o draw_compare ~/work/fortran/subroutine/webda_tools/draw_compare.f -L/Applications/scisoft/i386/Packages/sm/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 -L/Users/apple/utils/sm_mac/scisoft/lib -laquaterm
      character H*1,input*19,master*6,ps*33,c(19)*1,output*19,num_ref*4
      integer nvi,nv,nbv,nub1,nub2,nvvi
      real del_y,x,y
      open(24,file='draw_compare.input')
      do i=1,50
       read(24,'(t1,a1,t40,a6)',end=99) H,master
       if(H.ne.'#'.and.master.eq.'master') then
        backspace(24)
        read(24,'(t1,a19,t47,f5.2)') input,del_y
        call read_master(input)
       endif
      enddo
99    close(24)


      open(25,file='draw_compare.input')
      do i=1,50
       read(25,'(t1,a1,t40,a6)',end=900) H,master
       if(H.ne.'#'.and.master.ne.'master') then
        backspace(25)
        read(25,'(t1,a19,t20,6i3)') input,nvvi,nvi,nv,nbv,nub1,nub2
        call read_dat(input)
        read(input,'(t1,a4)') num_ref
        read(input,'(19a1)') (c(j),j=1,19)
        do j=1,19
         if(c(j).eq.'.') then
          goto 100
         endif
        enddo
100     write(c(j+1),'(a1)') 'p'
        write(c(j+2),'(a1)') 's'
        write(output,'(a19)') '                                    '
        write(output,'(19a1)') (c(k),k=1,j+2)
        write(ps,'(a13,a19)') 'postportfile ',output
        call sm_device (ps)
        call sm_graphics
        call sm_erase()
 
        x=-0.4+0.55*(4.5+0.4)
        y=-del_y+0.8*(del_y+del_y)
        call sm_location(7000,22000,25000,28000)
        call sm_limits(-0.4, 4.5, -del_y, del_y)
        call draw_compare(num_ref,5,1,0,2,0,0,nvvi,nvi,nv,nbv,nub1,nub2,
     &  x,y)
        call sm_ylabel('\\gDV')
        call sm_location(7000,22000,22000,25000)
        call sm_limits(-0.4, 4.5, -del_y, del_y)
        call draw_compare(num_ref,5,5,1,2,0,0,nvvi,nvi,nv,nbv,nub1,nub2,
     &  x,y)
        call sm_xlabel('V-I')
        call sm_ylabel('\\gD(V-I)')
  
        call sm_location(7000,22000,16000,19000)
        call sm_limits(-0.4, 2.7, -del_y, del_y)
        x=-0.4+0.55*(2.7+0.4)
        call draw_compare(num_ref,2,1,0,2,0,0,nvvi,nvi,nv,nbv,nub1,nub2,
     &  x,y)
        call sm_ylabel('\\gDV')
        call sm_location(7000,22000,13000,16000)
        call sm_limits(-0.4, 2.7, -del_y, del_y)
        call draw_compare(num_ref,2,2,0,2,0,0,nvvi,nvi,nv,nbv,nub1,nub2,
     &  x,y)
        call sm_ylabel('\\gD(B-V)')
        call sm_location(7000,22000,10000,13000)
        call sm_limits(-0.4, 2.7, -del_y, del_y)
        call draw_compare(num_ref,2,3,1,2,0,0,nvvi,nvi,nv,nbv,nub1,nub2,
     &  x,y)
        call sm_ylabel('\\gD(U-B)')
        call sm_xlabel('B-V')
  
        call sm_location(7000,22000,4000,7000)
        call sm_limits(-1., 2.3, -del_y, del_y)
        x=-1.+0.55*(2.3+1.)
        call draw_compare(num_ref,3,3,1,2,0,0,nvvi,nvi,nv,nbv,nub1,nub2,
     &  x,y)
        call sm_ylabel('\\gD(U-B)')
        call sm_xlabel('U-B')
 
        call sm_gflush()
        call sm_hardcopy
        call sm_alpha
       endif        
      enddo 
      

900   close(25)
      stop
      end

      subroutine draw_compare(num_ref,nx,ny,nb1,nb2,nb3,nb4,
     &nvvi,nvi,nv,nbv,nub1,nub2,xl,yl)
      parameter (ns=10000,nd=100)
      common/dat/no(ns),dv(ns,5),nstar
      common/master/no_m(ns),dv_m(ns,5),nmaster
      integer nvvi,nx,ny,nb1,nb2,nb3,nb4,nvi,nv,nbv,nub1,nub2
      integer num_delta,ndelta
      real dx(nd),dy(nd),real x1,x2,y1,y2,pp(1),f,ylocate,xlocate
      character num_ref*4,label*28
      real average,stddev,xl,yl
      call sm_box(nb1,nb2,nb3,nb4)
      if(nstar.lt.30) then
       f=0.9
      elseif(nstar.lt.60) then
       f=0.8
      elseif(nstar.lt.100) then
       f=0.7
      else
       f=0.4
      endif
      pp(1)=403.+f
      call sm_ptype(pp,1)
      average=0.
      stddev=0.
      ndelta=0
      do i=1,nstar
       if(dv(i,nx).lt.99.and.dv(i,ny).lt.99) then
        do j=1,nmaster
         if(no(i).eq.no_m(j).and.dv_m(j,nx).lt.99.and.dv_m(j,ny).lt.99)
     &   then
          call sm_points(dv_m(j,nx),dv_m(j,ny)-dv(i,ny),1)
          average=average+dv_m(j,ny)-dv(i,ny)
          ndelta=ndelta+1
c          print *,nx,ny,ndelta,' ',num_ref
          goto 100
         endif
        enddo
       endif
100   enddo
      average=average/real(ndelta)
      do i=1,nstar
       if(dv(i,nx).lt.99.and.dv(i,ny).lt.99) then
        do j=1,nmaster
         if(no(i).eq.no_m(j).and.dv_m(j,nx).lt.99.and.dv_m(j,ny).lt.99)
     &   then
          stddev=stddev+(average-(dv_m(j,ny)-dv(i,ny)))**2.
          goto 200
         endif
        enddo
       endif
200   enddo
      if(ndelta.ge.2) then
       stddev=sqrt(stddev/real(ndelta-1))
       write(label,'(a7,f5.2,a5,f4.2,a2,i4,a1)') '\\gD = ', average,
     &  ' \\g+', stddev,' (',ndelta,')'
       call sm_relocate(xl,yl)
       call sm_expand(0.5)
       call sm_label(label)
       call sm_expand(1.)
      endif
      call sm_relocate(-5.,0.)
      call sm_draw(5.,0.)
      call sm_ctype('red')
      if(ny.eq.5.and.nvi.ne.0) then
       call read_delta(num_ref,nvi,2,dx,dy,n_delta)
       call sm_conn(dx,dy,n_delta)
      endif
      if(ny.eq.1.and.nx.eq.2.and.nv.ne.0) then
       call read_delta(num_ref,nv,3,dx,dy,n_delta)
       call sm_conn(dx,dy,n_delta)
      endif
      if(ny.eq.1.and.nx.eq.5.and.nvvi.ne.0) then
       call read_delta(num_ref,nvvi,1,dx,dy,n_delta)
       call sm_conn(dx,dy,n_delta)
      endif
      if(ny.eq.2.and.nbv.ne.0) then
       call read_delta(num_ref,nbv,4,dx,dy,n_delta)
       call sm_conn(dx,dy,n_delta)
      endif
      if(ny.eq.3.and.nx.eq.2.and.nub1.ne.0) then
       call read_delta(num_ref,nub1,5,dx,dy,n_delta)
       call sm_conn(dx,dy,n_delta)
      endif
      if(ny.eq.3.and.nx.eq.3.and.nub2.ne.0) then
       call read_delta(num_ref,nub2,6,dx,dy,n_delta)
       call sm_conn(dx,dy,n_delta)
      endif
      call sm_ctype('black')
        
      return
      end
      

      subroutine read_delta(num_ref,n_d,num_delta,dx,dy,ndelta)
      parameter (nd=100)
      integer num_delta,ndelta,n_d
      real dx(nd),dy(nd)
      character input*8,H*1,num_ref*4
      if(n_d.eq.0) then
       n_delta=0
       goto 100 
      endif 
      if(num_delta.eq.1) then
       write(input,'(a4,a4)') num_ref,'.vvi'
      elseif(num_delta.eq.2) then
       write(input,'(a4,a4)') num_ref,'.vi '
      elseif(num_delta.eq.3) then
       write(input,'(a4,a4)') num_ref,'.v  '
      elseif(num_delta.eq.4) then
       write(input,'(a4,a4)') num_ref,'.bv ' 
      elseif(num_delta.eq.5) then
       write(input,'(a4,a4)') num_ref,'.ub '
      elseif(num_delta.eq.6) then
       write(input,'(a4,a4)') num_ref,'.ub2'
      endif
      open(21,file=input)
      ndelta=0
      do i=1,nd+10
       read(21,'(t1,a1)',end=99) H
       if(H.ne.'#') then
        backspace(21)
        ndelta=ndelta+1
        read(21,*) dx(ndelta),dy(ndelta)
       endif
      enddo
99    close(21)
100   return
      end


      subroutine read_dat(input)
      parameter (ns=10000)
      character input*19,H*1
      common/dat/no(ns),dv(ns,5),nstar
      open(21,file=input)
      nstar=0
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.ne.'#') then
        nstar=nstar+1
        backspace(21)
        read(21,*) no(nstar),(dv(nstar,j),j=1,5)
       endif
      enddo
99    close(21)
      return
      end

      subroutine read_master(input)
      parameter (ns=10000)
      character input*19,H*1
      common/master/no_m(ns),dv_m(ns,5),nmaster
      open(21,file=input)
      nmaster=0
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.ne.'#') then
        nmaster=nmaster+1
        backspace(21)
        read(21,*) no_m(nmaster),(dv_m(nmaster,j),j=1,5)
       endif
      enddo
99    close(21)
      return
      end
