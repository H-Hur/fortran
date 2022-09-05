c for only fortran77 to draw color-magnitude daigram.
c f77 -o 'cor_red' 'cor_red.f' -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11
c g77 -o cor_red ~/work/fortran/subroutine/webda_tools/cor_red.f -L/Applications/scisoft/i386/Packages/sm/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11 -L/Users/apple/utils/sm_mac/scisoft/lib -laquaterm

      real Av, dist, ebv, evi,Rv_fg, Rv_cl, ebv_fg
      character output*8,H*1,output2*9
      open(21,file='distance.par')
      read(21,'(t1,a1)') H
      read(21,*) dist, ebv, Rv_fg, Rv_cl, ebv_fg
      Av=ebv_fg*Rv_fg+(ebv-ebv_fg)*Rv_cl
      evi=ebv_fg*Rv_fg/2.45+(ebv-ebv_fg)*Rv_cl/2.45
      write(*,'(a14,f6.3)') 'mean E(V-I) = ', evi
      write(*,'(a10,f6.3)') 'mean Av = ', Av

      close(21)

      print *,'reading'
      call read_combine()
      call read_ZAMS()
      print *,'drawing U-B:B-V'

      call draw_UB_BV(ebv_fg,ebv)
      print *,'drawing Rv_cl'

      call draw_rvcl(Rv_fg, Rv_cl, ebv_fg,0.)
      print *,'drawing Rv_cross'

      call draw_rv_cross(Rv_fg, Rv_cl, ebv_fg,0.)
      print *,'drawing distance histogarm'
      call write_distance(Rv_fg, Rv_cl,dist,0.)
      print *,'cal_Rvcl'
      call cal_Rvcl()
      print *,'drawing cmd'
      write(output,'(a8)') 'V.B-V.ps'
      call draw_cmd(2,1,ebv,ebv_fg,evi,av,dist,output,-0.4,2.5,20.,3.)
      write(output,'(a8)') 'V.V-I.ps'
      call draw_cmd(4,1,ebv,ebv_fg,evi,av,dist,output,-0.4,4.5,20.,3.)
      write(output,'(a8)') 'V.U-B.ps'
      call draw_cmd(3,1,ebv,ebv_fg,evi,av,dist,output,-1.4,3.0,20.,3.)


      print *,'drawing cmd0'
      write(output2,'(a9)') 'Mv.BV0.ps'
      call draw_cmd0(2,1,dist,output2)
      write(output2,'(a9)') 'Mv.VI0.ps'
      call draw_cmd0(4,1,dist,output2)
     
      stop
      end

      subroutine cal_Rvcl()
      parameter (ns=100000)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      common/red/vM(ns),bv0(ns),vi0(ns),ub0(ns),Av(ns),nearly(ns)
      real*8 x(ns),y(ns),slopeout,yout,ebv_fg
      integer nleast
      nleast=0
      do i=1,nstar
       if(nearly(i).eq.1.and.p(i).ge.0.7.and.nobs(i,4).ne.0) then
        nleast=nleast+1
        x(nleast)=com(i,2)-bv0(i) 
        y(nleast)=com(i,4)-vi0(i)
       endif
      enddo     
      call cal_least(x,y,nleast,1,0.d0,0.d0,slopeout,yout)
      if(yout.ge.0.05) then
       call cal_least(x,y,nleast,2,0.d0,0.d0,slopeout,yout)
      endif
      ebv_fg=yout/(slopeout-3.1d0)
      write(*,'(a10,f6.3)') 'E(B-V)fg= ', ebv_fg
      write(*,'(a7,f6.3)') 'Rv_cl= ', slopeout*2.45d0
      print *,'n stars to compute Rv = ', nleast
      return
      end


      subroutine draw_cmd0(nx,ny,dist,output)
      parameter (ns=100000)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      common/red/vM(ns),bv0(ns),vi0(ns),ub0(ns),Av(ns),nearly(ns)
      character xlabel*9,ylabel*5,output*9,ps*22
      real pp(1),x,vM_r
      if(nx.eq.2) then
       write(xlabel,'(a9)') '(B-V)\\d0'
      elseif(nx.eq.3) then
       write(xlabel,'(a9)') '(U-B)\\d0'
      elseif(nx.eq.4) then
       write(xlabel,'(a9)') '(V-I)\\d0'
      endif
      if(ny.eq.1) then
       write(ylabel,'(a5)') 'M\\dV'
      endif
      write(ps,'(a13,a9)') 'postportfile ',output
      call sm_device (ps)
      call sm_graphics
      call sm_erase()
      call sm_expand(1.)
      call sm_location(3000,30000,3000,30000)
      call sm_limits(-1.,2.,7.,-8.)
      call sm_xlabel(xlabel)
      call sm_ylabel(ylabel)
      call sm_box(1,2,0,0)
      do i=1,nstar
       if(p(i).ge.0.7) then
        pp(1)=403.99
        call sm_ctype('blue')
        if(nearly(i).eq.1) then
         call sm_ctype('black')
        endif
        if(p(i).gt.1.) then
         call sm_ctype('green')
        endif
       elseif(p(i).ge.0.5) then
        pp(1)=403.99
        call sm_ctype('cyan')
       elseif(p(i).ge.0.2) then
        pp(1)=403.7
        call sm_ctype('magenta')
       elseif(p(i).ge.0.2) then
        pp(1)=403.7
        call sm_ctype('red')
       else
        pp(1)=400.1
        call sm_ctype('black')
       endif
       if(nearly(i).eq.-1) then
        pp(1)=41.+0.99
        call sm_ctype('black')
       endif
       call sm_ptype(pp,1)
       if(nobs(i,nx).ne.0.and.nobs(i,ny).ne.0.and.nearly(i).ne.0.
     & and.Av(i).ne.0.) then
        if(nx.eq.2) then
         x=bv0(i)
        elseif(nx.eq.3) then
         x=ub0(i)
        elseif(nx.eq.4) then
         x=vi0(i)
        endif
        vM_r=com(i,1)-Av(i)-dist
        call sm_points(x,vM_r,1)
       endif
      enddo

      call draw_zams(nx,ny,0.,0.,3.,0)
      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      return
      end

      subroutine draw_cmd(nx,ny,ebv,ebv_fg,
     &evi,av_m,dist,output,x1,x2,y1,y2)
      parameter (ns=100000)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      common/red/vM(ns),bv0(ns),vi0(ns),ub0(ns),Av(ns),nearly(ns)
      character xlabel*3,ylabel*1,output*8,ps*21
      real x1,x2,y1,y2,dist,pp(1),ebv_fg,av_fg
      if(nx.eq.2) then
       write(xlabel,'(a3)') 'B-V'
      elseif(nx.eq.3) then
       write(xlabel,'(a3)') 'U-B'
      elseif(nx.eq.4) then
       write(xlabel,'(a3)') 'V-I'
      endif
      if(ny.eq.1) then
       write(ylabel,'(a1)') 'V'
      endif
      write(ps,'(a13,a8)') 'postportfile ',output
      call sm_device (ps)
      call sm_graphics
      call sm_erase()
      call sm_expand(1.)
      call sm_location(3000,30000,3000,30000)
      call sm_limits(x1,x2,y1,y2)
      call sm_xlabel(xlabel)
      call sm_ylabel(ylabel)
      call sm_box(1,2,0,0)
      do i=1,nstar
       if(p(i).ge.0.7) then
        pp(1)=403.99
        call sm_ctype('blue')
        if(nearly(i).eq.1) then
         call sm_ctype('black')
        endif
        if(p(i).gt.1.) then
         call sm_ctype('green')
        endif
       elseif(p(i).ge.0.5) then
        pp(1)=403.99
        call sm_ctype('cyan')
       elseif(p(i).ge.0.2) then
        pp(1)=403.7
        call sm_ctype('magenta')
       elseif(p(i).ge.0.2) then
        pp(1)=403.7
        call sm_ctype('red')
       else
        pp(1)=400.1
        call sm_ctype('black')
       endif
       if(nearly(i).eq.-1) then
        pp(1)=41.+0.99
        call sm_ctype('black')
       endif
       call sm_ptype(pp,1)
       if(nobs(i,nx).ne.0.and.nobs(i,ny).ne.0) then
        call sm_points(com(i,nx),com(i,ny),1)
       endif
      enddo
      av_fg=ebv_fg*3.1
      if(nx.eq.2) then
       call draw_zams(nx,ny,ebv,av_m+dist,3.,0)
       call draw_zams(nx,ny,ebv_fg,av_fg+dist,3.,1)
      elseif(nx.eq.3) then
       call draw_zams(nx,ny,ebv*0.72,av_m+dist,3.,0)
       call draw_zams(nx,ny,ebv_fg*0.72,av_fg+dist,3.,1)
      elseif(nx.eq.4) then
       call draw_zams(nx,ny,evi,av_m+dist,3.,0)
       call draw_zams(nx,ny,ebv_fg*3.1/2.45,av_fg+dist,3.,1)
      endif
      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      return
      end

      subroutine read_combine()
      parameter (ns=100000)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      character H*1
      open(21,file='combine.dat')
      read(21,'(t1,a1)') H
      do i=1,ns
       read(21,*,end=99) no(i),(com(i,j),j=1,4),(nobs(i,j),j=1,4),p(i)
      enddo
99    nstar=i-1
      close(21)
      return
      end

      subroutine write_distance(Rv_fg,Rv_cl,dist_in,ebv_fg)
      parameter (ns=100000)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      common/red/vM(ns),bv0(ns),vi0(ns),ub0(ns),Av(ns),nearly(ns)
      real V0,Rv,ebv_cl,ebv,Rv_fg,Rv_cl,ebv_fg,dist,dist_mean50
      real hist(40),dist_mean,hist70(40),hist50(40),x1,x2,y1,y2,hmax
      real dist_in
      integer ne
      character label*16
      open(21,file='distance70.dat')
      write(21,'(t1,a80)') '#   No   Mv    V0     Av    dist    mp   E(B
     &-V) E(U-B) E(V-I)   nobs                                        '
      ne=0
      dist_mean=0.
1     format(t1,i6,8f7.3,4i3)
      do i=1,ns
       if(nearly(i).eq.1) then
        ne=ne+1
        ebv=com(i,2)-bv0(i)
        if(ebv_fg.ne.0.) then
         if(ebv_fg.lt.ebv) then
          ebv_cl=ebv-ebv_fg
          Av(i)=Rv_fg*ebv_fg + Rv_cl * ebv_cl
         else
          Av(i)=Rv_fg*ebv
         endif
        else
         Av(i)=Rv_cl*ebv
        endif  
         dist=com(i,1)-Av(i)-vM(i)
         dist_mean=dist_mean+dist
        write(21,1) no(i),vM(i),com(i,1)-Av(i),Av(i),dist,p(i),
     &  com(i,2)-bv0(i),com(i,3)-ub0(i),com(i,4)-vi0(i),
     &   (nobs(i,j),j=1,4)
       endif
      enddo
      dist_mean=dist_mean/real(ne)
      print *,'Mean distance modulus of mp>70 = ', dist_mean
      close(21)
      do i=1,40
       hist(i)=dist_mean+0.2*real(i)-2.
       hist70(i)=0.
       hist50(i)=0.
      enddo
      
      open(22,file='distance50.dat')
      write(22,'(t1,a80)') '#   No   Mv    V0     Av    dist    mp    (B
     &-V)0 (U-B)0 (V-I)0  nobs                                        '
      ne=0
      dist_mean50=0.
      do i=1,ns
       if(nearly(i).eq.1.or.nearly(i).eq.2) then
        ne=ne+1
        ebv=com(i,2)-bv0(i)
        if(ebv_fg.lt.ebv) then
         ebv_cl=ebv-ebv_fg
         Av(i)=Rv_fg*ebv_fg + Rv_cl * ebv_cl
        else
         Av(i)=Rv_fg*ebv
        endif
        dist=com(i,1)-Av(i)-vM(i)
        dist_mean50=dist_mean50+dist
        write(22,1) no(i),vM(i),com(i,1)-Av(i),Av(i),dist,p(i),bv0(i),
     &  ub0(i),vi0(i), (nobs(i,j),j=1,4)
       endif
       do j=1,19
        if(dist.gt.hist(j).and.dist.le.hist(j+1)) then
         if(nearly(i).eq.1) then
          hist70(j)=hist70(j)+1.
          hist50(j)=hist50(j)+1.
         elseif(nearly(i).eq.2) then
          hist50(j)=hist50(j)+1.
         endif
        endif
       enddo
      enddo
      dist_mean50=dist_mean50/real(ne)
      print *,'Mean distance modulus of mp>50 = ', dist_mean50
      close(22)
      hmax=0.
      do i=1,20
       hist(i)=hist(i)+0.1
       if(hist50(i).gt.hmax) then
        hmax=hist50(i)
       endif
       if(hist(i)+0.1.lt.dist_in) then!.and.hist(i+1).gt.dist_in) then
        dist_h=hist50(i+1)
       endif
      enddo
      x1=hist(1)-0.1
      x2=hist(20)+0.1
      y1=0.
      y2=hmax*1.3
      if((x1.eq.0.and.x2.eq.0).or.(y1.eq.0.and.y2.eq.0)) then
       goto 900
      endif
      call sm_device ('postportfile distance.ps')
      call sm_graphics
      call sm_erase()
      call sm_ltype(0)
      call sm_limits(x1,x2,y1,y2)
      call sm_ctype('black')
      call sm_expand(1.)
      call sm_location(3000,25000,15000,26000)
      call sm_box(1,2,0,0)
      call sm_xlabel('V\\d0-M\\dv')
      call sm_ylabel('N')
      call sm_ltype(2)
      call sm_histogram(hist,hist50,40)
      call sm_ltype(0)     
      call sm_histogram(hist,hist70,40)
      call sm_relocate(dist_in,dist_h)
      call sm_draw(dist_in,dist_h+0.15*(y2-y1))
      write(label,'(a12,f4.1)') 'V\\d0-M\\dV=',dist_in
      call sm_relocate(dist_in-0.05*(x2-x1),dist_h+0.17*(y2-y1))
      call sm_expand(0.7)
      call sm_label(label)
      call sm_ptype(33.99,1)
      call sm_expand(3.)
      call sm_angle(180.)
      call sm_points(dist_in,dist_h+0.03*(y2-y1),1)
      call sm_angle(0.)
      call sm_expand(1.)
      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
900   return
      end

      subroutine draw_rv_cross(Rv_fg,Rv_cl,ebv_fg,del_ebv_fg)
      parameter (ns=100000)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      common/red/vM(ns),bv0(ns),vi0(ns),ub0(ns),Av(ns),nearly(ns)
      real pp(1),evi,ebv,Rv_fg,Rv_cl,ebv_fg,del_ebv_fg,evi_fg,del_evi_fg
      real Rv_cl2,xlocate,ylocate
      character label*20
      call sm_device ('postportfile Rv_cross.ps')
      call sm_graphics
      call sm_erase()
      call sm_ltype(0)
      call sm_limits(0.,2.2,0.,3.5)
      call sm_ctype('black')
      call sm_expand(1.)
      call sm_location(3000,25000,9000,25000)
      call sm_box(1,2,0,0)
      call sm_xlabel('E(B-V)')
      call sm_ylabel('E(V-I)')
      evi_fg=Rv_fg*ebv_fg/2.45
      del_evi_fg=Rv_fg*del_ebv_fg/2.45
      do i=1,nstar
       if(nearly(i).ne.0.and.nobs(i,4).ne.0) then
        evi=com(i,4)-vi0(i)
        ebv=com(i,2)-bv0(i)
        Rv_cl2=2.45*(evi-evi_fg-del_evi_fg)/(ebv-ebv_fg-del_ebv_fg)
        if(p(i).ge.0.7) then
         pp(1)=403.99
         call sm_ctype('blue')
         if(nearly(i).eq.1) then
          call sm_ctype('black')
         endif
         if(p(i).gt.1.) then
          call sm_ctype('green')
         endif
        elseif(p(i).ge.0.5) then
         pp(1)=403.99
         call sm_ctype('cyan')
        elseif(p(i).ge.0.2) then
         pp(1)=403.7
         call sm_ctype('magenta')
        elseif(p(i).ge.0.2) then
         pp(1)=403.7
         call sm_ctype('red')
        else
         pp(1)=400.3
         call sm_ctype('black')
        endif
        if(nearly(i).eq.-1) then
         pp(1)=41.+0.99
         call sm_ctype('black')
        endif
        call sm_ptype(pp,1)
        if(p(i).ge.0.5) then
         call sm_points(ebv,evi,1)
        endif
       endif
      enddo

      if(Rv_cl.ne.Rv_fg) then
       call sm_ctype('red')
       call sm_relocate(0.,0.)
       call sm_draw(2.2,2.7837)
       call sm_relocate(1.6,1.8)
       call sm_expand(0.9)
       write(label,'(a14,f3.1,a3)') 'R\\dV\\d,\\df\\dg=',Rv_fg,'      '
       call sm_label(label)
       call sm_expand(1.)
      endif

      ylocate=(Rv_cl/2.45)*(2.2-ebv_fg)+evi_fg
      xlocate=2.2
      call sm_relocate(ebv_fg,evi_fg)
      call sm_ctype('black')
      call sm_draw(xlocate,ylocate)
      call sm_ctype('black')

      ylocate=(Rv_cl/2.45)*(2.2-ebv_fg-del_ebv_fg)+evi_fg+del_evi_fg
      xlocate=2.2
      call sm_ltype(3)
      call sm_relocate(ebv_fg+del_ebv_fg,evi_fg+del_evi_fg)
      call sm_ctype('black')
      call sm_draw(xlocate,ylocate)
      call sm_ctype('black')
      call sm_ltype(0)

      call sm_relocate(ebv_fg,evi_fg)
      call sm_draw(ebv_fg,-0.2)
      call sm_relocate(ebv_fg-0.3,-0.13)
      call sm_expand(0.9)
      write(label,'(a15,f4.2)') 'E(B-V)\\df\\dg=',ebv_fg
      if(ebv_fg.ne.0.) then
       call sm_label(label)
      endif
      call sm_expand(1.)
      call sm_ptype(30.99,1)
      call sm_points(ebv_fg,evi_fg-0.02,1)

      write(label,'(a14,f3.1,a3)') 'R\\dV\\d,\\dc\\dl=',Rv_cl,'      '
      call sm_relocate(1.15,3.)
      call sm_expand(0.9)
      call sm_ctype('black')
      call sm_label(label)
 
      call sm_expand(1.)
      call sm_ctype('black')
      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      return
      end

      subroutine draw_rvcl(Rv_fg,Rv_cl,ebv_fg,del_ebv_fg)
      parameter (ns=100000)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      common/red/vM(ns),bv0(ns),vi0(ns),ub0(ns),Av(ns),nearly(ns)
      real pp(1),evi,ebv,Rv_fg,Rv_cl,ebv_fg,del_ebv_fg,evi_fg,del_evi_fg
      character label*21
      call sm_device ('postportfile Rv_cl.ps')
      call sm_graphics
      call sm_erase()
      call sm_ltype(0)
      call sm_limits(0.,2.,0.,4.)
      call sm_ctype('black')
      call sm_expand(1.)
      call sm_location(3000,25000,9000,25000)
      call sm_box(1,2,0,0)
      call sm_xlabel('E(B-V)\\dc\\dl')
      call sm_ylabel('E(V-I)\\dc\\dl/E(B-V)\\dc\\dl')
      evi_fg=Rv_fg*ebv_fg/2.45
      del_evi_fg=Rv_fg*del_ebv_fg/2.45
      do i=1,nstar
       if(nearly(i).ne.0.and.nobs(i,4).ne.0.and.p(i).ge.0.5) then
        evi=com(i,4)-vi0(i)
        ebv=com(i,2)-bv0(i)
        Rv_cl2=2.45*(evi-evi_fg-del_evi_fg)/(ebv-ebv_fg-del_ebv_fg)
        if(p(i).ge.0.7) then
         pp(1)=403.99
         call sm_ctype('blue')
         if(nearly(i).eq.1) then
          call sm_ctype('black')
         endif
         if(p(i).gt.1.) then
          call sm_ctype('green')
         endif
        elseif(p(i).ge.0.5) then
         pp(1)=403.99
         call sm_ctype('cyan')
        elseif(p(i).ge.0.2) then
         pp(1)=403.7
         call sm_ctype('magenta')
        elseif(p(i).ge.0.2) then
         pp(1)=403.7
         call sm_ctype('red')
        else
         pp(1)=400.3
         call sm_ctype('black')
        endif
        if(nearly(i).eq.-1) then
         pp(1)=41.+0.99
         call sm_ctype('black')
        endif
        call sm_ptype(pp,1)
        call sm_points(ebv-ebv_fg,(evi-evi_fg)/(ebv-ebv_fg),1)
       endif
      enddo
      call sm_ctype('black')
      call sm_relocate(0.,Rv_fg/2.45)
      call sm_ctype('red')
      call sm_draw(2.,Rv_fg/2.45)
      call sm_relocate(0.,Rv_cl/2.45)
      call sm_ctype('black')
      call sm_draw(2.,Rv_cl/2.45)
      call sm_ctype('black')
      call sm_expand(0.9)
      call sm_relocate(1.3,Rv_cl/2.45+0.1)
      call sm_ctype('black')
      write(label,'(a18,f3.1)') 'R\\dV\\d,\\dc\\dl=', rv_cl
      call sm_label(label)
      if(Rv_cl.ne.Rv_fg) then
       call sm_ctype('red')
       call sm_relocate(1.6,Rv_fg/2.45+0.1)
       call sm_label('R\\dV\\d,\\df\\dg=3.1')
       call sm_expand(1.)
       call sm_ctype('black')
      endif

      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      return
      end


      subroutine draw_UB_BV(ebv_fg,ebv_in)
      parameter (ns=100000)
      common/combine/no(ns),com(ns,4),nobs(ns,4),p(ns),nstar
      common/red/vM(ns),bv0(ns),vi0(ns),ub0(ns),Av(ns),nearly(ns)
      common/zams/zam(41,4)
      real pp(1),ebv,eub,q,ubz,zbv(41),zub(41),ebv_fg,ebv_in
      integer nebv,no_ex(ns),nex,no_in(ns),nin
      character label*18,H*1
      do i=1,41
       zbv(i)=zam(i,2) 
       zub(i)=zam(i,3)
      enddo
      open(25,file='except.dat')
      read(25,'(t1,a1)',end=9) H
9     nex=0
      do i=1,ns
       read(25,*,end=10) no_ex(i) 
      enddo
10    nex=i-1
11    close(25)

      open(25,file='include.dat')
      read(25,'(t1,a1)',end=19) H
19     nin=0
      do i=1,ns
       read(25,*,end=20) no_in(i)
      enddo
20    nin=i-1
21    close(25)
      call sm_device ('postlandfile U-B_B-V.ps')
      call sm_graphics
      call sm_erase()
      call sm_limits(-0.5,2.,2.,-1.2)
      call sm_ctype('black')
      call sm_expand(1.)
      call sm_location(3000,30000,3000,30000)
      call sm_box(1,2,0,0)
      call sm_xlabel('B-V')
      call sm_ylabel('U-B')
      open(21,file='reddening.dat')
      write(21,'(t1,a60)') '#   No    V      B-V     U-B     V-I    nobs
     &        mp                                                     '
      ebv=0.
      nebv=0
1     format(t1,i6,4f8.3,4i3,f6.2)
      do i=1,nstar
       nearly(i)=0 
       do j=1,nex
        if(no(i).eq.no_ex(j)) then
         call sm_ctype('black')
         pp(1)=41.+0.9
         call sm_lweight(1.5)
         call sm_ptype(pp,1)
         call sm_expand(1.5)
         call sm_points(com(i,2),com(i,3),1)
         call sm_lweight(0.5)
         nearly(i)=-1
         goto 100
        endif
       enddo
       do j=1,nin
        if(no(i).eq.no_in(j).and.p(i).lt.0.) then
         p(i)=2.
        endif
       enddo
       if(nobs(i,2).ne.0.and.nobs(i,3).ne.0) then
        call lint(com(i,2),zbv,zub,ubz,41,41)
        q=com(i,3)-0.72*com(i,2)
        if((q.le.-0.5.or.(com(i,3).lt.ubz-0.1.and.com(i,2).lt.0.5))) 
     &  then
         nearly(i)=3
         call cor_red(com(i,2),com(i,3),5,bv0(i),ub0(i),vi0(i))
         call cor_vi(bv0(i),vi0(i))
         call cor_Mv(bv0(i),vM(i))
        endif
        if(p(i).ge.0.7) then
         if((q.le.-0.5.or.(com(i,3).lt.ubz-0.1.and.com(i,2).lt.0.5)))
     &   then
          nearly(i)=1
          ebv=ebv+(com(i,2)-bv0(i))
          nebv=nebv+1
          write(21,1) no(i),(com(i,j),j=1,4)
     &    ,(nobs(i,j),j=1,4),p(i)
         endif
         call sm_ctype('blue')
         pp(1)=403.+0.7
         call sm_expand(1.5)
         if(nearly(i).eq.1) then
          call sm_ctype('black')
         endif
         if(p(i).gt.1..and.nearly(i).eq.1) then
          call sm_ctype('green')
          pp(1)=43.+0.7
         endif
        elseif(p(i).ge.0.5) then
         call sm_ctype('cyan')
         pp(1)=403.+0.7
         call sm_expand(1.5)
         if((q.le.-0.5.or.(com(i,3).lt.ubz-0.1.and.com(i,2).lt.0.5)))
     &   then
          nearly(i)=2
         endif
        elseif(p(i).ge.0.2) then
         call sm_ctype('magenta')
         pp(1)=403.+0.7
        elseif(p(i).ge.0.) then
         call sm_ctype('red')
         pp(1)=403.+0.7
        else
         pp(1)=400.+0.3
         call sm_ctype('black')
        endif
        call sm_ptype(pp,1)
        call sm_points(com(i,2),com(i,3),1)
       endif
100    call sm_expand(1.)

      enddo
      ebv=ebv/real(nebv) 
      call draw_zams(2,3,0.,0.,3.,0)
      if(ebv.ne.0.)then
       eub=0.72*ebv
       call draw_zams(2,3,ebv,eub,2.,1)
       write(label,'(a14,f4.2)') 'mean E(B-V) = ',ebv
       call sm_relocate(1.5,-0.6)
       call sm_label(label)
       print *,'mean E(B-V) = ',ebv
      endif
      if(ebv_fg.gt.0.) then
       eub=0.72*ebv
       call draw_zams(2,3,ebv_fg,eub,2.,1)
      endif
      call sm_ctype('blue')
      call draw_zams(2,3,ebv_in,ebv_in*0.72,2.,1)
      call sm_ctype('black')


      close(21)

900   call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      return
      end

      subroutine draw_zams(nx,ny,dx,dy,weight,ltype)
      common /zams/zam(41,4)
      integer nx,ny,ltype
      real dx,dy,x(41),y(41),weight
      call sm_lweight(weight)
      call sm_ltype(ltype)
      call sm_ctype('black')
      do i=1,41
       x(i)=zam(i,nx)+dx
       y(i)=zam(i,ny)+dy
      enddo
      call sm_conn(x,y,41)
      call sm_lweight(0.5)
      call sm_ltype(0)
      return
      end

      subroutine read_ZAMS()
      character H*1
      common /zams/zam(41,4)
      open(31,file='Sung2013.ZAMS.cc')
      read(31,'(a1)') H
      do i=1,41
       read(31,*,end=100) zam(i,1),zam(i,4),zam(i,2),zam(i,3)
      enddo
100   nzams=i-1
      close(31)
      return
      end

c      subroutine cor_red(sBV,sUB,unredx,unredy)
c      common /zams/zam(41,4)
c      real fx,fy,dy,unredx,unredy,sBV,sUB,slope,eBV,eUB,Rv
c      real zUB(41),zBV(41)
c      do i=1,41
c       zBV(i)=zam(i,2)
c       zUB(i)=zam(i,3)
c      enddo
c      slope=0.72
c      Rv=3.1
c      call lint(sUB,zUB,zBV,fx,41,41)
c      dx=fx-sBV
c      call lint(sBV,zBV,zUB,fy,41,41)
c      dy=sUB-fy
c      unredx=sBV
c      unredy=sUB
c      do j=1,1000
c       call lint(unredy,zUB,zBV,fx,41,41)
c       dx=fx-unredx
c       dy=dx*slope +dx*dx*0.025
c       unredx=unredx+dx
c       unredy=unredy+dy
c        if(abs(dx).lt.0.00001.and.abs(dy).lt.0.00001) then
c        goto 100
c       endif
c100    eBV=sBV-unredx
c      enddo
c      return
c      end
      subroutine cor_red(bv,ub,lclass,bv0,ub0,vi0)
      common /zams/zam(41,4)
c      common /zams/zMv(41),zVI(41),zBV(41),zUB(41),zhvi(18),zrha(18)
      common/CC_relation/bvref(46),ub_ms(46),ub_gt(46),ub_bg(46),
     &      ub_1b(46),ub_1a(46),vi_ms(46),vi_gt(46),ri_ms(46),
     &      ref_vi(30),ref_ms(30),ref_gt(30)
      real bv,ub,bv0,ub0,vi0,s,ss,bv1,ub1,ebv_sum,zUB(41),zBV(41)
      integer lclass
      do i=1,41
       zBV(i)=zam(i,2)
       zUB(i)=zam(i,3)
      enddo
      s=0.72
      ss=0.025
      ebv_sum=0.
      bv1=bv
      ub1=ub
100   if(lclass.eq.5.or.lclass.eq.6) call lint(ub1,zUB,zBV,bv0,41,41)
      if(lclass.eq.4) call lint(ub1,zUB,zBV,bv0,41,41)
      if(lclass.eq.3) call lint(ub1,ub_gt,bvref,bv0,46,46)
      if(lclass.eq.2) call lint(ub1,ub_bg,bvref,bv0,46,46)
      if(lclass.eq.1) call lint(ub1,ub_1b,bvref,bv0,46,46)
      if(lclass.eq.0) call lint(ub1,ub_1a,bvref,bv0,46,46)
      ebv_sum=ebv_sum+(bv1-bv0)
      ub0=ub-s*ebv_sum-ss*ebv_sum*ebv_sum
      bv0=bv-ebv_sum
      if((ub0.le.zub(1).or.bv0.le.zbv(1))) goto 500
      if(abs(ub0-ub1).gt.0.0001.or.bv1-bv0.gt.0.) then
       ub1=ub0
       bv1=bv0
       goto 100
      endif
500   if(lclass.eq.5.or.lclass.eq.6) call lint(bv0,zBV,zVI,vi0,41,41)
      if(lclass.eq.4) call lint(bv0,zBV,zVI,vi0,41,41)
      if(lclass.eq.3) call lint(bv0,bvref,vi_gt,vi0,46,46)
      if(lclass.le.2) call lint(bv0,zBV,zVI,vi0,41,41)
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

      subroutine cor_Mv(bv0,Mv0)
      common /zams/zam(41,4)
      real zMv(41),zBV(41)
      real bv0,Mv0
      do i=1,41
       zBV(i)=zam(i,2)
       zMv(i)=zam(i,1)
      enddo
      call lint(bv0,zBV,zMv,Mv0,41,41)
      return
      end

      subroutine cor_vi(bv0,vi0)
      common /zams/zam(41,4)
      real zVI(41),zBV(41)
      real bv0,vi0
      do i=1,41
       zBV(i)=zam(i,2)
       zVI(i)=zam(i,4)
      enddo
      call lint(bv0,zBV,zVI,vi0,41,41)
      return
      end

      subroutine cal_least(x,y,n,ntype,slopein,yin,slopeout,yout)
      parameter (ns=100000)
      real*8 x(ns),y(ns),sx,sy,sxx,sxy,slopein,yin,slopeout,yout
      real*8 sxslopein
      integer n,ntype
      sx=0.d0
      sy=0.d0
      sxx=0.d0
      sxy=0.d0
      sxslopein=0.d0
      do i=1,n
       sx=sx+x(i)
       sxx=sxx+x(i)**2.d0
       sxy=sxy+x(i)*y(i)
       sy=sy+y(i)
       sxslopein=sxslopein+(x(i)-slopein)**2.d0
      enddo
      if(ntype.eq.1) then
       slopeout=(real(n)*sxy-sx*sy)/(real(n)*sxx-sx**2.d0)
       yout=(sxx*sy-sx*sxy)/(real(n)*sxx-sx**2.d0)
      elseif(ntype.eq.2) then
       slopeout=sy/sx
       yout=0.d0
      elseif(ntype.eq.3) then
       slopeout=(sxy-yin*sx)/sxx
       yout=yin
      elseif(ntype.eq.4) then
       slopeout=(sxy-slopein*sy-yin*sx+real(n)*slopein*yin)/sxslopein
       yout=yin-slopeout*slopein
      elseif(ntype.eq.5) then
       yout=(sy-slopein*sx)/real(n)
       slopeout=slopein
      endif
      return
      end
      block data int_col
      common/CC_relation/bvref(46),ub_ms(46),ub_gt(46),ub_bg(46),
     &      ub_1b(46),ub_1a(46),vi_ms(46),vi_gt(46),ri_ms(46),
     &      ref_vi(30),ref_ms(30),ref_gt(30)
        data bvref/-.325,-0.30,-.275,-0.25,-.225,-0.20,-.175,-0.15,
     &             -.125,-0.10,-.075,-0.05,-.025, 0.00, 0.05, 0.10,
     &              0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
     &              0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90,
     &              0.95, 1.00, 1.05, 1.10, 1.15, 1.20, 1.25, 1.30,
     &              1.35, 1.40, 1.45, 1.50, 1.55, 1.60/
        data ub_ms/-1.18,-1.06,-0.99,-0.89,-0.79,-0.70,-0.60,-0.50,
     &             -0.40,-0.30,-0.20,-0.11,-0.05, 0.00,0.055,0.085,
     &             0.095,0.085, 0.06,0.035,0.005,-0.01,-0.02, 0.00,
     &              0.04, 0.09, 0.16, 0.23, 0.32, 0.41, 0.53, 0.65,
     &              0.76, 0.86, 0.95, 1.04, 1.08, 1.13, 1.17, 1.20,
     &              1.22, 1.22, 1.21, 1.17, 1.15, 1.19/
        data ub_gt/-1.20,-1.10,-1.00,-0.90,-0.79,-0.68,-0.57,-0.47,
     &             -0.37,-0.28,-0.20,-0.11,-.055, 0.00,0.065, 0.11,
     &              0.12, 0.11,0.095,0.075, 0.07, 0.07,0.075, 0.09,
     &              0.12, 0.16, 0.21, 0.28, 0.35, 0.44, 0.53, 0.62,
     &              0.72, 0.84, 0.93, 1.04, 1.14, 1.24, 1.34, 1.44,
     &              1.54, 1.64, 1.73, 1.81, 1.87, 1.89/
        data ub_bg/-1.19,-1.12,-1.06,-1.00,-0.94,-0.87,-0.78,-0.68,
     &             -0.57,-0.47,-0.38,-0.29,-0.21,-0.15,-0.03, 0.05,
     &              0.10, 0.12, 0.13, 0.14,0.155, 0.17,0.185, 0.21,
     &              0.23, 0.25, 0.28, 0.31, 0.35, 0.41, 0.50, 0.62,
     &              0.72, 0.81,0.905, 1.00, 1.10, 1.20, 1.29, 1.38,
     &              1.48, 1.58, 1.67, 1.77, 1.86, 1.96/
        data ub_1b/-1.18,-1.15,-1.12,-1.08,-1.04,-1.00,-0.95,-0.88,
     &             -0.80,-0.74,-0.66,-0.58,-0.48,-0.35,-0.15, 0.01,
     &              0.08, 0.16, 0.23, 0.27, 0.30,0.325, 0.35,0.375,
     &              0.40, 0.42, 0.43, 0.45, 0.47, 0.49, 0.55, 0.63,
     &              0.72, 0.81, 0.90, 0.98, 1.06, 1.15, 1.24, 1.33,
     &              1.42, 1.51, 1.60, 1.68, 1.77, 1.86/
        data ub_1a/-1.18,-1.15,-1.12,-1.09,-1.05,-1.01,-0.97,-0.92,
     &             -0.86,-0.81,-0.76,-0.71,-0.65,-0.59,-0.33,-0.06,
     &              0.06, 0.15, 0.23,0.275, 0.30, 0.33, 0.35,0.375,
     &              0.40, 0.42, 0.43, 0.45, 0.47, 0.49, 0.55, 0.63,
     &              0.72, 0.81, 0.90, 0.98, 1.06, 1.15, 1.24, 1.33,
     &              1.42, 1.51, 1.60, 1.68, 1.77, 1.86/
        data vi_ms/-.355,-0.33,-.305,-0.28,-0.25,-0.22,-0.19,-0.16,
     &             -0.13,-0.11,-0.08,-.056,-0.03,-.003, .054, 0.11,
     &              0.17, 0.23,0.294,0.353, 0.41, 0.47,0.525,0.575,
     &              0.62, 0.66, 0.70, 0.74, 0.79, 0.84, 0.89, 0.95,
     &             1.013, 1.08, 1.15, 1.22,1.293,1.375,1.464,1.555,
     &              1.66,1.775, 1.95, 2.25, 2.43, 2.60/
        data vi_gt/-0.34,-0.32,-.304,-0.27,-0.23,-0.20,-0.17,-0.14,
     &             -0.11,-0.08,-.053,-.026, 0.00, .027, 0.08, 0.13,
     &              0.19,0.244, 0.30,0.356, 0.41,0.467, 0.52, 0.57,
     &              0.62, 0.67, 0.72, 0.76, 0.80, 0.83, 0.86, 0.89,
     &              .923, 0.96, 1.00,1.045, 1.09, 1.14,1.195,1.253,
     &             1.314,1.386,1.473, 1.57, 1.68, 1.80/
        data ri_ms/-0.18,-0.18,-0.17,-0.15,-0.12,-0.11,-0.10,-0.08,
     &             -0.06,-.055,-.035,-0.02,-0.01, .005, 0.04, 0.06,
     &              .085, 0.12, 0.15, 0.18, 0.21, .235, 0.26, 0.29,
     &              0.31, 0.32, 0.34, .365, 0.38, 0.40, .425, .445,
     &              0.47, 0.50, .525, 0.56, 0.58, .635, 0.70, 0.75,
     &              0.82, 0.90, 1.03, 1.24, 1.35, 1.45/
      end

