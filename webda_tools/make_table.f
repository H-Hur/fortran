c f95 -o make_table make_table.f
c gfortran -o make_table ~/work/fortran/subroutine/webda_tools/make_table.f
      character input*16, clname*20, prname*20, ubvpename*20
      character vripename*20, vriccdname*20, ubvccdname*20, viname*20
      character ra_c*20, dec_c*20,c(18)*1,output*16, compare*18
      integer master
 
      write(input,'(a16)') 'make_table.input'
c      print *,'input file name=?'
c      read(*,'(a16)') input

      call read_input(input,clname,prname,ra_c,dec_c,ubvpename,
     &ubvccdname,vripename,vriccdname,viname,master)

      call read_proper(prname)
      call read_ubv_peo(ubvpename)
      call read_ubv_ccd(ubvccdname)
      call read_vri_peo(vripename)
      call read_vri_ccd(vriccdname)
      call read_vi_ccd(viname)
      call make_id_ref(master)


      write(compare,'(a18)') '                  '
      read(input,'(18a1)') (c(i),i=1,16)
      do i=1,16
       if(c(i).eq.'.') then
        goto 100
       endif
      enddo
100   write(c(i+1),'(a1)') 'c'
      write(c(i+2),'(a1)') 'o'
      write(c(i+3),'(a1)') 'm'
      write(c(i+4),'(a1)') 'p'
      write(c(i+5),'(a1)') 'a'
      write(c(i+6),'(a1)') 'r'
      write(c(i+7),'(a1)') 'e'
      write(compare,'(18a1)') (c(j),j=1,i+7)

      write(output,'(a16)') '                '
      read(input,'(16a1)') (c(i),i=1,16)
      do i=1,16
       if(c(i).eq.'.') then
       goto 200
       endif
      enddo
200   write(c(i+1),'(a1)') 't'
      write(c(i+2),'(a1)') 'a'
      write(c(i+3),'(a1)') 'b'
      write(c(i+4),'(a1)') 'l'
      write(c(i+5),'(a1)') 'e'
      write(output,'(18a1)') (c(j),j=1,i+5)

      call make_compare(master)

      stop
      end

      subroutine make_compare(master)
      parameter (ns=10000,nr=30)
      common/vri_ccd/no_vri_ccd(ns),no_vri_ccd_ref(ns),v_vri_ccd(ns),
     & vr_ccd(ns),vi_ccd(ns),nvri_ccd
      character vripename*20,H*1,input*20
      common/vri_pe/no_vri_pe(ns),no_vri_pe_ref(ns),v_vri_pe(ns),
     & vr_pe(ns),vi_pe(ns),nvri_pe
      common/ubv_ccd/no_ubv_ccd(ns),no_ubv_ccd_ref(ns),v_ccd(ns),
     & bv_ccd(ns),ub_ccd(ns),nubv_ccd
      common/ubv_pe/no_ubv_pe(ns),no_ubv_pe_ref(ns),v_pe(ns),
     & bv_pe(ns),ub_pe(ns),nubv_pe
      common/vi_ccd/no_vi_ccd(ns),no_vi_ccd_ref(ns),v_vi_ccd(ns),
     &vi_vi_ccd(ns),nvi_ccd
      common/proper/no_pr(ns),no_pr_ref(ns),prob(ns),npro
      common/id_ref/n_ref(nr),nref
      character output*16
      real vm(ns),bvm(ns),ubm(ns),vrm(ns),vim(ns)
      integer nc1,nc2,nom(ns),nm,master
      integer nvri_pe_used(ns),nvri_ccd_used(ns),nvi_used(ns)
      open(22,file='draw_compare.input')
      write(22,'(t1,a60)') '# input data files to draw compares.       
     &                                                     '
      write(22,'(t1,a60)') '# input file      vvi vi vbv bv ub ub2 maste
     &r  del_y                                              '
      do j=1,nref
       if(n_ref(j).lt.10) then
        write(output,'(a3,i1,a4)') '000',n_ref(j),'.dat'
       elseif(n_ref(j).lt.100) then
        write(output,'(a2,i2,a4)') '00',n_ref(j),'.dat'
       elseif(n_ref(j).lt.1000) then
        write(output,'(a1,i3,a4)') '0',n_ref(j),'.dat'
       else
        write(output,'(i4,a4)') n_ref(j),'.dat'
       endif
       if(n_ref(j).ne.master) then
        write(22,'(t1,a16,t20,6I3)') output,0,0,0,0,0,0
        open(21,file=output)
       write(21,'(t1,a80)') '#   NO     V     (B-V)   (U-B)   (V-R)   (V
     &-I)                                                          '
       else
        write(22,'(t1,a16,t20,6I3,t40,a6,t47,f5.2)') output,0,0,0,0,0,0,
     &  'master',0.65

        open(21,file=output)
       write(21,'(t1,a80)') '#   NO     V     (B-V)   (U-B)   (V-R)   (V
     &-I)                                                          '
       endif
1      format(t1,i6,5f8.3,5i3)
       do i=1,nubv_pe
        if(no_ubv_pe_ref(i).eq.n_ref(j)) then
         call find_same(n_ref(j),no_ubv_pe(i),0,no_ubv_pe(i),
     &   no_ubv_pe(i),no_ubv_pe(i),nt,nj,vi)
         if(nt.eq.1) then
          nvri_pe_used(nj)=1
         elseif(nt.eq.2) then
          nvri_ccd_used(nj)=1
         elseif(nt.eq.3) then
          nvi_used(nj)=1
         endif
         write(21,1)no_ubv_pe(i),v_pe(i),bv_pe(i),ub_pe(i),99.999,vi
c         print *,n_ref(j),(no_ubv_pe(i)),v_pe(i),bv_pe(i),ub_pe(i),
c     &   nt,nj,vi
        endif
       enddo
       do i=1,nubv_ccd
        if(no_ubv_ccd_ref(i).eq.n_ref(j)) then
         call find_same(n_ref(j),0,no_ubv_ccd(i),no_ubv_ccd(i),
     &   no_ubv_ccd(i),no_ubv_ccd(i),nt,nj,vi)
         if(nt.eq.1) then
          nvri_pe_used(nj)=1
         elseif(nt.eq.2) then
          nvri_ccd_used(nj)=1
         elseif(nt.eq.3) then
          nvi_used(nj)=1
         endif
         write(21,1)no_ubv_ccd(i),v_ccd(i),bv_ccd(i),ub_ccd(i),
     &   99.999,vi
c         print *,n_ref(j),(no_ubv_ccd(i)),v_ccd(i),bv_ccd(i),ub_ccd(i),
c     &   nt,nj,vi

         endif
       enddo

       do i=1,nvri_pe
        if(no_vri_pe_ref(i).eq.n_ref(j).and.nvri_pe_used(i).ne.1) then
         write(21,1)no_vri_pe(i),v_vri_pe(i),99.999,99.999,vr_pe(i),
     &   vi_pe(i)
        endif
       enddo
       do i=1,nvri_ccd
        if(no_vri_ccd_ref(i).eq.n_ref(j).and.nvri_ccd_used(i).ne.1) then
         write(21,1)no_vri_ccd(i),v_vri_ccd(i),99.999,99.999,
     &   vr_ccd(i),vi_ccd(i)
        endif
       enddo
       do i=1,nvi_ccd
        if(no_vi_ccd_ref(i).eq.n_ref(j).and.nvi_used(i).ne.1) then
         write(21,1)no_vi_ccd(i),v_vi_ccd(i),99.999,99.999,
     &   99.999,vi_vi_ccd(i)
        endif
       enddo

       close(21)
      enddo
      close(22)

      return
      end

      subroutine find_same(nref,n_ubv_pe,n_ubv_ccd,n_vri_pe,
     &n_vri_ccd,n_vi,m,nj,vi)
      parameter (ns=10000)      
      common/vri_ccd/no_vri_ccd(ns),no_vri_ccd_ref(ns),v_vri_ccd(ns),
     & vr_ccd(ns),vi_ccd(ns),nvri_ccd
      character vripename*20,H*1,input*20
      common/vri_pe/no_vri_pe(ns),no_vri_pe_ref(ns),v_vri_pe(ns),
     & vr_pe(ns),vi_pe(ns),nvri_pe
      common/ubv_ccd/no_ubv_ccd(ns),no_ubv_ccd_ref(ns),v_ccd(ns),
     & bv_ccd(ns),ub_ccd(ns),nubv_ccd
      common/ubv_pe/no_ubv_pe(ns),no_ubv_pe_ref(ns),v_pe(ns),
     & bv_pe(ns),ub_pe(ns),nubv_pe
      common/vi_ccd/no_vi_ccd(ns),no_vi_ccd_ref(ns),v_vi_ccd(ns),
     &vi_vi_ccd(ns),nvi_ccd
      integer m,nref,nj
      real vi
      m=0
      vi=99.999
      if(nubv_pe.eq.0) then
       goto 31
      endif
      do i=1,nubv_pe
       if(nref.eq.no_ubv_pe_ref(i)) then
        if(n_vri_pe.eq.0) then
         goto 30
        endif
        do j=1,nvri_pe
         if(n_ubv_pe.eq.no_vri_pe(j).and.abs(v_pe(i)-v_vri_pe(j))
     &   .lt.0.0001.and.nref.eq.no_vri_pe_ref(j)) then
          m=1
          nj=j
          vi=vi_pe(j)
          goto 100
         endif
        enddo
       endif
30    l=l
      enddo
31    if(nubv_ccd.eq.0) then
       goto 61
      endif
      do i=1,nubv_ccd
       if(nref.eq.no_ubv_ccd_ref(i)) then
        if(nvri_ccd.eq.0) then
         goto 50
        endif
        do j=1,nvri_ccd
         if(n_ubv_ccd.eq.no_vri_ccd(j).and.
     &   abs(v_ccd(i)-v_vri_ccd(j)).lt.0.0001
     &   .and.nref.eq.no_vri_ccd_ref(j)) then
          m=2
          nj=j
          vi=vi_ccd(j)
          goto 100
         endif
        enddo
50      if(nvi_ccd.eq.0) then
         goto 60
        endif
        do j=1,nvi_ccd
         if(n_ubv_ccd.eq.no_vi_ccd(j).and.
     &   abs(v_ccd(i)-v_vi_ccd(j)).lt.0.0001
     &   .and.nref.eq.no_vi_ccd_ref(j)) then
          m=3
          nj=j
          vi=vi_vi_ccd(j)
          goto 100
         endif
        enddo
       endif
60    l=l
      enddo

61    vi=99.999
100   return
      end

      subroutine make_id_ref(master)
      parameter (ns=10000,nr=30)
      common/vri_ccd/no_vri_ccd(ns),no_vri_ccd_ref(ns),v_vri_ccd(ns),
     & vr_ccd(ns),vi_ccd(ns),nvri_ccd
      character vripename*20,H*1,input*20
      common/vri_pe/no_vri_pe(ns),no_vri_pe_ref(ns),v_vri_pe(ns),
     & vr_pe(ns),vi_pe(ns),nvri_pe
      common/ubv_ccd/no_ubv_ccd(ns),no_ubv_ccd_ref(ns),v_ccd(ns),
     & bv_ccd(ns),ub_ccd(ns),nubv_ccd
      common/ubv_pe/no_ubv_pe(ns),no_ubv_pe_ref(ns),v_pe(ns),
     & bv_pe(ns),ub_pe(ns),nubv_pe
      common/vi_ccd/no_vi_ccd(ns),no_vi_ccd_ref(ns),v_vi_ccd(ns),
     &vi_vi_ccd(ns),nvi_ccd
      common/proper/no_pr(ns),no_pr_ref(ns),prob(ns),npro
      common/id_ref/n_ref(nr),nref
      integer nsame,master
      nref=1
      n_ref(1)=master
      do i=1,nubv_pe
       nsame=0
       do j=1,nref
        if(n_ref(j).eq.no_ubv_pe_ref(i)) then
         nsame=nsame+1
        endif
       enddo
       if(nsame.eq.0) then
        nref=nref+1
        n_ref(nref)=no_ubv_pe_ref(i)
       endif
      enddo
      do i=1,nubv_ccd
       nsame=0
       do j=1,nref
        if(n_ref(j).eq.no_ubv_ccd_ref(i)) then
         nsame=nsame+1
        endif
       enddo
       if(nsame.eq.0) then
        nref=nref+1
        n_ref(nref)=no_ubv_ccd_ref(i)
       endif
      enddo
      do i=1,nvri_pe
       nsame=0
       do j=1,nref
        if(n_ref(j).eq.no_vri_pe_ref(i)) then
         nsame=nsame+1
        endif
       enddo
       if(nsame.eq.0) then
        nref=nref+1
        n_ref(nref)=no_vri_pe_ref(i)
       endif
      enddo
      do i=1,nvri_ccd
       nsame=0
       do j=1,nref
        if(n_ref(j).eq.no_vri_ccd_ref(i)) then
         nsame=nsame+1
        endif
       enddo
       if(nsame.eq.0) then
        nref=nref+1
        n_ref(nref)=no_vri_ccd_ref(i)
       endif
      enddo
      do i=1,nvi_ccd
       nsame=0
       do j=1,nref
        if(n_ref(j).eq.no_vi_ccd_ref(i)) then
         nsame=nsame+1
        endif
       enddo
       if(nsame.eq.0) then
        nref=nref+1
        n_ref(nref)=no_vi_ccd_ref(i)
       endif
      enddo
      return
      end  

      subroutine read_vri_ccd(vriccdname)
      parameter (ns=10000)
      character vriccdname*20,H*1,input*20
      common/vri_ccd/no_vri_ccd(ns),no_vri_ccd_ref(ns),v_vri_ccd(ns),
     & vr_ccd(ns),vi_ccd(ns),nvri_ccd
      real ri
      real rno_n,ref_n,v_n,bv_n,ub_n
      call c_name(vriccdname,input)
      if(input.eq.'none') then
       nvri_ccd=0
       goto 900
      endif
      open(21,file=input)
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.eq.'-') then
        goto 100
       endif
      enddo
99    print *,'Reading error the ',vriccdname
100   do i=1,ns
101    read(21,*,end=199) no_vri_ccd(i),no_vri_ccd_ref(i),v_vri_ccd(i),
     & vr_ccd(i),ri
       backspace(21)
       read(21,*) rno_n,ref_n
       if(rno_n.eq.no_vri_ccd(i).and.ref_n.eq.no_vri_ccd_ref(i)) then
        vi_ccd(i)=vr_ccd(i)+ri
        if(rno_n.eq.int(ri))  then
         vi_ccd(i)=99.999
        endif
       else
        backspace(21)
        vi_ccd(i)=99.999     
       endif
      enddo
199   nvri_ccd=i-1
900   close(21)
      return
      end

      subroutine read_vri_peo(vripename)
      parameter (ns=10000)
      character vripename*20,H*1,input*20
      common/vri_pe/no_vri_pe(ns),no_vri_pe_ref(ns),v_vri_pe(ns),
     & vr_pe(ns),vi_pe(ns),nvri_pe
      real ri
      real rno_n,ref_n,v_n,bv_n,ub_n
      call c_name(vripename,input)
      if(input.eq.'none') then
       nvri_pe=0
       goto 900
      endif
      open(21,file=input)
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.eq.'-') then
        goto 100
       endif
      enddo
99    print *,'Reading error the ',vripename
100   do i=1,ns
101    read(21,*,end=199) no_vri_pe(i),no_vri_pe_ref(i),v_vri_pe(i),
     & vr_pe(i),ri
       vi_pe(i)=vr_pe(i)+ri
       backspace(21)
       read(21,*) rno_n,ref_n
       if(rno_n.eq.no_vri_pe(i).and.ref_n.eq.no_vri_pe_ref(i)) then
        vi_pe(i)=vr_pe(i)+ri
        if(rno_n.eq.int(ri))  then
         vi_pe(i)=99.999 
        endif
       else
        backspace(21)
        vi_pe(i)=99.999
       endif
      enddo
199   nvri_pe=i-1
900   close(21)
      return
      end

      subroutine read_ubv_ccd(ubvccdname)
      parameter (ns=10000)
      character ubvccdname*20,H*1,input*20
      common/ubv_ccd/no_ubv_ccd(ns),no_ubv_ccd_ref(ns),v_ccd(ns),
     & bv_ccd(ns),ub_ccd(ns),nubv_ccd
      real rno_n,ref_n,v_n,bv_n,ub_n
      call c_name(ubvccdname,input)
      if(input.eq.'none') then
       nubv_ccd=0
       goto 900
      endif
      open(21,file=input)
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.eq.'-') then
        goto 100
       endif
      enddo
99    print *,'Reading error the ',ubvccdname
100   do i=1,ns
101    read(21,*,end=199) no_ubv_ccd(i),no_ubv_ccd_ref(i),v_ccd(i),
     & bv_ccd(i),ub_ccd(i)
       backspace(21)
       read(21,*) rno_n,ref_n
       if(rno_n.eq.no_ubv_ccd(i).and.ref_n.eq.no_ubv_ccd_ref(i)) then
        if(rno_n.eq.int(ub_ccd(i)))
     &  then
         ub_ccd(i)=99.999
        endif
       else
        backspace(21)
        ub_ccd(i)=99.999
       endif
      enddo
199   nubv_ccd=i-1
900   close(21)
      return
      end

      subroutine read_ubv_peo(ubvpename)
      parameter (ns=10000)
      character ubvpename*20,H*1,input*20
      common/ubv_pe/no_ubv_pe(ns),no_ubv_pe_ref(ns),v_pe(ns),
     & bv_pe(ns),ub_pe(ns),nubv_pe
      real rno_n,ref_n,v_n,bv_n,ub_n
      call c_name(ubvpename,input)
      if(input.eq.'none') then
       nubv_pe=0
       goto 900
      endif
      open(21,file=input)
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.eq.'-') then
        goto 100
       endif
      enddo
99    print *,'Reading error the ',ubvpename
100   do i=1,ns
101    read(21,*,end=199) no_ubv_pe(i),no_ubv_pe_ref(i),v_pe(i),
     & bv_pe(i),ub_pe(i)
c       print *,no_ubv_pe(i),no_ubv_pe_ref(i),v_pe(i),
c     & bv_pe(i),ub_pe(i)

       backspace(21)
       read(21,*) rno_n,ref_n
       if(rno_n.eq.no_ubv_pe(i).and.ref_n.eq.no_ubv_pe_ref(i)) then
        if(rno_n.eq.int(ub_pe(i))) then
         ub_pe(i)=99.999
        endif
       else
        backspace(21)
        ub_pe(i)=99.999
       endif
      enddo
199   nubv_pe=i-1
900   close(21)
      return
      end

      subroutine read_vi_ccd(viname)
      parameter (ns=10000)
      character viname*20,H*1,input*20
      common/vi_ccd/no_vi_ccd(ns),no_vi_ccd_ref(ns),v_vi_ccd(ns),
     &vi_vi_ccd(ns),nvi_ccd
      real rno_n,ref_n,v_n,bv_n,ub_n
      call c_name(viname,input)
      if(input.eq.'none') then
       nvi_ccd=0
       goto 900
      endif
      open(21,file=input)
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.eq.'-') then
        goto 100
       endif
      enddo
99    print *,'Reading error the ',viname
100   do i=1,ns
101    read(21,*,end=199) no_vi_ccd(i),no_vi_ccd_ref(i),v_vi_ccd(i),
     & vi_vi_ccd(i)
       if(no_vi_ccd(i).eq.no_vi_ccd(i-1).and.no_vi_ccd_ref(i).eq.
     & no_vi_ccd_ref(i-1)) then
        goto 101
       endif
      enddo
199   nvi_ccd=i-1
900   close(21)
      return
      end

      subroutine read_proper(prname)
      parameter (ns=10000)
      character prname*20,H*1,input*20
      common/proper/no_pr(ns),no_pr_ref(ns),prob(ns),npro
      call c_name(prname,input)
      if(input.eq.'none') then
       npro=0
       goto 900
      endif
      open(21,file=input)
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.eq.'-') then
        goto 100
       endif
      enddo
99    print *,'Reading error the ',prname
100   do i=1,ns
       read(21,*,end=199) no_pr(i),no_pr_ref(i),prob(i)
      enddo
199   npro=i-1
900   close(21)
      return
      end                

      subroutine read_input(input,clname,prname,ra_c,dec_c,ubvpename,
     &ubvccdname,vripename,vriccdname,viname,master)
      character input*16, clname*20, prname*20, ubvpename*20,viname*20
      character vripename*20, vriccdname*20, ubvccdname*20,HH*20
      character line(41)*1, H*1, sub*20, tar*20, ra_c*20, dec_c*20
      integer master
      open(21,file=input)
      do i=1,100
       read(21,'(t1,a1)',end=99) H
       if(H.eq.'#') then
        backspace(21)
        read(21,'(t3,41a1)') (line(j),j=1,41) 
        do j=1,41
         if(line(j).eq.':') then
          goto 10
         endif
        enddo
10      write(sub,'(a20)') '                    '  
        write(sub,'(20a1)') (line(k),k=1,j-1)
        if(sub.eq.'name') then 
         write(clname,'(a20)') '                    '
         write(clname,'(20a1)') (line(k),k=j+1,j+20)
        elseif(sub.eq.'proper motion data') then
         write(prname,'(a20)') '                    '
         write(prname,'(20a1)') (line(k),k=j+1,j+20)
        elseif(sub.eq.'UBV PE data') then
         write(ubvpename,'(a20)') '                    '
         write(ubvpename,'(20a1)') (line(k),k=j+1,j+20)
        elseif(sub.eq.'VRI PE data') then
         write(vripename,'(a20)') '                    '
         write(vripename,'(20a1)') (line(k),k=j+1,j+20)
        elseif(sub.eq.'UBV CCD data') then
         write(ubvccdname,'(a20)') '                    '
         write(ubvccdname,'(20a1)') (line(k),k=j+1,j+20)
        elseif(sub.eq.'VRI CCD data') then
         write(vriccdname,'(a20)') '                    '
         write(vriccdname,'(20a1)') (line(k),k=j+1,j+20)
        elseif(sub.eq.'VI CCD data') then
         write(viname,'(a20)') '                    '
         write(viname,'(20a1)') (line(k),k=j+1,j+20)
        elseif(sub.eq.'RA') then
         write(ra_c,'(a20)') '                    '
         write(ra_c,'(20a1)') (line(k),k=j+1,j+20)
        elseif(sub.eq.'DEC') then
         write(dec_c,'(a20)') '                    '
         write(dec_c,'(20a1)') (line(k),k=j+1,j+20)
        elseif(sub.eq.'master') then
         write(HH,'(a20)') '                    '
         write(HH,'(20a1)') (line(k),k=j+1,j+20)
         read(HH,*) master
        endif
       endif
      enddo
99    close(21)
      return
      end

      subroutine c_name(input,output)
      character input*20,c(20)*1,output*20
      write(output,'(a20)') '                    '
      read(input,'(20a1)') (c(i),i=1,20)
      do i=1,19
       if(c(i).eq.' '.and.c(i+1).ne.' ') then
        goto 100
       endif
      enddo
100   write(output,'(20a1)') (c(j),j=i+1,20)
      return
      end
