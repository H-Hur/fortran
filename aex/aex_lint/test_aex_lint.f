c     f77 -o test_aex_lint test_aex_lint.f -L/usr/local/lib -L/usr/X11R6/lib -lplotsub -ldevices -lutils -lX11  
c     21,22,23,24,31
      character input*16,output*32

c      write(input,'(a16)') 'coeff_c1.input                            ' 
c      call read_coeff(input)
c      write(input,'(a16)') 'V_1_new.prn                                '
c      call read_stsv(input)
c      write(input,'(a16)') 'I_1_new.prn                                '
c      call read_stsi(input)
c      write(input,'(a16)') 'B_1_new.prn                                '
c      call read_stsb(input)
c      write(output,'(a32)') 'postportfile chip1.ps                     '
c      call sts_standard(output)

c      write(input,'(a16)') 'coeff_c2.input                            ' 
c      call read_coeff(input)
c      write(input,'(a16)') 'V_2_new.prn                                '
c      call read_stsv(input)
c      write(input,'(a16)') 'I_2_new.prn                                '
c      call read_stsi(input)
c      write(input,'(a16)') 'B_2_new.prn                                '
c      call read_stsb(input)
c      write(output,'(a32)') 'postportfile chip2.ps                     '
c      call sts_standard(output)

c      write(input,'(a16)') 'coeff_c3.input                            ' 
c      call read_coeff(input)
c      write(input,'(a16)') 'V_3_new.prn                                '
c      call read_stsv(input)
c      write(input,'(a16)') 'I_3_new.prn                                '
c      call read_stsi(input)
c      write(input,'(a16)') 'B_3_new.prn                                '
c      call read_stsb(input)
c      write(output,'(a32)') 'postportfile chip3.ps                     '
c      call sts_standard(output)

c      write(input,'(a16)') 'coeff_c4.input                            ' 
c      call read_coeff(input)
c      write(input,'(a16)') 'V_4_new.prn                                '
c      call read_stsv(input)
c      write(input,'(a16)') 'I_4_new.prn                                '
c      call read_stsi(input)
c      write(input,'(a16)') 'B_4_new.prn                                '
c      call read_stsb(input)
c      write(output,'(a32)') 'postportfile chip4.ps                     '
c      call sts_standard(output)

c      write(input,'(a16)') 'coeff_c5.input                            ' 
c      call read_coeff(input)
c      write(input,'(a16)') 'V_5_new.prn                                '
c      call read_stsv(input)
c      write(input,'(a16)') 'I_5_new.prn                                '
c      call read_stsi(input)
c      write(input,'(a16)') 'B_5_new.prn                                '
c      call read_stsb(input)
c      write(output,'(a32)') 'postportfile chip5.ps                     '
c      call sts_standard(output)

      write(input,'(a16)') 'coeff_c6.input                            ' 
      call read_coeff(input)
      write(input,'(a16)') 'V_6_new.prn                                '
      call read_stsv(input)
      write(input,'(a16)') 'I_6_new.prn                                '
      call read_stsi(input)
      write(input,'(a16)') 'B_6_new.prn                                '
      call read_stsb(input)
      write(output,'(a32)') 'postportfile chip6.ps                     '
      call sts_standard(output)

c      write(input,'(a16)') 'coeff_c7.input                            ' 
c      call read_coeff(input)
c      write(input,'(a16)') 'V_7_new.prn                                '
c      call read_stsv(input)
c      write(input,'(a16)') 'I_7_new.prn                                '
c      call read_stsi(input)
c      write(input,'(a16)') 'B_7_new.prn                                '
c      call read_stsb(input)
c      write(output,'(a32)') 'postportfile chip7.ps                     '
c      call sts_standard(output)

c      write(input,'(a16)') 'coeff_c8.input                            ' 
c      call read_coeff(input)
c      write(input,'(a16)') 'V_8_new.prn                                '
c      call read_stsv(input)
c      write(input,'(a16)') 'I_8_new.prn                                '
c      call read_stsi(input)
c      write(input,'(a16)') 'B_8_new.prn                                '
c      call read_stsb(input)
c      write(output,'(a32)') 'postportfile chip8.ps                     '
c      call sts_standard(output)



      stop
      end
c     Filter 1 = U, 2 = B, 3 = V, 4 = R, 5 = I,6 = Ha
c     Color number  1 = U-B, 2 = B-V, 3 = V-R, 4 = V-I, 5 = R-I  
      subroutine sts_standard(output)
      parameter(ns=100000)
      common/sts_i/ci,star_i,xc_i(ns),yc_i(ns),ci_ins(ns),cierr(ns),
     &X_i(ns),t_i(ns),ni
      common/sts_v/cv,star_v,xc_v(ns),yc_v(ns),v_ins(ns),verr(ns),
     &X_v(ns),t_v(ns),sv(ns),sbv(ns),svi(ns),sri(ns),nsv(ns),nsbv(ns),
     &nsvi(ns),nsri(ns),nv
      common/sts_b/cb,star_b,xc_b(ns),yc_b(ns),b_ins(ns),berr(ns),
     &X_b(ns),t_b(ns),nb
      character star_b(ns)*12,cb(ns)*1,output*32,star_vi(ns)*2
      character star_v(ns)*12,cv(ns)*1,star_i(ns)*12,ci(ns)*1
      character ocv(ns)*1,oci(ns)*1,ocb(ns)*1,ocvi(ns)*1,ocbv(ns)*1
      character star_bv(ns)*12
      integer nvi,nb
      real old_V,old_I,old_B,Ovi(ns),Ovier(ns),Obv(ns),Obver(ns)
      real Ov(ns),Over(ns),Oi(ns),Oier(ns),Ob(ns),Ober(ns)
      real osb(ns),osbv(ns),osv(ns),osvi(ns),osi(ns),yyy
      real ovx(ns),ovy(ns),oix(ns),oiy(ns),obx(ns),oby(ns)
      real wti,wtv,wtb,wtvi,wtbv,pv,pi,pvi,pb,pbv,siz
      call sm_device (output)
      call sm_graphics
      call sm_erase()
      call sm_expand(0.7)
      yyy=0.11 
      siz=0.9
      call sm_location(3000,15000,3000,7000)
      call sm_limits(0.,2048.,-yyy,yyy)
      call sm_box(1,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_xlabel('Xccd')
      call sm_ylabel('\\gDB\\\\dCTIO-S\\\\u')
      call sm_location(3000,15000,7000,11000)
      call sm_box(0,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_ylabel('\\gDV\\\\dCTIO-S\\\\u')
      call sm_location(3000,15000,11000,15000)
      call sm_box(0,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_ylabel('\\gDI\\\\dCTIO-S\\\\u')

      call sm_location(3000,15000,18000,22000)
      call sm_limits(0.,4096.,-yyy,yyy)
      call sm_box(1,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_xlabel('Yccd')
      call sm_ylabel('\\gDB\\\\dCTIO-S\\\\u')
      call sm_location(3000,15000,22000,26000)
      call sm_box(0,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_ylabel('\\gDV\\\\dCTIO-S\\\\u')
      call sm_location(3000,15000,26000,30000)
      call sm_box(0,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_ylabel('\\gDI\\\\dCTIO-S\\\\u')


      call sm_location(19000,31000,3000,7000)
      call sm_limits(-0.4,3.,-yyy,yyy)
      call sm_box(1,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_xlabel('(B-V)\\dS')
      call sm_ylabel('\\gDB\\\\dCTIO-S\\\\u')
      call sm_location(19000,31000,7000,11000)
      call sm_box(0,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_ylabel('\\gDV\\\\dCTIO-S\\\\u')
      call sm_location(19000,31000,11000,15000)
      call sm_box(0,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_ylabel('\\gD(B-V)\\\\dCTIO-S\\\\u')

      call sm_location(19000,31000,18000,22000)
      call sm_limits(-0.4,4.,-yyy,yyy)
      call sm_box(1,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_xlabel('(V-I)\\dS')
      call sm_ylabel('\\gDV\\\\dCTIO-S\\\\u')
      call sm_location(19000,31000,22000,26000)
      call sm_box(0,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_ylabel('\\gDI\\\\dCTIO-S\\\\u')
      call sm_location(19000,31000,26000,30000)
      call sm_box(0,2,0,0)
      call sm_relocate(-100.,0.)
      call sm_draw(10000.,0.)
      call sm_ylabel('\\gD(V-I)\\\\dCTIO-S\\\\u')


      call sm_expand(1.5)
      nvi=0
      do i=1,ni
       do j=1,nv
        if(star_i(i).eq.star_v(j)) then
         if(svi(j).eq.0) goto 100
         if(v_ins(j).le.11..or.ci_ins(i).le.11.) goto 100
         nvi=nvi+1
         write(star_vi(nvi),'(a12)') star_i(i)
         write(ocv(nvi),'(a1)') cv(j)
         write(oci(nvi),'(a1)') ci(i)
         write(ocvi(nvi),'(a1)') ' '
         if(ci(i).eq.'c'.or.cv(j).eq.'c') write(ocvi(nvi),'(a1)') 'c'
         if(ci(i).eq.'b'.or.cv(j).eq.'b') write(ocvi(nvi),'(a1)') 'b'
         if(ci(i).eq.'s'.or.cv(j).eq.'s') write(ocvi(nvi),'(a1)') 's'
         if(ci(i).eq.'r'.or.cv(j).eq.'r') write(ocvi(nvi),'(a1)') 'r'
         if(ci(i).eq.'e'.or.cv(j).eq.'e') write(ocvi(nvi),'(a1)') 'e'
         OVer(nvi)=verr(j)
         Oier(nvi)=cierr(i)
         Ovier(nvi)=sqrt(cierr(i)*cierr(i)+verr(j)*verr(j))
         Ovi(nvi)=v_ins(j)-ci_ins(i)
         Ov(nvi)=v_ins(j)
         Oi(nvi)=ci_ins(i)
         osv(nvi)=sv(j)
         osvi(nvi)=svi(j)
         osi(nvi)=sv(j)-svi(j)
         ovx(nvi)=xc_v(j)
         ovy(nvi)=yc_v(j)
         oix(nvi)=xc_i(i)
         oiy(nvi)=yc_i(i)
30       old_V=Ov(nvi)
         old_I=Oi(nvi)
         call standard(3,4,v_ins(j),Ovi(nvi),0.,0.,0.,Ovi(nvi),0.,X_v(j)
     &   ,t_v(j),xc_v(j),yc_v(j),Ov(nvi))
         call standard(5,4,ci_ins(i),Ovi(nvi),0.,0.,0.,Ovi(nvi),0.,
     &   X_i(i), t_i(i),xc_i(i),yc_i(i),Oi(nvi))
         Ovi(nvi)=Ov(nvi)-Oi(nvi)
         if(abs(Ov(nvi)-old_V).gt.0.0001) goto 30
         if(abs(Oi(nvi)-old_I).gt.0.0001) goto 30
         wtv=0.015/OVer(nvi)/OVer(nvi)
         wti=0.015/Oier(nvi)/Oier(nvi)
         wtvi=0.015/Ovier(nvi)/Ovier(nvi)
         wtvi=min(wtvi,1./0.015)*0.015
         wtvi=max(wtvi,0.15)
         wtv=min(wtv,1./0.015)*0.015
         wtv=max(wtv,0.15)
         wti=min(wti,1./0.015)*0.015
         wti=max(wti,0.15)
         pv=403.+wtv*siz
         pi=403.+wti*siz
         pvi=403.+wtvi*siz
         write(ocvi(nvi),'(a1)') ' '
         if(cv(j).eq.' '.and.ci(i).ne.' ') write(ocvi(nvi),'(a1)') ci(i)
         if(cv(j).ne.' '.and.ci(i).eq.' ') write(ocvi(nvi),'(a1)') cv(j)
         if(cv(j).ne.' '.and.ci(i).ne.' '.and.cv(j).eq.ci(i)) 
     &   write(ocvi(nvi),'(a1)') ci(i)
         if(cv(j).ne.' '.and.ci(i).ne.' '.and.cv(j).ne.ci(i)) 
     &   write(ocvi(nvi),'(a1)') 'F'
         if(cv(j).ne.' ') pv=403.+wtv*siz*0.5
         if(ocvi(nvi).ne.' ') pvi=403.+wtvi*siz*0.5
         if(ci(i).ne.' ') pi=403.+wti*siz*0.5

         call sm_limits(0.,2048.,-yyy,yyy)     !X : Xc
         call sm_location(3000,15000,7000,11000) ! V
         call sm_ctype('yellow')
         if(Over(nvi).le.0.015) call sm_ctype('black')
         if(Over(nvi).gt.0.015) pv=403.+wtv*siz*0.5
         if(cv(j).eq.'c') call sm_ctype('blue')
         if(cv(j).eq.'b') call sm_ctype('green')
         if(cv(j).eq.'s') call sm_ctype('magenta')
         if(cv(j).eq.'r') call sm_ctype('cyan')
         if(cv(j).eq.'e') call sm_ctype('red')
         call sm_ptype(pv,1)
c         call sm_points(ovx(nvi),-osv(nvi)+Ov(nvi),1)
         call sm_location(3000,15000,11000,15000)! I 
         call sm_ctype('yellow')
         if(Oier(nvi).le.0.015) call sm_ctype('black')
         if(Oier(nvi).gt.0.015) pi=403.+wti*siz*0.5
         if(ci(i).eq.'c') call sm_ctype('blue')
         if(ci(i).eq.'b') call sm_ctype('green')
         if(ci(i).eq.'s') call sm_ctype('magenta')
         if(ci(i).eq.'r') call sm_ctype('cyan')
         if(ci(i).eq.'e') call sm_ctype('red')
         call sm_ptype(pi,1)
         call sm_points(oix(nvi),-osi(nvi)+Oi(nvi),1)

         call sm_limits(0.,4096.,-yyy,yyy)    ! X :Yc
         call sm_location(3000,15000,22000,26000) !V
         if(Over(nvi).le.0.015) call sm_ctype('black')
         if(Over(nvi).gt.0.015) pi=403.+wtvi*siz*0.5
         if(cv(j).eq.'c') call sm_ctype('blue')
         if(cv(j).eq.'b') call sm_ctype('green')
         if(cv(j).eq.'s') call sm_ctype('magenta')
         if(cv(j).eq.'r') call sm_ctype('cyan')
         if(cv(j).eq.'e') call sm_ctype('red')
         call sm_ptype(pv,1)
c         call sm_points(ovy(nvi),-osv(nvi)+Ov(nvi),1)
         call sm_location(3000,15000,26000,30000) !I
         call sm_ctype('yellow')
         if(Oier(nvi).le.0.015) call sm_ctype('black')
         if(Oier(nvi).gt.0.015) pi=403.+wti*siz*0.5
         if(ci(i).eq.'c') call sm_ctype('blue')
         if(ci(i).eq.'b') call sm_ctype('green')
         if(ci(i).eq.'s') call sm_ctype('magenta')
         if(ci(i).eq.'r') call sm_ctype('cyan')
         if(ci(i).eq.'e') call sm_ctype('red')
         call sm_ptype(pi,1)
         call sm_points(oiy(nvi),-osi(nvi)+Oi(nvi),1)
 
         call sm_limits(-0.4,4.,-yyy,yyy)          !X : V-I
         call sm_location(19000,31000,18000,22000) !V
         call sm_ctype('yellow')
         if(Over(nvi).le.0.015) call sm_ctype('black')
         if(Over(nvi).gt.0.015) pv=403.+wtv*siz*0.5
         if(cv(j).eq.'c') call sm_ctype('blue')
         if(cv(j).eq.'b') call sm_ctype('green')
         if(cv(j).eq.'s') call sm_ctype('magenta')
         if(cv(j).eq.'r') call sm_ctype('cyan')
         if(cv(j).eq.'e') call sm_ctype('red')
         call sm_ptype(pv,1)
         call sm_points(osvi(nvi),-osv(nvi)+Ov(nvi),1)
         call sm_location(19000,31000,22000,26000) !I
         call sm_ctype('yellow')
         if(Oier(nvi).le.0.015) call sm_ctype('black')
         if(Oier(nvi).gt.0.015) pi=403.+wti*siz*0.5
         if(ci(i).eq.'c') call sm_ctype('blue')
         if(ci(i).eq.'b') call sm_ctype('green')
         if(ci(i).eq.'s') call sm_ctype('magenta')
         if(ci(i).eq.'r') call sm_ctype('cyan')
         if(ci(i).eq.'e') call sm_ctype('red')
         call sm_ptype(pi,1)
         call sm_points(osvi(nvi),-osi(nvi)+Oi(nvi),1)
         call sm_location(19000,31000,26000,30000) !V-I
         call sm_ctype('yellow')
         if(Ovier(nvi).le.0.015) call sm_ctype('black')
         if(Ovier(nvi).gt.0.015) pvi=403.+wtvi*siz*0.5
         if(ocvi(nvi).eq.'c') call sm_ctype('blue')
         if(ocvi(nvi).eq.'b') call sm_ctype('green')
         if(ocvi(nvi).eq.'s') call sm_ctype('magenta')
         if(ocvi(nvi).eq.'r') call sm_ctype('cyan')
         if(ocvi(nvi).eq.'e') call sm_ctype('red')
         call sm_ptype(pvi,1)
         if(ocvi(nvi).eq.'F') call sm_ptype(400.+wtvi*siz*0.5*0.7,1)
         call sm_points(osvi(nvi),-osvi(nvi)+Ovi(nvi),1)
         goto 100
        endif
       enddo
100   enddo

      nbv=0
      do i=1,nb
       do j=1,nv
        if(star_b(i).eq.star_v(j)) then
         if(sbv(j).eq.0) goto 200
         if(v_ins(j).le.11..or.b_ins(i).le.11.) goto 200
         nbv=nbv+1
         write(star_bv(nbv),'(a12)') star_b(i)
         write(ocv(nbv),'(a1)') cv(j)
         write(ocb(nbv),'(a1)') cb(i)
         write(ocbv(nbv),'(a1)') ' '
         OVer(nbv)=verr(j)
         Ober(nbv)=berr(i)
         Obver(nbv)=sqrt(berr(i)*berr(i)+verr(j)*verr(j))
         Obv(nbv)=b_ins(i)-v_ins(j)
         Ov(nbv)=v_ins(j)
         Ob(nbv)=b_ins(i)
         osv(nbv)=sv(j)
         osbv(nbv)=sbv(j)
         osb(nbv)=sbv(j)+sv(j)
         ovx(nbv)=xc_v(j)
         ovy(nbv)=yc_v(j)
         obx(nbv)=xc_b(i)
         oby(nbv)=yc_b(i)
130      old_V=Ov(nbv)
         old_B=Ob(nbv)
         call standard(3,2,v_ins(j),Obv(nbv),0.,Obv(nbv),0.,0.,0.,X_v(j)
     &   ,t_v(j),xc_v(j),yc_v(j),Ov(nbv))
         call standard(2,2,b_ins(i),Obv(nbv),0.,Obv(nbv),0.,0.,0.,
     &   X_b(i), t_b(i),xc_b(i),yc_b(i),Ob(nbv))
         Obv(nbv)=Ob(nbv)-Ov(nbv)
         if(abs(Ov(nbv)-old_V).gt.0.0001) goto 130
         if(abs(Ob(nbv)-old_B).gt.0.0001) goto 130
         wtv=0.015/OVer(bvi)/OVer(nbv)
         wtb=0.015/Ober(nbv)/Ober(nbv)
         wtbv=0.015/Obver(nbv)/Obver(nbv)
         wtbv=min(wtbv,1./0.015)*0.015
         wtbv=max(wtbv,0.15)
         wtv=min(wtv,1./0.015)*0.015
         wtv=max(wtv,0.15)
         wtb=min(wtb,1./0.015)*0.015
         wtb=max(wtb,0.15)
         pv=403.+wtv*siz
         pb=403.+wtb*siz
         pbv=403.+wtbv*siz
         write(ocbv(nbv),'(a1)') ' '
         if(cv(j).eq.' '.and.cb(i).ne.' ') write(ocbv(nbv),'(a1)') cb(i)
         if(cv(j).ne.' '.and.cb(i).eq.' ') write(ocbv(nbv),'(a1)') cv(j)
         if(cv(j).ne.' '.and.cb(i).ne.' '.and.cv(j).eq.cb(i))
     &   write(ocbv(nbv),'(a1)') cb(i)
         if(cv(j).ne.' '.and.cb(i).ne.' '.and.cv(j).ne.cb(i))
     &   write(ocbv(nbv),'(a1)') 'F'
         if(cv(j).ne.' ') pv=403.+wtv*siz*0.5
         if(ocbv(nbv).ne.' ') pbv=403.+wtbv*siz*0.5
         if(cb(i).ne.' ') pb=403.+wtb*siz*0.5
         call sm_ctype('yellow')
         if(Ober(nbv).le.0.015) call sm_ctype('black')
         if(Ober(nbv).gt.0.015) pb=403.+wtb*siz*0.5
         if(cb(i).eq.'c') call sm_ctype('blue')
         if(cb(i).eq.'b') call sm_ctype('green')
         if(cb(i).eq.'s') call sm_ctype('magenta')
         if(cb(i).eq.'r') call sm_ctype('cyan')
         if(cb(i).eq.'e') call sm_ctype('red')
         call sm_limits(0.,2048.,-yyy,yyy)     !X : Xc
         call sm_location(3000,15000,3000,7000) ! B
         call sm_ptype(pb,1)
         call sm_points(obx(nbv),-osb(nbv)+Ob(nbv),1)
         call sm_location(3000,15000,7000,11000)! V 
         call sm_ctype('yellow')
         if(Over(nbv).le.0.015) call sm_ctype('black')
         if(Over(nbv).gt.0.015) pv=403.+wtv*siz*0.5
         if(cv(j).eq.'c') call sm_ctype('blue')
         if(cv(j).eq.'b') call sm_ctype('green')
         if(cv(j).eq.'s') call sm_ctype('magenta')
         if(cv(j).eq.'r') call sm_ctype('cyan')
         if(cv(j).eq.'e') call sm_ctype('red')
         call sm_ptype(pv,1)
         call sm_points(ovx(nbv),-osv(nbv)+Ov(nbv),1)

         call sm_limits(0.,4096.,-yyy,yyy)    ! X :Yc
         call sm_location(3000,15000,18000,22000) !B
         call sm_ctype('yellow')
         if(Ober(nbv).le.0.015) call sm_ctype('black')
         if(Ober(nbv).gt.0.015) pb=403.+wtb*siz*0.5
         if(cb(i).eq.'c') call sm_ctype('blue')
         if(cb(i).eq.'b') call sm_ctype('green')
         if(cb(i).eq.'s') call sm_ctype('magenta')
         if(cb(i).eq.'r') call sm_ctype('cyan')
         if(cb(i).eq.'e') call sm_ctype('red')
         call sm_ptype(pb,1)
         call sm_points(oby(nbv),-osb(nbv)+Ob(nbv),1)
         call sm_location(3000,15000,22000,26000) !V
         call sm_ctype('yellow')
         if(Over(nbv).le.0.015) call sm_ctype('black')
         if(Over(nbv).gt.0.015) pv=403.+wtv*siz*0.5
         if(cv(j).eq.'c') call sm_ctype('blue')
         if(cv(j).eq.'b') call sm_ctype('green')
         if(cv(j).eq.'s') call sm_ctype('magenta')
         if(cv(j).eq.'r') call sm_ctype('cyan')
         if(cv(j).eq.'e') call sm_ctype('red')
         call sm_ptype(pv,1)
         call sm_points(ovy(nbv),-osv(nbv)+Ov(nbv),1)

         call sm_limits(-0.4,3.,-yyy,yyy)          !X : B-V
         call sm_location(19000,31000,3000,7000) !B
         call sm_ctype('yellow')
         if(Ober(nbv).le.0.015) call sm_ctype('black')
         if(Ober(nbv).gt.0.015) pb=403.+wtb*siz*0.5
         if(cb(i).eq.'c') call sm_ctype('blue')
         if(cb(i).eq.'b') call sm_ctype('green')
         if(cb(i).eq.'s') call sm_ctype('magenta')
         if(cb(i).eq.'r') call sm_ctype('cyan')
         if(cb(i).eq.'e') call sm_ctype('red')
         call sm_ptype(pb,1)
         call sm_points(osbv(nbv),-osb(nbv)+Ob(nbv),1)
         call sm_location(19000,31000,7000,11000) !V
         call sm_ctype('yellow')
         if(Over(nbv).le.0.015) call sm_ctype('black')
         if(Over(nbv).gt.0.015) pv=403.+wtv*siz*0.5
         if(cv(j).eq.'c') call sm_ctype('blue')
         if(cv(j).eq.'b') call sm_ctype('green')
         if(cv(j).eq.'s') call sm_ctype('magenta')
         if(cv(j).eq.'r') call sm_ctype('cyan')
         if(cv(j).eq.'e') call sm_ctype('red')
         call sm_ptype(pv,1)
         call sm_points(osbv(nbv),-osv(nbv)+Ov(nbv),1)
         call sm_location(19000,31000,11000,15000) !B-V
         call sm_ctype('yellow')
         if(Obver(nbv).le.0.015) call sm_ctype('black')
         if(Obver(nbv).gt.0.015) pbv=403.+wtbv*siz*0.5
         if(ocbv(nbv).eq.'c') call sm_ctype('blue')
         if(ocbv(nbv).eq.'b') call sm_ctype('green')
         if(ocbv(nbv).eq.'s') call sm_ctype('magenta')
         if(ocbv(nbv).eq.'r') call sm_ctype('cyan')
         if(ocbv(nbv).eq.'e') call sm_ctype('red')
         call sm_ptype(pbv,1)
         if(ocbv(nbv).eq.'F') call sm_ptype(400.+wtbv*siz*0.5*0.7,1)
         call sm_points(osbv(nbv),-osbv(nbv)+Obv(nbv),1)
         goto 200
        endif
       enddo
200   enddo
900   call sm_ctype('black')
      call sm_expand(1.)
      call sm_gflush()
      call sm_hardcopy
      call sm_alpha
      return
      end

c     nfil = filter number
c     nc = color number
c     Color number  1 = U-B, 2 = B-V, 3 = V-R, 4 = V-I, 5 = R-I  
c     Other number 100 = Ut   101=xc,      102=yc,    103=rad,   104=theta
      subroutine standard(nfil,nc,cmag1,col,ub,bv,vr,vi,ri,X,UT,xc,yc,
     &omag)
      parameter(nf=20,ns=100000)
      common/coeff_const/num_filter(nf),num_color(nf),aeK1(nf),aeK2(nf)
     &,zero(nf),xcenter(nf),ycenter(nf),color_dp,ut_dp,x_dp,y_dp,r_dp,
     &theta_dp,add1_dp,add2_dp,add3_dp,nconst
      character color_dp(nf)*16,ut_dp(nf)*16,x_dp(nf)*16,y_dp(nf)*16
     &,r_dp(nf)*16,theta_dp(nf)*16,add1_dp(nf)*16,add2_dp(nf)*16
     &,add3_dp(nf)*16
      integer nfil,nc,nt
      real cmag1,omag,r,th,ub,bv,vr,vi,ri,X,UT,xc,yc,pi,col
      do i=1,nf
       if(nfil.eq.num_filter(i).and.nc.eq.num_color(i)) then
        nt=i
        goto 100
       endif
100   enddo
      r=sqrt((xc-xcenter(nt))**2+(yc-ycenter(nt))**2)
      pi=3.141593
      if(yc.ge.CY) then
       th=asin(y/r)*180./pi
      else
       th=asin(-y/r)*180./pi
      endif
      call cal_trans(color_dp(nt),ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_col)
      call cal_trans(ut_dp(nt),ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_ut)
      call cal_trans(x_dp(nt),ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_xc)
      call cal_trans(y_dp(nt),ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_yc)
      call cal_trans(r_dp(nt),ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_rad)
      call cal_trans(theta_dp(nt),ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_th)
      call cal_trans(add1_dp(nt),ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_a1)
      call cal_trans(add2_dp(nt),ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_a2)
      call cal_trans(add3_dp(nt),ub,bv,vr,vi,ri,X,UT,xc,yc,r,th,dZ_a3)
      omag=cmag1-(aeK1(nt)-aeK2(nt)*col)*X+zero(nt)
     &+dZ_col+dZ_ut+dZ_xc+dZ_yc+dZ_rad+dZ_th+dZ_a1+dZ_a2+dZ_a3
c      if(xc.eq. 1779.284) then
c       print *,cmag1,aeK1(nt),aeK2(nt),col,X,zero(nt),dZ_col,dZ_ut,dZ_xc
c     & ,dZ_yc,dZ_rad,dZ_th,dZ_a1,dZ_a2,dZ_a3
c      endif
      return
      end
 
c     Color number  1 = U-B, 2 = B-V, 3 = V-R, 4 = V-I, 5 = R-I  
c     Other number 100 = Ut   101=xc,      102=yc,    103=rad,   104=theta
      subroutine cal_trans(input,ub,bv,vr,vi,ri,Air,UT,xc,yc,r,th,yyy)
      character input*16,H*1
      real x(100),y(100),ub,bv,vr,vi,ri,Air,UT,xc,yc,r,th,xxx,yyy
      integer n,npoints,n_depend
      open(31,file=input)
      read(31,*) npoints
      yyy=0. 
      if(npoints.eq.0) goto 900   
      read(31,*) n_depend
      n=0
      do i=1,100
       read(31,'(t1,a1)',end=99) H
       if(H.ne.'#') then
        backspace(31)
        n=n+1
        read(31,*) x(n),y(n) 
       endif
      enddo
99    if(n.ne.npoints) print *,input,' has different points!!!'
      close(31)
      if(n_depend.eq.1) xxx=ub
      if(n_depend.eq.2) xxx=bv
      if(n_depend.eq.3) xxx=vr
      if(n_depend.eq.4) xxx=vi
      if(n_depend.eq.5) xxx=ri
      if(n_depend.eq.100) xxx=UT
      if(n_depend.eq.101) xxx=xc
      if(n_depend.eq.102) xxx=yc
      if(n_depend.eq.103) xxx=r
      if(n_depend.eq.104) xxx=theta
      call lint(xxx,x,y,yyy,n,n)
900   return
      end

      subroutine read_stsi(input)
      parameter(ns=100000)
      character H*1
      common/sts_i/ci,star_i,xc_i(ns),yc_i(ns),ci_ins(ns),cierr(ns),
     &X_i(ns),t_i(ns),ni
      character star_i(ns)*12,ci(ns)*1,line*110,input*16
      open(23,file=input)
      read(23,'(t1,a1)') H
      do i=1,ns
       read(23,'(t1,a1,t3,a12,t13,a110)',end=99) ci(i),star_i(i),line
       read(line,*) xc_i(i),yc_i(i),ci_ins(i),cierr(i),x_i(i),t_i(i)
      enddo
99    close(23)
      ni=i-1
      return
      end

      subroutine read_stsv(input)
      parameter(ns=100000)
      character H*1
      common/sts_v/cv,star_v,xc_v(ns),yc_v(ns),v_ins(ns),verr(ns),
     &X_v(ns),t_v(ns),sv(ns),sbv(ns),svi(ns),sri(ns),nsv(ns),nsbv(ns),
     &nsvi(ns),nsri(ns),nv
      character star_v(ns)*12,cv(ns)*1,line*110,input*16
      open(22,file=input)
      read(22,'(t1,a1)') H
      do i=1,ns
       read(22,'(t1,a1,t3,a12,t13,a110)',end=99) cv(i),star_v(i),line
       read(line,*) xc_v(i),yc_v(i),v_ins(i),verr(i),x_v(i),t_v(i),
     & sv(i),sbv(i),svi(i),sri(i),nsv(i),nsbv(i),nsvi(i),nsri(i)     
      enddo
99    close(22)
      nv=i-1
      return
      end

      subroutine read_stsb(input)
      parameter(ns=100000)
      character H*1
      common/sts_b/cb,star_b,xc_b(ns),yc_b(ns),b_ins(ns),berr(ns),
     &X_b(ns),t_b(ns),nb
      character star_b(ns)*12,cb(ns)*1,line*110,input*16
      open(24,file=input)
      read(24,'(t1,a1)') H
      do i=1,ns
       read(24,'(t1,a1,t3,a12,t13,a110)',end=99) cb(i),star_b(i),line
       read(line,*) xc_b(i),yc_b(i),b_ins(i),berr(i),x_b(i),t_b(i)
      enddo
99    close(24)
      nb=i-1
      return
      end

      subroutine read_coeff(input)
      parameter(nf=20,ns=100000)
      character H*1,input*16
      common/coeff_const/num_filter(nf),num_color(nf),aeK1(nf),aeK2(nf)
     &,zero(nf),xcenter(nf),ycenter(nf),color_dp,ut_dp,x_dp,y_dp,r_dp,
     &theta_dp,add1_dp,add2_dp,add3_dp,nconst
      character color_dp(nf)*16,ut_dp(nf)*16,x_dp(nf)*16,y_dp(nf)*16
     &,r_dp(nf)*16,theta_dp(nf)*16,add1_dp(nf)*16,add2_dp(nf)*16
     &,add3_dp(nf)*16
      open(21,file=input)
      nconst=0 
      do i=1,ns
       read(21,'(t1,a1)',end=99) H
       if(H.ne.'#') then
        backspace(21)
        nconst=nconst+1
        read(21,*) num_filter(nconst),num_color(nconst),aeK1(nconst),
     &  aeK2(nconst),zero(nconst),xcenter(nconst),ycenter(nconst),
     &  color_dp(nconst),ut_dp(nconst),x_dp(nconst),y_dp(nconst),
     &  r_dp(nconst),theta_dp(nconst),add1_dp(nconst),add2_dp(nconst),
     &  add3_dp(nconst)
       endif
      enddo
99    close(21)
      return
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
12     y=yt1(i)
       goto 14
10    continue
11    y=yt1(k)+(yt1(k+1)-yt1(k))*((x-xt1(k))/(xt1(k+1)-xt1(k)))
14    return
      end
