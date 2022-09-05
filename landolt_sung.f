        subroutine SSO_Landolt(rmg,err,nobs,filter,nx,nflt,nf,nfilter,
     &                   smg,color,ecolr,nobs_c,ncol,col_name,xair,nsmg)
c
c       *** updated on 15/04/2005 using tr_Landolt
        real rmg(nx,nflt),smg(nx,nflt),color(nx,nflt),ecolr(nx,nflt),
     &       err(nx,nflt),ref_bv(21),ref_u0(21),xair(nflt)
        integer nobs(nx,nflt),nobs_c(nx,nflt),ncol
        character*1 filter(nflt),f1(5),f2(5)
        character*3 col_name(nflt)
        data ref_bv/-.35,-.33,-0.3,-0.2,-0.1, 0.0, 0.1, 0.2, 0.3, 0.4,
     *     0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0/
        data ref_u0/-.132,-.132,-.126,-.074,-.030,0.002,0.024,0.026,
     *              0.010,-.017,-.060,-.104,-.132,-.132,-.132,-.132,
     *              -.132,-.132,-.132,-.132,-.132/
15      format(a1,1x,a1)
25      format(a3)
        do i = 1, nfilter
           print*,'Filter ',i,'=',filter(i)
           if ( filter(i) . eq . 'V' ) then
              nv = i
           elseif ( filter(i) . eq . 'I' ) then
              ni = i
           elseif ( filter(i) . eq . 'B' ) then
              nb = i
           elseif ( filter(i) . eq . 'U' ) then
              nu = i
           elseif ( filter(i) . eq . 'R' ) then
              nr = i
           elseif ( filter(i) . eq . 'C' ) then
              nc = i
           elseif ( filter(i) . eq . 'H' ) then
              nh = i
           endif
        enddo
        print*,'                       '
        print*,'What is the number of colors?'
        read(*,*)ncol
        do 110 j = 1, ncol
           print*,'Print the name of color ',j
           read(*,'(a3)') col_name(j)
           if ( col_name(j) .eq. 'V-I' ) then
              print*,'Print the zero point for V wrt (V-I)!'
              read(*,*)zv2
              print*,'Print the zero point for I!'
              read(*,*)zi0
              nvi = j
              do 20 i = 1, nf
                 if ( nobs(i,nv) .eq. 0 .or. nobs(i,ni) .eq. 0 ) goto 20
                 vi = rmg( i, nv ) - rmg( i, ni )
                 nnn = 1
10               if ( vi .ge. 0.8 ) then
                    corr_i = zi0 + 0.0361 * ( vi - 0.8 )
                 else
                    corr_i = zi0 + 0.0
                 endif
                 corr_v = zv2 + 0.0647 * vi
                 smg( i, ni ) = rmg( i, ni ) + corr_i
                 smg( i, nv ) = rmg( i, nv ) + corr_v
                 color( i, j ) = smg( i, nv ) - smg( i, ni )
                 delta = abs( color( i, j ) - vi )
                 if ( delta .gt. 0.001 .and. nnn .le. 10 ) then
                    nnn = nnn + 1

                    vi = color( i, j )
                    goto 10
                 endif
                 ecolr( i, j ) = sqrt( err( i, nv ) * err( i, nv )
     &                          + err( i, ni ) * err( i, ni ) )
                 nobs_c( i, j ) = nmin( nobs( i, nv ),nobs( i, ni ) )
20            enddo
           elseif ( col_name( j ) . eq . 'B-V' ) then
              print*,'Print the zero point for V wrt (B-V)!'
              read(*,*)zv1
              print*,'Print the zero point for B!'
              read(*,*)zb0
              nbv = j
              do 40 i = 1, nf
                 if ( nobs(i,nv) .eq. 0 .or. nobs(i,nb) .eq. 0 ) goto 40
                 bv = rmg( i, nb ) - rmg( i, nv )
30               smg( i, nv ) = rmg( i, nv ) + 0.0693* bv + zv1
                 corr_b = zb0 - 0.1006 * bv + 0.0314 * bv * xair( nb )
                 smg( i, nb ) = rmg( i, nb ) + corr_b
                 color( i, j ) = smg( i, nb ) - smg( i, nv )
                 delta = abs( color( i, j ) - bv )
                 if ( delta . gt . 0.001 ) then
                    bv = color(i,j)
                    goto 30
                 endif
                 ecolr( i, j ) = sqrt( err( i, nv ) * err( i, nv )
     &                            + err( i, nb ) * err( i, nb ) )
                 nobs_c( i, j ) = nmin( nobs( i, nv ),nobs( i, nb ) )
40            enddo
           elseif ( col_name( j ) . eq . 'U-B' ) then
              print*,'Print the zero point for U!'
              read(*,*)zu0
              print*,'Print the mean reddening of the field?'
              read(*,*)ebv_fg
              nub = j
              do 60 i = 1, nf
                 if ( nobs(i,nu) .eq. 0 .or. nobs(i,nb) .eq. 0
     &                .or. nobs(i,nv) .eq. 0 ) goto 60
                 nnn = 1
                 bv0 = color( i, nbv ) - ebv_fg
                 if ( bv0 .lt. -0.3 ) bv0 = -0.3
                 call lint(bv0,ref_bv,ref_u0,del_u0,21,21)
                 smg( i, nu ) = rmg( i, nu ) + del_u0 + zu0
                 ub = smg( i, nu ) - smg( i, nb )
50               qqq = ub - 0.72 * color( i, nbv )
                 if ( qqq .lt. -0.05 .and. ub .lt. -0.05 ) then
                    bv0 = 0.332 * qqq
                 else
                    bv0 = color( i, nbv ) - ebv_fg
                 endif
                 call lint(bv0,ref_bv,ref_u0,del_u0,21,21)
                 corr_u = zu0 + del_u0 + (0.1023 + 0.0129*xair(nu)) * ub
                 smg( i, nu ) = rmg( i, nu ) + corr_u
                 color( i, j ) = smg( i, nu ) - smg( i, nb )
                 delta = abs( color( i, j ) - ub )
                 if ( delta .gt. 0.001 .and. nnn .lt. 10 ) then
                    ub = color( i, j )
                    nnn = nnn + 1
                    goto 50
                 endif
                 ecolr( i, j ) = sqrt( err( i, nu ) * err( i, nu )
     &                          + err( i, nb ) * err( i, nb ) )
                 nobs_c( i, j ) = nmin( nobs( i, nb ),nobs( i, nu ) )
60            enddo
           elseif ( col_name( j ) . eq . 'R-H' ) then
              print*,'Print zero point for R!'
              read(*,*)zr2
              print*,'Print zero point for H_alpha!'
              read(*,*)zh0
              nha = j
              do 70 i = 1, nf
                 if ( nobs(i,nr) .eq. 0 .or. nobs(i,nh) .eq. 0 ) goto 70
                 color( i, j ) = rmg(i,nr) - rmg(i,nh) + zr2 - zh0
                 ecolr( i, j ) = sqrt( err( i, nr ) * err( i, nr )
     &                            + err( i, nh ) * err( i, nh ) )
                 nobs_c( i, j ) = nmin( nobs( i, nr ),nobs( i, nh ) )
70            enddo
           elseif ( col_name( j ) . eq . 'R-I' ) then
              print*,'Print zero point for R!'
              read(*,*)zr1
              if ( zi0 .eq. 0.0 ) then
                 print*,'Print zero point for I!'
                 read(*,*)zi0
                 nnni = 1
              endif
              nri = j
              do 90 i = 1, nf
                 if ( nobs(i,nr) .eq. 0 .or. nobs(i,ni) .eq. 0 ) goto 90
                 nnn = 1
                 if ( nnni .eq. 1 ) then
                    smg( i, ni ) = rmg( i, ni ) + zi0
                 endif
                 ri = rmg( i, nr ) + zr1 - smg( i, ni )
                 smg( i, nr ) = rmg( i, nr ) + zr1
80               if ( nnni .eq. 1 ) then
                    if ( ri .le. 0.4 ) then
                       corr_i = zi0
                    else
                       corr_i = zi0 + 0.0735 * ri
                    endif
                    smg( i, ni ) = rmg( i, ni ) + corr_i
                 endif
                 color( i, j ) = smg( i, nr ) - smg( i, ni )
                 delta = abs( color( i, j ) - ri )
                 if ( delta .gt. 0.001 .and. nnn .lt. 10 ) then
                    nnn = nnn + 1
                    ri = color( i, j )
                    goto 80
                 endif
                 ecolr( i, j ) = sqrt( err( i, nr ) * err( i, nr )
     &                            + err( i, ni ) * err( i, ni ) )
                 nobs_c( i, j ) = nmin( nobs( i, nr ),nobs( i, ni ) )
90            enddo
           elseif ( col_name( j ) . eq . 'H-C' ) then
              if ( nc .ne. 0 ) then
                 print*,'Print zero point for C!'
                 read(*,*)zc0
              else
                 zc0 = 0.0
              endif
              print*,'Print zero point for H_alpha!'
              read(*,*)zh0
              nha = j
              do 100 i = 1, nf
                 if(nc.eq.0.and.nobs(i,ni).ne.0.and.nobs(i,nv).ne.0)then
                    continuum = ( smg(i,ni) + smg(i,nv) ) / 2.
                    err_cont  = ecolr(i,nvi)
                    nobs_cont = nmin( nobs(i,nv),nobs(i,ni) )
                 elseif ( nc .eq. 0 ) then
                    nobs_cont = nobs( i, nh )
                 else
                    continuum = rmg(i,nc) - zc0
                    err_cont  = err(i,nc)
                    nobs_cont = nobs(i,nc)
                 endif
                 if ( nobs_cont .eq. 0 .or. nobs(i,nh) .eq. 0 ) goto 100

                 color( i, j ) = rmg(i,nh) - continuum + zh0
                 ecolr( i, j ) = sqrt( err_cont * err_cont
     &                            + err( i, nh ) * err( i, nh ) )
                 nobs_c( i, j ) = nmin( nobs_cont, nobs( i, nh ) )
100           enddo
           endif
110     enddo
c
c       calculate standard magnitudes
c
        print*,'-----------------------------------------------------'
        print*,'                 Standard Magnitude'
        print*,'-----------------------------------------------------'
        do i = 1, nfilter
           print*,'Filter ',i,'=',filter(i)
        enddo
        print*,' select filter number for standard magnitude!'
        read(*,*)nsmg
        if ( filter(nsmg) .eq. 'V' ) then
           do i = 1, nf
              if ( nobs(i,nv) .ne. 0 ) then
                 if ( nobs_c(i,nvi) .eq. 0 ) then
                    wt_vi = 0.0
                 elseif ( ecolr(i,nvi) .ge. 0.1 ) then
                    wt_vi = 0.1
                 else
                    wt_vi = 1.0
                 endif
                 if ( nobs_c(i,nbv) .eq. 0 ) then
                    wt_bv = 0.0
                 elseif ( ecolr(i,nbv) .ge. 0.1 ) then
                    wt_bv = 0.1
                 else
                    wt_bv = 1.0
                 endif
                 sm1 = rmg( i, nv ) + 0.0693 * color( i, nbv ) + zv1
                 sm2 = rmg( i, nv ) + 0.0647 * color( i, nvi ) + zv2
                 if ( (wt_bv+wt_vi) .ne. 0.0 ) then
                 smg( i, nv ) = ( sm1*wt_bv + sm2*wt_vi )/(wt_bv+wt_vi)
                 else
                 smg( i, nv ) = rmg( i, nv ) + ( zv1 + zv2 ) / 2.
                 endif
              else
                 smg( i, nv ) = 99.999
              endif
           enddo
        elseif ( filter(nsmg) .eq. 'I' ) then
           do i = 1, nf
              if ( nobs( i, ni ) .ne. 0 ) then
                 if ( nobs_c(i,nvi) .eq. 0 ) then
                    wt_vi = 0.0
                 elseif ( ecolr(i,nvi) .ge. 0.1 ) then
                    wt_vi = 0.1
                 else
                    wt_vi = 1.0
                 endif
                 if ( nobs_c(i,nri) .eq. 0 ) then
                    wt_ri = 0.0
                 elseif ( ecolr(i,nri) .ge. 0.1 ) then
                    wt_ri = 0.1
                 else
                    wt_ri = 1.0
                 endif
                 if ( color(i,nvi) .gt. 0.8 ) then
                    cor_i1 = 0.0361 * ( color(i,nvi) - 0.8 )
                 else
                    cor_i1 = 0.0
                 endif
3                 if ( color(i,nri) .gt. 0.4 ) then
                    cor_i2 = 0.0735 * ( color(i,nri) - 0.4 )
                 else
                    cor_i2 = 0.0
                 endif
                 sm1 = rmg( i, ni ) + cor_i2 + zi0
                 sm2 = rmg( i, ni ) + cor_i1 + zi0
                 if ( wt_ri + wt_vi .ne. 0.0 ) then
                 smg( i, ni ) = ( sm1*wt_ri + sm2*wt_vi )/(wt_ri+wt_vi)
                 else
                 smg( i, ni ) = rmg( i, ni ) + zi0
                 endif
              else
                 smg( i, ni ) = 99.999
              endif
           enddo
        endif
        return
        end

