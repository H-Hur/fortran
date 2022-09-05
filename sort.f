      subroutine tndsort(x,n)
      real x(100,2),x1,x2
      integer n
      do i=1,n
       do j=1,n-1
        if(x(j+1,1).lt.x(j,1)) then
         x1=x(j,1)
         x2=x(j,2)
         x(j,1)=x(j+1,1)
         x(j,2)=x(j+1,2)
         x(j+1,1)=x1
         x(j+1,2)=x2
        endif
       enddo
      enddo
      return
      end

      subroutine datasort(val,ch,nstar)
      integer nstar
      real val(nstar),ref
      character ch(nstar)*250,chref*250
      do i=1,nstar
       do j=1,nstar-1
        if(val(j+1).lt.val(j)) then
         ref=val(j)
         write(chref,'(a250)') ch(j)
         val(j)=val(j+1)
         write(ch(j),'(a250)') ch(j+1)
         val(j+1)=ref
         write(ch(j+1),'(a250)') chref
        endif
       enddo
      enddo
      return
      end

      subroutine sort(ns,E)
      real C,E(ns)
      integer i,j
      do i=1,ns
       do j=1,ns-1
        if (E(j+1).lt.E(j)) then
         C=E(j)
         E(j)=E(j+1)
         E(j+1)=C
        endif
       enddo
      enddo
      return
      end

      subroutine sort2(ns,E1,E2)
      real C1,C2,E1(ns),E2(ns)
      integer i,j
      do i=1,ns
       do j=1,ns-1
        if (E1(j+1).lt.E1(j)) then
         C1=E1(j)
         C2=E2(j)
         E1(j)=E1(j+1)
         E2(j)=E2(j+1)
         E1(j+1)=C1
         E2(j+1)=C2
        endif
       enddo
      enddo
      return
      end

