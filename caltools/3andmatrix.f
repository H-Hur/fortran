       real*8 t31(3,3),res1(3,1),t32(3,3),res2(3,1),eq(3),r
       open(1,file='3t.dat')
       do i=1,3
       read(1,*) t31(i,1),t31(i,2),t31(i,3),res1(i,1)
       enddo
       call recall(t31,t32,res1,res2)
       call cr2(t32,res2,1)
       call mateq(t32,r)
       call recall(t31,t32,res1,res2)
       print *,r
       call cr2(t32,res2,2)
       call mateq(t32,r)
       call recall(t31,t32,res1,res2)
       print *,r
       call cr2(t32,res2,3)
       call mateq(t32,r)
       call recall(t31,t32,res1,res2)
       print *,r

       close(1)
       stop
       end




c 행렬m2=m1,s2=s1가 되게 하는 subroutine. 
       subroutine recall(m1,m2,s1,s2)
       real*8 m1(3,3),m2(3,3),s1(3,1),s2(3,1)
       do i=1,3
        s2(i,1)=s1(i,1)
        do j=1,3
        m2(i,j)=m1(i,j)
        enddo
       enddo
       return
       end

c 행렬과 행렬식을 print.
       subroutine printing(m2,s1)
       real*8 m2(3,3),s1(3,1),meq
       call mateq(m2,meq)
       print *,'|',m2(1,1),m2(1,2),m2(1,3),'|  ',s1(1,1)
       print *,'|',m2(2,1),m2(2,2),m2(2,3),'| =',s1(2,1)
       print *,'|',m2(3,1),m2(3,2),m2(3,3),'|  ',s1(3,1)
       print *,'meq=',meq
       print *,''
       return
       end

c 행렬식을 계산.
       subroutine mateq(m,eq)
       real*8 m(3,3),eq
       eq=m(1,1)*m(2,2)*m(3,3)+m(1,2)*m(2,3)*m(3,1)+m(1,3)*m(2,1)*m(3,2)
     *-m(1,2)*m(2,1)*m(3,3)-m(1,1)*m(2,3)*m(3,2)-m(1,3)*m(2,2)*m(3,1)
       return
       end

c 행렬식에 크래머공식을 적용(사용후 반드시 마지막에 recall루트를 넣을것.)
       subroutine cr2(matrix,change,l)
       integer i,l
       real*8 matrix(3,3),change(3,1),slot(3,1)
       do i=1,3
        slot(i,1)=matrix(i,l)
        matrix(i,l)=change(i,1)
        change(i,1)=slot(i,1)
       enddo
       return
       end

