c 행렬m2=m1,s2=s1가 되게 하는 subroutine. 
       subroutine recall(m1,m2,s1,s2)
       real*8 m1(2,2),m2(2,2),s1(2,1),s2(2,1)
       do i=1,2
        s2(i,1)=s1(i,1)
        do j=1,2
        m2(i,j)=m1(i,j)
        enddo
       enddo
       return
       end

c 행렬과 행렬식을 print.
       subroutine printing(m2,s1)
       real*8 m2(2,2),s1(2,1),meq
       call mateq(m2,meq)
       print *,'|',m2(1,1),m2(1,2),'| =',s1(1,1)
       print *,'|',m2(2,1),m2(2,2),'| =',s1(2,1)
       print *,'meq=',meq
       print *,''
       return
       end

c 행렬식을 계산.
       subroutine mateq(matrix,eq)
       real*8 matrix(2,2),eq
       eq=matrix(1,1)*matrix(2,2)-matrix(1,2)*matrix(2,1)
       return
       end

c 행렬식에 크래머공식을 적용(사용후 반드시 마지막에 recall루트를 넣을것.)
       subroutine cr2(matrix,change,l)
       integer i,l
       real*8 matrix(2,2),change(2,1),slot(2,1)
       do i=1,2
        slot(i,1)=matrix(i,l)
        matrix(i,l)=change(i,1)
        change(i,1)=slot(i,1)
       enddo
       return
       end

