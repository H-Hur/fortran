       real*8 a,mean,S(9999),sv,std,S2(9999)
       integer m
      print *,'Input numbers to compute mean & STDDEV.'
      print *,'  '
      print *,'Input the 99999 to compute inputed numbers.'
      print *,'  '
      print *,'Input the -99999 to end.'   
       
       do j=1,9999
        do i=1,9999
        S(i)=0.d0 
        S2(i)=0.d0
        enddo
       a=0.d0
       sv=0.d0
       m=0


        do i=1,9999
        read (*,*) S(i)
         
         if (S(i).ge.99990.) then
         goto 100
         elseif (S(i).lt.-99990.) then
         goto 1000
         endif
        a=a+S(i)
        S2(i)=S(i) 
        m=m+1
        enddo !i
100    mean=a/real(m)

        do k=1,m
        sv=sv+(mean-S2(k))**2.d0
        enddo !k
       std=sqrt(sv/real(m-1))
        
       print *,'Nnumber=',m,' ','mean=',mean,' ','STDDEV=',std
       enddo !j
1000   i=i
       stop
       end
                 
