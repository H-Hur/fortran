      subroutine cid_coo2(x1,y1,n1,x2,y2,n2,rad,
     &np1,np2,num1to2,num2to1)
      parameter(ns=200000)
      integer n1,n2,np1(n1),np2(n2),num1to2(ns,10),num2to1(ns,10)
      real x1(n1),y1(n1),x2(n2),y2(n2),rad,dist
      print *,n1,n2,'@@@'
      do j=1,n2
       np2(j)=0
      enddo
      print *,n1,n2,'@@@'
      do i=1,n1
       np1(i)=0
       do j=1,n2
        dist=sqrt((x1(i)-x2(j))**2+(y1(i)-y2(j))**2)
        if(dist.le.rad) then
         np1(i)=np1(i)+1
         np2(j)=np2(j)+1
         num1to2(i,np1(i))=j
         num2to1(j,np2(j))=i
        endif
       enddo
      enddo
      print *,n1,n2,'###'
      return
      end

      subroutine cid_coo_is_db(x1,y1,e1,n1,x2,y2,e2,n2,
     &rad_min,rad_max,th,
     &num1to2,num2to1,num_single,
     &nsingle,ndouble1,ndouble2)
      integer n1,n2,num1to2(n1),num2to1(n2),nsingle,ndouble1,ndouble2
      integer np1(n1),np2(n2),num_single(n1)
      real*8 x1(n1),y1(n1),x2(n2),y2(n2),dist,close1(n1),close2(n2)
      real*8 e1(n1),e2(n2),rad_min,rad_max,r_t,th
      do j=1,n2
       num2to1(j)=0
       np2(j)=0
       close2(j)=rad_max
      enddo
      do i=1,n1
       num1to2(i)=0
       np1(i)=0
       close1(i)=rad_max
       do j=1,n2
        dist=dsqrt((x1(i)-x2(j))**2+(y1(i)-y2(j))**2)
        r_t=th*dsqrt(e1(i)**2+e2(j)**2)
        if(r_t.gt.rad_max) r_t=rad_max
        if(r_t.gt.rad_min) r_t=rad_min
        if(dist.le.r_t) then
         np1(i)=np1(i)+1
         np2(j)=np2(j)+1
         if(np1(i).eq.1.and.np2(j).eq.1) then
          close1(i)=dist
          close2(j)=dist
          num1to2(i)=j
          num2to1(j)=i
         elseif(np1(i).eq.1.and.np2(j).ge.2.and.dist.lt.close2(j)) then
          close1(i)=dist
          close2(j)=dist
          num1to2(num2to1(j))=0
          num1to2(i)=j
          num2to1(j)=i
         elseif(np1(i).ge.2.and.np2(j).eq.1.and.dist.lt.close1(i)) then
          close1(i)=dist
          close2(j)=dist
          num2to1(num1to2(i))=0
          num2to1(j)=i
          num1to2(i)=j
         elseif(np1(i).ge.2.and.np2(j).ge.2.and.
     &   dist.lt.close1(i).and.dist.lt.close2(j)) then
          close1(i)=dist
          close2(j)=dist
          num1to2(num2to1(j))=0
          num2to1(num1to2(i))=0
          num2to1(j)=i
          num1to2(i)=j
         endif
        endif
       enddo
      enddo
      nsingle=0
      ndouble1=0
      ndouble2=0
      do i=1,n1
       num_single(i)=0
       if(np1(i).eq.1.and.np2(num1to2(i)).eq.1) nsingle=nsingle+1
       if(np1(i).eq.1.and.np2(num1to2(i)).eq.1) num_single(i)=1
       if(np1(i).ge.2) ndouble1=ndouble1+1
      enddo
      do j=1,n2
       if(np2(j).ge.2) ndouble2=ndouble2+1
      enddo
      return
      end


      subroutine cid_coo_db(x1,y1,n1,x2,y2,n2,rad,
     &num1to2,num2to1,num_single,
     &nsingle,ndouble1,ndouble2)
      integer n1,n2,num1to2(n1),num2to1(n2),nsingle,ndouble1,ndouble2
      integer np1(n1),np2(n2),num_single(n1)
      real*8 x1(n1),y1(n1),x2(n2),y2(n2),rad,dist,close1(n1),close2(n2)
      do j=1,n2
       num2to1(j)=0
       np2(j)=0
       close2(j)=rad
      enddo
      do i=1,n1
       num1to2(i)=0
       np1(i)=0
       close1(i)=rad
       do j=1,n2
        dist=dsqrt((x1(i)-x2(j))**2+(y1(i)-y2(j))**2)
        if(dist.le.rad) then
         np1(i)=np1(i)+1
         np2(j)=np2(j)+1
         if(np1(i).eq.1.and.np2(j).eq.1) then
          close1(i)=dist
          close2(j)=dist
          num1to2(i)=j  
          num2to1(j)=i  
         elseif(np1(i).eq.1.and.np2(j).ge.2.and.dist.lt.close2(j)) then
          close1(i)=dist
          close2(j)=dist
          num1to2(num2to1(j))=0
          num1to2(i)=j
          num2to1(j)=i
         elseif(np1(i).ge.2.and.np2(j).eq.1.and.dist.lt.close1(i)) then
          close1(i)=dist
          close2(j)=dist
          num2to1(num1to2(i))=0
          num2to1(j)=i
          num1to2(i)=j
         elseif(np1(i).ge.2.and.np2(j).ge.2.and.
     &   dist.lt.close1(i).and.dist.lt.close2(j)) then
          close1(i)=dist
          close2(j)=dist
          num1to2(num2to1(j))=0
          num2to1(num1to2(i))=0
          num2to1(j)=i
          num1to2(i)=j      
         endif
        endif
       enddo
      enddo
      nsingle=0
      ndouble1=0
      ndouble2=0
      do i=1,n1
       num_single(i)=0
       if(np1(i).eq.1.and.np2(num1to2(i)).eq.1) nsingle=nsingle+1
       if(np1(i).eq.1.and.np2(num1to2(i)).eq.1) num_single(i)=1
       if(np1(i).ge.2) ndouble1=ndouble1+1       
      enddo
      do j=1,n2
       if(np2(j).ge.2) ndouble2=ndouble2+1       
      enddo
      return
      end
      
      subroutine cid_coo(x1,y1,n1,x2,y2,n2,rad,
     &num1to2,num2to1,num_single,
     &nsingle,ndouble1,ndouble2)
      integer n1,n2,num1to2(n1),num2to1(n2),nsingle,ndouble1,ndouble2
      integer np1(n1),np2(n2),num_single(n1)
      real x1(n1),y1(n1),x2(n2),y2(n2),rad,dist,close1(n1),close2(n2)
      do j=1,n2
       num2to1(j)=0
       np2(j)=0
       close2(j)=rad
      enddo
      do i=1,n1
       num1to2(i)=0
       np1(i)=0
       close1(i)=rad
       do j=1,n2
        dist=sqrt((x1(i)-x2(j))**2+(y1(i)-y2(j))**2)
        if(dist.le.rad) then
         np1(i)=np1(i)+1
         np2(j)=np2(j)+1
         if(np1(i).eq.1.and.np2(j).eq.1) then
          close1(i)=dist
          close2(j)=dist
          num1to2(i)=j  
          num2to1(j)=i  
         elseif(np1(i).eq.1.and.np2(j).ge.2.and.dist.lt.close2(j)) then
          close1(i)=dist
          close2(j)=dist
          num1to2(num2to1(j))=0
          num1to2(i)=j
          num2to1(j)=i
         elseif(np1(i).ge.2.and.np2(j).eq.1.and.dist.lt.close1(i)) then
          close1(i)=dist
          close2(j)=dist
          num2to1(num1to2(i))=0
          num2to1(j)=i
          num1to2(i)=j
         elseif(np1(i).ge.2.and.np2(j).ge.2.and.
     &   dist.lt.close1(i).and.dist.lt.close2(j)) then
          close1(i)=dist
          close2(j)=dist
          num1to2(num2to1(j))=0
          num2to1(num1to2(i))=0
          num2to1(j)=i
          num1to2(i)=j      
         endif
        endif
       enddo
      enddo
      nsingle=0
      ndouble1=0
      ndouble2=0
      do i=1,n1
       num_single(i)=0
       if(np1(i).eq.1.and.np2(num1to2(i)).eq.1) nsingle=nsingle+1
       if(np1(i).eq.1.and.np2(num1to2(i)).eq.1) num_single(i)=1
       if(np1(i).ge.2) ndouble1=ndouble1+1       
      enddo
      do j=1,n2
       if(np2(j).ge.2) ndouble2=ndouble2+1       
      enddo
      return
      end
     
      subroutine cal_draddec(ra,dec,cra,cdec,dra,ddec)
      real*8 ra,dec,cra,cdec,dra,ddec,pi
      pi=3.14159265358979d0
      dra=(ra-cra)*60.d0*dcos(dec*pi/180.d0)
      ddec=(dec-cdec)*60.d0
      return
      end 
 
