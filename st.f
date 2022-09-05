!    for fortran95
!    to transform observation times on the output file from aperture photomerty !    by txdump task.

      integer h,m
      real*8 s,time
      character line1*30,pe*9,input*16,output*16

      print *,'input st file=?'
      read(*,'(a16)') input
      print *,'output st file=?'
      read(*,'(a16)') output
      open(21,file=input)
      open(22,file=output)
  
100   format(t1,a30,t32,i2,t35,i2,t38,f6.3,t45,a9)
101   format(t1,a30,t32,f12.9,t47,a9)
      do i=1,1000
       read(21,100,end=200) line1,h,m,s,pe
       time=real(h)+real(m)/60.d0+s/3600.d0
       write(22,101) line1,time,pe
      enddo
200   i=i
      stop
      en

