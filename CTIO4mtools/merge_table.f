c     f95 -o merge_table merge_table.f
      character input*32,output*32,H*1
      open(21,file='merge_table.input')
      read(21,'(t1,a1)') 
      read(21,'(t1,a32)') output
      open(22,file=output)
      write(22,'(a135)') '# Star             Xc        Yc      Mag      
     &Merr  Airmass UT             V       U-B     B-V      V-I    R-I  
     &nV nUB nBV nVI nRI                                               '
      do i=1,1000
       read(21,'(t1,a1,t1,a32)',end=99) H,input
       if(H.ne.'#') call write_table(input,22)
      enddo
99    close(22)
      close(21)
      return
      end 
      subroutine write_table(input,nnn)
      character input*32,H*1,line*135
      open(23,file=input)
      do i=1,100000
       read(23,'(t1,a1,t1,a135)',end=99) H,line 
c       print *,line
       if(H.ne.'#') write(nnn,'(t1,a135)') line
      enddo
99    close(23)
      return
      end
