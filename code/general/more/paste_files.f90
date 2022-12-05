         program cutc
         real dat(8192*144*40)
         integer i, nsubf, nt, nsta, srec

         nsta = 40                 !stations
         nt   = 8192               !time samples
         nsubf= 144                !subfaults per file
         srec = nt*nsta*nsubf      !size of record


         open(10,file='src1/SIGMA_XX1_C1',status='old',access='DIRECT',recl=srec*4) 
         open(11,file='src2/SIGMA_XX2_C1',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_XX_C1',status='unknown',access='DIRECT',recl=srec*4)
        
         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)
     

         open(10,file='src1/SIGMA_YY1_C1',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_YY2_C1',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_YY_C1',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)


         open(10,file='src1/SIGMA_ZZ1_C1',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_ZZ2_C1',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_ZZ_C1',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)

         open(10,file='src1/SIGMA_XY1_C1',status='old',access='DIRECT',recl=srec*4) 
         open(11,file='src2/SIGMA_XY2_C1',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_XY_C1',status='unknown',access='DIRECT',recl=srec*4)
        
         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)
     

         open(10,file='src1/SIGMA_XZ1_C1',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_XZ2_C1',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_XZ_C1',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)


         open(10,file='src1/SIGMA_YZ1_C1',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_YZ2_C1',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_YZ_C1',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)











         open(10,file='src1/SIGMA_XX1_C2',status='old',access='DIRECT',recl=srec*4) 
         open(11,file='src2/SIGMA_XX2_C2',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_XX_C2',status='unknown',access='DIRECT',recl=srec*4)
        
         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)
     

         open(10,file='src1/SIGMA_YY1_C2',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_YY2_C2',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_YY_C2',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)


         open(10,file='src1/SIGMA_ZZ1_C2',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_ZZ2_C2',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_ZZ_C2',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)

         open(10,file='src1/SIGMA_XY1_C2',status='old',access='DIRECT',recl=srec*4) 
         open(11,file='src2/SIGMA_XY2_C2',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_XY_C2',status='unknown',access='DIRECT',recl=srec*4)
        
         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)
     

         open(10,file='src1/SIGMA_XZ1_C2',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_XZ2_C2',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_XZ_C2',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)


         open(10,file='src1/SIGMA_YZ1_C2',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_YZ2_C2',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_YZ_C2',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)











         open(10,file='src1/SIGMA_XX1_C3',status='old',access='DIRECT',recl=srec*4) 
         open(11,file='src2/SIGMA_XX2_C3',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_XX_C3',status='unknown',access='DIRECT',recl=srec*4)
        
         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)
     

         open(10,file='src1/SIGMA_YY1_C3',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_YY2_C3',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_YY_C3',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)


         open(10,file='src1/SIGMA_ZZ1_C3',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_ZZ2_C3',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_ZZ_C3',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)

         open(10,file='src1/SIGMA_XY1_C3',status='old',access='DIRECT',recl=srec*4) 
         open(11,file='src2/SIGMA_XY2_C3',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_XY_C3',status='unknown',access='DIRECT',recl=srec*4)
        
         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)
     

         open(10,file='src1/SIGMA_XZ1_C3',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_XZ2_C3',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_XZ_C3',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)


         open(10,file='src1/SIGMA_YZ1_C3',status='old',access='DIRECT',recl=srec*4)       
         open(11,file='src2/SIGMA_YZ2_C3',status='old',access='DIRECT',recl=srec*4)
         open(12,file='SIGMA_YZ_C3',status='unknown',access='DIRECT',recl=srec*4)

         read(10,rec=1) dat
         close(10)
         write(12,rec=1) dat
         read(11,rec=1) dat
         close(11)
         write(12,rec=2) dat
         close(10)






         endprogram cutc
