  do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
     call TRN(n,x,fcost,grad,optim,FLAG)
     if(FLAG.eq.'GRAD') then
        call Rosenbrock(x,fcost,grad)        
     elseif(FLAG.eq.'HESS') then        
        call Rosenbrock_Hess(x,optim%d,optim%Hd)
     endif
  enddo
  


