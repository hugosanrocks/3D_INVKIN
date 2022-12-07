  do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))     
     call PSTD(n,x,fcost,grad,optim,FLAG)
     if(FLAG.eq.'GRAD') then        
        !compute cost and gradient at point x
        call Rosenbrock(x,fcost,grad)
     elseif(FLAG.eq.'NSTE') then
        write(*,*) x(:) ! or save it into disk or... 
     endif
  enddo
  


