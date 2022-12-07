  do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
     call LBFGS(n,x,fcost,grad,optim,FLAG)
     if(FLAG.eq.'GRAD') then        
        !compute cost and gradient at point x
        call Rosenbrock(x,fcost,grad)        
     endif
  enddo
