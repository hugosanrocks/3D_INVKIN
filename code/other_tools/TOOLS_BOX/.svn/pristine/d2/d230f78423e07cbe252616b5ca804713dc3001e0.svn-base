  do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
     call PLBFGS(n,x,fcost,grad,grad_preco,optim,FLAG)
     if(FLAG.eq.'GRAD') then        
        call Rosenbrock(x,fcost,grad)        
     endif
     if(FLAG.eq.'PREC') then
        !apply preconditioning to optim%q_plb
        !if nothing is done, PLBFGS is equivalent to LBFGS
        optim%q_plb(:)=optim%q_plb(:)
     endif
  enddo
