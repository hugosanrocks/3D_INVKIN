  do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))     
     call PSTD(n,x,fcost,grad,grad_preco,optim,FLAG)
     if(FLAG.eq.'GRAD') then        
        !compute cost and gradient at point x
        call Rosenbrock(x,fcost,grad)
        !no preconditioning
        grad_preco(:)=grad(:)
     endif
  enddo
  


