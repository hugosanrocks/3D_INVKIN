  do while ((FLAG.ne.'CONV').and.(FLAG.ne.'FAIL'))
     call PTRN(n,x,fcost,grad,grad_preco,optim,FLAG)
     if(FLAG.eq.'GRAD') then        
        call Rosenbrock(x,fcost,grad)        
        grad_preco(:)=grad(:)
     elseif(FLAG.eq.'HESS') then
        call Rosenbrock_Hess(x,optim%d,optim%Hd)
     elseif(FLAG.eq.'PREC') then
        optim%residual_preco(:)=optim%residual(:)
     endif
  enddo
  
