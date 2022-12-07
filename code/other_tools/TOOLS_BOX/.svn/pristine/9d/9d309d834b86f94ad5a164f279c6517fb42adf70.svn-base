!*********************************************!
!*    SEISCOPE OPTIMIZATION TOOLBOX          *!
!*********************************************!
! These routine encode the Rosenbrock         !
! and its derivatives for basic for tests of  !
! the minimization algorithms                 !
!---------------------------------------------!
! The Rosenbrock function is defined by       !
! f(x1,x2)  = (1-x1)**2+100.*(x2-x1**2)**2    !
!---------------------------------------------!
! The gradient of the Rosenbrock function is  !
! dx1f(x1,x2)= -2(x1-1)-400x1*(x2-x1**2)      !
! dx2f(x1,x2)= 200(x2-x1**2)                  !
!---------------------------------------------!
! The Hessian operator of the Rosenbrok       !
! function is                                 !
! dx1x1f(x1,x2)=-2-400(x2-x1**2)+800x1**2     !
! dx1x2f(x1,x2)=-400x1                        !
! dx2x2f(x1,x2)=200                           !
!---------------------------------------------!



!---------------------------------------------!
!  The routine Rosenbrock returns             !
!  f(x1,x2) in fcost                          !
!  (dx1f(x1,x2),dx2f(x1,x2)) in grad          !
!  for input parameter (x1,x2) in x           !
!---------------------------------------------!
subroutine rosenbrock(x,fcost,grad)
  
  !IN
  real,dimension(2) :: x
  !IN/OUT
  real :: fcost
  real,dimension(2) :: grad
 
  fcost=(1-x(1))**2+100.*(x(2)-x(1)**2)**2
  grad(1)=2.*(x(1)-1)-400.*x(1)*(x(2)-x(1)**2)
  grad(2)=200.*(x(2)-x(1)**2)

end subroutine rosenbrock


!---------------------------------------------!
!  The routine Rosenbrock_Hess returns        !
!  Hessian-vector product H(x)d in output Hd  !
!  for input parameters x and d               !
!  H is the Hessian matrix                    !
!  x=(x1,x2), d=(d1,d2) are two vector of R^2 !
!---------------------------------------------!
subroutine rosenbrock_hess(x,d,Hd)
  
  !IN
  real,dimension(2) :: x,d
  !IN/OUT
  real,dimension(2) :: Hd
  
  Hd(1)=(1200.*x(1)**2-400.*x(2)+2.)*d(1)&
       -400*x(1)*d(2)
  Hd(2)=-400.*x(1)*d(1)+&
       200.*d(2)
  
end subroutine rosenbrock_hess


