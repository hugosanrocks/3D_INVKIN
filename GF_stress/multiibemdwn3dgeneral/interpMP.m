function yi=interpMP(y,nrate)

n1  =length(y);
x0=linspace(0,1,n1);
x1=linspace(0,1,nrate*n1);
yi = interp1(x0,y,x1,'spline');