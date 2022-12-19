function [x,dr]=rec_malla_dx(drmin,xi,xf,n,cont)

x1   	= linspace(xi,xf,n);
x1tmp   = x1-(cont.xa+cont.a);
z1      = eq_contour(x1tmp,cont)+cont.za;

r       = zeros(n,1);
dr      = zeros(n,1);
r(2:n)	= cumsum(sqrt(diff(x1).^2+diff(z1).^2));
r(1)  	= 0;
i      	= 2:(n-1);
dr(n)	= (r(n)-r(n-1))/2;
dr(i)  	= (r(i+1)-r(i-1))/2;
dr(1) 	= (r(2)-r(1))/2;

x       = x1;
x0      = x;
dx      = abs(xf-xi)/n;
% if (max(dr)/dx>1e3)
%     toto=2;
% end
while (max(dr)>drmin) && dx>1e-3
    ind   	= find(dr>drmin);
    %separacion en partes contiguas
    ind2    = find(diff(ind)>1);
    n2      = length(ind2);
    indi            = zeros(n2+1,1);
    indf            = zeros(n2+1,1);
    indi(2:(n2+1))  = ind(ind2+1)-1;
    indi(1)         = ind(1);
    indf(n2+1)      = min(ind(end)+1,n);
    indf(1:n2)      = ind(ind2)+1;
    
    for i=1:n2+1
        xi  	= x0(indi(i));
        xf      = x0(indf(i));
        n     	= 4*(indf(i)-indi(i)+1);
        [x1,dr] = rec_malla_dx(drmin,xi,xf,n,cont);
        x       = [x,x1];
    end
end