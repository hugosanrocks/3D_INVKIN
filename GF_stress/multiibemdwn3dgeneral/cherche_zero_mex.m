function sol=cherche_zero_mex(x,y,ind)

% sol=zeros(
y   = y./max(abs(y));%pour eviter pb numerique
if ind==1
    %interp lin
    a   = (y(2)-y(1))/(x(2)-x(1));
    b   = y(2)-a*x(2);
    sol= -b/a;
else
    fp  = (y(ind+1)-y(ind-1))./(x(ind+1)-x(ind-1));
    fpp = 2*(y(ind+1)-y(ind)-fp.*(x(ind+1)-x(ind)))./(x(ind+1)-x(ind)).^2;
    a   = fpp/2;
    b   = fp-x(ind)*fpp;
    c   = y(ind)-fp*x(ind)+1/2*fpp*x(ind)^2;
    b   = b*sign(a);
    c   = c*sign(a);
    a   = a*sign(a);
    delta=b^2-4*a*c;
    
    sol1=(-b+sqrt(delta))/(2*a);
    sol2=(-b-sqrt(delta))/(2*a);
    
    sol    = sol2.*((sol2)>x(ind))+sol1.*((sol1)<=x(ind+1));
    if (sol2<x(ind) || sol1>x(ind+1)) || (sol2>x(ind) && sol1<x(ind+1)) || a==0 || delta<0
        a       = (y(ind(indpb)+1)-y(ind(indpb)))./(x(ind(indpb)+1)-x(ind(indpb)));
        b       = y(ind(indpb)+1)-a.*x(ind(indpb)+1);
        sol     =-b./a;
    end
end
sol     = real(sol);