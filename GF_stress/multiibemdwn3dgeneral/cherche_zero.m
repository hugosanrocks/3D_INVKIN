function [sol,sollin]=cherche_zero(x,y,ind)
%ne peut s utiliser que quand tous les indices se suivent

y   = y./max(abs(y));%pour eviter pb numerique

ind2    = ind;
%pour determiner quel point ind ou ind+1 s utilisera pour l
%interpolation, on recherche d abord de quel point le 0 lineaire est le
%plus proche
%interp lin
a       = (y(ind2+1)-y(ind2))./(x(ind2+1)-x(ind2));
b       =  y(ind2+1)-a.*x(ind2+1);
sollin  = -b./a;

%verificacion de la cercania del zero lin con los puntos
test    = (sollin-x(ind2))>(x(ind2+1)-sollin);
ind3    = ind2+test;
if ind3(end)==length(y)
    ind3(end)=ind2(end);
end
if ind3(1)==1
    ind3(1)=2;
end
%verificacion problema de paso no constante
indpb   = find(abs(1-abs(x(ind3+1)-x(ind3))./abs(x(ind3)-x(ind3-1)))>1e-6);
ind3(indpb)=ind3(indpb)-test(indpb)+(test(indpb)==0);
if ind3(end)==length(y)
    ind3(end)=ind2(end);
end

%ojo por lo menos entre ind3-1 y ind3 +1 paso constante
fp      = (y(ind3+1)-y(ind3-1))./(x(ind3+1)-x(ind3-1));
fpp     = (y(ind3+1)+y(ind3-1)-2*y(ind3)) ./(x(ind3+1)-x(ind3)).^2;%2*(y(ind3+1)-y(ind3)-fp.*(x(ind3+1)-x(ind3)))./(x(ind3+1)-x(ind3)).^2;
a       = fpp/2;
b       = fp-x(ind3).*fpp;
c       = y(ind3)-fp.*x(ind3)+1/2*fpp.*x(ind3).^2;
b       = b.*sign(a);
c       = c.*sign(a);
a       = a.*sign(a);
delta   = b.^2-4*a.*c;

sol1    = (-b+sqrt(complex(delta)))./(2*a);
sol2    = (-b-sqrt(complex(delta)))./(2*a);

sol3    = sol2.*((sol2)>x(ind2))+sol1.*((sol1)<=x(ind2+1));
indpb   = find(((sol2>x(ind2))+(sol1<x(ind2+1)))==0);
indpb   = [indpb find(((sol2>x(ind2))+(sol1<x(ind2+1)))==2)];
indpb   = [indpb find(a==0)];
indpb   = [indpb find(delta<0)];
indpb   = unique(indpb);
a       = (y(ind2(indpb)+1)-y(ind2(indpb)))./(x(ind2(indpb)+1)-x(ind2(indpb)));
b       = y(ind2(indpb)+1)-a.*x(ind2(indpb)+1);
sol3(indpb)	= -b./a;
if ~isempty(indpb)
    sol3(indpb)	= sol3(indpb).*((sol3(indpb))>x(ind2(indpb))).*((sol3(indpb))<=x(ind2(indpb)+1));
end

sol     = real(sol3);