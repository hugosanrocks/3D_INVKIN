function [indtri,listc] = suit_cont(indtri,dfcont,listc,coord,drref)

ncont   = size(dfcont,1);
j       = 1:ncont;
for k=length(listc):-1:1
    j(j==abs(listc(k)))=[];
end
if isempty(j)
    return
end
n       = length(j);
x1j = coord.x(dfcont(j,1));
z1j = coord.z(dfcont(j,1));
x2j = coord.x(dfcont(j,2));
z2j = coord.z(dfcont(j,2));

xi	= coord.x(indtri(1));
zi  = coord.z(indtri(1));
xf  = coord.x(indtri(end));
zf  = coord.z(indtri(end));

dr(     1:n )= sqrt((xi-x1j).^2+(zi-z1j).^2);
dr(  n+(1:n))= sqrt((xi-x2j).^2+(zi-z2j).^2);
dr(2*n+(1:n))= sqrt((xf-x1j).^2+(zf-z1j).^2);
dr(3*n+(1:n))= sqrt((xf-x2j).^2+(zf-z2j).^2);
[~,j0] = min(dr);
if dr(j0)>5*drref
    %es un contorno aparte cual es possible solo en caso de medio 1
    listc=[listc,0,length(indtri),0];
end    
if j0>3*n
    j1=j(j0-3*n);
    indtri=[indtri,(dfcont(j1,2):-1:dfcont(j1,1))];
    listc=[listc,-j1];
elseif j0>2*n
    j1=j(j0-2*n);
    indtri=[indtri,(dfcont(j1,1):dfcont(j1,2))];
    listc=[listc,j1];
elseif j0>n
    j1=j(j0-n);
    indtri=[(dfcont(j1,1):dfcont(j1,2)),indtri];
    listc=[j1,listc];
else
    j1=j(j0);
    indtri=[(dfcont(j1,2):-1:dfcont(j1,1)),indtri];
    listc=[-j1,listc];
end

[indtri,listc] = suit_cont(indtri,dfcont,listc,coord,drref);