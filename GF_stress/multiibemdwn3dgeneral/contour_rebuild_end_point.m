%situer si le dernier point est dehors ou dedans d' un autre
%milieu
x2=x(end);
z2=z(end);
m0=inclusiontest(x2,z2,para,1);
if m0<2 || m0==m%m0=m
    %pour corriger certaine erreur d inclusion
    %le dernier point n'appartient qu au milieu m
    m2=m;
    c2=c;
else
    %on change le dernier point pour le point de gauche du
    %contour du milieu m0
    m2=m0;
    x2=cont(m1,1).xa;
    z2=cont(m1,1).za;
    if z1<cont(m,1).za;
        c2=1;
    else
        c2=2;
    end
end