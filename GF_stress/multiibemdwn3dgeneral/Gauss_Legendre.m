function gaussL=Gauss_Legendre(n)
%fonction donnant les positions des points xi et les poids wi pour l integration de Gauss-Legendre

%definition
%Px(i,j) : coefficient j du polynôme de Legendre de degré i-1.
%la definitions des coefficients est celle de matlab
%Qx      : dérivée du polynôme de Legendre Px

%ordre et coefficients des polynômes
Px=zeros(n+2,n+2);
Px(1,n+2)=1;%P0=1
Px(2,n+1)=1;%P1=x
for i=2:n+1
    j=i-1;
    for k=(n+1-i+2):n+2
        Px(i+1,k-1)=Px(i+1,k-1)+(2*j+1)/(j+1)*Px(i,k);
        Px(i+1,k  )=-j/(j+1)*Px(i-1,k);
    end
end
Pn=Px(n+1,:);
xi=sort(roots(Pn));
xi=xi.';

Qx=zeros(n+1,n+2);
for i=1:n+1
    for k=2:n+2
        Qx(i,k)=(n+2-(k-1))*Px(i,k-1);
    end
end
Pn1=polyval(Px(n+2,:),xi);
Qn =polyval(Qx(n+1,:),xi);

wi=-2./((n+1)*Qn.*Pn1);

gaussL.ngau  = n;
gaussL.xgau  = xi;
gaussL.wgau  = wi;
