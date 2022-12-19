function Hr=besselh_R_Im3(ord,type,k,x,z,dr)
% [Ht,Hr,Hi]=besselh_R_Im3(ord,type,k,x,z)

if sign(imag(k))==-1 && type==1
    k=-real(k)-1i*imag(k);
elseif sign(imag(k))==1 && type==2
    k=real(k)-1i*imag(k);
end

%hay que rotar el sistema de referencia en el sistema local del elemento
%Hr corresponde a las ondas homogeneas
%Hi corresponde a las ondas heterogeneas por la cual hay que buscar el
%limite de la integracion
kr      = abs(real(k));
Hr      = integral(@(kx) DComp_hankel_OP(kx,ord,type,k,x,z,dr),-kr ,kr,'AbsTol',1e-3);

% % limite de l integration a chercher
% % ref
% kx      = 1.01*kr;
% cOPref  = DComp_hankel_OP(kx,ord,type,k,x,z);
% cOP     = cOPref;
% j=1;
% while abs(cOP)>abs(1e-5*cOPref) && j<10
%     kx  = kx*10;
%     cOP = DComp_hankel_OP(kx,ord,type,k,x,z);
%     j=j+1;
% end
% kxf     = kx;
% 
% Hi      = integral(@(kx) DComp_hankel_OP(kx,ord,type,k,x,z),kr,kxf,'RelTol',1e-5) + ...
%     integral(@(kx) DComp_hankel_OP(kx,ord,type,k,x,z),-kxf,-kr,'RelTol',1e-5);
% 
% Ht      = Hr+Hi;