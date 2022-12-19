function T=Tij_PSV(ks,kp,rij,g,C,vn)

%version MP 
d   = eye(2);
T   = zeros(2,2,length(rij));

h2P = besselh(2,2,kp*rij);
h2S = besselh(2,2,ks*rij);

B   = h2P/C(1,1)-h2S/C(6,6);

DP  = kp*rij.*besselh(1,2,kp*rij);
DS  = ks*rij.*besselh(1,2,ks*rij);

C0  = DP/C(1,1)-DS/C(6,6);

F1  = B+C(1,2)/(2*C(6,6)*C(1,1))*DP;
F2  = B+1/(2*C(6,6))*DS;
F3  = C0-4*B;

gknk=g(1,:).*vn(1,:)+g(2,:).*vn(2,:);
for i=1:2
    for j=1:2
        T(i,j,:)=1i*C(6,6)./(2*rij).*(...
            F1.* g(j,:).*vn(i,:)+...
            F2.*(g(i,:).*vn(j,:)+gknk*d(i,j))+...
            F3.* g(i,:).*g(j,:).*gknk);
    end
end

% %version FJSS
% g12=g(1,:).^2;
% g32=g(2,:).^2;
% 
% gg13=2*(3*g12-g32);
% gg31=2*(3*g32-g12);
% 
% 
% fac = 2.0*C(6,6)/C(1,1);
% fac2= 1i/4./rij;
% fac4= fac2*fac/2.*DP;
% e1  = fac4.*g(1,:);
% e3  = fac4.*g(2,:);
% 
% B   = h2P/C(1,1)*C(6,6)-h2S;
% 
% e131= fac2.*g(2,:).*( fac*g12.*DP + (g32-g12).*DS - gg13.*B);
% e133= fac2.*g(1,:).*( fac*g32.*DP + (g12-g32).*DS - gg31.*B);
% e111= fac2.*g(1,:).*( fac*g12.*DP + 2.0* g32 .*DS + gg31.*B);
% e331= fac2.*g(1,:).*( fac*g32.*DP - 2.0* g32 .*DS - gg31.*B);
% e113= fac2.*g(2,:).*( fac*g12.*DP - 2.0* g12 .*DS - gg13.*B);
% e333= fac2.*g(2,:).*( fac*g32.*DP + 2.0* g12 .*DS + gg13.*B);
% 
% fac3= C(1,2)/C(6,6);
% t(1,1,:)= fac3*e1.*vn(1)+e111.*vn(1)+e131.*vn(2);%t11
% t(2,1,:)= fac3*e1.*vn(2)+e131.*vn(1)+e331.*vn(2);%t31
% t(1,2,:)= fac3*e3.*vn(1)+e113.*vn(1)+e133.*vn(2);%t13
% t(2,2,:)= fac3*e3.*vn(2)+e133.*vn(1)+e333.*vn(2);%t33

% fac3= C(1,2)/C(6,6);
% t(1,1,:)= fac3*e1.*vn(1)+e111.*vn(1)+e131.*vn(2);%t11
% t(2,1,:)= fac3*e1.*vn(2)+e131.*vn(1)+e331.*vn(2);%t31
% t(1,2,:)= fac3*e3.*vn(1)+e113.*vn(1)+e133.*vn(2);%t13
% t(2,2,:)= fac3*e3.*vn(2)+e133.*vn(1)+e333.*vn(2);%t33

end