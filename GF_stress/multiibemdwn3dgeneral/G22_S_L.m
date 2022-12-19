function G22r=G22_S_L(ks,x,z,rij,mu,dr)
% [G22r,G22i]=G22_S_L(ks,x,z,rij,mu)

% h0P = besselh(0,2,kp*rij);
% h0S = besselh(0,2,ks*rij);
% h2P = besselh(2,2,kp*rij);
% h2S = besselh(2,2,ks*rij);
n   = length(rij);
% h0St= zeros(1,n);
h0Sr= zeros(1,n);
% h0Si= zeros(1,n);

% for k=1:n
%     if x(k)==0
%         [h0St(k),h0Sr(k),h0Si(k)]=besselh_R_Im2(0,2,ks,0,rij(k));
%     elseif abs(z(k))<1e-2*abs(x(k))
%         z1=-1e-2*abs(x(k))*sign(z(k));
%         [h0St(k),h0Sr(k),h0Si(k)]=besselh_R_Im3(0,2,ks,x(k),z1);
%     else
%         [h0St(k),h0Sr(k),h0Si(k)]=besselh_R_Im3(0,2,ks,x(k),z(k));
%     end
% end

for k=1:n
    if x(k)==0
        h0Sr(k)=besselh_R_Im2(0,2,ks,0,rij(k),dr(k));
%     elseif abs(z(k))<1e-2*abs(x(k))
%         z1=-1e-2*abs(x(k))*sign(z(k));
%         h0Sr(k)=besselh_R_Im3(0,2,ks,x(k),z1,dr(k));
    else
        h0Sr(k)=besselh_R_Im3(0,2,ks,x(k),z(k),dr(k));
    end
end


G22r=1/(4*1i*mu).*h0Sr;
% G22i=1/(4*1i*mu).*h0Si;


