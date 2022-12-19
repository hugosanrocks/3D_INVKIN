function G=Gij_3D(ks,kp,rij,g,C,n)

ba      = kp/ks;%beta/alpha
kpr     = kp*rij;
ksr     = ks*rij;
kprm1   = 1./kpr;
ksrm1   = 1./ksr;

f1=ba^2*(1-2i*kprm1-2*kprm1.^2).*exp(-1i*kpr)+...
    (2i*ksrm1+2*ksrm1.^2).*exp(-1i*ksr);

f2=ba^2*(1i*kprm1+kprm1.^2).*exp(-1i*kpr)+...
    (1-1i*ksrm1-ksrm1.^2).*exp(-1i*ksr);

G   = zeros(3,3,n);
d   = eye(3);

for i=1:3
    for j=i:3
        G(i,j,:)=(f2*d(i,j)+(f1-f2).*g(i,:).*g(j,:))./(4*pi*C(6,6)*rij);
        G(j,i,:)=G(i,j,:);
    end
end