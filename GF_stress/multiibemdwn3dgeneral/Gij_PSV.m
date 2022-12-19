function G=Gij_PSV(ks,kp,rij,g,C,n)

h0P = besselh(0,2,kp*rij);
h0S = besselh(0,2,ks*rij);
h2P = besselh(2,2,kp*rij);
h2S = besselh(2,2,ks*rij);

A   = h0P/C(1,1)+h0S/C(6,6);
B   = h2P/C(1,1)-h2S/C(6,6);

G   = zeros(2,2,n);
d   = eye(2);

i=1;
for j=1:2
    G(i,j,:)=1/(1i*8)*(d(i,j)*A-(2*g(i,:).*g(j,:)-d(i,j)).*B);
end

i=2;j=2;
G(i,j,:)=1/(1i*8)*(d(i,j)*A-(2*g(i,:).*g(j,:)-d(i,j)).*B);
G(2,1,:)=G(1,2,:);

end