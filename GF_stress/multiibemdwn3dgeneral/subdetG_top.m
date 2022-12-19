function det0=subdetG_top(k2,kS,ln,cP,sP,cS,sS,k1P,k1S,C66,det0,ic)

k22     = k2.^2;
k23     = k22.*k2; 
ln2     = ln.^2;
k1PS    = k1P.*k1S;
cPS     = cP.*cS;
sPS     = sP.*sS;
cPsS    = cP.*sS;
sPcS    = sP.*cS;
k22m1   = 1./k22;
uP      = cP.^2 + sP.^2;
uS      = cS.^2 + sS.^2;
lnk2m12 = (ln./k2).^2;
kS2     = kS.^2;
k1S2    = k1S.^2;

det0(:,ic,1,1)=k22.*ln.*(...
    -2*(uP + uS)+ ...
    (4*k22./ln+ln.*k22m1).*cPS ...
    -(ln./k1PS+4*k1PS./ln).*sPS);

det0(:,ic,2,1)=-k22.*kS2*C66.*(4*k1S.*cPsS+lnk2m12./k1P.*sPcS);

det0(:,ic,3,1)=1i*C66*k23.*ln.*(...
    +4*(uP) ...
    +2*(uS).*ln.*k22m1 ...
    +(8*k1PS./ln +ln2./k1PS.*k22m1).*sPS ...
    -2*(3*k22-k1S2).*k22m1.*cPS ...
    );

det0(:,ic,4,1)=det0(:,ic,3,1);

det0(:,ic,5,1)=C66*k22.*kS2.*(4*k1P.*sPcS+lnk2m12./k1S.*cPsS);

det0(:,ic,6,1)=-C66^2.*k22.*ln2.*(...
    4*(uP+uS-2*cPS) ...
    +(ln2./(k1PS.*k22)+16*(k2./ln).^2.*k1PS).*sPS);