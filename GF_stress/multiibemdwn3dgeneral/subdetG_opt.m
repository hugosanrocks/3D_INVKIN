function det0=subdetG_opt(k2,kS,ln,cP,sP,cS,sS,k1P,k1S,C66,det0,ic)

k22     = k2.^2;
k23     = k22.*k2; 
ln2     = ln.^2;
k1PS    = k1P.*k1S;
k1PSm1  = 1./k1PS;
lnk1    = ln./k1PS;
cPS     = cP.*cS;
sPS     = sP.*sS;
cPsS    = cP.*sS;
sPcS    = sP.*cS;
k22m1   = 1./k22;
uP      = cP.^2 + sP.^2;
uS      = cS.^2 + sS.^2;
lnk2m12 = (ln./k2).^2;
kS2     = kS.^2;
k1Sm1   = 1./k1S;
k1Pm1   = 1./k1P;
k1S2    = k1S.^2;

det0(:,ic,1,1)=k22.*ln.*(...
    -2*(uP + uS)+ ...
    (4*k22./ln+ln.*k22m1).*cPS ...
    -(lnk1+4*k1PS./ln).*sPS);

det0(:,ic,1,2)=k22.*kS2/C66.*(cPsS.*k1Sm1+k1P.*k22m1.*sPcS);

det0(:,ic,1,3)=1i/C66*k23.*(...
    +2*(uP-cPS) ...
    +(uS-cPS).*ln.*k22m1 ...
    +(lnk1+2*k1PS.*k22m1).*sPS);

det0(:,ic,1,4)=1i/C66*k23.*(...
    +2*(uS-cPS) ...
    +(uP-cPS).*ln.*k22m1 ...
    +(lnk1+2*k1PS.*k22m1).*sPS);

det0(:,ic,1,5)=-k22.*kS2/C66.*(sPcS.*k1Pm1+k1S.*k22m1.*sS.*cP);

det0(:,ic,1,6)=-k22/C66^2.*(...
    +uP+uS-2*cPS ...
    +(k22.*k1PSm1+k1PS.*k22m1).*sPS);

det0(:,ic,2,1)=-k22.*kS2*C66.*(4*k1S.*cPsS+lnk2m12.*k1Pm1.*sPcS);

det0(:,ic,2,2)=k22.*( lnk2m12 + 4*k1S2 ).*cPS;

det0(:,ic,2,3)=1i*k2.*kS2.*(ln.*k1Pm1.*sPcS+2*k1S.*cPsS);

det0(:,ic,2,4)=det0(:,ic,2,3);

det0(:,ic,2,5)=k22.*k1S.*k1Pm1.*(4*k1S2 +lnk2m12).*sPS;

det0(:,ic,2,6)=det0(:,ic,1,5);

det0(:,ic,3,1)=1i*C66*k23.*ln.*(...
    +4*(uP) ...
    +2*(uS).*ln.*k22m1 ...
    +(8*k1PS./ln +ln2.*k1PSm1.*k22m1).*sPS ...
    -2*(3*k22-k1S2).*k22m1.*cPS ...
    );

det0(:,ic,3,2)=-1i*k2.*kS2.*(2*k1P.*sPcS+ln.*k1Sm1.*cPsS);

det0(:,ic,3,3)=k22.*(...
    +lnk2m12.*(uS) ...
    +4*k22.*(uP) ...
    -4*ln.*cPS ...
    +(ln2.*k1PSm1+4*k1PS ).*sPS);

det0(:,ic,3,4)=k22.*ln.*(...
    +2*(uP+uS-2*cPS) ...
    +(4*k1PS./ln +lnk1).*sPS ...
    );

det0(:,ic,3,5)=det0(:,ic,2,3);

det0(:,ic,3,6)=det0(:,ic,1,3);

det0(:,ic,4,1)=det0(:,ic,3,1);
det0(:,ic,4,2)=det0(:,ic,3,2);
det0(:,ic,4,3)=det0(:,ic,3,4);
det0(:,ic,4,4)=det0(:,ic,3,3);
det0(:,ic,4,5)=det0(:,ic,2,3);
det0(:,ic,4,6)=det0(:,ic,1,3);


det0(:,ic,5,1)=C66*k22.*kS2.*(4*k1P.*sPcS+lnk2m12.*k1Sm1.*cPsS);

det0(:,ic,5,2)=k22.*k1P.*k1Sm1.*sPS.*(4*k1S2 +lnk2m12);

det0(:,ic,5,3)=det0(:,ic,3,2);

det0(:,ic,5,4)=det0(:,ic,3,2);
det0(:,ic,5,5)=det0(:,ic,2,2);
det0(:,ic,5,6)=det0(:,ic,1,2);

det0(:,ic,6,1)=-C66^2.*k22.*ln2.*(...
    4*(uP+uS-2*cPS) ...
    +(ln2./(k1PS.*k22)+16*(k2./ln).^2.*k1PS).*sPS);

det0(:,ic,6,2)=det0(:,ic,5,1);
det0(:,ic,6,3)=det0(:,ic,3,1);
det0(:,ic,6,4)=det0(:,ic,3,1);
det0(:,ic,6,5)=det0(:,ic,2,1);
det0(:,ic,6,6)=det0(:,ic,1,1);