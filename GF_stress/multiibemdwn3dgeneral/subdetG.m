function det0=subdetG(k2,kS,ln,cP,sP,cS,sS,k1P,k1S,C66,det0,ic)

det0(ic,1,1,:)=k2.^2.*ln.*(...
    -2*(cP.^2 + sP.^2 + cS.^2 + sS.^2)+ ...
    (4*k2.^2./ln+ln./k2.^2).*cP.*cS ...
    -(ln./(k1P.*k1S)+4*k1P.*k1S./ln).*sP.*sS);

det0(ic,1,2,:)=k2.^2.*kS.^2/C66.*(cP.*sS./k1S+k1P./k2.^2.*sP.*cS);

det0(ic,1,3,:)=1i/C66*k2.^3.*(...
    +2*(cP.^2+sP.^2-cP.*cS) ...
    +(cS.^2+sS.^2-cP.*cS).*ln./k2.^2 ...
    +(ln./(k1P.*k1S)+2*k1P.*k1S./k2.^2).*sP.*sS);

det0(ic,1,4,:)=1i/C66*k2.^3.*(...
    +2*(sS.^2+cS.^2-cP.*cS) ...
    +(cP.^2+sP.^2-cP.*cS).*ln./k2.^2 ...
    +(ln./(k1P.*k1S)+2*k1P.*k1S./k2.^2).*sP.*sS);

det0(ic,1,5,:)=-k2.^2.*kS.^2/C66.*(sP.*cS./k1P+k1S./k2.^2.*sS.*cP);

det0(ic,1,6,:)=-k2.^2/C66^2.*(...
    +cP.^2+sP.^2+cS.^2+sS.^2-2*cP.*cS ...
    +(k2.^2./(k1P.*k1S)+k1P.*k1S./k2.^2).*sP.*sS);

det0(ic,2,1,:)=-k2.^2.*kS.^2*C66.*(4*k1S.*cP.*sS+(ln./k2).^2./k1P.*sP.*cS);

det0(ic,2,2,:)=k2.^2.*( (ln./k2).^2 + 4*k1S.^2 ).*cP.*cS;

det0(ic,2,3,:)=1i*k2.*kS.^2.*(ln./k1P.*sP.*cS+2*k1S.*cP.*sS);

det0(ic,2,4,:)=det0(ic,2,3,:);

det0(ic,2,5,:)=k2.^2.*k1S./k1P.*(4*k1S.^2 +(ln./k2).^2).*sP.*sS;

det0(ic,2,6,:)=det0(ic,1,5,:);

det0(ic,3,1,:)=1i*C66*k2.^3.*ln.*(...
    +4*(cP.^2+sP.^2) ...
    +2*(cS.^2+sS.^2).*ln./k2.^2 ...
    +(8*k1P.*k1S./ln +ln.^2./(k1P.*k1S)./k2.^2).*sP.*sS ...
    -2*(3*k2.^2-k1S.^2)./k2.^2.*cP.*cS ...
    );

det0(ic,3,2,:)=-1i*k2.*kS.^2.*(2*k1P.*sP.*cS+ln./k1S.*cP.*sS);

det0(ic,3,3,:)=k2.^2.*(...
    +(ln./k2).^2.*(sS.^2+cS.^2) ...
    +4*k2.^2.*(cP.^2+sP.^2) ...
    -4*ln.*cP.*cS ...
    +(ln.^2./(k1P.*k1S)+4*k1P.*k1S ).*sP.*sS);

det0(ic,3,4,:)=k2.^2.*ln.*(...
    +2*(sP.^2+cP.^2+sS.^2+cS.^2-2*cP.*cS) ...
    +(4*k1P.*k1S./ln +ln./(k1P.*k1S)).*sP.*sS ...
    );

det0(ic,3,5,:)=det0(ic,2,3,:);

det0(ic,3,6,:)=det0(ic,1,3,:);

%segun wathelet no cheque analiticalmente pero con las graficas si salen
det0(ic,4,1,:)=det0(ic,3,1,:);%a ver
det0(ic,4,2,:)=det0(ic,3,2,:);%a ver
det0(ic,4,3,:)=det0(ic,3,4,:);%a ver
det0(ic,4,4,:)=det0(ic,3,3,:);%a ver
det0(ic,4,5,:)=det0(ic,2,3,:);%a ver
det0(ic,4,6,:)=det0(ic,1,3,:);%a ver


det0(ic,5,1,:)=C66*k2.^2.*kS.^2.*(4*k1P.*sP.*cS+(ln./k2).^2./k1S.*cP.*sS);

det0(ic,5,2,:)=k2.^2.*k1P./k1S.*sP.*sS.*(4*k1S.^2 +(ln./k2).^2);

det0(ic,5,3,:)=det0(ic,3,2,:);

det0(ic,5,4,:)=det0(ic,3,2,:);
det0(ic,5,5,:)=det0(ic,2,2,:);%a ver
det0(ic,5,6,:)=det0(ic,1,2,:);%a ver

det0(ic,6,1,:)=-C66^2.*k2.^2.*ln.^2.*(...
    4*(cP.^2+sP.^2+cS.^2+sS.^2-2*cP.*cS) ...
    +(ln.^2./(k1P.*k1S.*k2.^2)+16*(k2./ln).^2.*k1P.*k1S).*sP.*sS);

det0(ic,6,2,:)=det0(ic,5,1,:);%a ver
det0(ic,6,3,:)=det0(ic,3,1,:);%a ver
det0(ic,6,4,:)=det0(ic,3,1,:);%a ver
det0(ic,6,5,:)=det0(ic,2,1,:);%a ver
det0(ic,6,6,:)=det0(ic,1,1,:);%a ver

