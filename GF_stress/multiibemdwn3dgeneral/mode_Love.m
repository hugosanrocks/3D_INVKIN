function det0=mode_Love(para,DWN)

ncapas      = para.nsubmed;
sub         = para.sub;
k2          = DWN.k2;
nk2       	= length(k2(1,:));
if nk2==1
    nk2=length(DWN.omegac);
end
p           = complex(1,1);
k1          = zeros(nk2,ncapas,'like',p);


for ic=1:ncapas
    ksi 	= DWN.omegac/sub(ic).bet;
    k1(:,ic)= sqrt(k2.^2-ksi.^2); %k1S%%%%%%% attention c est a l envers (*1i)
    k1(k1(:,ic)==0,ic)=1e-10;
end
% k1 =k1.*(sign(imag(k1)).*(imag(k1)~=0)+(imag(k1)==0));

%%%%%%%%%%%%%%%%
% calcul des G %
%%%%%%%%%%%%%%%%
% G   = zeros(2,2,nk2);

ic=1;

% sub     = sub(ic);
S12     = k1(:,ic)*sub(ic).C(6,6); %sigma_zy sigma12
kh      = k1(:,ic)*sub(ic).h;

% sh = sinh(kh);
% ch = cosh(kh);

sh = zeros(nk2,1,'like',p);
ch = zeros(nk2,1,'like',p);

indi        = (imag(kh)==0);
sh(~indi,1) = sinh(kh(~indi));
ch(~indi,1) = cosh(kh(~indi));
sh(indi,1)  = 0.5*(1-exp(-2*kh(indi)));
ch(indi,1)  = 0.5*(1+exp(-2*kh(indi)));
% sh      = 0.5*(1-exp(-2*kh));
% ch      = 0.5*(1+exp(-2*kh));
% G(1,1,:)= cosh(k1(:,ic)*MAT.h);
% G(1,2,:)= sh./S12;
% G(2,1,:)= sh.*S12;
% G(2,2,:)= G(1,1,:);
% G1=G;

P11 = ch;
P21 = sh.*S12;
if ncapas>2
    P12 = sh./S12;
    P22 = P11;
else
    %for compilation
    P12 = 0;
    P22 = 0;
end

for ic=2:ncapas-1
%     MAT     = sub(ic);
    S12     = k1(:,ic)*sub(ic).C(6,6); %sigma_zy sigma12
    kh      = k1(:,ic)*sub(ic).h;
    
    indi        = (imag(kh)==0);
    sh(~indi,1) = sinh(kh(~indi));
    ch(~indi,1) = cosh(kh(~indi));
    sh(indi)    = 0.5*(1-exp(-2*kh(indi)));
    ch(indi)    = 0.5*(1+exp(-2*kh(indi)));
    
    

%         G(1,1,:)= cosh(k1(:,ic)*MAT.h);
%         G(1,2,:)= sh./S12;
%         G(2,1,:)= sh.*S12;
%         G(2,2,:)= G(1,1,:);
    G11 = ch;
    G12 = sh./S12;
    G21 = sh.*S12;
    G22 = G11;
    
    T11 = G11.*P11+G12.*P21;
    T21 = G21.*P11+G22.*P21;
    P11 = T11;
    P21 = T21;
    if ic~=(ncapas-1)
        T12 = G11.*P12+G12.*P22;
        T22 = G21.*P12+G22.*P22;
        P12 = T12;
        P22 = T22;
    end
    
%         for i=1:nk2
%             G1(:,:,i) = G(:,:,i)*G1(:,:,i);
%         end
    %a developper analytiquement
end
% Tm1         = zeros(2,2,nk2);
ic          = ncapas;
% MAT         = sub(ic);
% Tm1(1,1,:)  = 1;
% Tm1(1,2,:)  =-1./k1(:,ic)/MAT.Ci(6,6);
% Tm1(2,1,:)  = 1;
% Tm1(2,2,:)  =-Tm1(1,2,:);

% det0 = squeeze(Tm1(2,1,:).*G1(1,1,:)+Tm1(2,2,:).*G1(2,1,:));


% det0 = squeeze(G1(1,1,:))+squeeze(1./k1(:,ic)/MAT.Ci(6,6)).*squeeze(G1(2,1,:));
det0 = P11+1./k1(:,ic)/sub(ic).C(6,6).*P21;
