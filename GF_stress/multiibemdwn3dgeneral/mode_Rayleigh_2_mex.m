function det0=mode_Rayleigh_2_mex(DWN,ncapas,h,vS,vP,C66)

%se calcula el determinante con el metodo de Dunkin cf. Wathelet
%con subdeterminate, la clave es que los subdeterminante estan calculados
%explicitamente, lo que reduce errores numericos
p           = complex(1,1);

k2          = DWN.k2;
nk2       	= length(k2(1,:));
if nk2==1
    nk2=length(DWN.omegac);
end

pos=[...
    1 2
    1 3
    1 4
    2 3
    2 4
    3 4];
det0        = zeros(nk2,ncapas,6,6,'like',p);

k1S          = zeros(nk2,ncapas,'like',p);
k1P          = zeros(nk2,ncapas,'like',p);
for ic=1:ncapas
    ksi         = DWN.omegac/vS(ic);
    kpi         = DWN.omegac/vP(ic);
    k1S(:,ic)   = sqrt(complex(ksi.^2-k2.^2)); %k1S
    k1P(:,ic)   = sqrt(complex(kpi.^2-k2.^2)); %k1P
    k1S(k1S(:,ic)==0,ic)=1e-10;
    k1P(k1P(:,ic)==0,ic)=1e-10;
    k1S(:,ic) = k1S(:,ic).*(-sign(imag(k1S(:,ic))).*(imag(k1S(:,ic))~=0)+(imag(k1S(:,ic))==0));
    k1P(:,ic) = k1P(:,ic).*(-sign(imag(k1P(:,ic))).*(imag(k1P(:,ic))~=0)+(imag(k1P(:,ic))==0));
end

%%%%%%%%%%%%%%%%
% calcul des G %
%%%%%%%%%%%%%%%%
ic=1;

k1Ph        = k1P(:,ic)*h(ic);
k1Sh        = k1S(:,ic)*h(ic);

sP = zeros(nk2,1,'like',p);
cP = zeros(nk2,1,'like',p);
sS = zeros(nk2,1,'like',p);
cS = zeros(nk2,1,'like',p);


indi        = (imag(k1P(:,ic))==0); % a fus avec au dessus k1 sign
sP(indi,1)  = sin(k1Ph(indi));
cP(indi,1)  = cos(k1Ph(indi));
sP(~indi,1) = (1-exp(-2i*k1Ph(~indi)))/2i;
cP(~indi,1) = (1+exp(-2i*k1Ph(~indi)))/2;

sS(indi,1)  = sin(k1Sh(indi));
cS(indi,1)  = cos(k1Sh(indi));
sS(~indi,1) = (exp(1i*k1Sh(~indi)-1i*k1Ph(~indi))-exp(-1i*k1Sh(~indi)-1i*k1Ph(~indi)))/2i;
cS(~indi,1) = (exp(1i*k1Sh(~indi)-1i*k1Ph(~indi))+exp(-1i*k1Sh(~indi)-1i*k1Ph(~indi)))/2 ;

% sP0(:,1)	= sin(k1Ph(:,ic));
% cP0(:,1)	= cos(k1Ph(:,ic));
% sS0(:,1)	= sin(k1Sh(:,ic));
% cS0(:,1)	= cos(k1Sh(:,ic));

k2          = k2.';
ln          = k2.^2-k1S(:,ic).^2;
C66ic         = C66(ic);

% G           = zeros(4,2,nk2,'like',p);
% G(1,1,:)    = 2*k2.^2.*cP-ln.*cS;
% G(1,2,:)    = 1i*k2./k1P(:,ic).*(ln.*sP+2*k1P(:,ic).*k1S(:,ic).*sS);
% G(2,1,:)    =-1i*k2./k1S(:,ic).*(2*k1P(:,ic).*k1S(:,ic).*sP+ln.*sS);
% G(2,2,:)    =-ln.*cP+2*k2.^2.*cS;
% G(3,1,:)    = 2i*C66ic*k2.*ln.*(cP-cS);
% G(3,2,:)    =-C66ic./k1P(:,ic).*(ln.^2.*sP+4.*k2.^2.*k1P(:,ic).*k1S(:,ic).*sS);
% G(4,1,:)    =-C66ic./k1S(:,ic).*(4*k2.^2.*k1P(:,ic).*k1S(:,ic).*sP+ln.^2.*sS);
% G(4,2,:)    = G(3,1,:);

% for jic=1:6
%     a=pos(jic,1);
%     b=pos(jic,2);
%     det0(ic,jic,1,:)=subdet(a,b,1,2,G);
% end

kS	= DWN.omegac.'/vS(ic);
det0= subdetG_top(k2,kS,ln,cP,sP,cS,sS,k1P(:,ic),k1S(:,ic),C66ic,det0,ic);


G               = zeros(4,4,nk2,'like',p);
for ic=2:ncapas-1
    k1Ph        = k1P(:,ic)*h(ic);
    k1Sh        = k1S(:,ic)*h(ic);
    
    indi        = (imag(k1P(:,ic))==0); % a fus avec au dessus k1 sign
    sP(indi,1)  = sin(k1Ph(indi));
    cP(indi,1)  = cos(k1Ph(indi));
    sP(~indi,1) = (1-exp(-2i*k1Ph(~indi)))/2i;
    cP(~indi,1) = (1+exp(-2i*k1Ph(~indi)))/2;
    
    sS(indi,1)  = sin(k1Sh(indi));
    cS(indi,1)  = cos(k1Sh(indi));
    sS(~indi,1) = (exp(1i*k1Sh(~indi)-1i*k1Ph(~indi))-exp(-1i*k1Sh(~indi)-1i*k1Ph(~indi)))/2i;
    cS(~indi,1) = (exp(1i*k1Sh(~indi)-1i*k1Ph(~indi))+exp(-1i*k1Sh(~indi)-1i*k1Ph(~indi)))/2 ;
    
    ln          = (k2.^2-k1S(:,ic).^2);
    C66ic         = C66(ic);
    
    
%     G(1,1,:)    = 2*k2.^2.*cP-ln.*cS;
%     G(1,2,:)    = 1i*k2./k1P(:,ic).*(ln.*sP+2*k1P(:,ic).*k1S(:,ic).*sS);
%     G(1,3,:)    = 1i*k2/C66ic.*(cP-cS);
%     G(1,4,:)    = 1./(C66ic*k1P(:,ic)).*(k2.^2.*sP+k1P(:,ic).*k1S(:,ic).*sS);
%     
%     G(2,1,:)    =-1i*k2./k1S(:,ic).*(2*k1P(:,ic).*k1S(:,ic).*sP+ln.*sS);
%     G(2,2,:)    =-ln.*cP+2*k2.^2.*cS;
%     G(2,3,:)    = 1./(C66ic*k1S(:,ic)).*(k1P(:,ic).*k1S(:,ic).*sP+k2.^2.*sS);
%     G(2,4,:)    = G(1,3,:);
%     
%     G(3,1,:)    = 2i*C66ic*k2.*ln.*(cP-cS);
%     G(3,2,:)    =-C66ic./k1P(:,ic).*(ln.^2.*sP+4.*k2.^2.*k1P(:,ic).*k1S(:,ic).*sS);
%     G(3,3,:)    = G(2,2,:);
%     G(3,4,:)    = G(1,2,:);
%     
%     G(4,1,:)    =-C66ic./k1S(:,ic).*(4*k2.^2.*k1P(:,ic).*k1S(:,ic).*sP+ln.^2.*sS);
%     G(4,2,:)    = G(3,1,:);
%     G(4,3,:)    = G(2,1,:);
%     G(4,4,:)    = G(1,1,:);
    
    
    %     for jic=1:6
    %         a=pos(jic,1);
    %         b=pos(jic,2);
    %         for jic2=1:6
    %             c=pos(jic2,1);
    %             d=pos(jic2,2);
    %             det0(ic,jic,jic2,:)=squeeze(subdet(a,b,c,d,G));
    %         end
    %     end
    kS	= DWN.omegac.'/vS(ic);
%     det0= subdetG(k2,kS,ln,cP,sP,cS,sS,k1P(:,ic),k1S(:,ic),C66ic,det0,ic);
    det0= subdetG_opt(k2,kS,ln,cP,sP,cS,sS,k1P(:,ic),k1S(:,ic),C66ic,det0,ic);
end
% clear G

%semi-espacio
ic          = ncapas;

Tm1         = zeros(2,4,nk2,'like',p);

ln          = k2.^2-k1S(:,ic).^2;
C66ic         = C66(ic);

Tm1(1,1,:)  = 2i*C66ic*k2.*k1P(:,ic).*k1S(:,ic);
Tm1(1,2,:)  = 1i*C66ic*ln.*k1S(:,ic);
Tm1(1,3,:)  =-k1S(:,ic).*k1P(:,ic);
Tm1(1,4,:)  = k2.*k1S(:,ic);

Tm1(2,1,:)  =-1i*C66ic*ln.*k1P(:,ic);
Tm1(2,2,:)  = Tm1(1,1,:);
Tm1(2,3,:)  = k2.*k1P(:,ic);
Tm1(2,4,:)  =-Tm1(1,3,:);

coeffTm       = 1./(k1P(:,ic).*k1S(:,ic));

for jic=1:6
    a=pos(jic,1);
    b=pos(jic,2);
    det0(:,ic,1,jic)=subdet(1,2,a,b,Tm1);
end
% clear Tm1

% det     = ones(6^(ncapas-1),nk2);
% j       = zeros(ncapas-1,1);
% det     = subdetab(det0,j,det,1,ncapas);
% % indpb   = sign(sum(det==0,1));
% 
% det     = sum(det,1);


det 	= subdetab2(det0,ncapas,nk2);
det  	= sum(det,1);

% if nk2>100 && sum(indpb)>0
%     indpb=(1:nk2).*indpb;
%     indpb(indpb==0)=[];
%     det(indpb)=1/5*(det(indpb-2)+det(indpb-1)+det(indpb)+det(indpb+1)+det(indpb+2));
% end
det0  = det.'.*coeffTm.^2;