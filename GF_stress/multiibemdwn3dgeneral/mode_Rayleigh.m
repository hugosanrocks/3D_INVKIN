function det0=mode_Rayleigh(para,DWN)

%se calcula el determinante con el metodo de Dunkin cf. Wathelet
%con subdeterminate
%problema los subdeterminantes estan calculados numericamente y conllevan a
%errores numericas, la mejora esta en "mode_Rayleigh_2"

ncapas      = para.nsubmed;
k2          = DWN.k2;
nk2       	= length(k2(1,:));
sub         = para.sub;
if nk2==1
    nk2=length(DWN.omegac);
end
k1S          = zeros(nk2,ncapas);
k1P          = zeros(nk2,ncapas);

for ic=1:ncapas
    ksi         = DWN.omegac/sub(ic).bet;
    kpi         = DWN.omegac/sub(ic).alpha;
    k1S(:,ic)   = sqrt(ksi.^2-k2.^2); %k1S
    k1P(:,ic)   = sqrt(kpi.^2-k2.^2); %k1P
    k1S(k1S(:,ic)==0,ic)=1e-10;
    k1P(k1P(:,ic)==0,ic)=1e-10;
    k1S(:,ic) = k1S(:,ic).*(-sign(imag(k1S(:,ic))).*(imag(k1S(:,ic))~=0)+(imag(k1S(:,ic))==0));
    k1P(:,ic) = k1P(:,ic).*(-sign(imag(k1P(:,ic))).*(imag(k1P(:,ic))~=0)+(imag(k1P(:,ic))==0));
end

%%%%%%%%%%%%%%%%
% calcul des G %
%%%%%%%%%%%%%%%%
% G   = zeros(2,2,nk2);

ic=1;

k1Ph        = k1P(:,ic)*sub(ic).h;
k1Sh        = k1S(:,ic)*sub(ic).h;

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
C66         = sub(ic).C(6,6);

G           = zeros(4,2,nk2);

G(1,1,:)    = 2*k2.^2.*cP-ln.*cS;
G(1,2,:)    = 1i*k2./k1P(:,ic).*(ln.*sP+2*k1P(:,ic).*k1S(:,ic).*sS);
% G(1,3,:)    = 1i*k2/C66.*(cP-cS);
% G(1,4,:)    = 1./(C66*k1P(:,ic)).*(k2.^2.*sP+k1P(:,ic).*k1S(:,ic).*sS);

G(2,1,:)    =-1i*k2./k1S(:,ic).*(2*k1P(:,ic).*k1S(:,ic).*sP+ln.*sS);
G(2,2,:)    =-ln.*cP+2*k2.^2.*cS;
% G(2,3,:)    = 1./(C66*k1S(:,ic)).*(k1P(:,ic).*k1S(:,ic).*sP+k2.^2.*sS);
% G(2,4,:)    = G(1,3,:);

G(3,1,:)    = 2i*C66*k2.*ln.*(cP-cS);
G(3,2,:)    =-C66./k1P(:,ic).*(ln.^2.*sP+4.*k2.^2.*k1P(:,ic).*k1S(:,ic).*sS);
% G(3,3,:)    = G(2,2,:);
% G(3,4,:)    = G(1,2,:);

G(4,1,:)    =-C66./k1S(:,ic).*(4*k2.^2.*k1P(:,ic).*k1S(:,ic).*sP+ln.^2.*sS);
G(4,2,:)    = G(3,1,:);
% G(4,3,:)    = G(2,1,:);
% G(4,4,:)    = G(1,1,:);
% G           = G/max(max(max(abs(G))));

G1          = G;
for ic=2:ncapas-1
    k1Ph        = k1P(:,ic)*sub(ic).h;
    k1Sh        = k1S(:,ic)*sub(ic).h;
    
    indi        = (imag(k1P(:,ic))==0); % a fus avec au dessus k1 sign
    sP(indi,1)  = sin(k1Ph(indi));
    cP(indi,1)  = cos(k1Ph(indi));
    sP(~indi,1) = (1-exp(-2i*k1Ph(~indi)))/2i;
    cP(~indi,1) = (1+exp(-2i*k1Ph(~indi)))/2;
    
    sS(indi,1)  = sin(k1Sh(indi));
    cS(indi,1)  = cos(k1Sh(indi));
    sS(~indi,1) = (exp(1i*k1Sh(~indi)-1i*k1Ph(~indi))-exp(-1i*k1Sh(~indi)-1i*k1Ph(~indi)))/2i;
    cS(~indi,1) = (exp(1i*k1Sh(~indi)-1i*k1Ph(~indi))+exp(-1i*k1Sh(~indi)-1i*k1Ph(~indi)))/2 ;
    
    % sP0(:,1)	= sin(k1Ph);
    % cP0(:,1)	= cos(k1Ph);
    % sS0(:,1)	= sin(k1Sh);
    % cS0(:,1)	= cos(k1Sh);
    
    ln          = k2.^2-k1S(:,ic).^2;
    C66         = sub(ic).C(6,6);
    
    G           = zeros(4,4,nk2);
    
    G(1,1,:)    = 2*k2.^2.*cP-ln.*cS;
    G(1,2,:)    = 1i*k2./k1P(:,ic).*(ln.*sP+2*k1P(:,ic).*k1S(:,ic).*sS);
    G(1,3,:)    = 1i*k2/C66.*(cP-cS);
    G(1,4,:)    = 1./(C66*k1P(:,ic)).*(k2.^2.*sP+k1P(:,ic).*k1S(:,ic).*sS);
    
    G(2,1,:)    =-1i*k2./k1S(:,ic).*(2*k1P(:,ic).*k1S(:,ic).*sP+ln.*sS);
    G(2,2,:)    =-ln.*cP+2*k2.^2.*cS;
    G(2,3,:)    = 1./(C66*k1S(:,ic)).*(k1P(:,ic).*k1S(:,ic).*sP+k2.^2.*sS);
    G(2,4,:)    = G(1,3,:);
    
    G(3,1,:)    = 2i*C66*k2.*ln.*(cP-cS);
    G(3,2,:)    =-C66./k1P(:,ic).*(ln.^2.*sP+4.*k2.^2.*k1P(:,ic).*k1S(:,ic).*sS);
    G(3,3,:)    = G(2,2,:);
    G(3,4,:)    = G(1,2,:);
    
    G(4,1,:)    =-C66./k1S(:,ic).*(4*k2.^2.*k1P(:,ic).*k1S(:,ic).*sP+ln.^2.*sS);
    G(4,2,:)    = G(3,1,:);
    G(4,3,:)    = G(2,1,:);
    G(4,4,:)    = G(1,1,:);
%     G           = G/max(max(max(abs(G))));
    for i=1:nk2
        G1(:,:,i) = G(:,:,i)*G1(:,:,i);
    end
    %a developper analytiquement
end

%semi-espacio
ic          = ncapas;

% Tm1         = zeros(4,4,nk2);
Tm1         = zeros(2,4,nk2);

ln          = k2.^2-k1S(:,ic).^2;
C66         = sub(ic).C(6,6);

Tm1(1,1,:)  = 2i*C66*k2.*k1P(:,ic).*k1S(:,ic);
Tm1(1,2,:)  = 1i*C66*ln.*k1S(:,ic);
Tm1(1,3,:)  =-k1S(:,ic).*k1P(:,ic);
Tm1(1,4,:)  = k2.*k1S(:,ic);

Tm1(2,1,:)  =-1i*C66*ln.*k1P(:,ic);
Tm1(2,2,:)  = Tm1(1,1,:);
Tm1(2,3,:)  = k2.*k1P(:,ic);
Tm1(2,4,:)  =-Tm1(1,3,:);

% Tm1(3,1,:)  = Tm1(1,1,:);
% Tm1(3,2,:)  =-Tm1(1,2,:);
% Tm1(3,3,:)  = Tm1(1,3,:);
% Tm1(3,4,:)  =-Tm1(1,4,:);
% 
% Tm1(4,1,:)  =-Tm1(2,1,:);
% Tm1(4,2,:)  = Tm1(2,2,:);
% Tm1(4,3,:)  =-Tm1(2,3,:);
% Tm1(4,4,:)  = Tm1(2,4,:);

coeffTm       = 1./(k1P(:,ic).*k1S(:,ic));

% G1f(2,2,:) = squeeze(Tm1(2,1,:).*G1(1,2,:)+Tm1(2,2,:).*G1(2,2,:)+Tm1(2,3,:).*G1(3,2,:)+Tm1(2,4,:).*G1(4,2,:));
% G1f(1,1,:) = squeeze(Tm1(1,1,:).*G1(1,1,:)+Tm1(1,2,:).*G1(2,1,:)+Tm1(1,3,:).*G1(3,1,:)+Tm1(1,4,:).*G1(4,1,:));
% G1f(1,2,:) = squeeze(Tm1(1,1,:).*G1(1,2,:)+Tm1(1,2,:).*G1(2,2,:)+Tm1(1,3,:).*G1(3,2,:)+Tm1(1,4,:).*G1(4,2,:));
% G1f(2,1,:) = squeeze(Tm1(2,1,:).*G1(1,1,:)+Tm1(2,2,:).*G1(2,1,:)+Tm1(2,3,:).*G1(3,1,:)+Tm1(2,4,:).*G1(4,1,:));
% 
% % for i=1:nk2
% %     G1(:,:,i) = Tm1(:,:,i)*G1(:,:,i)*coeffTm(i);
% % end
% 
% det0        = squeeze(G1f(1,1,:).*G1f(2,2,:)-G1f(1,2,:).*G1f(2,1,:)).*coeffTm.^2;

% det0        = coeffTm.^2.*squeeze(...
%     +(Tm1(1,1,:).*G1(1,1,:)+Tm1(1,2,:).*G1(2,1,:)+Tm1(1,3,:).*G1(3,1,:)+Tm1(1,4,:).*G1(4,1,:)).*...
%     +(Tm1(2,1,:).*G1(1,2,:)+Tm1(2,2,:).*G1(2,2,:)+Tm1(2,3,:).*G1(3,2,:)+Tm1(2,4,:).*G1(4,2,:)) ...
%     -(Tm1(1,1,:).*G1(1,2,:)+Tm1(1,2,:).*G1(2,2,:)+Tm1(1,3,:).*G1(3,2,:)+Tm1(1,4,:).*G1(4,2,:)).*...
%     +(Tm1(2,1,:).*G1(1,1,:)+Tm1(2,2,:).*G1(2,1,:)+Tm1(2,3,:).*G1(3,1,:)+Tm1(2,4,:).*G1(4,1,:)));

% det0        = coeffTm.^2.*squeeze(...
%     +Tm1(1,1,:).*Tm1(2,1,:).*G1(1,1,:).*G1(1,2,:) ...
%     -Tm1(1,1,:).*Tm1(2,1,:).*G1(1,2,:).*G1(1,1,:) ...
%     +Tm1(1,1,:).*Tm1(2,2,:).*G1(1,1,:).*G1(2,2,:) ...
%     -Tm1(1,1,:).*Tm1(2,2,:).*G1(1,2,:).*G1(2,1,:) ...
%     +Tm1(1,1,:).*Tm1(2,3,:).*G1(1,1,:).*G1(3,2,:) ...
%     -Tm1(1,1,:).*Tm1(2,3,:).*G1(1,2,:).*G1(3,1,:) ...
%     +Tm1(1,1,:).*Tm1(2,4,:).*G1(1,1,:).*G1(4,2,:) ...
%     -Tm1(1,1,:).*Tm1(2,4,:).*G1(1,2,:).*G1(4,1,:) ...
%     +Tm1(1,2,:).*Tm1(2,1,:).*G1(2,1,:).*G1(1,2,:) ...
%     -Tm1(1,2,:).*Tm1(2,1,:).*G1(2,2,:).*G1(1,1,:) ...
%     +Tm1(1,2,:).*Tm1(2,2,:).*G1(2,1,:).*G1(2,2,:) ...
%     -Tm1(1,2,:).*Tm1(2,2,:).*G1(2,2,:).*G1(2,1,:) ...
%     +Tm1(1,2,:).*Tm1(2,3,:).*G1(2,1,:).*G1(3,2,:) ...
%     -Tm1(1,2,:).*Tm1(2,3,:).*G1(2,2,:).*G1(3,1,:) ...
%     +Tm1(1,2,:).*Tm1(2,4,:).*G1(2,1,:).*G1(4,2,:) ...
%     -Tm1(1,2,:).*Tm1(2,4,:).*G1(2,2,:).*G1(4,1,:) ...
%     +Tm1(1,3,:).*Tm1(2,1,:).*G1(3,1,:).*G1(1,2,:) ...
%     -Tm1(1,3,:).*Tm1(2,1,:).*G1(3,2,:).*G1(1,1,:) ...
%     +Tm1(1,3,:).*Tm1(2,2,:).*G1(3,1,:).*G1(2,2,:) ...
%     -Tm1(1,3,:).*Tm1(2,2,:).*G1(3,2,:).*G1(2,1,:) ...
%     +Tm1(1,3,:).*Tm1(2,3,:).*G1(3,1,:).*G1(3,2,:) ...
%     -Tm1(1,3,:).*Tm1(2,3,:).*G1(3,2,:).*G1(3,1,:) ...
%     +Tm1(1,3,:).*Tm1(2,4,:).*G1(3,1,:).*G1(4,2,:) ...
%     -Tm1(1,3,:).*Tm1(2,4,:).*G1(3,2,:).*G1(4,1,:) ...
%     +Tm1(1,4,:).*Tm1(2,1,:).*G1(4,1,:).*G1(1,2,:) ...
%     -Tm1(1,4,:).*Tm1(2,1,:).*G1(4,2,:).*G1(1,1,:) ...
%     +Tm1(1,4,:).*Tm1(2,2,:).*G1(4,1,:).*G1(2,2,:) ...
%     -Tm1(1,4,:).*Tm1(2,2,:).*G1(4,2,:).*G1(2,1,:) ...
%     +Tm1(1,4,:).*Tm1(2,3,:).*G1(4,1,:).*G1(3,2,:) ...
%     -Tm1(1,4,:).*Tm1(2,3,:).*G1(4,2,:).*G1(3,1,:) ...
%     +Tm1(1,4,:).*Tm1(2,4,:).*G1(4,1,:).*G1(4,2,:) ...
%     -Tm1(1,4,:).*Tm1(2,4,:).*G1(4,2,:).*G1(4,1,:) ...
% );

det0  = coeffTm.^2.*squeeze(...
    +(Tm1(1,1,:).*Tm1(2,2,:)-Tm1(1,2,:).*Tm1(2,1,:)).*(G1(1,1,:).*G1(2,2,:)-G1(1,2,:).*G1(2,1,:)) ...
    +(Tm1(1,1,:).*Tm1(2,3,:)-Tm1(1,3,:).*Tm1(2,1,:)).*(G1(1,1,:).*G1(3,2,:)-G1(1,2,:).*G1(3,1,:)) ...
    +(Tm1(1,1,:).*Tm1(2,4,:)-Tm1(1,4,:).*Tm1(2,1,:)).*(G1(1,1,:).*G1(4,2,:)-G1(1,2,:).*G1(4,1,:)) ...
    +(Tm1(1,2,:).*Tm1(2,3,:)-Tm1(1,3,:).*Tm1(2,2,:)).*(G1(2,1,:).*G1(3,2,:)-G1(2,2,:).*G1(3,1,:)) ...
    +(Tm1(1,2,:).*Tm1(2,4,:)-Tm1(1,4,:).*Tm1(2,2,:)).*(G1(2,1,:).*G1(4,2,:)-G1(2,2,:).*G1(4,1,:)) ...
    +(Tm1(1,3,:).*Tm1(2,4,:)-Tm1(1,4,:).*Tm1(2,3,:)).*(G1(3,1,:).*G1(4,2,:)-G1(3,2,:).*G1(4,1,:)) ...
    );
