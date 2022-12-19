function [Ua,para]=calculU_fint_ncouche_semi_espace_iso(omega,para,MAT)
% METHODE : @ la resolution de l eq de dispersion donne les valeurs possibles de k1
%           @ calcul de la solution homogene sous forme exp
%           @ calcul de la contrainte issue de la solution homogene sous forme sinus
%           @ calcul de la solution particuliere
%           @ calcul de la contrainte issue de la solution particuliere
%           @ calcul des coeff de la sol h grace aux conditions limites
tic
para.ncouche=4;

k2          = para.k2;
Kpos2       = k2;
Kneg2       = k2(2:(para.nbptkx/2));


para.fint(1,:)= [1 0 1];
para.fint(2,:)= [0 1 1];
para.zs     = [5.5 5.5 5.5];%0:.01:2;
nfint       = length(para.zs);%para.nfint;%nombre de sources
para.zrec   = [0.5 4   5.8 8];
para.xrec   = [2.5 2.5 2.5 2.5];
nrec        = length(para.zrec);

para.X20_G=.2;
FWHM=para.X20_G/sqrt(8*log(2));
mxspec=exp(-(Kpos2*FWHM).^2/2);
para.Tps_G=3e-2;
FWHM=para.Tps_G/sqrt(8*log(2));
iom=0.05;
omg=((1:nw)-1)*2*pi*dfe-1i*iom;
mtspec=exp(-(omg*FWHM).^2/2);%.*exp(-complex(0,1)*omg*FWHM*10);%exp(-complex(0,1)*omg*FWHM*10) pour rendre le signal causal et l empecher d exploser aux tps long
para.zero_pad_w=5;

U1XW =zeros(nw,nrec,nfint);
U1KW =zeros(1,para.nbptkx);
U2XW =zeros(nw,nrec,nfint);
U2KW =zeros(1,para.nbptkx);
U1KW0=zeros(nw,para.nbptkx/2+1);
S11XW=zeros(nw,nrec,nfint);
S11KW=zeros(1,para.nbptkx);
S22XW=zeros(nw,nrec,nfint);
S22KW=zeros(1,para.nbptkx);
S12XW=zeros(nw,nrec,nfint);
S12KW=zeros(1,para.nbptkx);


clear lent
% for i=1:para.ncouche
% alpha(1)=1500;
% beta(1)=800;
% rho(1)=1000;
% end
% function [Ua,para]=calculU_fint_ncouche_semi_espace(k2,omega,para,MAT)
% METHODE : @ la resolution de l eq de dispersion donne les valeurs possibles de k1
%           @ calcul de la solution homogene sous forme exp
%           @ calcul de la contrainte issue de la solution homogene sous forme sinus
%           @ calcul de la solution particuliere
%           @ calcul de la contrainte issue de la solution particuliere
%           @ calcul des coeff de la sol h grace aux conditions limites


n=length(k2(1,:));

u1= zeros(2,n,para.ncouche);
u2= zeros(2,n,para.ncouche);
k1= zeros(2,n,para.ncouche);

A0      =zeros(2,n);
C0      =zeros(2,n);
A       =zeros(4,n);
Ampli   =zeros(4*(para.ncouche-1)+2,4*(para.ncouche-1)+2,n);
Fsource =zeros(4*(para.ncouche-1)+2,n);
Solution=zeros(4*(para.ncouche-1)+2,n);
detA    =zeros(n,1);
Gr      =zeros(2,2);

omegac=omega-1i*iom;
viscosite;

%*********************************************************************%
%  calcul des polarisations et nombre d'onde de la solution homogene  %
%*********************************************************************%

for ic=1:para.ncouche
    %*********************%
    %**     Milieu i    **%
    %*********************%
    
    
    % Determination des k1
    A(1,:)= MATnc(ic).MAT.rho*omegac^2-MATnc(ic).MAT.C(6,6)*k2.*k2;
    A(4,:)= MATnc(ic).MAT.rho*omegac^2-MATnc(ic).MAT.C(2,2)*k2.*k2;
    A(2,:)=(MATnc(ic).MAT.C(1,2)+MATnc(ic).MAT.C(6,6))*k2;
    A(3,:)=(MATnc(ic).MAT.C(1,2)+MATnc(ic).MAT.C(6,6))*k2;
    
    [tmpk1]     = sol_k1(MATnc(ic).MAT.C,A);
    
    
    MAT.C = MAT.TC+1i*omgcpl_plus*MAT.TN;
    
    k1(:,1,ic) =sqrt((w/alpha(ic))^2-k2.^2);
    k1(:,2,ic) =sqrt((w/beta(ic) )^2-k2.^2);
    
    if ic==para.ncouche
        %onde evanescente derniere couche, cf. convention exp(+ikx)
        k1(:,:,ic) =tmpk1.*sign(imag(tmpk1)).*(imag(tmpk1)~=0) + ...
            tmpk1.*(imag(tmpk1)==0);
    else
        k1(:,:,ic)  =tmpk1;
    end
    
    % et des vecteurs propres homogene : vecteur solution homogene I a une cst pres
    [tmpu1,tmpu2]=sol_homogene( k1(:,:,ic),MATnc(ic).MAT.C,A,n);
    u1(:,:,ic)=tmpu1;%solution pour k1+
    u2(:,:,ic)=tmpu2;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul de la solution homogene %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Ampli   =zeros(4*(para.ncouche-1)+2,4*(para.ncouche-1)+2,n);

for ic=1
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    
    for i=1:2
        A0(i,:)= -1i*(k1(i,:,ic).*MATnc(ic).MAT.C(1,1).*u1(i,:,ic) +k2.*MATnc(ic).MAT.C(1,2).*u2(i,:,ic)); %sigma_11
        C0(i,:)= -1i*MATnc(ic).MAT.C(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic)); %sigma_11
    end
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    
    mic	= (-1)^ic;
    
    %SIGMA11 EN 0
    Ampli(1,1,:)=mic*A0(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude longi sens +
    Ampli(1,2,:)=mic*A0(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude trans sens +
    if para.ncouche>1
        Ampli(1,3,:)=-mic*A0(1,:); %.*exp( 1i.*k1(1,:,ic).*0);		%Amplitude longi sens -
        Ampli(1,4,:)=-mic*A0(2,:); %.*exp( 1i.*k1(2,:,ic).*0);		%Amplitude trans sens -
    end
    %SIGMA12 EN 0
    Ampli(2,1,:)=mic* C0(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude longi sens +
    Ampli(2,2,:)=mic* C0(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude trans sens +
    if para.ncouche>1
        Ampli(2,3,:)=mic* C0(1,:); %.*exp( 1i.*k1(1,:,ic).*0);       %Amplitude longi sens -
        Ampli(2,4,:)=mic* C0(2,:); %.*exp( 1i.*k1(2,:,ic).*0);	    %Amplitude trans sens -
        
        exp1=exp(-1i.*k1(1,:,ic).*MATnc(ic).h);
        exp2=exp(-1i.*k1(2,:,ic).*MATnc(ic).h);
        exp3=exp( 1i.*k1(1,:,ic).*MATnc(ic).h);
        exp4=exp( 1i.*k1(2,:,ic).*MATnc(ic).h);
        
        %SIGMA11 EN h(ic)
        Ampli(3,1,:)=mic* A0(1,:).*exp1;		%Amplitude longi sens +
        Ampli(3,2,:)=mic* A0(2,:).*exp2;		%Amplitude trans sens +
        Ampli(3,3,:)=mic*-A0(1,:).*exp3;		%Amplitude longi sens -
        Ampli(3,4,:)=mic*-A0(2,:).*exp4;		%Amplitude trans sens -
        
        %SIGMA12 EN h(ic)
        Ampli(4,1,:)=mic* C0(1,:).*exp1;		%Amplitude longi sens +
        Ampli(4,2,:)=mic* C0(2,:).*exp2;		%Amplitude trans sens +
        Ampli(4,3,:)=mic* C0(1,:).*exp3;        %Amplitude longi sens -
        Ampli(4,4,:)=mic* C0(2,:).*exp4;	    %Amplitude trans sens -
        
        %U1 EN h(ic)
        Ampli(5,1,:)=mic* u1(1,:,ic).*exp1;		%Amplitude longi sens +
        Ampli(5,2,:)=mic* u1(2,:,ic).*exp2;		%Amplitude trans sens +
        Ampli(5,3,:)=mic* u1(1,:,ic).*exp3;		%Amplitude longi sens -
        Ampli(5,4,:)=mic* u1(2,:,ic).*exp4;		%Amplitude trans sens -
        
        %U2 EN h(ic)
        Ampli(6,1,:)=mic* u2(1,:,ic).*exp1;		%Amplitude longi sens +
        Ampli(6,2,:)=mic* u2(2,:,ic).*exp2;		%Amplitude trans sens +
        Ampli(6,3,:)=mic*-u2(1,:,ic).*exp3;     %Amplitude longi sens -
        Ampli(6,4,:)=mic*-u2(2,:,ic).*exp4;	    %Amplitude trans sens -
    end
end

for ic=2:(para.ncouche-1)
    %***********************************************************%
    %   calcul de la contrainte issue de la solution homogene   %
    %***********************************************************%
    
    for i=1:2
        A0(i,:)= -1i*(k1(i,:,ic).*MATnc(ic).MAT.C(1,1).*u1(i,:,ic) +k2.*MATnc(ic).MAT.C(1,2).*u2(i,:,ic));
        C0(i,:)= -1i*MATnc(ic).MAT.C(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic));
    end
    
    %//////////////////////////////////////////////////////////////%
    %// ecriture sous forme exponentielle (normal p/r au calcul) //%
    %//////////////////////////////////////////////////////////////%
    exp1= exp(-1i.*k1(1,:,ic).*MATnc(ic).h);
    exp2= exp(-1i.*k1(2,:,ic).*MATnc(ic).h);
    exp3= exp( 1i.*k1(1,:,ic).*MATnc(ic).h);
    exp4= exp( 1i.*k1(2,:,ic).*MATnc(ic).h);
    mic	= (-1)^ic;
    
    %SIGMA11 EN 0
    Ampli((ic-2)*4+2+1,(ic-1)*4+1,:)=mic* A0(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude longi sens +
    Ampli((ic-2)*4+2+1,(ic-1)*4+2,:)=mic* A0(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude trans sens +
    Ampli((ic-2)*4+2+1,(ic-1)*4+3,:)=mic*-A0(1,:); %.*exp( 1i.*k1(1,:,ic).*0);		%Amplitude longi sens -
    Ampli((ic-2)*4+2+1,(ic-1)*4+4,:)=mic*-A0(2,:); %.*exp( 1i.*k1(2,:,ic).*0);		%Amplitude trans sens -
    
    %SIGMA12 EN 0
    Ampli((ic-2)*4+2+2,(ic-1)*4+1,:)=mic* C0(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude longi sens +
    Ampli((ic-2)*4+2+2,(ic-1)*4+2,:)=mic* C0(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude trans sens +
    Ampli((ic-2)*4+2+2,(ic-1)*4+3,:)=mic* C0(1,:); %.*exp( 1i.*k1(1,:,ic).*0);      	%Amplitude longi sens -
    Ampli((ic-2)*4+2+2,(ic-1)*4+4,:)=mic* C0(2,:); %.*exp( 1i.*k1(2,:,ic).*0);	    %Amplitude trans sens -
    
    %U1 EN 0
    Ampli((ic-2)*4+2+3,(ic-1)*4+1,:)=mic* u1(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude longi sens +
    Ampli((ic-2)*4+2+3,(ic-1)*4+2,:)=mic* u1(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude trans sens +
    Ampli((ic-2)*4+2+3,(ic-1)*4+3,:)=mic* u1(1,:,ic); %.*exp( 1i.*k1(1,:,ic).*0);	%Amplitude longi sens -
    Ampli((ic-2)*4+2+3,(ic-1)*4+4,:)=mic* u1(2,:,ic); %.*exp( 1i.*k1(2,:,ic).*0);	%Amplitude trans sens -
    
    %U2 EN 0
    Ampli((ic-2)*4+2+4,(ic-1)*4+1,:)=mic* u2(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude longi sens +
    Ampli((ic-2)*4+2+4,(ic-1)*4+2,:)=mic* u2(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude trans sens +
    Ampli((ic-2)*4+2+4,(ic-1)*4+3,:)=mic*-u2(1,:,ic); %.*exp( 1i.*k1(1,:,ic).*0);	%Amplitude longi sens -
    Ampli((ic-2)*4+2+4,(ic-1)*4+4,:)=mic*-u2(2,:,ic); %.*exp( 1i.*k1(2,:,ic).*0);	%Amplitude trans sens -
    
    %SIGMA11 EN h(ic)
    Ampli((ic-2)*4+2+5,(ic-1)*4+1,:)=mic* A0(1,:).*exp1;		%Amplitude longi sens +
    Ampli((ic-2)*4+2+5,(ic-1)*4+2,:)=mic* A0(2,:).*exp2;		%Amplitude trans sens +
    Ampli((ic-2)*4+2+5,(ic-1)*4+3,:)=mic*-A0(1,:).*exp3;		%Amplitude longi sens -
    Ampli((ic-2)*4+2+5,(ic-1)*4+4,:)=mic*-A0(2,:).*exp4;		%Amplitude trans sens -
    
    %SIGMA12 EN h(ic)
    Ampli((ic-2)*4+2+6,(ic-1)*4+1,:)=mic* C0(1,:).*exp1;		%Amplitude longi sens +
    Ampli((ic-2)*4+2+6,(ic-1)*4+2,:)=mic* C0(2,:).*exp2;		%Amplitude trans sens +
    Ampli((ic-2)*4+2+6,(ic-1)*4+3,:)=mic* C0(1,:).*exp3;        %Amplitude longi sens -
    Ampli((ic-2)*4+2+6,(ic-1)*4+4,:)=mic* C0(2,:).*exp4;	    %Amplitude trans sens -
    
    %U1 EN h(ic)
    Ampli((ic-2)*4+2+7,(ic-1)*4+1,:)=mic* u1(1,:,ic).*exp1;		%Amplitude longi sens +
    Ampli((ic-2)*4+2+7,(ic-1)*4+2,:)=mic* u1(2,:,ic).*exp2;		%Amplitude trans sens +
    Ampli((ic-2)*4+2+7,(ic-1)*4+3,:)=mic* u1(1,:,ic).*exp3;		%Amplitude longi sens -
    Ampli((ic-2)*4+2+7,(ic-1)*4+4,:)=mic* u1(2,:,ic).*exp4;		%Amplitude trans sens -
    
    %U2 EN h(ic)
    Ampli((ic-2)*4+2+8,(ic-1)*4+1,:)=mic* u2(1,:,ic).*exp1;		%Amplitude longi sens +
    Ampli((ic-2)*4+2+8,(ic-1)*4+2,:)=mic* u2(2,:,ic).*exp2;		%Amplitude trans sens +
    Ampli((ic-2)*4+2+8,(ic-1)*4+3,:)=mic*-u2(1,:,ic).*exp3;   	%Amplitude longi sens -
    Ampli((ic-2)*4+2+8,(ic-1)*4+4,:)=mic*-u2(2,:,ic).*exp4;	    %Amplitude trans sens -
end

if para.ncouche>1
    for ic=para.ncouche
        %***********************************************************%
        %   calcul de la contrainte issue de la solution homogene   %
        %***********************************************************%
        for i=1:2
            A0(i,:)= -1i*(k1(i,:,ic).*MATnc(ic).MAT.C(1,1).*u1(i,:,ic) +k2.*MATnc(ic).MAT.C(1,2).*u2(i,:,ic));
            C0(i,:)= -1i*MATnc(ic).MAT.C(6,6)*( k1(i,:,ic).*u2(i,:,ic) +k2.*u1(i,:,ic));
        end
        mic	= (-1)^ic;
        
        %//////////////////////////////////////////////////////////////%
        %// ecriture sous forme exponentielle (normal p/r au calcul) //%
        %//////////////////////////////////////////////////////////////%
        
        %SIGMA11 EN 0
        Ampli((ic-2)*4+2+1,(ic-1)*4+1,:)=mic*-A0(1,:); %.*exp(1i.*k1(1,:,ic).*0);		%Amplitude longi sens +
        Ampli((ic-2)*4+2+1,(ic-1)*4+2,:)=mic*-A0(2,:); %.*exp(1i.*k1(2,:,ic).*0);		%Amplitude trans sens +
        %         Ampli((ic-2)*4+2+1,(ic-1)*4+3,:)=mic*-A0(1,:); %.*exp(-1i.*k1(1,:,ic).*0);		%Amplitude longi sens -
        %         Ampli((ic-2)*4+2+1,(ic-1)*4+4,:)=mic*-A0(2,:); %.*exp(-1i.*k1(2,:,ic).*0);		%Amplitude trans sens -
        
        %SIGMA12 EN 0
        Ampli((ic-2)*4+2+2,(ic-1)*4+1,:)=mic* C0(1,:); %.*exp( 1i.*k1(1,:,ic).*0);		%Amplitude longi sens +
        Ampli((ic-2)*4+2+2,(ic-1)*4+2,:)=mic* C0(2,:); %.*exp( 1i.*k1(2,:,ic).*0);		%Amplitude trans sens +
        %         Ampli((ic-2)*4+2+2,(ic-1)*4+3,:)=mic* C0(1,:); %.*exp(-1i.*k1(1,:,ic).*0);       %Amplitude longi sens -
        %         Ampli((ic-2)*4+2+2,(ic-1)*4+4,:)=mic* C0(2,:); %.*exp(-1i.*k1(2,:,ic).*0);	    %Amplitude trans sens -
        
        %U1 EN 0
        Ampli((ic-2)*4+2+3,(ic-1)*4+1,:)=mic* u1(1,:,ic); %.*exp( 1i.*k1(1,:,ic).*0);  	%Amplitude longi sens +
        Ampli((ic-2)*4+2+3,(ic-1)*4+2,:)=mic* u1(2,:,ic); %.*exp( 1i.*k1(2,:,ic).*0);     %Amplitude trans sens +
        %         Ampli((ic-2)*4+2+3,(ic-1)*4+3,:)=mic* u1(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude longi sens -
        %         Ampli((ic-2)*4+2+3,(ic-1)*4+4,:)=mic* u1(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude trans sens -
        
        %U2 EN 0
        Ampli((ic-2)*4+2+4,(ic-1)*4+1,:)=mic*-u2(1,:,ic); %.*exp( 1i.*k1(1,:,ic).*0);	%Amplitude longi sens +
        Ampli((ic-2)*4+2+4,(ic-1)*4+2,:)=mic*-u2(2,:,ic); %.*exp( 1i.*k1(2,:,ic).*0);	%Amplitude trans sens +
        %         Ampli((ic-2)*4+2+4,(ic-1)*4+3,:)=mic*-u2(1,:,ic); %.*exp(-1i.*k1(1,:,ic).*0);	%Amplitude longi sens -
        %         Ampli((ic-2)*4+2+4,(ic-1)*4+4,:)=mic*-u2(2,:,ic); %.*exp(-1i.*k1(2,:,ic).*0);	%Amplitude trans sens -
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%
% inversion de la matrice %
%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:n
    %phase de conditionnement
    detA(i)=det(Ampli(:,:,i))^(1/4);
    Ampli(:,:,i)=Ampli(:,:,i)/detA(i);
    Ampli(:,:,i)=inv(Ampli(:,:,i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% attention, pour calculer les differentes composantes
% il est nécessaire de subdiviser le calcul en fonction
% des parties paire et impaire (k2)
% cf G11, G22 paires et G12 impaire
% ce qui peut être fait en traitant séparement fint1 et fint2

for icf=1:nfint
    
    %identification de la couche dans laquelle se trouve la force
    ics =1;
    hc  =MATnc(1).h;
    while para.zs(icf)>hc && ics<para.ncouche
        ics=ics+1;
        hc=hc+MATnc(ics).h;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % identification des ondes %
    % k1(1,:) =k1L             %
    % k1(2,:) =k1T             %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    lent(1,:)=abs(k1(1,:,ics).^2+k2.^2);
    lent(2,:)=abs(k1(2,:,ics).^2+k2.^2);
    
    [~,ind]=sort(lent,1);
    tmp2=k1(:,:,ics);
    k1f=zeros(2,n);
    for i=1:n
        k1f(1,i)=tmp2(ind(1,i),i);%sqrt(omegac^2/( MATnc(ics).MAT.C(1,1)/MATnc(ics).MAT.rho)-k2.^2);
        k1f(2,i)=tmp2(ind(2,i),i);%sqrt(omegac^2/( MATnc(ics).MAT.C(6,6)/MATnc(ics).MAT.rho)-k2.^2);
    end
    k1f=-k1f.*(sign(imag(k1f)).*(sign(imag(k1f))~=0)+(sign(imag(k1f))==0));
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul du vecteur source et resolution des coefficients %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    % le calcul du vecteur source se fait en utilisant
    % la decomposition en ondes planes
    % cependant l expression du champ incident se reintroduit a la fin
    % sous la forme des fonctions de Hankel
    % de maniere a converger plus rapidement
    if ics>1
        htop=0;
        for inch=1:ics-1
            htop=htop+MATnc(inch).h;
        end
    else
        htop=0;
    end
    hbot=0;
    for inch=1:ics
        hbot=hbot+MATnc(inch).h;
    end
    
    zt      = htop-para.zs(icf);%x=z dans l epaisseur
    exp1t   = exp(-1i.*k1f(1,:).*abs(zt));
    exp2t   = exp(-1i.*k1f(2,:).*abs(zt));
    
    zb      = hbot-para.zs(icf);%x=z dans l epaisseur
    exp1b   = exp(-1i.*k1f(1,:).*abs(zb));
    exp2b   = exp(-1i.*k1f(2,:).*abs(zb));
    
    mic	=-(-1)^ics;
    
    %%%%%%%%%%%%%%%%%
    % interface sup %
    %%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul des deplacements %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    fac = 1/(4*pi*MATnc(ics).MAT.rho*omegac^2);
    %pair en k2
    G110=-1i*fac*( k1f(1,:).*exp1t         + k2.^2 .*exp2t./k1f(2,:));
    %pair en k2
    G220=-1i*fac*( k2.^2  .*exp1t./k1f(1,:)+k1f(2,:).*exp2t);
    %impair en k2
    G120=-k2*sign(zt)*1i*fac.*(exp1t-exp2t);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul des derivees utiles au calcul des contraintes normales %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % derivee p/x1 en X (x1=z)
    %pair
    G1110=-sign(zt)*fac*(k1f(1,:).^2.*exp1t+k2.^2     .*exp2t);
    %pair
    G2210=-sign(zt)*fac*(k2.^2     .*exp1t+k1f(2,:).^2.*exp2t);
    %impair
    G1210=-k2*fac.*(+k1f(1,:).*exp1t-k1f(2,:).*exp2t);
    
    % derivee p/x2 en X (x2=x)
    %impair
    G1120=-k2.*fac.*( k1f(1,:).*exp1t+ k2.^2./k1f(2,:) .*exp2t);
    %impair
    G2220=-k2.*fac.*(k2.^2   .*exp1t./k1f(1,:)+k1f(2,:).*exp2t);
    %pair
    G1220=-k2.^2*sign(zt).*fac.*(+exp1t-exp2t);
    
    if ics~=para.ncouche
        %pair en k2
        G11h=-1i*fac*( k1f(1,:).*exp1b + k2.^2 .*exp2b./k1f(2,:));
        %pair en k2
        G22h=-1i*fac*( k2.^2  .*exp1b./k1f(1,:)+k1f(2,:).*exp2b);
        %impair en k2
        G12h=-k2*sign(zb)*1i*fac.*(exp1b-exp2b);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % calcul des derivees utiles au calcul des contraintes normales %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % derivee p/x1 en X (x1=z)
        %pair
        G111h=-sign(zb)*fac*(k1f(1,:).^2.*exp1b+k2.^2     .*exp2b);
        %pair
        G221h=-sign(zb)*fac*(k2.^2     .*exp1b+k1f(2,:).^2.*exp2b);
        %impair
        G121h=-k2*fac.*(+k1f(1,:).*exp1b-k1f(2,:).*exp2b);
        
        % derivee p/x2 en X (x2=x)
        %impair
        G112h=-k2.*fac.*( k1f(1,:).*exp1b+ k2.^2./k1f(2,:) .*exp2b);
        %impair
        G222h=-k2.*fac.*(k2.^2   .*exp1b./k1f(1,:)+k1f(2,:).*exp2b);
        %pair
        G122h=-k2.^2*sign(zb).*fac.*(+exp1b-exp2b);
    end
    
    %indice de los terminos en el vector fuente
    if ics==1
        is110=1;%SIGMA11 EN 0
        is120=2;%SIGMA12 EN 0
        if para.ncouche>1
            is11h=3;%SIGMA11 EN h(ics)
            is12h=4;%SIGMA12 EN h(ics)
            iu1h =5;%U1 EN h(ics)
            iu2h =6;%U2 EN h(ics)
        end
    elseif ics>1
        is110=(ics-2)*4+2+1;%SIGMA11 EN 0
        is120=(ics-2)*4+2+2;%SIGMA12 EN 0
        iu10=(ics-2)*4+2+3;%U1 EN 0
        iu20=(ics-2)*4+2+4;%U2 EN 0
        if ics<=(para.ncouche-1)
            is11h=(ics-2)*4+2+5;%SIGMA11 EN h(ics)
            is12h=(ics-2)*4+2+6;%SIGMA12 EN h(ics)
            iu1h =(ics-2)*4+2+7;%U1 EN h(ics)
            iu2h =(ics-2)*4+2+8;%U2 EN h(ics)
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % calcul des conditions aux interface dues aux forces internes %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    MAT=MATnc(ics).MAT;
    for jfint=1:2
        if jfint==1
            fint1=para.fint(1,icf);
            fint2=0;
        else
            fint1=0;
            fint2=para.fint(2,icf);
        end
        if fint1~=0 || fint2~=0
            
            Fsource(is110,:)=MAT.C(1,1)*(fint1*G1110+fint2.*G1210) + MAT.C(1,2).*(fint2.*G2220+fint1*G1220);  %sigma11 en 0
            Fsource(is120,:)=MAT.C(6,6)*(fint1*(G1120+G1210)+fint2.*(G1220+G2210));                           %sigma12 en 0
            if ics>1
                Fsource(iu10,:)=fint1 *G110+fint2.*G120;%U1 EN 0
                Fsource(iu20,:)=fint1 *G120+fint2.*G220;%U2 EN 0
            end
            if (para.ncouche>1 && ics==1) || (ics<=(para.ncouche-1))
                Fsource(is11h,:)=MAT.C(1,1)*(fint1*G111h+fint2.*G121h) + MAT.C(1,2).*(fint2.*G222h+fint1*G122h);  %sigma11 en h
                Fsource(is12h,:)=MAT.C(6,6)*(fint1*(G112h+G121h)+fint2.*(G122h+G221h));                           %sigma12 en h
                Fsource(iu1h,:) =fint1*G11h+fint2.*G12h;%U1 EN h
                Fsource(iu2h,:) =fint1*G12h+fint2.*G22h;%U2 EN h
            end
            Fsource=Fsource*mic;
            
            %Resolution du systeme des conditions aux limites
            for i=1:n
                Fsource(:,i)    = Fsource(:,i)/detA(i);
                Solution(:,i)   = squeeze(Ampli(:,:,i))*squeeze(Fsource(:,i));
            end
            %                 Sol(:,:,jfint,icf) = Solution;
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % calcul des champs au niveau du recepteur %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            for ir=1:nrec
                %Ecriture de la solution en fonction de la position en x1 du
                %recepteur
                zr=para.zrec(ir);
                icr=1;
                while zr>MATnc(icr).h && icr<para.ncouche
                    zr  = zr-MATnc(icr).h;
                    icr = icr+1;
                end
                
                exp1=exp(-1i.*k1(1,:,icr).*zr);
                exp2=exp(-1i.*k1(2,:,icr).*zr);
                exp3=exp( 1i.*k1(1,:,icr).*zr);
                exp4=exp( 1i.*k1(2,:,icr).*zr);
                
                %calcul des deplacements diffractes
                if icr<para.ncouche
                    U1KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*u1(1,:,icr).* exp1 ...
                        +Solution(2+(icr-1)*4,:).*u1(2,:,icr).* exp2 ...
                        +Solution(3+(icr-1)*4,:).*u1(1,:,icr).* exp3 ...
                        +Solution(4+(icr-1)*4,:).*u1(2,:,icr).* exp4 ;
                else
                    U1KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*u1(1,:,icr).* exp3 ...
                        +Solution(2+(icr-1)*4,:).*u1(2,:,icr).* exp4 ;
                end
                tmp=find(isnan(U1KW));
                U1KW(tmp)=0; %#ok<FNDSB>
                
                if icr<para.ncouche
                    U2KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*u2(1,:,icr).* exp1 ...
                        +Solution(2+(icr-1)*4,:).*u2(2,:,icr).* exp2 ...
                        -Solution(3+(icr-1)*4,:).*u2(1,:,icr).* exp3 ...
                        -Solution(4+(icr-1)*4,:).*u2(2,:,icr).* exp4;
                else
                    U2KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*u2(1,:,icr).* exp3 ...
                        +Solution(2+(icr-1)*4,:).*u2(2,:,icr).* exp4 ;
                end
                tmp=find(isnan(U2KW));
                U2KW(tmp)=0; %#ok<FNDSB>
                
                
                %calcul des contraintes diffractees
                for i=1:2
                    A0(i,:)= -1i*(k1(i,:,icr).*MATnc(icr).MAT.C(1,1).*u1(i,:,icr) +k2.*MATnc(icr).MAT.C(1,2).*u2(i,:,icr));
                end
                if icr<para.ncouche
                    S11KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*A0(1,:).* exp1 ...
                        +Solution(2+(icr-1)*4,:).*A0(2,:).* exp2 ...
                        -Solution(3+(icr-1)*4,:).*A0(1,:).* exp3 ...
                        -Solution(4+(icr-1)*4,:).*A0(2,:).* exp4;
                else
                    S11KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*A0(1,:).* exp3 ...
                        +Solution(2+(icr-1)*4,:).*A0(2,:).* exp4;
                end
                tmp=find(isnan(S11KW));
                S11KW(tmp)=0; %#ok<FNDSB>
                
                for i=1:2
                    A0(i,:)= -1i*(k1(i,:,icr).*MATnc(icr).MAT.C(1,2).*u1(i,:,icr) +k2.*MATnc(icr).MAT.C(2,2).*u2(i,:,icr));
                end
                if icr<para.ncouche
                    S22KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*A0(1,:).* exp1 ...
                        +Solution(2+(icr-1)*4,:).*A0(2,:).* exp2 ...
                        -Solution(3+(icr-1)*4,:).*A0(1,:).* exp3 ...
                        -Solution(4+(icr-1)*4,:).*A0(2,:).* exp4;
                else
                    S22KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*A0(1,:).* exp3 ...
                        +Solution(2+(icr-1)*4,:).*A0(2,:).* exp4;
                end
                tmp=find(isnan(S22KW));
                S22KW(tmp)=0; %#ok<FNDSB>
                
                for i=1:2
                    C0(i,:)= -1i*MATnc(icr).MAT.C(6,6)*( k1(i,:,icr).*u2(i,:,icr) +k2.*u1(i,:,icr));
                end
                if icr<para.ncouche
                    S12KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*C0(1,:).* exp1 ...
                        +Solution(2+(icr-1)*4,:).*C0(2,:).* exp2 ...
                        +Solution(3+(icr-1)*4,:).*C0(1,:).* exp3 ...
                        +Solution(4+(icr-1)*4,:).*C0(2,:).* exp4;
                else
                    S12KW(1:n) = ...
                        +Solution(1+(icr-1)*4,:).*C0(1,:).* exp3 ...
                        +Solution(2+(icr-1)*4,:).*C0(2,:).* exp4;
                end
                tmp=find(isnan(S12KW));
                S12KW(tmp)=0; %#ok<FNDSB>
                
                %rajouts des champs incidents
                if icr==ics
                    zrs=para.zrec(ir)-para.zs(icf);
                    exp1t= exp(-1i.*k1f(1,:).*abs(zrs));
                    exp2t= exp(-1i.*k1f(2,:).*abs(zrs));
                    
                    fac=1/(4*pi*MATnc(icr).MAT.rho*omegac^2);
                    
                    %pair en k2
                    G11r=-1i*fac*( k1f(1,:)       .*exp1t + k2.^2./k1f(2,:) .*exp2t);
                    %pair en k2
                    G22r=-1i*fac*( k2.^2./k1f(1,:).*exp1t + k1f(2,:)        .*exp2t);
                    %impair en k2
                    G12r=-k2*sign(zrs)*1i*fac.*(exp1t-exp2t);
                    
                    %calcul des deplacements incidents
                    U1KW(1:n) = U1KW(1:n) + fint1*G11r + fint2.*G12r;
                    U2KW(1:n) = U2KW(1:n) + fint2*G22r + fint1.*G12r;
                    
                    %%%derivee p/x1 en X (x1=z)
                    %pair
                    G111=-sign(zrs)*fac*( k1f(1,:).^2.*exp1t + k2.^2      .*exp2t);
                    %pair
                    G221=-sign(zrs)*fac*( k2.^2      .*exp1t + k1f(2,:).^2.*exp2t);
                    %impair
                    G121=-k2*fac.*(+k1f(1,:).*exp1t - k1f(2,:).*exp2t);
                    
                    % derivee p/x2 en X (x2=x)
                    %impair
                    G112=-k2.*fac.*( k1f(1,:)        .*exp1t + k2.^2./k1f(2,:).*exp2t);
                    %impair
                    G222=-k2.*fac.*( k2.^2  .*exp1t./k1f(1,:)+ k1f(2,:).*exp2t);
                    %pair
                    G122=-k2.^2*sign(zrs).*fac.*( +exp1t - exp2t);
                    
                    %calcul des tractions t1j
                    S11KW(1:n) = S11KW(1:n) + MATnc(icr).MAT.C(1,1)*(fint1*G111+fint2.*G121) + MATnc(icr).MAT.C(1,2).*(fint2.*G222+fint1*G122);
                    S22KW(1:n) = S22KW(1:n) + MATnc(icr).MAT.C(2,2)*(fint1*G122+fint2.*G222) + MATnc(icr).MAT.C(1,2).*(fint1.*G111+fint2*G121);
                    S12KW(1:n) = S12KW(1:n) + MATnc(icr).MAT.C(6,6)*(fint1*(G112+G121)+fint2.*(G122+G221));
                end
                
                %correction spectre
                U1KW(1:n)=U1KW(1:n).*mxspec;
                if jfint==1
                    U1KW(para.nbptkx:-1:(para.nbptkx/2+2)) = U1KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                else
                    U1KW(para.nbptkx:-1:(para.nbptkx/2+2)) =-U1KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                end
                U1KW(1:para.nbptkx/2+1) = U1KW(1:para.nbptkx/2+1).*exp(-1i*para.xrec(ir)*Kpos2);
                U1KW0(ifq,:)=U1KW(1:n);
                tmp=DK*2*pi*para.nbptkx*ifft(U1KW,para.nbptkx);
                U1XW(ifq,ir,icf)=U1XW(ifq,ir,icf)+tmp(1);
                
                %correction spectre
                U2KW(1:n)=U2KW(1:n).*mxspec;
                if jfint==1
                    U2KW(para.nbptkx:-1:(para.nbptkx/2+2)) =-U2KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                else
                    U2KW(para.nbptkx:-1:(para.nbptkx/2+2)) = U2KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                end
                U2KW(1:para.nbptkx/2+1) = U2KW(1:para.nbptkx/2+1).*exp(-1i*para.xrec(ir)*Kpos2);
                tmp=DK*2*pi*para.nbptkx*ifft(U2KW,para.nbptkx);
                U2XW(ifq,ir,icf)=U2XW(ifq,ir,icf)+tmp(1);
                
                %correction spectre
                S11KW(1:n)=S11KW(1:n).*mxspec;
                if jfint==1
                    S11KW(para.nbptkx:-1:(para.nbptkx/2+2)) = S11KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                else
                    S11KW(para.nbptkx:-1:(para.nbptkx/2+2)) =-S11KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                end
                S11KW(1:para.nbptkx/2+1) = S11KW(1:para.nbptkx/2+1).*exp(-1i*para.xrec(ir)*Kpos2);
                tmp=DK*2*pi*para.nbptkx*ifft(S11KW,para.nbptkx);
                S11XW(ifq,ir,icf)=S11XW(ifq,ir,icf)+tmp(1);
                
                %correction spectre
                S22KW(1:n)=S22KW(1:n).*mxspec;
                if jfint==1
                    S22KW(para.nbptkx:-1:(para.nbptkx/2+2)) = S22KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                else
                    S22KW(para.nbptkx:-1:(para.nbptkx/2+2)) =-S22KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                end
                S22KW(1:para.nbptkx/2+1) = S22KW(1:para.nbptkx/2+1).*exp(-1i*para.xrec(ir)*Kpos2);
                tmp=DK*2*pi*para.nbptkx*ifft(S22KW,para.nbptkx);
                S22XW(ifq,ir,icf)=S22XW(ifq,ir,icf)+tmp(1);
                
                %correction spectre
                S12KW(1:n)=S12KW(1:n).*mxspec;
                if jfint==1
                    S12KW(para.nbptkx:-1:(para.nbptkx/2+2)) =-S12KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                else
                    S12KW(para.nbptkx:-1:(para.nbptkx/2+2)) = S12KW(2:(para.nbptkx/2)).*exp( 1i*para.xrec(ir)*Kneg2);
                end
                S12KW(1:para.nbptkx/2+1) = S12KW(1:para.nbptkx/2+1).*exp(-1i*para.xrec(ir)*Kpos2);
                tmp=DK*2*pi*para.nbptkx*ifft(S12KW,para.nbptkx);
                S12XW(ifq,ir,icf)=S12XW(ifq,ir,icf)+tmp(1);
                
            end
        end
        
        %             clear g
        %             rij     = abs(sqrt((para.zs-para.zrec).^2+para.xrec.^2));
        %             g(1,:)  = (para.zrec-para.zs)./rij;
        %             g(2,:)  = para.xrec./rij;
        %         solucion analitica para el campo incidente
        %         alpha= sqrt(MATnc(ics).MAT.C(1,1)/MATnc(ics).MAT.rho);
        %         beta = sqrt(MATnc(ics).MAT.C(6,6)/MATnc(ics).MAT.rho);
        %         kp   = omg/alpha;
        %         ks   = omg/beta;
        %
        %         h0P = besselh(0,2,kp*rij);
        %         h0S = besselh(0,2,ks*rij);
        %         h2P = besselh(2,2,kp*rij);
        %         h2S = besselh(2,2,ks*rij);
        %
        %         AA   = h0P/MATnc(ics).MAT.C(1,1)+h0S/MATnc(ics).MAT.C(6,6);
        %         BB   = h2P/MATnc(ics).MAT.C(1,1)-h2S/MATnc(ics).MAT.C(6,6);
        %
        %         d   = eye(2);
        %         i=1;
        %         for j=1:2
        %             Gr(i,j)=1/(1i*8)*(d(i,j)*AA-(2*g(i).*g(j)-d(i,j)).*BB);
        %         end
        %
        %         i=2;j=2;
        %         Gr(i,j)=1/(1i*8)*(d(i,j)*AA-(2*g(i).*g(j)-d(i,j)).*BB);
        %
        %         U1XW(ifq)=U1XW(ifq)+para.fint1*Gr(1,1)+para.fint2.*Gr(1,2);
        
        %             U_2=para.fint2.*G22+para.fint1*G12;
        
        
    end
end