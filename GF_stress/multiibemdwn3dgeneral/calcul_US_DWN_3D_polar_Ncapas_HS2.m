function [UXW_fx,SXW_fx,UXW_fy,SXW_fy,UXW_fz,SXW_fz]=calcul_US_DWN_3D_polar_Ncapas_HS2(para,rec,salu,sals,coordf,fij,DWN)
% METHODE : @ recuperacion de kz,kr y de las matrices de las soluciones
%             homogenea
%           @ calculo del vector fuente
%           @ calculo de la solucion
%           zr=zr0(izr0)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% calcul des conditions aux interface dues aux forces internes, %
%    calcul du vecteur source et resolution des coefficients    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% el vector fuente se calcula usando una descomposicion en ondas
% cilindricas y corresponde a la evaluacion del campo incidente al
% nivel de la interface superior y inferior de los estratos en los cuales
% pertenecen las fuentes
%
% en esta version 2 se intenta acelerar los calculos
% las fuentes de misma posicion en z son procesadas al mismo tiempo
% no se considera integracion de las fuentes en campo cercano para el campo
% difractado (DWN), solamente para el campo incidente
% asi se considera para el calculo del campo difractado solamente una
% posicion zs, y varias posciones xs, ys
% al final del programa se le agrega el campo incidente,
% asi que coordf sigue teniendo la informacion de las normales y
% superificie de cada fuente

ncapas      = para.nsubmed;
kr          = DWN.kr;
nkr       	= length(kr);
kz          = DWN.kz;
A_DWN       = DWN.A_DWN; %matrices ya inversada
B_DWN       = DWN.B_DWN;

gaussian       = para.gaussian;
if para.dim == 3
gaussex     = para.gaussex;
end
MAT         = para.reg(1).sub;
DK          = DWN.dkr;

xs          = coordf.xs;
ys          = coordf.ys;
zs          = coordf.zs;
nxs       	= length(xs);

xr          = rec.xr;
yr          = rec.yr;
zr          = rec.zr;
nrec        = length(xr);

zr0         = rec.zr0;
izr0        = rec.izr0;
nrecz       = length(zr0);
zricrall    = rec.zricrall;
icrall      = rec.icrall;

%campo de desplazamientos en funcion de cada componente de la fuerza
UXW_fx      = zeros(3,nrec,nxs);
UXW_fy      = zeros(3,nrec,nxs);
UXW_fz      = zeros(3,nrec,nxs);

SXW_fx    	= zeros(3,3,nrec,nxs);
SXW_fy    	= zeros(3,3,nrec,nxs);
SXW_fz    	= zeros(3,3,nrec,nxs);

Ur          = zeros(3,nkr);
Uz          = zeros(3,nkr);
Szr         = zeros(3,nkr);
Szz         = zeros(3,nkr);
Srr         = zeros(3,nkr);

FsourcePSV  = zeros(4*(ncapas-1)+2,nkr);
SolutionPSV = zeros(4*(ncapas-1)+2,nkr);
FsourceSH   = zeros(2*(ncapas-1)+1,nkr);
SolutionSH  = zeros(2*(ncapas-1)+1,nkr);

% para no perder tiempo a calcular contribuciones que contribuen poco
if sum(sum(fij)) ~= size(fij,1)*size(fij,2)
for ixs=1:nxs
    maxf=max(abs(fij(ixs,:)));
    for j=1:3
        if abs(fij(ixs,j))<1e-3*maxf
            fij(ixs,j)=0;
        end
    end
end
end
%str = toc; disp(['calcu_0 ',num2str(str)])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         calculo del vector fuente          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%identification des couches dans lesquelles se trouvent les forces et
%les coordonnees z auxquels doivent etre evaluees les conditions aux
%interface
ics = 1;
htop= 0;
hbot= 0;
hc 	= MAT(1).h;
hbot= hbot+hc;
while zs(1)>hc && ics<ncapas
    ics	= ics+1;
    htop= hc;
    hc 	= hc+MAT(ics).h;
    hbot= hc;
end

%posicion de las fuentes para la resolucion de las amplitudes difractadas
coordfi.xs	= 0;
coordfi.ys	= 0;
coordfi.zs	= zs(1);

%posiciones de las interfases (no son los receptores!), no mas se ocupa
%para resolver el sistema
coordr.z(1) = htop;
if ics~=ncapas
    coordr.z(2)=hbot;
end
coordr.x    = 0*coordr.z;
coordr.y    = 0*coordr.z;

%signe oppose a celui de la matrice A utile a l'arrangement des termes
%incident et diffracte
mic	=-(-1)^ics;

%indice de los terminos en el vector fuente
if ics==1
    iszz0=1;%SIGMAzz EN 0
    iszr0=2;%SIGMAzr EN 0
    if ncapas>1
        iszzh=3;%SIGMAzz EN h(ics)
        iszrh=4;%SIGMAzr EN h(ics)
        iuzh =5;%Uz EN h(ics)
        iurh =6;%Ur EN h(ics)
    end
elseif ics>1
    iszz0=(ics-2)*4+2+1;%SIGMAzz EN 0
    iszr0=(ics-2)*4+2+2;%SIGMAzr EN 0
    iuz0 =(ics-2)*4+2+3;%Uz EN 0
    iur0 =(ics-2)*4+2+4;%Ur EN 0
    if ics<=(ncapas-1)
        iszzh=(ics-2)*4+2+5;%SIGMAzz EN h(ics)
        iszrh=(ics-2)*4+2+6;%SIGMAzr EN h(ics)
        iuzh =(ics-2)*4+2+7;%Uz EN h(ics)
        iurh =(ics-2)*4+2+8;%Ur EN h(ics)
    end
end

%calculo del vector fuente
%el numero 2 al final de los argumentos de la funcion es para indicar
%que no se hace la suma, ni tampoco se multiplica por los coeficientes
%respectivos a las fuerzas horizontal y verticales, y que se restringe
%el calculo para solamente la coordenadas polares
[~,Gp,Spz,Spx]=Gij_3D_DWN_polar(para,DWN,coordr,coordfi,ics,2);

%str = toc; disp(['calcu_1 ',num2str(str)])


%%%%%%di
% fx %
%%%%%%
if max(abs(fij(:,1)))~=0 || max(abs(fij(:,2)))~=0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % resolucion de las amplitudes difractadas para una fuerza horizontal %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %3 eme indice de Spx ou Gp correspond au recepteur (1 : interface du
    %haut;2 : interface du bas)
    %4 eme source
    %5 eme kr
    FsourcePSV(iszz0,:)    = Spx(3,3,1,1,:);      %sigmazz en 0
    FsourcePSV(iszr0,:)    = Spx(3,1,1,1,:);      %sigmazr en 0
    if ics>1
        FsourcePSV(iuz0,:) = Gp(3,1,1,1,:);   	%Uz en 0
        FsourcePSV(iur0,:) = Gp(1,1,1,1,:);   	%Ur en 0
    end
    if ics<=(ncapas-1)
        FsourcePSV(iszzh,:)= Spx(3,3,2,1,:);      %sigmazz en h
        FsourcePSV(iszrh,:)= Spx(3,1,2,1,:);      %sigmazr en h
        FsourcePSV(iuzh,:) = Gp(3,1,2,1,:);       %Uz en h
        FsourcePSV(iurh,:) = Gp(1,1,2,1,:);       %Ur en h
    end
    FsourcePSV= FsourcePSV*mic;
    
    FsourceSH(iszr0/2,:)    = Spx(2,2,1,1,:);   	%sigmazr en 0 %falsos indices
    if ics>1
        FsourceSH(iur0/2,:) = Gp(2,2,1,1,:);   	%Ur en 0 %falsos indices
    end
    if ics<=(ncapas-1)
        FsourceSH(iszrh/2,:)= Spx(2,2,2,1,:);     %sigmazr en h %falsos indices
        FsourceSH(iurh/2,:) = Gp(2,2,2,1,:);      %Ur en h %falsos indices
    end
    FsourceSH= FsourceSH*mic;
    
    %Resolution du systeme des conditions aux limites
    % for i=1:nkr
    %    SolutionPSV(:,i)   = squeeze(A_DWN(:,:,i))*squeeze(FsourcePSV(:,i));
    % end
    %attention dependance de A pour force horizontal
    for j=1:4*(ncapas-1)+2
        SolutionPSV(j,:)= sum(squeeze(A_DWN(j,:,:)).*FsourcePSV,1);
    end
    if ncapas==1
        SolutionSH(1,:)	= (FsourceSH.*squeeze(B_DWN(1,:,:)).');
    else
        for j=1:2*(ncapas-1)+1
            SolutionSH(j,:)	= sum(squeeze(B_DWN(j,:,:)).*FsourceSH,1);
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evaluacion de los campos difractados en los receptores %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%str = toc; disp(['calcu_2 ',num2str(str)])
    for irz=1:nrecz %ciclo sobre las profundidades de los receptores
        
        ixr     = find(izr0==irz);
        nrecx   = length(ixr);
        
        %Ecriture de la solution en fonction de la position du
        %recepteur et de la distance ds le plan (r,theta) entre le
        %recepteur et la source
         
        zricr   = zricrall(irz); % profundidad relativa a la interface de la capa del receptor
        icr     = icrall(irz);   % couche du recepteur et profondeur relative
        ksi 	= para.reg(1).sub(icr).ksi;
        kpi 	= para.reg(1).sub(icr).kpi;
        exp1    = exp(-1i.*kz(1,:,icr).* zricr);
        exp2    = exp(-1i.*kz(2,:,icr).* zricr);
        exp3    = exp( 1i.*kz(1,:,icr).*(zricr-MAT(icr).h));
        exp4    = exp( 1i.*kz(2,:,icr).*(zricr-MAT(icr).h));
        
        if max(salu(ixr))==1
            %coeficientes de los desplazamientos en funcion de los
            %potenciales
            Uz(1,:) =-1i*kz(2,:,icr);                                   %Amplitude Phi sens z +
            Uz(2,:) = kr.^2;                                            %Amplitude Psi z +
            Uz(3,:) = 0;                                                %Amplitude Chi z +
            
            Ur(1,:) = 1;                                                %Amplitude Phi z +
            Ur(2,:) =-1i*kz(1,:,icr);                                   %Amplitude Psi z +
            Ur(3,:) = 1;                                                %Amplitude Chi z +
            
            %Ut=Ur;
            if icr<ncapas
                %coeficientes de los desplazamientos en funcion de los
                %potenciales
                facUPSV=...
                    +Ur(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Ur(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    +Ur(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4 ...
                    -Ur(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
                facUSH=...
                    +Ur(3,:).*SolutionSH(1+(icr-1)*2,:).*exp1 ...
                    +Ur(3,:).*SolutionSH(2+(icr-1)*2,:).*exp3;
                facUz=...
                    +Uz(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Uz(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    -Uz(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4 ...
                    +Uz(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
            else
                facUPSV=...
                    +Ur(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Ur(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
                facUSH=...
                    +Ur(3,:).*SolutionSH(1+(icr-1)*2,:).*exp1;
                facUz=...
                    +Uz(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Uz(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
            end
        end
        
        if max(sals(ixr))==1
            %coeficientes de los esfuerzos en funcion de los
            %potenciales
            xi = kz(1,:,icr).^2-kr.^2;
            
            Szz(1,:)=-xi;                                           %Amplitude Phi sens z +
            Szz(2,:)=-2i*kz(1,:,icr).*kr.^2;                      	%Amplitude Psi z +
            Szz(3,:)= 0;                                         	%Amplitude Chi z +
            
            Szr(1,:)=-2i*kz(2,:,icr);                            	%Amplitude Phi z +
            Szr(2,:)=-xi;                                        	%Amplitude Psi z +
            Szr(3,:)=-1i*kz(1,:,icr);                             	%Amplitude Chi z +
            
            Srr(1,:)= 1;                                            %Amplitude Phi sens z +
            Srr(2,:)=-1i*kz(1,:,icr);                               %Amplitude Psi z +
            Srr(3,:)= 1;                                         	%Amplitude Chi z +
            
            %Szt = Szr;
            %Srt = Srr;
            %Stt = Srr;
            
            if icr<ncapas
                facSzrPSV=...
                    +Szr(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Szr(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    -Szr(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4 ...
                    +Szr(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
                facSzrSH=...
                    +Szr(3,:).*SolutionSH(1+(icr-1)*2,:).*exp1 ...
                    -Szr(3,:).*SolutionSH(2+(icr-1)*2,:).*exp3;
                facSzz= ...
                    +Szz(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Szz(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    +Szz(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4 ...
                    -Szz(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
                
                facSrrP= ...
                    +Srr(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Srr(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4;
                facSrrSV=...
                    +Srr(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    -Srr(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
                facSrrSH=...
                    +Srr(3,:).*SolutionSH(1+(icr-1)*2,:).*exp1 ...
                    +Srr(3,:).*SolutionSH(2+(icr-1)*2,:).*exp3;
            else
                facSzrPSV=...
                    +Szr(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Szr(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
                facSzrSH=...
                    +Szr(3,:).*SolutionSH(1+(icr-1)*2,:).*exp1;
                facSzz= ...
                    +Szz(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Szz(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
                
                facSrrP= ...
                    +Srr(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2;
                facSrrSV=...
                    +Srr(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
                facSrrSH=...
                    +Srr(3,:).*SolutionSH(1+(icr-1)*2,:).*exp1;
            end
        end
        
        for ixs=1:nxs %ciclo sobre las posiciones radiales de las fuentes
            if max(fij(ixs,1))~=0 || max(fij(ixs,2))~=0
                for irx=1:nrecx%ciclo sobre las posiciones radiales de los receptores
                    ir = ixr(irx);
                    
                    Xfs     = xr(ir)-xs(ixs);
                    Yfs     = yr(ir)-ys(ixs);
                    r       = sqrt(Xfs^2+Yfs^2);
                    
                    J0      = besselj(0,kr*r);
                    J1      = besselj(1,kr*r);
                    dJ1     = kr.*J0-J1/r;%kr.*(J0-J1./(kr*r));
                    d2J1    =-kr/r.*J0+J1.*(2/r^2-kr.^2);
                    if r<1e-6
                        ct  = 1/sqrt(2);
                        st  = 1/sqrt(2);
                    else
                        ct  = Xfs/r;
                        st  = Yfs/r;
                    end
                    
                    % calculo de los desplazamientos difractados en caso de que
                    % esten pedidos
                    if salu(ir)==1
                        Urx = dJ1 .*facUPSV+J1/r.*facUSH;
                        Utx =-J1/r.*facUPSV-dJ1 .*facUSH;
                        Uzx =   J1.*facUz;
                        
                        Urx(isnan(Urx))=0;
                        Utx(isnan(Utx))=0;
                        Uzx(isnan(Uzx))=0;
                        
                        Urx	= sum(DK.*Urx);
                        Utx = sum(DK.*Utx);
                        Uzx = sum(DK.*Uzx);
                        
                        UXW_fx(1,ir,ixs) = fij(ixs,1)*(Urx*ct^2-Utx*st^2);          %Uxx
                        UXW_fx(2,ir,ixs) = fij(ixs,1)*(Urx+Utx)*ct*st;              %Uyx
                        UXW_fx(3,ir,ixs) = fij(ixs,1)* Uzx*ct;                      %Uzx
                        
                        UXW_fy(1,ir,ixs) = fij(ixs,2)*(Urx+Utx)*ct*st;              %Uxy=Uyx, G12=G21
                        UXW_fy(2,ir,ixs) = fij(ixs,2)*(Urx*st^2-Utx*ct^2);          %Uyy=Uxx, th= -(pi/2-th), ct->st,st>-ct
                        UXW_fy(3,ir,ixs) = fij(ixs,2)* Uzx*st;
                    end
                    
                    % calculo de los esfuerzos difractados en caso de que
                    % esten pedidos
                    if sals(ir)==1
                        Szrx = dJ1 .*facSzrPSV + J1/r.*facSzrSH;
                        Sztx = J1/r.*facSzrPSV + dJ1 .*facSzrSH;
                        Szzx = J1  .*facSzz;
                        Srrx = 2*((d2J1-(ksi^2/2-kpi^2)*J1)         .*facSrrP + d2J1          .*facSrrSV + (dJ1/r-J1/r^2).*facSrrSH);
                        Sttx = 2*( (-(ksi^2/2-kpi^2+1/r^2)*J1+dJ1/r).*facSrrP + (dJ1/r-J1/r^2).*facSrrSV + (J1/r^2-dJ1/r).*facSrrSH);
                        Srtx = 2*(dJ1/r-J1/r^2).*(facSrrP+facSrrSV)                               + (d2J1-dJ1/r+J1/(r^2)).*facSrrSH;
                        
                        Szrx(isnan(Szrx))=0;
                        Sztx(isnan(Sztx))=0;
                        Szzx(isnan(Szzx))=0;
                        Srrx(isnan(Srrx))=0;
                        Sttx(isnan(Sttx))=0;
                        Srtx(isnan(Srtx))=0;
                        
                        Szrx = sum(DK.*Szrx);
                        Sztx = sum(DK.*Sztx);
                        Szzx = sum(DK.*Szzx);
                        Srrx = sum(DK.*Srrx);
                        Sttx = sum(DK.*Sttx);
                        Srtx = sum(DK.*Srtx);
                        
                        % Sp=[Srrx Srtx Srzx; Strx Sttx Stzx; Szrx Sztx Szzx];
                        Sp=[ct*Srrx -st*Srtx ct*Szrx; -st*Srtx ct*Sttx -st*Sztx; ct*Szrx -st*Sztx ct*Szzx];
                        T =[ ct    -st    0  ;   st    ct    0  ;   0     0     1  ];
                        Sc= T*Sp*(T.');
                        %Sc=[SxxxW SxyxW SxzxW; SyxxW SyyxW SyzxW; SzxxW SzyxW SzzxW];
                        SXW_fx(:,:,ir,ixs) = fij(ixs,1)*MAT(icr).Ci(6,6)*Sc;
                        
                        % theta de la fuerza solamente: th= -(pi/2-th), ct->st,st>-ct
                        Sp=[st*Srrx ct*Srtx st*Szrx; ct*Srtx st*Sttx ct*Sztx; st*Szrx ct*Sztx st*Szzx];
                        T =[ ct    -st    0  ;   st    ct    0  ;   0     0     1  ];
                        Sc= T*Sp*(T.');
                        SXW_fy(:,:,ir,ixs) = fij(ixs,2)*MAT(icr).Ci(6,6)*Sc;
                    end
                end
            end
        end
        
%str = toc; disp(['calcu_3 :',num2str(irz),' ',num2str(str)])
    end
end


%%%%%%
% fz %
%%%%%%
if max(abs(fij(:,3)))~=0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % resolucion de las amplitudes difractadas para una fuerza vertical %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %3 eme indice de Spz ou Gp correspond au recepteur (1 : interface du
    %haut;2 : interface du bas)
    %4 eme source
    %5 eme kr
    FsourcePSV(iszz0,:)    = Spz(3,3,1,1,:);	%sigmazz en 0
    FsourcePSV(iszr0,:)    = Spz(3,1,1,1,:);	%sigmazr en 0
    if ics>1
        FsourcePSV(iuz0,:) = Gp(3,3,1,1,:);   	%Uz en 0
        FsourcePSV(iur0,:) = Gp(1,3,1,1,:);   	%Ur en 0
    end
    if ics<=(ncapas-1)
        FsourcePSV(iszzh,:)= Spz(3,3,2,1,:); %sigmazz en h
        FsourcePSV(iszrh,:)= Spz(3,1,2,1,:); %sigmazr en h
        FsourcePSV(iuzh,:) = Gp(3,3,2,1,:); 	%Uz en h
        FsourcePSV(iurh,:) = Gp(1,3,2,1,:); 	%Ur en h
    end
    FsourcePSV= FsourcePSV*mic;
    
    %Resolution du systeme des conditions aux limites
    %         for i=1:nkr
    %             SolutionPSV(:,i)   = squeeze(A_DWN(:,:,i))*squeeze(FsourcePSV(:,i));
    %         end
    %attention dependance de A pour force horizontal
    for j=1:4*(ncapas-1)+2
        SolutionPSV(j,:)   = sum(squeeze(A_DWN(j,:,:)).*FsourcePSV,1);
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % evaluacion de los campos difractados en los receptores %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%str = toc; disp(['calcu_4 ',num2str(str)])
    for irz=1:nrecz %ciclo sobre las profundidades de los receptores
        
        ixr     = find(izr0==irz);
        nrecx   = length(ixr);
        %Ecriture de la solution en fonction de la position du
        %recepteur et de la distance ds le plan (r,theta) entre le
        %recepteur et la source
        zricr   = zricrall(irz); % profundidad relativa a la interface de la capa del receptor
        icr     = icrall(irz);   %couche du recepteur et profondeur relative
        ksi     = para.reg(1).sub(icr).ksi;
        kpi     = para.reg(1).sub(icr).kpi;
        exp1    = exp(-1i.*kz(1,:,icr).* zricr);
        exp2    = exp(-1i.*kz(2,:,icr).* zricr);
        exp3    = exp( 1i.*kz(1,:,icr).*(zricr-MAT(icr).h));
        exp4    = exp( 1i.*kz(2,:,icr).*(zricr-MAT(icr).h));
        
        if max(salu(ixr))==1
            %calcul des coefficients des deplacements diffractes

            Uz(1,:) =-1i*kz(2,:,icr);                           	%Amplitude Phi sens z +
            Uz(2,:) = kr.^2;                                      	%Amplitude Psi z +
            
            Ur(1,:) = 1;                                            %Amplitude Phi z +
            Ur(2,:) =-1i*kz(1,:,icr);                              	%Amplitude Psi z +
            
            if icr<ncapas
                farUrz= ...
                    +Ur(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Ur(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    +Ur(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4 ...
                    -Ur(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
                facUzz= ...
                    +Uz(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Uz(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    -Uz(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4 ...
                    +Uz(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
            else
                farUrz= ...
                    +Ur(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Ur(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
                facUzz= ...
                    +Uz(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Uz(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
            end
        end
        
        if max(sals(ixr))==1
            %calcul des coefficients des contraintes diffractees
            xi      = kz(1,:,icr).^2-kr.^2;
            
            Szz(1,:)=-xi;                                           %Amplitude Phi sens z +
            Szz(2,:)=-2i*kz(1,:,icr).*kr.^2;                      	%Amplitude Psi z +
            
            Szr(1,:)=-2i*kz(2,:,icr);                            	%Amplitude Phi z +
            Szr(2,:)=-xi;                                        	%Amplitude Psi z +
            
            Srr(1,:)= 1;                                            %Amplitude Phi sens z +
            Srr(2,:)=-1i*kz(1,:,icr);                               %Amplitude Psi z +
            
            %Stt = Srr;
            
            if icr<ncapas
                facSzr= ...
                    +Szr(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Szr(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    -Szr(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4 ...
                    +Szr(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
                facSzz= ...
                    +Szz(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Szz(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    +Szz(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4 ...
                    -Szz(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
                facSrrP= ...
                    +Srr(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Srr(1,:).*SolutionPSV(3+(icr-1)*4,:).*exp4;
                facSrrSV= ...
                    +Srr(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1 ...
                    -Srr(2,:).*SolutionPSV(4+(icr-1)*4,:).*exp3;
            else
                facSzr= ...
                    +Szr(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Szr(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
                facSzz= ...
                    +Szz(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2 ...
                    +Szz(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
                facSrrP= ...
                    +Srr(1,:).*SolutionPSV(1+(icr-1)*4,:).*exp2;
                facSrrSV= ...
                    +Srr(2,:).*SolutionPSV(2+(icr-1)*4,:).*exp1;
            end
        end
        
        for ixs=1:nxs %ciclo sobre las posiciones radiales de las fuentes
            if fij(ixs,3)~=0 
                for irx=1:nrecx%ciclo sobre las posiciones radiales de los receptores
                    ir = ixr(irx);
                    
                    Xfs = xr(ir)-xs(ixs);
                    Yfs = yr(ir)-ys(ixs);
                    r   = sqrt(Xfs^2+Yfs^2);
                    J0 	= besselj(0,kr*r);
                    J1 	= besselj(1,kr*r);
                    dJ1	= kr.*J0-J1/r;%kr.*(J0-J1./(kr*r));
                    dJ0 =-kr.*J1;
                    d2J0=-kr.*dJ1;
                    
                    if r<1e-6
                        ct  = 1/sqrt(2);
                        st  = 1/sqrt(2);
                    else
                        ct  = Xfs/r;
                        st  = Yfs/r;
                    end
                    
                    if salu(ir)==1
                        Urz = dJ0.*farUrz;
                        Uzz = J0 .*facUzz;
                        
                        Urz(isnan(Urz))=0;
                        Uzz(isnan(Uzz))=0;
                        %disp([sum(DK) length(DK)])
%                         DKacum = zeros(1,length(DK));
%                         for idkacum = 1:length(DK)
%                         DKacum(idkacum) = sum(DK(1:idkacum));
%                         end
%                         figure(424); set(gcf,'Name','Uzz'); hold on
%                         scatter3(DKacum, ...
%                           ones(1,length(DK))*real(DWN.omegac),...
%                           abs(Uzz),...
%                           [],abs(Uzz),'.');
%                         figure(425); set(gcf,'Name','Uxz'); hold on
%                         scatter3(DKacum, ...
%                           ones(1,length(DK))*real(DWN.omegac),...
%                           abs(Urz*ct),[],...
%                           abs(Urz*ct),'.');
%                         figure(426); set(gcf,'Name','Uyz'); hold on
%                         scatter3(DKacum, ...
%                           ones(1,length(DK))*real(DWN.omegac),...
%                           abs(Urz*st),...
%                           [],abs(Urz*st),'.');
                        
                        
                        Urz	= sum(DK.*Urz);
                        Uzz = sum(DK.*Uzz);
                        
                        UXW_fz(1,ir,ixs) = fij(ixs,3)*Urz*ct;                   %UxzW
                        UXW_fz(2,ir,ixs) = fij(ixs,3)*Urz*st;                   %UyzW
                        UXW_fz(3,ir,ixs) = fij(ixs,3)*Uzz;                      %UzzW
                    end
                    
                    if sals(ir)==1
                        Szzz = J0 .*facSzz;
                        Szrz = dJ0.*facSzr;
                        Sztz = 0*Szrz;
                        Srtz = 0*Szrz;
                        Srrz = 2*( (d2J0-(ksi^2/2-kpi^2)*J0)  .*facSrrP + d2J0   .*facSrrSV);
                        Sttz = 2*( (-(ksi^2/2-kpi^2)*J0+dJ0/r).*facSrrP + (dJ0/r).*facSrrSV);
                        
                        Szrz(isnan(Szrz))=0;
                        Sztz(isnan(Sztz))=0;
                        Szzz(isnan(Szzz))=0;
                        Srrz(isnan(Srrz))=0;
                        Sttz(isnan(Sttz))=0;
                        Srtz(isnan(Srtz))=0;
                        
%                         DKacum = zeros(1,length(DK));
%                         for idkacum = 1:length(DK)
%                         DKacum(idkacum) = sum(DK(1:idkacum));
%                         end
%                         figure(427); set(gcf,'Name','Szrz'); hold on
%                         scatter3(DKacum, ...
%                           ones(1,length(DK))*real(DWN.omegac),...
%                           abs(Szrz),...
%                           [],abs(Szrz),'.');
%                         figure(428); set(gcf,'Name','Sztz'); hold on
%                         scatter3(DKacum, ...
%                           ones(1,length(DK))*real(DWN.omegac),...
%                           abs(Sztz),[],...
%                           abs(Sztz),'.');
%                         figure(429); set(gcf,'Name','Szzz'); hold on
%                         scatter3(DKacum, ...
%                           ones(1,length(DK))*real(DWN.omegac),...
%                           abs(Szzz),...
%                           [],abs(Szzz),'.');
%                         figure(430); set(gcf,'Name','Srrz'); hold on
%                         scatter3(DKacum, ...
%                           ones(1,length(DK))*real(DWN.omegac),...
%                           abs(Srrz),...
%                           [],abs(Srrz),'.');
%                         figure(431); set(gcf,'Name','Sttz'); hold on
%                         scatter3(DKacum, ...
%                           ones(1,length(DK))*real(DWN.omegac),...
%                           abs(Sttz),[],...
%                           abs(Sttz),'.');
%                         figure(432); set(gcf,'Name','Srtz'); hold on
%                         scatter3(DKacum, ...
%                           ones(1,length(DK))*real(DWN.omegac),...
%                           abs(Srtz),...
%                           [],abs(Srtz),'.');
                        
                        Szrz = sum(DK.*Szrz);
                        Sztz = sum(DK.*Sztz);
                        Szzz = sum(DK.*Szzz);
                        Srrz = sum(DK.*Srrz);
                        Sttz = sum(DK.*Sttz);
                        Srtz = sum(DK.*Srtz);
                        
                        %Sp=[Srrz Srtz Srzz; Strz Sttz Stzz; Szrz Sztz Szzz];
                        Sp=[Srrz Srtz Szrz; Srtz Sttz Sztz; Szrz Sztz Szzz];
                        T =[ ct   -st   0 ;   st   ct   0 ;   0    0    1 ];
                        Sc= T*Sp*T.'; % a coordendas cartesianas
                        
                        SXW_fz(:,:,ir,ixs) = fij(ixs,3)*MAT(icr).Ci(6,6)*Sc;
                    end
                end
            end
        end  
%str = toc; disp(['calcu_5 :',num2str(irz),' ',num2str(str)])
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rajout du champ incident %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ir=1:nrec
    
    %Ecriture de la solution en fonction de la position du recepteur
    icr = icrall(izr0(ir)); % couche du recepteur et profondeur relative
    
    % solamente se agrega el campo incidente si el receptor y la fuente
    % pertencen al mismo estrato:
    if icr==ics
        
        ksi 	= para.reg(1).sub(ics).ksi;
        kpi 	= para.reg(1).sub(ics).kpi;
        Ci      = MAT(ics).Ci;
        
        for ixs=1:nxs%a changer pour version vectorielle
            
            xij     = xr(ir) - xs(ixs);
            yij     = yr(ir) - ys(ixs);
            zij     = zr(ir) - zs(ixs);
            rij     = sqrt(xij.^2+yij.^2+zij.^2);
            g(1,1)  = xij./rij;
            g(2,1)  = yij./rij;
            g(3,1)  = zij./rij;
            
            if salu(ir)==1
                if isfield(coordf,'vnx')
                    dr      = coordf.dr(ixs);
                    if rij<=para.npplo/2*dr %campo cercano
                        coordf0.x       = coordf.xs(ixs);
                        coordf0.y       = coordf.ys(ixs);
                        coordf0.z       = coordf.zs(ixs);
                        coordf0.vnx     = coordf.vnx(ixs);
                        coordf0.vny     = coordf.vny(ixs);
                        coordf0.vnz     = coordf.vnz(ixs);
                        if para.dim ~= 4
                        coordf0.x0      = coordf.x0(ixs);
                        coordf0.th      = coordf.th(ixs);
                        coordf0.drxz    = coordf.drxz(ixs);
                        coordf0.nbptc   = coordf.nbptc(ixs);
                        else % para.dim == 4
                        coordf0.CubPtsex = coordf.CubPtsex(:,:,ixs);
                        coordf0.CubPts   = coordf.CubPts  (:,:,ixs);
                        end
                        
                        if para.dim == 4 % 3Dgeneral
%                           if rij==0
%                             ex = true; %gaussex
%                             Gij0 = Gij_3D_r_smallGEN(coordf0,[xr(ir),yr(ir),zr(ir)],1,ksi,kpi,para,Ci,ex);
%                           else
                            ex = false;
                            Gij0 = Gij_3D_r_smallGEN(coordf0,[xr(ir),yr(ir),zr(ir)],1,ksi,kpi,para,Ci,ex);
%                           end
                        else %3D axisimétrico
                          if rij==0
                            Gij0 = Gij_3D_r_small(coordf0,xr(ir),yr(ir),zr(ir),1,ksi,kpi,gaussex,Ci);
                          else
                            Gij0 = Gij_3D_r_small(coordf0,xr(ir),yr(ir),zr(ir),1,ksi,kpi,gaussian,Ci);
                          end
                        end
                    else % suficientemente lejos
                        Gij0 = Gij_3D(ksi,kpi,rij,g,Ci,1);
                    end
                else
                    Gij0 = Gij_3D(ksi,kpi,rij,g,Ci,1);
                end
                UXW_fx(:,ir,ixs) = UXW_fx(:,ir,ixs) + fij(ixs,1)*Gij0(:,1);
                UXW_fy(:,ir,ixs) = UXW_fy(:,ir,ixs) + fij(ixs,2)*Gij0(:,2);
                UXW_fz(:,ir,ixs) = UXW_fz(:,ir,ixs) + fij(ixs,3)*Gij0(:,3);
            end
            
            if sals(ir)==1
                if isfield(coordf,'vnx')
                    dr      = coordf.dr(ixs);
                    if rij<=para.npplo/2*dr
                        coordf0.x       = coordf.xs(ixs);
                        coordf0.y       = coordf.ys(ixs);
                        coordf0.z       = coordf.zs(ixs);
                        coordf0.vnx     = coordf.vnx(ixs);
                        coordf0.vny     = coordf.vny(ixs);
                        coordf0.vnz     = coordf.vnz(ixs);
                        if para.dim ~= 4
                        coordf0.x0      = coordf.x0(ixs);
                        coordf0.th      = coordf.th(ixs);
                        coordf0.drxz    = coordf.drxz(ixs);
                        coordf0.nbptc   = coordf.nbptc(ixs);
                        else % para.dim == 4
                        coordf0.CubPtsex = coordf.CubPtsex(:,:,ixs);
                        coordf0.CubPts   = coordf.CubPts  (:,:,ixs);
                        end
                        
                        if para.dim == 4 % 3Dgeneral
%                           if rij==0
%                             ex = true; %gaussex
%                             [S_fx,S_fy,S_fz] = S_3D_r_smallGEN(coordf0,[xr(ir),yr(ir),zr(ir)],1,ksi,kpi,para,Ci,ex);
%                           else
                            ex = false;
                            [S_fx,S_fy,S_fz] = S_3D_r_smallGEN(coordf0,[xr(ir),yr(ir),zr(ir)],1,ksi,kpi,para,Ci,ex);
%                           end
                        else %3D axisimétrico
                        if rij==0
                            [S_fx,S_fy,S_fz] = S_3D_r_small(coordf0,1,xr(ir),yr(ir),zr(ir),ksi,kpi,fij(ixs,:),gaussex);
                        else
                            [S_fx,S_fy,S_fz] = S_3D_r_small(coordf0,1,xr(ir),yr(ir),zr(ir),ksi,kpi,fij(ixs,:),gaussian);
                        end
                        end
                    else
                        [S_fx,S_fy,S_fz] = S_3D(rij,g,ksi,kpi,fij(ixs,:));
                    end
                else
                    [S_fx,S_fy,S_fz] = S_3D(rij,g,ksi,kpi,fij(ixs,:));
                end
                SXW_fx(:,:,ir,ixs) = SXW_fx(:,:,ir,ixs) + fij(ixs,1)*S_fx;
                SXW_fy(:,:,ir,ixs) = SXW_fy(:,:,ir,ixs) + fij(ixs,2)*S_fy;
                SXW_fz(:,:,ir,ixs) = SXW_fz(:,:,ir,ixs) + fij(ixs,3)*S_fz;
            end
        end
    end
    
%str = toc; disp(['calcu_6 :', num2str(ir), ' ' ,num2str(str)])
end
end