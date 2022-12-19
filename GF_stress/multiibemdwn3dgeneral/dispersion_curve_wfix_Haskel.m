function dispersion_curve_wfix_Haskel(para)

% funcion que busqua a establecer las curvas de dispersion por un medio
% estratificado sobre un semi espacio
% se busca por cada w(wfix) 
% los ceros de la fase de un sub-determinante de la matriz de haskel 
% en un gran intervalo busquando un cambio de signo
% metodo mas rapido que con la matriz global
% no hay falso ceros con las asintotas intermedias y eso permite no perder
% ceros

tic
pol     = para.pol;
para.sub= para.reg(1).sub;
beta    = zeros(para.nsubmed,1);
for ms=1:para.nsubmed
    para.sub(ms).C(6,6)  = para.sub(ms).rho*para.sub(ms).bet^2;
    beta(ms)  = para.sub(ms).bet;
end
vmin        = min(beta)*(1+1e-6);
[vmax,imax] = max(beta*(1-1e-6));
vdeep       = beta(end)*(1-1e-6);

wmax    = 2*pi*para.fmax;



%numero de frecuencias a buscar los k2 soluciones
nf      = 100; 
dw      = 2*pi*para.fmax/nf;

if vmax~=vdeep
    n0=para.nsubmed;
    para.nsubmed=imax;
    [k2c,w0c]= dispersion_curve_k_critik_Haskel(para,wmax,vmax);
    nmode  	= length(k2c);
    para.nsubmed=n0;
    k2sol2  	= zeros(nf,nmode);
else
    [k2c,w0c]= dispersion_curve_k_critik_Haskel(para,wmax,vdeep);
    nmode  	= length(k2c);
end
k2sol  	= zeros(nf,nmode);

ikmax	= zeros(1,nmode);

nk    	= 1e2*nmode;
indk 	= 1:(nk-1);

for j=1:nf
    wj          = j*dw + 0.0001*dw*(j==0);
    DWN         = struct('omegac',wj);
    DWN.k2    	= linspace(wj/vdeep,wj/vmin*1.2,nk);
   
    if  pol==1
        tmp=mode_Love(para,DWN);
    else
        tmp=mode_Rayleigh_2(para,DWN);
    end
    tmp(tmp==inf)=0;
    tmp1    = angle(tmp)-1e-6;
    indd1   = logical(tmp1(1:(nk-1)).*tmp1(2:nk)<0);
    indd1   = indk(indd1);
    if ~isempty(indd1)
        k20   	= cherche_zero(DWN.k2,real(tmp).',indd1).';
        nind    = length(k20);
        if nind>nmode
            k20=k20(1:nmode);
        end
        k2sol(j,:)= [flipud(k20);zeros(nmode-nind,1)];
    end
    
    if vdeep~=vmax
        %attention ces modes peuvent ne pas se propager longtemps car si il
        %atteignent l interface de la couche inférieur, ils s'attenuent
        %rapidement
        DWN.k2    	= linspace(wj/vmax,wj/(vdeep*(1+2e-6)),nk);
        para.nsubmed=imax;
        if  pol==1
            tmp=mode_Love(para,DWN);
        else
            tmp=mode_Rayleigh_2(para,DWN);
        end
        tmp(tmp==inf)=0;
        tmp1    = angle(tmp)-1e-6;
        indd1   = logical(tmp1(1:(nk-1)).*tmp1(2:nk)<=0);
        indd1   = indk(indd1);
        if ~isempty(indd1)
            k20   	= cherche_zero(DWN.k2,real(tmp).',indd1).';
            nind    = length(k20);
            if nind>nmode
                k20=k20(1:nmode);
            end
            k2sol2(j,:)= [flipud(k20);zeros(nmode-nind,1)];
        end
        para.nsubmed=n0;
    end
end
toc

j   = 1:nf;
wj	= j*dw + 0.0001*dw*(j==0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica w=f(k)                            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(205);hold on;
for j=1:nmode
    indi    = find(k2sol(:,j)~=0,1,'first');
    indf    = find(k2sol(:,j)~=0,1,'last');
    
    if ~isempty(indf)
        if vdeep~=vmax
            %partes entrecortadas
            indi1    = find(k2sol(indi:indf,j)==0,1,'first')-1;
            indf1    = find(k2sol(indi:indf,j)==0,1,'last')+1;
            ikmax(j)= indf-indi+1;
            if isempty(indi1)
                plot(k2sol(indi:indf,j),wj(indi:indf),'')
            else
                plot([k2c(j);k2sol(indi:indi1,j)],[w0c(j);wj(indi:indi1).'],'')
                plot(k2sol(indf1:indf,j),wj(indf1:indf),'')
            end
        else
            ikmax(j)= indf-indi+1;
            plot([k2c(j);k2sol(indi:indf,j)],[w0c(j);wj(indi:indf).'],'.')
        end
    end
end
k2sol3=k2sol2;
k2sol2=k2sol3;
for j=1:nmode
    plot(k2sol2(:,j),wj.','.r')
end
if vdeep~=vmax
    for j=1:nmode
        indi    = find(k2sol2(:,j)~=0,1,'first');
        indf    = find(diff(k2sol2(:,j))<0,1,'first')+1;
        if indf>indi+3
            if ~isempty(indi)
                for jj=nmode:-1:j+1
                    k2sol2(indi:end,jj)=k2sol2(indi:end,jj-1);
                    k2sol2(indi:end,jj-1)=0;
                end
            end
        else
            k2sol2(1:indf,j)=0;
            indf    = find(diff(k2sol2(:,j))<0,1,'first')+1;
            if ~isempty(indf)
                for jj=nmode:-1:j+1
                    k2sol2(indf:end,jj)=k2sol2(indf:end,jj-1);
                    k2sol2(indf:end,jj-1)=0;
                end
            end
        end
    end
    
    for j=1:nmode
        indi    = find(k2sol2(:,j)~=0,1,'first');
        indf    = find(k2sol2(:,j)~=0,1,'last');
        if ~isempty(indf)
            plot(k2sol2(indi:indf,j),wj(indi:indf),'r')
        end
    end
end
xlabel('k_x');
ylabel('\omega');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               grafica asintotas                         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kmax= wmax/min(beta);
k2	= linspace(0,kmax,1e3);
for ic=1:para.nsubmed
    C   = para.sub(ic).C;
    rho = para.sub(ic).rho;
    vT  = sqrt(C(6,6)/rho);
    if  para.pol==1
        w=k2*vT;
        i0=find(w>wmax,1);
        if isempty(i0)
            i0=length(k2);
        end
        plot(k2(1:i0),w(1:i0),'k');
    else
        vL  = sqrt(C(1,1)/rho);
        w   = k2*vT;
        i0  = find(w>wmax,1);
        if isempty(i0)
            i0=length(k2);
        end
        plot(k2(1:i0),w(1:i0),'k');
        w   = k2*vL;
        i0  = find(w>wmax,1);
        if isempty(i0)
            i0=length(k2);
        end
        plot(k2(1:i0),w(1:i0),'r')
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vp=w/k =f(w)                        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure(206);hold on;

for imode=1:nmode
    indi    = find(k2sol(:,imode)~=0,1,'first');
    indf    = find(k2sol(:,imode)~=0,1,'last');
    if vdeep~=vmax
        indi1    = find(k2sol(indi:indf,imode)==0,1,'first')-1;
        indf1    = find(k2sol(indi:indf,imode)==0,1,'last')+1;
        if isempty(indi1)
            vp      = wj(indi:indf)./k2sol(indi:indf,imode).';
            plot(wj(indi:indf)/2/pi,vp,'k')
        else
            vp = wj(indi:indi1)./k2sol(indi:indi1,imode).';
            plot(wj(indi:indi1)/2/pi,vp,'k')
            vp = wj(indf1:indf)./k2sol(indf1:indf,imode).';
            plot(wj(indf1:indf)/2/pi,vp,'k')
        end
    else
        vp      = wj(indi:indf)./k2sol(indi:indf,imode).';
        plot(wj(indi:indf)/2/pi,vp,'k')
    end
end

if vdeep~=vmax
    for j=1:nmode
        indi    = find(k2sol2(:,j)~=0,1,'first');
        indf    = find(k2sol2(:,j)~=0,1,'last');
        if ~isempty(indi)
            vp      = wj(indi:indf)./k2sol2(indi:indf,j).';
            plot(wj(indi:indf)/2/pi,vp,'r')
        else
            break
        end
    end
end

xlabel('frecuencia');
ylabel('Vp');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    cambio a grafica vg=dw/dk =f(w)                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vg=zeros(nmode,max(ikmax)-1);
f1=zeros(nmode,max(ikmax)-1);

for imode=1:nmode
    indi    = find(k2sol(:,imode)~=0,1,'first');
    indf    = find(k2sol(:,imode)~=0,1,'last');
    if vdeep~=vmax
        indi1    = find(k2sol(indi:indf,imode)==0,1,'first')-1;
        indf1    = find(k2sol(indi:indf,imode)==0,1,'last')+1;
        if isempty(indi1)
            vg(imode,1:indf-indi)=diff(wj(indi:indf))./diff(k2sol(indi:indf,imode)).';
            f1(imode,1:indf-indi)=(wj(indi:(indf-1))+wj((indi+1):indf))/2/2/pi;
            ikmax(imode)=indf-indi;
        else
            vg(imode,1:indi1-indi)=diff(wj(indi:indi1))./diff(k2sol(indi:indi1,imode)).';
            f1(imode,1:indi1-indi)=(wj(indi:(indi1-1))+wj((indi+1):indi1))/2/2/pi;
            n0=indi1-indi;
            vg(imode,n0+(1:indf-indf1))=diff(wj(indf1:indf))./diff(k2sol(indf1:indf,imode)).';
            f1(imode,n0+(1:indf-indf1))=(wj(indf1:(indf-1))+wj((indf1+1):indf))/2/2/pi;
            ikmax(imode)=n0+indf-indf1;
        end
    else
        vg(imode,1:indf-indi)=diff(wj(indi:indf))./diff(k2sol(indi:indf,imode)).';
        f1(imode,1:indf-indi)=(wj(indi:(indf-1))+wj((indi+1):indf))/2/2/pi;
        ikmax(imode)=indf-indi;
    end
end
figure(207);hold on;
for imode=1:nmode
    plot( f1(imode,1:ikmax(imode)),vg(imode,1:ikmax(imode)),'k')
end
xlabel('frecuencia');
ylabel('Vg');

