function checktime_disp(para,wi)

% funcion que busqua a establecer las curvas de dispersion de velocidad de
% fase
% por un medio estratificado sobre un semi espacio
% se busca por cada w(wfix)
% los ceros de la fase de un sub-determinante de la matriz de haskel
% en un gran intervalo busquando un cambio de signo
% metodo mas rapido que con la matriz global
% no hay falso ceros
%
% diferencia con dispersion_curve_wfix_Haskel_4inv:
% no mas las velocidad de fase

para.sub= para.reg(1).sub;
pol     = para.pol;

beta    = zeros(para.nsubmed,1);
alpha   = zeros(para.nsubmed,1);
for ms=1:para.nsubmed
    para.sub(ms).C(6,6)  = para.sub(ms).rho*para.sub(ms).bet^2;
    beta(ms)  = para.sub(ms).bet;
    alpha(ms) = para.sub(ms).alpha;
end

vmin    = min(beta);
vmax    = max(beta);

vmax    = (1-1e-8)*vmax;

wi      = unique(wi);
wmax    = max(wi);
nf      = length(wi);
%&precalculo lista indices
if para.pol==2
    jj=ind4subdetab(para.nsubmed);
    para.jj=jj;
end
[k2c,w0c]   = dispersion_curve_k_critik_Haskel(para,wmax,vmax);
nmode       = length(k2c);
k01         = zeros(nf,nmode);
for j=50%1:nf
    % parfor j=1:nf
    wj          = wi(j);
    if wj==0
        k20     = 0;
        %automaticamente hecho porque k20=0
        % nmodej  = 1;
        % k01(j,:)= [flipud(k20);zeros(nmode-nmodej,1)].';
        continue;
    end
    nmodej      = sum(w0c<=wj);
    DWN0        = struct('omegac',wj);
    
    %limites en k dados por asintotas
    ki          = wj/vmax;
    kf          = wj/vmin*1.3;
    
    j0=0;
    for nk1         = 1:1:20
        nk1
        j0=j0+1;
        
        tic
        
        for lp=1:50%average
            nk    	= nmodej*nk1;%min(floor(max(dk0/dwi*1e3,nmodej*3)),5e3);
            DWN0.k2 = linspace(ki,kf,nk);
            nk      = length(DWN0.k2);
            indk    = 1:(nk-1);
            
            if  pol==1
                tmp=mode_Love(para,DWN0);
            else
                tmp=mode_Rayleigh_2(para,DWN0);
            end
            tmp(tmp==inf)=0;
            tmp1    = angle(tmp)-pi/2;
            indd1   = logical(tmp1(1:(nk-1)).*tmp1(2:nk)<=0);
            indd1   = indk(indd1);
            nind    = length(indd1);
        end
        elapsedTime(j0,1)= toc;
        elapsedTime(j0,2)= nk;
    end
    figure;plot( elapsedTime(:,2), elapsedTime(:,1))
end
