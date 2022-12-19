function [k0f,indf,nmodej,det0]=solve_modos_falt_u2d(nf,j,wi,k01,DWN0,indd1,nmodej,para,ki)
%busqueda de los modos faltantes
nmodef  = nmodej-length(indd1);
k0f     = zeros(nmodef,1);
det0    = zeros(nmodef*2,1);
if nmodef>2
    disp('jamais teste')
end
if j<nf-3
    nmodejp1= sum(k01(j+1,:)~=0);
    if nmodejp1==nmodej
        %guess basado en la extrapolacion
        k0g     = interp1(wi(j+1:j+3),k01(j+1:j+3,1:nmodejp1),wi(j),'spline','extrap');
        k0c     = fliplr(DWN0.k2(indd1));
        
        %identificacion de los modos faltantes
        %a priori los modos que faltan son pares muy cercanas
        ii  = 1;
        i   = 1;
        indf= zeros(nmodef,1);
        iii = 0;
        while i<length(k0c)
            if abs(k0g(ii)-k0c(i))<abs(k0g(ii+1)-k0c(i))
                ii=ii+1;
            else
                %indice faltante i
                iii         = iii+1;
                indf(iii)   = i;
            end
            i=i+1;
        end
        indf(iii+1:nmodej-length(indd1))=ii+1:nmodejp1;
        
        for i=1:2:nmodef
            k21    = k0g(indf(i));
            k22    = k0g(indf(i+1));
            flagpbsain=0;
            kic     = k21-2*(k22-k21);
            if kic<ki
                flagpbsain=1;
            end
            kfc     = k22+2*(k22-k21);
            DWN1.omegac    = DWN0.omegac;
            DWN1.k2	= linspace(kic,kfc,1e3);
            nk      = length(DWN1.k2);
            indk    = 1:(nk-1);
            
            if  para.pol==1
                det = mode_Love(para,DWN1);
            else
                det = mode_Rayleigh_2(para,DWN1);
            end
            det(det==inf)=0;
            detf    = angle(det)-pi/2;
            indd1   = logical(detf(1:(nk-1)).*detf(2:nk)<=0);
            indd1   = indk(indd1);
            nind    = length(indd1);
            if nind==1 && flagpbsain==1
                %se debe borrar el primer modo
                nmodej  = nmodej-1;
                k0f(i)  = DWN1.k2(indd1);
                indf(i+1)= 0;
                det0(2*(i-1)+1:2*(i-1)+2)=det(indd1:indd1+1);
%                 indf(i)=
            else
                k0f(i:i+1)=DWN1.k2(indd1);
            end
        end
    end
    
else
    disp('a faire')
end
k0f(k0f==0)=[];
indf(indf==0)=[];
det0(det0==0)=[];
indf=nmodej+1-indf;