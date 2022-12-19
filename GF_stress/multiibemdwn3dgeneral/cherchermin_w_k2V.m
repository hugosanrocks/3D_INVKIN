function w0=cherchermin_w_k2V(y,para,DWN,v,np)
%para buscar los puntos cerca de zero que no cambian de signo
%primero se busca los cambios de signo de la derivada que sean - -> +
%despues se refina la buscqueda alrededor de estos puntos y se verifica que
%la extrapolacion de la parte antes de la curva cambie de signo
dy  =diff(y);
n   =length(dy);
indd =find((dy(1:n-1).*dy(2:n)<=0));
indd=indd(dy(indd)<0);
DWN0=DWN;
n0=length(indd);
w0=zeros(n0,1);
dw=DWN.omegac(2)-DWN.omegac(1);
for i=n0:-1:1
    ind=indd(i)+1;
    k2i         = DWN0.k2(max(ind-1,1));
    k2f         = DWN0.k2(min(ind+1,length(DWN0.omegac)));
    DWN.k2      = linspace(k2i,k2f,np);%deberia mejorarse con buscqueda bitcorr
    fj          = DWN.k2*v/2/pi;
    DWN.omegac	= 2*pi*fj;

    for ic=1:para.nsubmed
        para.reg(1).sub(ic).ksi  = DWN.omegac/para.reg(1).sub(ic).bet;
    end
    if  para.pol==1
        DWN     = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
    else
        for ic=1:para.nsubmed
            para.reg(1).sub(ic).kpi  = DWN.omegac/para.reg(1).sub(ic).alpha;
        end
        DWN     = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
    end
    tmp     = zeros(np,1);
    for k=1:np
        tmp(k)=abs(det(DWN.A_DWN(:,:,k)));
    end
    %check if 2 min are present in the same interval
    dy  =diff(y);
    n   =length(dy);
    indd =find((dy(1:n-1).*dy(2:n)<=0));
    indd=indd(dy(indd)<0);
    
    [~,ind1]=min(tmp);%figure;plot(tmp)
    
    ind1=ind1-1;%a veces el min esta sobre la pendiente positiva
    if ind1==0
        ind1=1;%a dichotomier // presence de 2 min
        nextw0=1;
    elseif ind1==1
        ind1=2;
        nextw0 = interp1(DWN.omegac(1:ind1),tmp(1:ind1),DWN.omegac(ind1)+3*dw,'pchip','extrap');
    else
        nextw0 = interp1(DWN.omegac(1:ind1),tmp(1:ind1),DWN.omegac(ind1)+3*dw,'pchip','extrap');
    end
    if sign(nextw0)==-1
        w0(i)=DWN.omegac(ind1);
    else
        w0(i)=[];
    end
    
end