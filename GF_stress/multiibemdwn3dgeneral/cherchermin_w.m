function k20=cherchermin_w(y,para,DWN,np)
%para buscar los puntos cerca de zero que no cambian de signo
%primero se busca los cambios de signo de la derivada que sean - -> +
%despues se refina la buscqueda alrededor de estos puntos y se verifica que
%la extrapolacion de la parte antes de la curva cambie de signo despues del
%cambio de signo

%busqueda de cambio de signo de la derivada
dy      = diff(y);
n       = length(dy);
indd1   = find((dy(1:n-1).*dy(2:n)<=0));
%selection de los que sean - -> +
indd1   = indd1(dy(indd1)<0);
n0      = length(indd1);

DWN0    = DWN;
dk      = DWN.k2(2)-DWN.k2(1);

if n0==0
    %cas ou le changement de sign initial est sur les bords de l'intervale
    %et que soit lamplitude au point final apres le changement de signe est
    %plus petit en valeur absolue que l' anterieur,
    %soit que le second est plus petit en val. abs. que le premier
    n0      =2;
    indd1   =[1 n];
    DWN0.k2(1)  =DWN0.k2(1)-dk;
    DWN0.k2(end)=DWN0.k2(end)+dk;
end
k20=zeros(n0,1);

%se verifica que son verdamente unos ceros
for i=n0:-1:1
    ind=indd1(i)+1;
    k2i      	= DWN0.k2(max(ind-1,1));
    k2f      	= DWN0.k2(min(ind+1,length(DWN0.k2)));
    DWN.k2   	= linspace(k2i,k2f,np);
    paratmp     = para;
    
    if  para.pol==1
        DWN     = calcul_A_DWN_SH_Ncapas_HS(paratmp,DWN);
    else
        DWN     = calcul_A_DWN_PSV_Ncapas_HS(paratmp,DWN);
    end
    
    tmp     = zeros(np,1);
    for k=1:np
        tmp(k)=abs(det(DWN.A_DWN(:,:,k)));
    end
    
    %     tmp=abs(determinant_vec_de_mat(DWN.A_DWN));
    
    %     %check if 2 min are present in the same interval
    %     dy  =diff(y);
    %     n   =length(dy);
    %     indd =find((dy(1:n-1).*dy(2:n)<=0));
    %     indd=indd(dy(indd)<0);
    
    [~,ind1]=min(tmp);%figure;plot(DWN.k2,tmp)
    
    ind1=ind1-1;%a veces el min esta sobre la pendiente positiva
    if ind1==0
        ind1=1;%a dichotomier // presence de 2 min
        nextw0=1;
    else
        if ind1==1
            nextw0 = interp1(DWN.k2(3:6),tmp(3:6),DWN.k2(3)-5*dk/np,'pchip','extrap');
        else
            nextw0 = interp1(DWN.k2(1:ind1),tmp(1:ind1),DWN.k2(ind1)+6*dk/np,'pchip','extrap');
        end
    end
    if sign(nextw0)==-1
        k20(i)=DWN.k2(ind1);
    else
        k20(i)=[];
    end
end