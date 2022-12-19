function detA=zerodeAk(para,DWN,k2)
%fonction utilisee dans la recherche de zero du determinant de la matrice
%du DWN 

DWN.k2  = k2;
if  para.pol==1
    DWN     = calcul_A_DWN_SH_Ncapas_HS(para,DWN);
else
    DWN     = calcul_A_DWN_PSV_Ncapas_HS(para,DWN);
end

nk2=length(DWN.k2);
detA=zeros(nk2,1);
for i=1:nk2;
    detA(i)=real(det(DWN.A_DWN(:,:,i)));
end