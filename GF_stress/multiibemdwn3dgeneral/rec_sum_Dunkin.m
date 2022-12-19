function det=rec_sum_Dunkin(nc,k0,ncapas,det0)
%equivalente a subdetab 1 2 3 4
%recursivo para mejorar tiempos de calculo
% if nc==ncapas
%     det = det0(:,ncapas,1,k0);
%marche pas on repasse toujours ds les meme boucles
if nc==ncapas-1
    det=0;
    for k=1:6
        det  = det+det0(:,nc,k,k0).*det0(:,ncapas,1,k);
    end
elseif nc==1
    det=0;
    for k=1:6
        prod = rec_sum_Dunkin(nc+1,k,ncapas,det0);
        det  = det+det0(:,nc,k,1).*prod;
    end
else
    det=0;
    for k=1:6
        prod = rec_sum_Dunkin(nc+1,k,ncapas,det0);
        det  = det+det0(:,nc,k,k0).*prod;
    end
end