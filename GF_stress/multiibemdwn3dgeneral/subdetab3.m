function det=subdetab3(det0,ncapas,nk2)
%calculo del determinante para todos los k2 de una vez
%det0(nk2,ncapas,6,6)
%cf routine_dev_subdet_txt
%cf Dunkin

det     = zeros(1,nk2)+0i;
j       = zeros(ncapas-1,1);
for k=1:(6^(ncapas-1))
    %encontrar indices
    for i=1:ncapas-1
        kk=k;
        for ii=1:i-1
            kk=kk-(j(ii)-1)*(6^(ncapas-1-ii));
        end
        j(i)=ceil(kk/(6^(ncapas-1-i)));
    end
    tmp=det0(:,1,j(1),1).*det0(:,ncapas,1,j(ncapas-1));
    for i=2:(ncapas-1)
        tmp=tmp.*det0(:,i,j(i),j(i-1));
    end
    det=det+tmp.';
end