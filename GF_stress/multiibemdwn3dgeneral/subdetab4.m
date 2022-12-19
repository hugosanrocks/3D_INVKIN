function det=subdetab4(det0,ncapas,nk2,jj)
%calculo del determinante para todos los k2 de una vez
%det0(nk2,ncapas,6,6)
%cf routine_dev_subdet_txt
%cf Dunkin
%las 6 possibilidades por cada capa
% 1 2
% 1 3
% 1 4
% 2 3
% 2 4
% 3 4
%se ocupa recurencias para acelera el calculo
det     = zeros(1,nk2);
for k=1:(6^(ncapas-1))
    j=jj(k,:);
    tmp=det0(:,1,j(1),1).*det0(:,ncapas,1,j(ncapas-1));
    for i=2:(ncapas-1)
        tmp=tmp.*det0(:,i,j(i),j(i-1));
    end
    det=det+tmp.';
end