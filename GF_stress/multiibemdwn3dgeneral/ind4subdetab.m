function jj=ind4subdetab(ncapas)
%se establece de una vez por todas la lista de los indices necesarios
%al calculo de dunkin:
% det(capa_1,1,a)*det(capa_2,a,b)*det(capa_3,b,c)*...*det(capa_n-1,y,z)*det(capa_n,z,1)
%con a,b,c ..,y,z perteneciendo a 
%las 6 possibilidades por cada capa
% 1 2
% 1 3
% 1 4
% 2 3
% 2 4
% 3 4

%init lista
jj  = zeros(6^(ncapas-1),ncapas-1);
j 	= zeros(ncapas-1,1);

%encontrar indices
for k=1:(6^(ncapas-1))
    for i=1:ncapas-1
        kk=k;
        for ii=1:i-1
            kk=kk-(j(ii)-1)*(6^(ncapas-1-ii));
        end
        j(i)=ceil(kk/(6^(ncapas-1-i)));
    end
    jj(k,:)=j;
end
