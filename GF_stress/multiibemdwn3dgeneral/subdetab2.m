function det=subdetab2(det0,ncapas,nk2)
det     = ones(6^(ncapas-1),nk2)+0i;
j       = zeros(ncapas-1,1);
% det0    = shiftdim(det0,3);
for k=1:(6^(ncapas-1))
    for i=1:ncapas-1
        kk=k;
        for ii=1:i-1
            kk=kk-(j(ii)-1)*(6^(ncapas-1-ii));
        end
        j(i)=ceil(kk/(6^(ncapas-1-i)));
    end
%     tmp=det0(1,j(1),1,:).*det0(ncapas,1,j(ncapas-1),:);
%     for i=2:(ncapas-1)
%         tmp=tmp.*det0(i,j(i),j(i-1),:);
%     end
%     det(k,:)=squeeze(tmp);
    
    tmp=det0(:,1,j(1),1).*det0(:,ncapas,1,j(ncapas-1));
    for i=2:(ncapas-1)
        tmp=tmp.*det0(:,i,j(i),j(i-1));
    end
    det(k,:)=tmp;
end