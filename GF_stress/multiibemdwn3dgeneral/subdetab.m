function det=subdetab(det0,j,det,ic,ncapas)

if ic~=ncapas
    for j1=1:6
        j(ic)=j1;
        det=subdetab(det0,j,det,ic+1,ncapas);
    end
else
    % ic==ncapas
    k=1;
    for i=1:(ncapas-1)
        k=k+6^(i-1)*(j(i)-1);
    end
    for i=1:ncapas
        if i==1
            det(k,:)=det(k,:).*squeeze(det0(:,i,j(i),1)).';
        elseif i==ncapas
            det(k,:)=det(k,:).*squeeze(det0(:,i,1,j(i-1))).';
        else
            det(k,:)=det(k,:).*squeeze(det0(:,i,j(i),j(i-1))).';
        end
    end
end