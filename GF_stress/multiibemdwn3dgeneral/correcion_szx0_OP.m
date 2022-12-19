function szx0c = correcion_szx0_OP(szx0,ns,~,para)
% almacenar sólo las variables solicitadas por el usuario
[~,nrec,ninc]  = size(szx0);
szx0c           = zeros(ns,nrec,ninc);
sal             = para.sortie;

for i=1:nrec
    j=1;
    if sal.Ut==1
        if sal.sxx==1
            szx0c(j,i,:)=             szx0(1,i,:);
            j=j+1;
        end
        if sal.sxz==1
            szx0c(j,i,:)=             szx0(2,i,:);
            j=j+1;
        end
        if sal.szz==1
            szx0c(j,i,:)=             szx0(2,i,:);
%             j=j+1;
        end
    end
end
end