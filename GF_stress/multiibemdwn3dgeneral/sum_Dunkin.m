function det=sum_Dunkin(nk2,ncapas,det0)
%equivalente a subdetab 1 2 3 4

if ncapas>2
    %je remonte des couches les plus eloignees au plus superficielle
    %et je stocke ds det
    prod    = zeros(nk2,6);
    prod2   = zeros(nk2,6);
    
    %2 ultimas capas
    nc      = ncapas-1;
    for k0=1:6
        det=0;
        for k=1:6
            det  = det+det0(:,nc,k,k0).*det0(:,ncapas,1,k);
        end
        %stock k0
        prod(:,k0)=det;
    end
    for nc=ncapas-2:-1:2
        for k0=1:6
            det=0;
            for k=1:6
                det  = det+det0(:,nc,k,k0).*prod(:,k);
            end
            prod2(:,k0)=det;
        end
        prod=prod2;
    end
    
    nc  =1;
    det =0;
    for k=1:6
        det  = det+det0(:,nc,k,1).*prod(:,k);
    end
else
    %2 ultimas capas
    nc      = ncapas-1;
    det=0;
    for k0=1:6
        for k=1:6
            det  = det+det0(:,nc,k,k0).*det0(:,ncapas,1,k);
        end
    end
end