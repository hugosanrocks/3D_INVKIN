function det=sum_Dunkin_opt(nk2,ncapas,det0)
%equivalente a subdetab 1 2 3 4
%aqui se simplifica algunos terminos

if ncapas>2
    %je remonte des couches les plus eloignees au plus superficielle
    %et je stocke ds det
    prod    = zeros(nk2,6);
    prod2   = zeros(nk2,6);
    
    %2 ultimas capas
    nc      = ncapas-1;
%     for k0=[1 2 3 5 6]
%         prod(:,k0)  = ...
%             det0(:,nc,1,k0).*det0(:,ncapas,1,1)+det0(:,nc,2,k0).*det0(:,ncapas,1,2)+(det0(:,nc,3,k0)+det0(:,nc,4,k0)).*det0(:,ncapas,1,3)+...
%             det0(:,nc,5,k0).*det0(:,ncapas,1,5)+det0(:,nc,6,k0).*det0(:,ncapas,1,6);
%     end
    prod(:,1)  = ...
        det0(:,nc,1,1).*det0(:,ncapas,1,1)+det0(:,nc,2,1).*det0(:,ncapas,1,2)+2*det0(:,nc,3,1).*det0(:,ncapas,1,3)+...
        det0(:,nc,5,1).*det0(:,ncapas,1,5)+det0(:,nc,6,1).*det0(:,ncapas,1,6);
    
    prod(:,2)  = ...
        det0(:,nc,1,2).*det0(:,ncapas,1,1)+det0(:,nc,2,2).*det0(:,ncapas,1,2)+2*det0(:,nc,3,2).*det0(:,ncapas,1,3)+...
        det0(:,nc,5,2).*det0(:,ncapas,1,5)+det0(:,nc,6,2).*det0(:,ncapas,1,6);
    
    prod(:,3)  = ...
        det0(:,nc,1,3).*det0(:,ncapas,1,1)+det0(:,nc,2,3).*det0(:,ncapas,1,2)+(det0(:,nc,3,3)+det0(:,nc,4,3)).*det0(:,ncapas,1,3)+...
        det0(:,nc,5,3).*det0(:,ncapas,1,5)+det0(:,nc,6,3).*det0(:,ncapas,1,6);
    
    prod(:,5)  = ...
        det0(:,nc,1,5).*det0(:,ncapas,1,1)+det0(:,nc,2,5).*det0(:,ncapas,1,2)+2*det0(:,nc,3,5).*det0(:,ncapas,1,3)+...
        det0(:,nc,5,5).*det0(:,ncapas,1,5)+det0(:,nc,6,5).*det0(:,ncapas,1,6);
    
    prod(:,6)  = ...
        det0(:,nc,1,6).*det0(:,ncapas,1,1)+det0(:,nc,2,6).*det0(:,ncapas,1,2)+2*det0(:,nc,3,6).*det0(:,ncapas,1,3)+...
        det0(:,nc,5,6).*det0(:,ncapas,1,5)+det0(:,nc,6,6).*det0(:,ncapas,1,6);

    prod(:,4)=prod(:,3);
    
    for nc=ncapas-2:-1:2
%         for k0=[1 2 3 5 6]
%             prod2(:,k0)  =det0(:,nc,1,k0).*prod(:,1)+det0(:,nc,2,k0).*prod(:,2)+(det0(:,nc,3,k0)+det0(:,nc,4,k0)).*prod(:,3)+det0(:,nc,5,k0).*prod(:,5)+det0(:,nc,6,k0).*prod(:,6);
%         end
%         prod=prod2;
        
        prod2(:,1)  = det0(:,nc,1,1).*prod(:,1)+det0(:,nc,2,1).*prod(:,2)+2*det0(:,nc,3,1).*prod(:,3)+det0(:,nc,5,1).*prod(:,5)+det0(:,nc,6,1).*prod(:,6);
        prod2(:,2)  = det0(:,nc,1,2).*prod(:,1)+det0(:,nc,2,2).*prod(:,2)+2*det0(:,nc,3,2).*prod(:,3)+det0(:,nc,5,2).*prod(:,5)+det0(:,nc,6,2).*prod(:,6);
        prod2(:,3)  = det0(:,nc,1,3).*prod(:,1)+det0(:,nc,2,3).*prod(:,2)+(det0(:,nc,3,3)+det0(:,nc,4,3)).*prod(:,3)+det0(:,nc,5,3).*prod(:,5)+det0(:,nc,6,3).*prod(:,6);
        prod2(:,5)  = det0(:,nc,1,5).*prod(:,1)+det0(:,nc,2,5).*prod(:,2)+2*det0(:,nc,3,5).*prod(:,3)+det0(:,nc,5,5).*prod(:,5)+det0(:,nc,6,5).*prod(:,6);
        prod2(:,6)  = det0(:,nc,1,6).*prod(:,1)+det0(:,nc,2,6).*prod(:,2)+2*det0(:,nc,3,6).*prod(:,3)+det0(:,nc,5,6).*prod(:,5)+det0(:,nc,6,6).*prod(:,6);
        prod        = prod2;
        prod(:,4)   = prod(:,3);
    end
    
    nc  =1;
    %     det =0;
    %     for k=1:6
    %         det  = det+det0(:,nc,k,1).*prod(:,k);
    %     end
    det  = det0(:,nc,1,1).*prod(:,1)+ det0(:,nc,2,1).*prod(:,2)+ 2*det0(:,nc,3,1).*prod(:,3)+ det0(:,nc,5,1).*prod(:,5)+ det0(:,nc,6,1).*prod(:,6);
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