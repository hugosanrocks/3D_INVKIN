function uxz0c=correcion_uxz0_OP(uxz0,ns,polOP,para)
%en el caso de las ondas planas, la descomposicion del campo incidente es
%trivial:

[~,nrec,ninc]  = size(uxz0);
uxz0c           = zeros(ns,nrec,ninc);
sal             = para.sortie;

for i=1:nrec
    j=1;
    if sal.Ut==1
        if sal.Ux==1
            uxz0c(j,i,:)=             uxz0(1,i,:);
            j=j+1;
        end
        if sal.Uz==1
            uxz0c(j,i,:)=             uxz0(2,i,:);
            j=j+1;
        end
    end
    if sal.UPh==1
        if sal.Ux==1
            uxz0c(j,i,:)= (polOP==1).'.*squeeze(uxz0(1,i,:));
            j=j+1;
        end
        if sal.Uz==1
            uxz0c(j,i,:)= (polOP==1).'.*squeeze(uxz0(2,i,:));
            j=j+1;
        end
    end
    if sal.USh==1
        if sal.Ux==1
            uxz0c(j,i,:)= (polOP==2).'.*squeeze(uxz0(1,i,:));
            j=j+1;
        end
        if sal.Uz==1
            uxz0c(j,i,:)= (polOP==2).'.*squeeze(uxz0(2,i,:));
            j=j+1;
        end
    end
    if sal.UIh==1
        if sal.Ux==1
            uxz0c(j,i,:)= (polOP==3).'.*squeeze(uxz0(1,i,:));
            j=j+1;
        end
        if sal.Uz==1
            uxz0c(j,i,:)= (polOP==3).'.*squeeze(uxz0(2,i,:));
        end
    end
    if sal.UPt==1
        if sal.Ux==1
            uxz0c(j,i,:)= (polOP==1).'.*squeeze(uxz0(1,i,:));
            j=j+1;
        end
        if sal.Uz==1
            uxz0c(j,i,:)= (polOP==1).'.*squeeze(uxz0(2,i,:));
            j=j+1;
        end
    end
    if sal.USt==1
        if sal.Ux==1
            uxz0c(j,i,:)= (polOP==2).'.*squeeze(uxz0(1,i,:));
            j=j+1;
        end
        if sal.Uz==1
            uxz0c(j,i,:)= (polOP==2).'.*squeeze(uxz0(2,i,:));
        end
    end
end