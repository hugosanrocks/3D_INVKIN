function solw=cherche_zero_dic(para,DWN1,det)
phase	= angle(det)-pi/2;
indd1   = find(phase(1:end-1).*phase(2:(end))<=0);
DWN     = DWN1;
solw    = indd1;
for i=1:length(indd1)
    indd0   = indd1(i);
    wi      = DWN.omegac(indd0);
    wf      = DWN.omegac(indd0+1);
    
    while (wf-wi)/3>1e-5
        [~,~,~,det,~,DWN1]  = checkchgmtsign_H_vec(wi,wf,para,DWN1,1e2);
        phase               = angle(det)-pi/2;
        indd0               = find(phase(1:end-1).*phase(2:(end))<=0);
        if length(indd0)>1
            solw=[];
            return;
        end
        wi                  = DWN1.omegac(indd0);
        wf                  = DWN1.omegac(indd0+1);
    end
    %pour etre sur que indd0 n est pas 1 ni 5
    [~,~,~,det,~,DWN1]  = checkchgmtsign_H_vec(wi,wf,para,DWN1,5);
    phase               = angle(det)-pi/2;
    if max(phase==0)
        solw(i)     = DWN1.omegac(phase==0);
    else
        indd0               = find(phase(1:end-1).*phase(2:(end))<=0);
        if indd0==1
            %interp lin
            y   = real(det(1:2));
            x   = DWN1.omegac(1:2);
            a   = (y(2)-y(1))/(x(2)-x(1));
            b   = y(2)-a*x(2);
            solw(i)             = -b/a;
        elseif indd0==5
            %interp lin
            y   = real(det(4:5));
            x   = DWN1.omegac(4:5);
            a   = (y(2)-y(1))/(x(2)-x(1));
            b   = y(2)-a*x(2);
            solw(i)             = -b/a;
        else
            %interp quad
            solw(i)             = cherche_zero(DWN1.omegac,real(det).',indd0);
        end
    end
end