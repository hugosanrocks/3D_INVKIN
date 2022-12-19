function paratmp=attenuation(para,fj)
%compute wave vector of L, T and Rayleigh waves and their attenuation
%  
paratmp = para;

for m=1:para.nmed
    if m==1 && para.geo(1)==3
        nsubmed = para.nsubmed;
    else
        nsubmed = 1;
    end
    
    for ms=1:nsubmed
        if m==1 && para.geo(1)==3
            paraM   = para.reg(1).sub(ms);
        else
            paraM   = para.reg(m);
        end
        
        C   = paraM.C;
        rho = paraM.rho;
        q   = paraM.qd;
        if paraM.tipoatts==0 || paraM.tipoatts==3
            Ci  = C;
            vL  = sqrt(C(1,1)/rho);
            vT  = sqrt(C(6,6)/rho);
        elseif paraM.tipoatts==1
            %Q
            Ci  = C/(1-1i/(2*q))^2;
            vL  = sqrt(Ci(1,1)/rho);
            vT  = sqrt(Ci(6,6)/rho);
        elseif paraM.tipoatts==2
            % Kelvin-Voigt-Model
            if length(fj)==1
                Ci  = C.*(1+1i*2*pi*fj*q);
                vL  = sqrt(Ci(1,1)/rho);
                vT  = sqrt(Ci(6,6)/rho);
            else
                Ci=zeros(6,6,length(fj));
                Ci(1,1,:)  = C(1,1).*(1+1i*2*pi*fj*q);
                Ci(2,2,:)  = C(2,2).*(1+1i*2*pi*fj*q);
                Ci(6,6,:)  = C(6,6).*(1+1i*2*pi*fj*q);
                vL  = squeeze(sqrt(Ci(1,1,:)/rho)).';
                vT  = squeeze(sqrt(Ci(6,6,:)/rho)).';
            end
        end
        
        if m==1 && para.geo(1)==3
            paratmp.reg(m).sub(ms).ksi = 2*pi*fj./vT;
            paratmp.reg(m).sub(ms).kpi = 2*pi*fj./vL;
            paratmp.reg(m).sub(ms).Ci  = Ci;
        else
            paratmp.reg(m).ksi = 2*pi*fj./vT;
            paratmp.reg(m).kpi = 2*pi*fj./vL;
            paratmp.reg(m).Ci  = Ci;    
        end
        
        if C(1,1)~=0
            x0              = Ci(6,6)/Ci(1,1);
            tmp             = roots([-1,8,8*(2*x0-3),16*(-x0+1)]);
            vR              = vT*sqrt(tmp(abs(tmp)<1));
            paratmp.reg(m).kri = 2*pi*fj/vR;
        end
    end
end