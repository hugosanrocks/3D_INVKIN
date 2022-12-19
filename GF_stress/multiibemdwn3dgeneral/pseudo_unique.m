function [zr1,izr1,zrref]=pseudo_unique(zr,para)
%function para encontrar las posiciones z en el espesor 
%las posiciones muy cercanas estan agrupadas
%zr=zr1(izr1)
zrref   = zr;
[zr0,~,izr0]  = unique(zr);%zr0=zr(izr)  zr=zr0(izr0)
izr1    = izr0;
zh      = 0;
zr1     = [];
j       = 1;
jzr1    = 1;
% if length(zr0)>3
%     for i = 1:para.nsubmed
%         hmin = para.reg(1).sub(i).h + 1e6*(i==para.nsubmed);
%         
%         dzmin  	= min(abs(para.reg(1).sub(i).bet/para.fj),hmin)/para.npplo/10;
%         if i == para.nsubmed
%             indzhm=length(zr0);
%         else
%             zh      = zh+para.reg(1).sub(i).h;
%             indzhm  = find(zr0<zh,1,'last');
%         end
%         inddzs  = find(diff(zr0(1:indzhm))<dzmin); %indice de los dz chiquito
%         while j<=indzhm
%             if sum(j==inddzs)==1
%                 tmp = mean(zr0(j:j+1));
%                 zr1                 = [zr1; tmp]; %#ok<*AGROW>
%                 zrref(zr==zr0(j))   = tmp;
%                 zrref(zr==zr0(j+1)) = tmp;
%                 izr1(izr1==jzr1+1)  = mean(izr1(izr1==jzr1));
%                 izr1(izr1>jzr1)     = izr1(izr1>jzr1)-1;
%                 jzr1= jzr1+1;
%                 j   = j+2;
%             else
%                 zr1	= [zr1; zr0(j)];
%                 j   = j+1;
%                 jzr1= jzr1+1;
%             end
%         end
%     end
% else
    zr1=zr0.';
% end