
%resout les pb lorsqu il y a une discontinuite
%par ex lorsque 2 contours sont tres proche
tmp=(abs(phi_fv(:,iinc)));
indtmp=find( ...
    ((tmp(2:(nbpt-2)))>2*(tmp(3:(nbpt-1)))) .* ...
    ((tmp(2:(nbpt-2)))>2*(tmp(1:(nbpt-3)))) ...
    )+1;
phi_fv(indtmp,iinc)=0.5*(phi_fv(indtmp+1,iinc)+phi_fv(indtmp-1,iinc));

%resout les pb lorsqu il y a un manque de discretisation qui se traduit
%par des oscillation rapide de phi_fv
%     for k=1%:5
%         tmp     = diff(real(phi_fv(:,iinc)));
%         indtmp  = find(tmp(1:(nbpt-2)).*tmp(2:(nbpt-1))<0)+1;
%         while indtmp(end)>=nbpt-2
%             indtmp=indtmp(1:(end-1));
%         end
%         indtmp  = indtmp(find(tmp(indtmp+1).*tmp(indtmp+2)<0));
%         tmp     = diff(imag(phi_fv(:,iinc)));
%         indtmp1 = find(tmp(1:(nbpt-2)).*tmp(2:(nbpt-1))<0)+1;
%         while indtmp1(end)>=nbpt-2
%             indtmp1=indtmp1(1:(end-1));
%         end
%         indtmp1 = indtmp1(find(tmp(indtmp1+1).*tmp(indtmp1+2)<0));
%         indtmp  = [indtmp; indtmp+1; indtmp1; indtmp1+1; ];
%         indtmp  = unique(indtmp);
%         if ~isempty(indtmp)
%             if indtmp(1)==1
%                 indtmp=indtmp(2:end);
%             end
%         end
%         phi_fv(indtmp,iinc)=0.25*(phi_fv(indtmp,iinc)+phi_fv(indtmp+1,iinc)+phi_fv(indtmp-1,iinc));
%     end

%     w = blackmanharris(nbptE);
%     w=circshift(w,floor(nbptE/2+1));
% 
%     sp=fft(real(phi_fv(1:nbptE,iinc)));
%     sp=sp.*w;
%     spr=real(ifft(sp));
% 
%     sp=fft(imag(phi_fv(1:nbptE,iinc)));
%     sp=sp.*w;
%     spi=real(ifft(sp));
%     phi_fv(1:nbptE,iinc)=spr+1i*spi;