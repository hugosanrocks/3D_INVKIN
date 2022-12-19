function [newphi,newcoord]=resample_phi_fv(coord,phi_fv,para,coord0)
%1) reorganizacion de los phi segun r
%2) interpolacion fina de phi
%3) remuestro segun fuerzas homogeneas

pplot=1;
if pplot==1
    figure(154);
    hold on;
    map=colormap;
end
newcoord= struct( ...
    'x'     ,[],	'z'     ,[], ...
    'vnx'   ,[],	'vnz'   ,[], ...
    'dr'    ,[],	'phi'  	,[]);

newphi  =[];

%los medios que hay que remuestrear (solo los de los receptores)
im       = unique(para.rec.medio).';
im(im==0)= [];
nnew	 = zeros(max(im),1);

for m=im
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % organizar phi segun r creciente %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    indtri  = coord.m(m).ind;
    xtmp	= coord.x(indtri);
    ztmp	= coord.z(indtri);
    drj   	= coord.dr(indtri);
    drref   = mean(drj);
    clear r
    r       = cumsum(drj)-drj(1)/2;
    n    	= length(r);
    phi1	= phi_fv(coord.phi(indtri,m)).';
    if pplot==1
        figure(90);hold on;
        for ikm=1:1:n
            plot(xtmp(ikm),ztmp(ikm),'marker','.','markersize',30,'MarkerEdgeColor',map(round(abs(phi1(ikm))/max(abs(phi1))*63+1),:))
        end
        plot(xtmp(1),ztmp(1),'*r');plot(xtmp(2),ztmp(2),'ok')
        figure(154);plot(r,abs(phi1),'r.');
    end
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % interpolacion fina de los phi %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     long            = r(n)+drj(n)/2;
%     nf              = 100*n;
%     rfin            = linspace(0,long,nf);
%     if m>1
%         %para cerrar el contorno
%         r(n+1)      = rfin(end);
%         phi1(n+1)   = phi1(1);
%     end
%     phifin          = interp1(r,phi1,rfin,'spline');
% %     plot(rfin,abs(phifin),'.','markersize',2);
%     i               = 2:(nf-1);
%     drv             = zeros(1,nf);
%     drv(i)       	= (rfin(i+1)-rfin(i-1))/2;
%     drv(1)        	= (rfin(2)-rfin(1))/2;
%     drv(nf)         = (rfin(nf)-rfin(nf-1))/2;
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % interpolacion fina de la malla %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%&%
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % construction de fuerzas iguales %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %fuerza promedia
%     df  = mean(abs(phi1(1:n)).*drj(1:n))/4;
%     j   = 0;
%     r2  = zeros(1,10*n);
%     phi2= zeros(10*n);
%     i   = 0;
%     while i <10*n && j+1<nf
%         i   = i+1;
%         j   = j+1;
%         j0  = j;
%         dfc = cumsum(abs(phifin(j0)).*drv(j0));%fuerza del tramito
%         drvc= drv(j0);
%         while dfc<df && j<nf && drvc<drref
%             j   = j+1;
%             dfc = sum(abs(phifin(j0:j)).*drv(j0:j));
%             drvc= sum(drv(j0:j));
%         end
%         r2(i)   = (rfin(j0)+rfin(j))/2;
%         phi2(i) = sum(phifin(j0:j).*drv(j0:j))/sum(drv(j0:j));
%     end
%     n               =i;
%     r2(n+1:end)     =[];
%     phi2(n+1:end)   =[];
%     newphi          =[newphi,phi2];
% %     plot(r2,abs(phi2),'*k');
%     
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     % expression en la malla inicial  %
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     
%     %no se conoce las coordenadas, normales en las posiciones r2 por lo
%     %cual se busca puntos en la malla coriente y se amplifica esta malla
%     %con interpolacion
%     lr      = r2*0;
%     indxz   = coord0.m(m).ind;
%     r0      = coord0.m(m).r;
%     i       = 1;
%     [~,imin]= min(abs(r0-r2(i)));
%     lr(i)   = imin;
%     i       = i+1;
%     j       = imin;
%     while i<n+1
%         [~,imin]= min(abs(r2(i)-r0(j:end)));
%         imin    = j-1+imin;
%         lr(i)   = imin;
%         i       = i+1;
%         j       = imin;
%     end
%     lxz     = indxz(lr);
%     nnew(m)	= n;
%     
%     xtmp            = coord0.x(lxz);
%     ztmp            = coord0.z(lxz);
%     r               = r0(lr);
%     i               = 2:n-1;
%     dr              = zeros(n,1);
%     dr(2:n-1)       = (r(i+1)-r(i-1))/2;
%     dr(1)           = (r(2)+r(1))/2;
%     dr(n)           = r0(end)-(r(n)+r(n-1))/2;
%     
%     newcoord.x      = [newcoord.x  ,xtmp];
%     newcoord.z      = [newcoord.z  ,ztmp];
%     newcoord.vnx    = [newcoord.vnx,coord0.vnx(lxz)];
%     newcoord.vnz    = [newcoord.vnz,coord0.vnz(lxz)];
%     newcoord.dr     = [newcoord.dr ,dr.'];
end
% newcoord.nbeq   = sum(nnew);
% newcoord.nbpt   = newcoord.nbeq;
% newcoord.phi	= zeros(sum(nnew),max(im));
% i=0;
% for m=im
%     newcoord.phi(i+(1:nnew(m)),m)    = i+(1:nnew(m));
%     newcoord.indm(m).ind         = false(newcoord.nbeq,1);
%     newcoord.indm(m).ind(i+(1:nnew(m)),1)= true(nnew(m),1);
%     i=i+nnew(m);
% end
% newphi=newphi.';