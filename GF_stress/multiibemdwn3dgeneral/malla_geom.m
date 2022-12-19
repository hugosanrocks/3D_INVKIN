function [coord,list]=malla_geom(para,fj)
% funcion que permite definir la posicion de los puntos de las fuentes
% virtuales, de las normales a los contornes y del los medios entre los
% puntos. Se ocupa coord para memorizar la informacion de los puntos de
% colocacion (2 medios colocados en 1 punto)
% tambien esta funcion da la posicion de los phis (densidad de fuerzas) que
% se occupa en la matriz y en el vector fuente

nbpt0  	= 0; %numero de puntos a esta frecuencia
cont1   = para.cont1;
medio(para.nmedf).ind =[];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% discretizacion de los contornos %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%---------------------------------------%
% coordenadas de los puntos y de las fv %
%---------------------------------------%

%para la initializacion
% l=zeros(para.nmed1,1);
l(para.nmed1).ind=0;
switch nargin
    case 2
        %ya estamos en el ciclo de frecuencias
        fj      = abs(fj);
        for m1=1:para.nmed1
            % indices de los puntos que cumplen la discretizacion por la frecuencia dada
            l(m1).ind	= malla_contorno_fast2(cont1(m1),para,fj);
            nbpt0      	= nbpt0+length(l(m1).ind);
        end
    case 1
        %estamos en el caso de la mas alta frecuencia
        %se necesita todos los puntos para despues obtener la continuidad
        %de los phis
        for m1=1:para.nmed1
            l(m1).ind	= 1:length(cont1(m1).vec.xc);
            nbpt0      	= nbpt0+length(l(m1).ind);
        end
end

% numeros de puntos del sistema
n1          = para.nmed;
coord       = struct( ...
    'x'     ,zeros(1,nbpt0),	'z'     ,zeros(1,nbpt0), ...
    'vnx'   ,zeros(1,nbpt0),	'vnz'   ,zeros(1,nbpt0), ...
    'dr'    ,zeros(1,nbpt0),    'cv'    ,zeros(nbpt0,1), ...
    'Xim' 	,zeros(nbpt0,n1),   'phi'  	,zeros(nbpt0,n1), ...
    'iectz'	,zeros(nbpt0,1),   'iecuz'	,zeros(nbpt0,1), ...
    'fl' 	,zeros(nbpt0,n1));
coord.nbpt  = nbpt0;

j   = 0;	%variable temporal para principios de coord de cada medio
jp  = 0;    %initializacion del indice de los phi

for m1=1:para.nmed1
    m   = cont1(m1).m;	% el medio original
    l1 	= l(m1).ind;
    if isempty(l1)
        continue
    end
    n   = length(l1);
    jj  = j+(1:n);
    j   = j+n;
    
    coord.x(jj)     = cont1(m1).vec.xc(l1);
    coord.z(jj)     = cont1(m1).vec.zc(l1);
    coord.vnx(jj)   = cont1(m1).vec.vnx(l1);
    coord.vnz(jj)   = cont1(m1).vec.vnz(l1);
    coord.cv(jj)	= cont1(m1).vec.cv(l1);
    clear r
    r               = cont1(m1).vec.r(l1);
    
%     %promedio de las posiciones y normales para mejorar un poco los
%     %resultados
%     switch nargin
%         case 2
%             l0  = round((l1(1:(end-1))+l1(2:end))/2);
%             l2  = [ 1;l0];
%             l3  = [ l0;length(cont1(m1).vec.xc)];
%             for k=1:n
%                 coord.x(jj(k))  = mean(cont1(m1).vec.xc(l2(k):l3(k)));
%                 coord.z(jj(k))	= mean(cont1(m1).vec.zc(l2(k):l3(k)));
%attention aux saut ds vn, le moyennage introduit des erreur
%                 coord.vnx(jj(k))= mean(cont1(m1).vec.vnx(l2(k):l3(k)));
%                 coord.vnz(jj(k))= mean(cont1(m1).vec.vnz(l2(k):l3(k)));
%                 coord.cv(jj(k)) = round(median(cont1(m1).vec.cv(l2(k):l3(k))));
%                 r(k)            = mean(cont1(m1).vec.r(l2(k):l3(k)));
%             end
%             norm        = sqrt(coord.vnx.^2+coord.vnz.^2);
%             coord.vnx   = coord.vnx./norm;
%             coord.vnz   = coord.vnz./norm;
%     end
    
    i = 2:n-1;
    if n==1
        dr              = cont1(m1).vec.r(end);
    else
        dr              = zeros(n,1);
        dr(2:n-1)       = (r(i+1)-r(i-1))/2;
        dr(1)           = (r(2)+r(1))/2;
        dr(n)           = cont1(m1).vec.r(end)-(r(n)+r(n-1))/2;
    end
    coord.dr(jj)    = dr;
    
    mext	= cont1(m1).mv; %medio exterior
    mint	= cont1(m1).mi;	%medio interior
    c       = coord.cv(jj); %posicion en el contorno
    
    %indexacion (ver si hay un punto o dos), medios en contacto, y normal
    if mext==0 || para.reg(mext).rho==0 %2015-01-16
        %superficie libre y hueco afuera
        coord.Xim(jj,mint)  = c;	%ademas de saber que el punto Xi pertenece al medio m, se conoce su posicion arriba/abajo en el contorno
        coord.phi(jj,mint)  = jp+(1:n); %posicion de los phi en el vector X del sistema: AX=B
        jp                  = jp+n;
    elseif mint==0 || para.reg(m).rho==0
        %para tomar en cuenta los huecos
        coord.Xim(jj,mext)  = (c==2)+2*(c==1);%inversion de la posicion
        coord.phi(jj,mext)  = jp+(1:n);
        jp                  = jp+n;
    else
        coord.Xim(jj,mext)  = (c==2)+2*(c==1);
        coord.phi(jj,mext)  = jp+(1:n);
        jp                  = jp+n;
        
        coord.Xim(jj,mint)	= c;
        coord.phi(jj,mint) 	= jp+(1:n);
        jp                  = jp+n;
    end
    
     
    if mint~=0
        medio(mint).ind	= [medio(mint).ind;[jj(1) jj(end)]];%principio y fin
        if para.reg(mint).bet==0 && para.reg(mint).rho~=0 
            %indicacion de presencia de fluido
            coord.fl(jj,mint)	= 1;
        end
    end
    if mext~=0
        medio(mext).ind	= [medio(mext).ind;[jj(1) jj(end)]];%principio y fin
        if  para.reg(mext).bet==0 && para.reg(mext).rho~=0 
            %indicacion de presencia de fluido
            coord.fl(jj,mext)    = 1;
        end
    end
end

%se censa los puntos que pertenecen al medio m o mas bien al sub medio mi
for mi=1:size(coord.Xim,2)
    coord.indm(mi).ind=(coord.Xim(:,mi)~=0);
end

% switch nargin
%     case 2
%         go=1;
%     otherwise
%         go=0;
% end
%
% for mi=1:size(coord.Xim,2)
%     if ~isempty(medio(mi).ind)
%         if go==0
%             drref           = mean(coord.dr);
%             indtri          = (medio(mi).ind(1,1):medio(mi).ind(1,2));
%             listc           = 1;
%             [indtri,listc]  = suit_cont(indtri,medio(mi).ind,listc,coord,drref);
%             list(mi).cont   = listc;
%         else
%             indtri=zeros(1,sum(medio(mi).ind(:,2)-medio(mi).ind(:,1)));
%             i1=1;
%             for i=1:length(para.listc(mi).cont)
%                 j   = para.listc(mi).cont(i);
%                 ja  = abs(j);
%                 n   = medio(mi).ind(ja,2)-medio(mi).ind(ja,1);
%                 i2  = i1+n;
%                 if sign(j)==1
%                     indtri(i1:i2)= medio(mi).ind(ja,1):   medio(mi).ind(ja,2);
%                 else
%                     indtri(i1:i2)= medio(mi).ind(ja,2):-1:medio(mi).ind(ja,1);
%                 end
%                 i1=i2+1;
%             end
%         end
%         coord.m(mi).ind = indtri;
% %         
% %         if max(list(mi).cont==0)==1
% %             ncontindpt  = sum(list(mi).cont==0)/2+1;
% %             k=find(list(mi).cont==0);
% %             for kk=1:k/2
% %                 intri
% %             
% %         else
% %             ncontindpt=1;
% %         end
% %         for k=1:ncontindpt
%             
%             xtmp    = coord.x(indtri);
%             ztmp    = coord.z(indtri);
%             
%             clear r
%             n               = length(indtri);
%             r(1)            = 0;
%             r(2:n)          = cumsum(sqrt(diff(xtmp).^2+diff(ztmp).^2));
%             coord.m(mi).r   = r+coord.dr(indtri(1))/2;
%             switch nargin
%                 case 1
%                     figure(1);plot(xtmp(1:10:end),ztmp(1:10:end),'c')
%                     toto=0;
%             end
% %         end
%     end
% end
list=[];

%se censa los puntos que no tienen frontera libre de manera a aplicar la
%continuidad de los desplazamientos
coord.ju        = cumsum(sum(coord.Xim~=0,2)==2);
coord.nbeq      = coord.ju(coord.nbpt)+coord.nbpt;
coord.nfl       = sum(sum(coord.fl));


if para.pol==2 && para.dim==1
    %en caso que hay presencia de fluido (P-SV)
    %los phi_z sont reacomodados:
    %phi_z(elmt==i)=phi_x(elmt==i)+nbeq-sum_j(fl(j,:)==1)
    %sino phi_z(elmt==i)=phi_x(elmt==i)+nbeq
    if max(max(coord.fl))>0
        phiz    = coord.phi.*(coord.fl==0);
        indphiz = reshape(phiz,coord.nbpt*para.nmedf,1);
        indphiz(indphiz==0)=[];
        logicphi= false(coord.nbeq,1);
        logicphi(indphiz) = true;
    else
        logicphi= true(coord.nbeq,1);
    end
    coord.logicphi=logicphi;
    coord.indphiz0=cumsum(logicphi).*logicphi;
    
    %indice de las ecuaciones tz y uz o segun tt y ut
    %los tz se toman en cuenta si no hay 2 fluidos en contacto ni si hay
    %fluido/libre
    coord.ectz  = (sum(coord.fl,2)==1).*(sum(coord.Xim~=0,2)==2)+(sum(coord.fl,2)==0);
    coord.iectz = coord.nbeq+cumsum(coord.ectz);
    
    %los uz se toman en cuenta si no hay fluidos y si la frontera no es
    %libre
    coord.iecuz=coord.iectz(coord.nbpt)+cumsum((sum(coord.Xim~=0,2)==2).*(sum(coord.fl,2)<1));
end