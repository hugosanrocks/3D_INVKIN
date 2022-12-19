function vec=malla_contorno_dx_new(nseg,cont,geo)
% funcion que permite construir una malla fina de puntos perteneciendo a un
% contorno y calcula tambien sus normales
% la malla se calcula con un factor fac mas fina que por la frequencia mas
% alta que se tendra que calcular

%initialisation

%xc,zc      : posicion centro
%vnx, vnz   : vector normal al contorno
%r          : longitud lineic

dr      = cont.long/nseg;
fac     = 20;
drmin   = dr/fac;

%tant que la derive dz/dx est inferieure a 5 il n y aura pas de pb
nseg2   = 5*nseg*fac;


[x,z]       = visu_courbe(geo,cont,nseg2);
r           = zeros(1,nseg2);
r(1)        = 0;
r(2:nseg2)	= cumsum(sqrt(diff(x).^2+diff(z).^2));
i           = 2:(nseg2-1);
dr          = zeros(1,nseg2);
dr(i)       = (r(i+1)-r(i-1))/2;
dr(1)       = (r(2)-r(1));
dr(nseg2)   = (r(nseg2)-r(nseg2-1));
x0          = x;

if max(dr)>drmin
    ind   	= find(dr>drmin);
    %separacion en partes contiguas
    ind2    = find(diff(ind)>1);
    n2      = length(ind2);
    indi(2:(n2+1))  = ind(ind2+1)-1;
    indi(1)         = ind(1);
    indf(n2+1)      = min(ind(end)+1,nseg2);
    indf(1:n2)      = ind(ind2)+1;
    
    for i=1:n2+1
        xi  	= x0(indi(i));
        xf      = x0(indf(i));
        n     	= 4*(indf(i)-indi(i)+1);
        [x1,dr] = rec_malla_dx(drmin,xi,xf,n,cont);
        x       = [x,x1];
    end
end
x       = sort(x);
x       = unique(x);

% if geo==1
    xtmp    = x-(cont.xa+cont.a);
    z   	= eq_contour(xtmp,cont)+cont.za;
% elseif geo==2
%     xtmp    = x;
%     z   	= eq_contour_P(xtmp,cont)+cont.za;
% elseif geo==3
%     contL=cont;
%     contL.a =contL.ba;
%     contL.ba=0;
%     if contL.a==0 || contL.h==0
%         z   = eq_contour_P(x,cont);
%         z   = z+cont.za-cont.h;
%     else
%         dxL = -x(1)-contL.a;
%         xL  = x+dxL;
%         xL  = xL(1:find(xL>contL.a));
%         zL  = eq_contour(xL,contL);
%         nL2 = round(length(xL)/2);
%         xL  = xL(1:nL2)-dxL;
%         zL  = zL(1:nL2)+contL.za;
%         
%         xR   = x(nL2+1:end);
%         zR   = eq_contour_P(xR,cont);
%         zR   = zR+cont.za;
%         
%         x=[xL xR];
%         z=[zL zR];
%     end
%     xtmp    = x;
% elseif geo==4
%     contR=cont;
%     contR.a =contR.ba;
%     contR.ba=0;
%     if contR.a==0 || contR.h==0
%         z   = eq_contour_P(x,cont);
%         z   = z+cont.za-cont.h;
%     else
%         dxR = -x(end)+contR.a;
%         xR  = x+dxR;
%         xR  = xR(find(xR>-contR.a):end);
%         zR  = eq_contour(xR,contR);
%         nR2 = round(length(xR)/2);
%         xR  = xR(nR2:end)-dxR;
%         zR  = zR(nR2:end)+contR.za;
%         
%         xL   = x(1:end-nR2-1);
%         zL   = eq_contour_P(xL,cont);
%         zL   = zL+cont.za;
%         
%         x=[xL xR];
%         z=[zL zR];
%     end
%     xtmp    = x;
% end
n       = length(x);
r(1)  	= 0;
r(2:n)	= cumsum(sqrt(diff(x).^2+diff(z).^2));
i      	= 2:(n-1);
dr(i)  	= (r(i+1)-r(i-1))/2;
dr(1) 	= (r(2)-r(1))/2;
dr(n)	= (r(n)-r(n-1))/2;
dzdx	= dzdx_contorno(xtmp.',cont);
norm	= sqrt(1.0+dzdx.^2);

%problemas con bordes verticales o con dz/dz demasiado grande
if max(dr)>drmin
    ind             = find(dr>drmin);
    %separacion en partes contiguas
    ind2            = find(diff(ind)>1);
    n2              = length(ind2);
    indi            = zeros(n2+1,1);
    indf            = zeros(n2+1,1);
    indi(2:(n2+1))  = ind(ind2+1)-1;
    indi(1)         = ind(1);
    indf(n2+1)      = min(ind(end)+1,n);
    indf(1:n2)      = ind(ind2)+1;
    
    for i=n2+1:-1:1
        xi  	= x(indi(i));
        xf      = x(indf(i));
        zi  	= z(indi(i));
        zf      = z(indf(i));
        r0      = sqrt((zi-zf)^2+(xi-xf)^2);
        n     	= ceil(r0/drmin);
        x1      = linspace(xi,xf,n);
        z1      = linspace(zi,zf,n);
        dzdx1	= (zf-zi)/(xf-xi)*ones(n,1);
        norm1	= sqrt(1.0+dzdx1.^2);
        x       = [x(1:indi(i)-1),x1,x(indf(i)+1:end)];
        z       = [z(1:indi(i)-1),z1,z(indf(i)+1:end)];
        dzdx	= [dzdx(1:indi(i)-1);dzdx1;dzdx(indf(i)+1:end)];
        norm	= [norm(1:indi(i)-1);norm1;norm(indf(i)+1:end)];
    end
    n       = length(x);
    r(1)  	= 0;
    r(2:n)	= cumsum(sqrt(diff(x).^2+diff(z).^2));
%     i      	= 2:(n-1);
%     dr(i)  	= (r(i+1)-r(i-1))/2;
%     dr(1) 	= (r(2)-r(1))/2;
%     dr(n)	= (r(n)-r(n-1))/2;
%     if geo==1
%         xtmp    = x-(cont.xa+cont.a);
%     else
%         xtmp    = x;
%     end
end



vec.xc  = x.';
vec.zc  = z.';
vec.r   = r;
vec.vnx	= dzdx./norm;
vec.vnz	=  -1 ./norm;%le vector apunta siempre hacia abajo