function vec=malla_contorno_dx(nseg,cont,geo)
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

while max(dr)>drmin
    ind   	= find(dr>drmin);
    ind2    = find(diff(ind)>1,1,'first');
    indi    = ind(1);%[ind(1) ind(ind2+1)];
    if isempty(ind2)
        indf = n;
    else
        indf = ind(ind2)+1;
    end
    
    x0i  	= x(indi);
    x0f     = x(indf);

    n2 	= 3*(indf-indi);
    dr  = 2*drmin;
    while max(dr)>drmin
        n2      = n2*2;
        x1   	= linspace(x0i,x0f,n2);
        x1tmp   = x1-(cont.xa+cont.a);
        z1      = eq_contour(x1tmp,cont)+cont.za;
        
        r(1)  	= 0;
        r(2:n2)	= cumsum(sqrt(diff(x1).^2+diff(z1).^2));
        i      	= 2:(n2-1);
        dr      = zeros(1,n2);
        dr(i)  	= (r(i+1)-r(i-1))/2;
        dr(1) 	= (r(2)-r(1))/2;
        dr(n2)	= (r(n2)-r(n2-1))/2;
    end
    x       = [x x1(2:n2-1)];
    x       = sort(x);
    x       = unique(x);
    xtmp    = x-(cont.xa+cont.a);
    z   	= eq_contour(xtmp,cont)+cont.za;
    n       = length(x);
    r(1)  	= 0;
    r(2:n)	= cumsum(sqrt(diff(x).^2+diff(z).^2));
    i      	= 2:(n-1);
    dr      = zeros(1,n);
    dr(i)  	= (r(i+1)-r(i-1))/2;
    dr(1) 	= (r(2)-r(1))/2;
    dr(n)	= (r(n)-r(n-1))/2;
end
x       = sort(x);
x       = unique(x);


xtmp    = x-(cont.xa+cont.a);
z   	= eq_contour(xtmp,cont)+cont.za;

n       = length(x);
r(1)  	= 0;
r(2:n)	= cumsum(sqrt(diff(x).^2+diff(z).^2));
dzdx	= dzdx_contorno(xtmp.',cont);
norm	= sqrt(1.0+dzdx.^2);

i      	= 2:(n-1);
dr(i)  	= (r(i+1)-r(i-1))/2;
dr(1) 	= (r(2)-r(1))/2;
dr(n)	= (r(n)-r(n-1))/2;

vec.xc  = x.';
vec.zc  = z.';
vec.r   = r;
vec.vnx	= dzdx./norm;
vec.vnz	=  -1 ./norm;%le vector apunta siempre hacia abajo