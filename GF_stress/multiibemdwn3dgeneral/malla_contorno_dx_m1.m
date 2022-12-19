function vec=malla_contorno_dx_m1(nseg,cont,geo)
% funcion que permite construir una malla fina de puntos perteneciendo al
% background, y calcula tambien sus normales
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

[x,z]       = visu_courbe_m1(geo,cont,nseg2);
r(1)        = 0;
r(2:nseg2)	= cumsum(sqrt(diff(x).^2+diff(z).^2));
i           = 2:(nseg2-1);
dr(i)       = (r(i+1)-r(i-1))/2;
dr(1)       = (r(2)-r(1))/2;
dr(nseg2)   = (r(nseg2)-r(nseg2-1))/2;

if max(dr)>drmin
    ind   	= find(dr>drmin);
    ind2    = find(diff(ind)>1);
    indi    = [ind(1) ind(ind2+1)];
    indf    = [ind(ind2)+1 max(ind(end)+1,nseg2)];
    x0      = x;
    n     	= length(indi);
    for j=1:n
        n2 	= 3*(indf(j)-indi(j));
        dr  = 2*drmin;
        while max(dr)>drmin
            n2      = n2*2;
            x1   	= linspace(x0(indi(j)),x0(indf(j)),n2);
            z1   	= eq_contour_m1(x1,cont,geo);
            
            r(1)  	= 0;
            r(2:n2)	= cumsum(sqrt(diff(x1).^2+diff(z1).^2));
            i      	= 2:(n2-1);
            dr(i)  	= (r(i+1)-r(i-1))/2;
            dr(1) 	= (r(2)-r(1))/2;
            dr(n2)	= (r(n2)-r(n2-1))/2;
        end
        x   = [x x1(2:n2-1)];
    end
end
x       = sort(x);
x       = unique(x);
z   	= eq_contour_m1(x,cont,geo);
n       = length(x);
r(1)  	= 0;
r(2:n)	= cumsum(sqrt(diff(x).^2+diff(z).^2));
dzdx    = dzdx_contorno_m1(x.',cont,geo);

norm	= sqrt(1.0+dzdx.^2);

vec.xc  = x.';
vec.zc  = z.';
vec.r   = r;
vec.vnx	= dzdx./norm;
vec.vnz	=  -1 ./norm;%le vecteur pointe toujours vers le bas

