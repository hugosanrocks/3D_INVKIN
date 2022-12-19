function coord3D=malla_geom_axi(coord,para)

% funcion que permite definir la posicion de los puntos de las fuentes
% virtuales, de las normales a los contornes y del los medios entre los
% puntos. Se ocupa coord para memorizar la informacion de los puntos de
% colocacion (2 medios colocados en 1 punto)
% tambien esta funcion da la posicion de los phis (densidad de fuerzas) que
% se occupa en la matriz y en el vector fuente
% los puntos son calculados en funcion de los coord 2D y transpuestos en 3D
% con axisimetria

nbpt2D  	= coord.nbpt; %numero de puntos en 2D a esta frecuencia

%l initialisation des vecteurs requiert le calcul du nombre de points en 3D
% calculo del numeros de puntos en 3D
nbpt3D=0;
for i=1:nbpt2D
    %numero de puntos associado a cada coord
    nbptc   = max(ceil(2*pi*abs(coord.x(i))/coord.dr(i)),4);
    %para verificiar la simetria, pero lleva a pb numericos
    while mod(nbptc,4)~=0
        nbptc=nbptc+1;
    end
    
    nbpt3D  = nbpt3D+nbptc;
end

n1          = para.nmed;
coord3D    	= struct( ...
    'x'     ,zeros(1,nbpt3D),   'y'     ,zeros(1,nbpt3D),	'z'     ,zeros(1,nbpt3D), ...
    'vnx'   ,zeros(1,nbpt3D),	'vny'   ,zeros(1,nbpt3D),   'vnz'   ,zeros(1,nbpt3D), ...
    'drxy'  ,zeros(1,nbpt3D),   'drxz'  ,zeros(1,nbpt3D),   'dA'    ,zeros(1,nbpt3D), ...
    'x0'    ,zeros(1,nbpt3D),   'th'    ,zeros(1,nbpt3D),   'dth'   ,zeros(1,nbpt3D), ...
    'mi'    ,zeros(nbpt3D,1),	'mv'    ,zeros(nbpt3D,1),   'cv'    ,zeros(nbpt3D,1), ...
    'Xim' 	,zeros(nbpt3D,n1),  'phi'  	,zeros(nbpt3D,n1),  'nbptc' ,zeros(1,nbpt3D));
coord3D.nbpt  = nbpt3D;

j   = 0;	%variable temporal para principios de coord de cada medio
jp  = 0;    %initializacion del indice de los phi

for i=1:nbpt2D
    nbptc   = max(ceil(2*pi*abs(coord.x(i))/coord.dr(i)),4);
    %para verificiar la simetria, pero lleva a pb numericos
    while mod(nbptc,4)~=0
        nbptc=nbptc+1;
    end
    
    jj              = j+(1:nbptc);
    j               = j+nbptc;
    
    dteta           = 2*pi/nbptc;
    coord3D.th(jj)  = (0:(nbptc-1))*dteta;
    coord3D.nbptc(jj)= nbptc;
    coord3D.dth(jj) = dteta;
    coord3D.x0(jj)  = coord.x(i);
    coord3D.x(jj)   = coord.x(i)*cos((0:(nbptc-1))*dteta);
    coord3D.y(jj)   = coord.x(i)*sin((0:(nbptc-1))*dteta);
    coord3D.z(jj)   = coord.z(i);
    coord3D.vnx(jj) = coord.vnx(i)*cos((0:(nbptc-1))*dteta);
    coord3D.vny(jj) = coord.vnx(i)*sin((0:(nbptc-1))*dteta);
    coord3D.vnz(jj) = coord.vnz(i);
    coord3D.cv(jj)	= coord.cv(i);
    coord3D.drxy(jj)= 2*pi*abs(coord.x(i))/nbptc;
    coord3D.drxz(jj)= coord.dr(i);
    %int 2*pi*r*dr/nbptc de -r0/2 a +r0/2=2*pi*x*r0
    %otra manera 2pi int (x*sqrt((dx/dt)^2+(dy/dt)^2))dt t=[t1,t2]
    coord3D.dA(jj)  = 2*pi*abs(coord.x(i))*coord.dr(i)/nbptc;
    for k=1:n1
        coord3D.Xim(jj,k)   = coord.Xim(i,k);	%ademas de saber que el punto Xi pertenece al medio m, se conoce su posicion arriba/abajo en el contorno
        if coord.phi(i,k)~=0
            coord3D.phi(jj,k)   = jp+(1:nbptc); %posicion de los phi en el vector X del sistema: AX=B
            jp              	= jp+nbptc;
        end
    end
end

%se censa los puntos que pertenecen al medio m o mas bien al sub medio mi
for mi=1:size(coord3D.Xim,2)
    coord3D.indm(mi).ind=(coord3D.Xim(:,mi)~=0);
end


%se censa los puntos que no tienen frontera libre de manera a aplicar la
%continuidad de los desplazamientos
coord3D.ju        = cumsum(sum(coord3D.Xim~=0,2)==2);
coord3D.nbeq      = coord3D.ju(coord3D.nbpt)+coord3D.nbpt;