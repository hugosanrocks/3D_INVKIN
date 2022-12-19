function cont1=nuevo_contorno(cont,m1,c1,xmin,xmax,cont1,im,fp,para)

x	= cont(m1,c1).vec.xc;
z	= cont(m1,c1).vec.zc;
vnx	= cont(m1,c1).vec.vnx;
vnz	= cont(m1,c1).vec.vnz;

ind     = find((x>=xmin).*(x<=xmax));
if fp==1
    ind = fliplr(ind);
end

x2	= x(ind);
z2	= z(ind);
vnx2= vnx(ind);
vnz2= vnz(ind);

if m1==cont1(im).m
    %el nuevo contorno es el mismo que antes
    %los medios en contacto son entonces con m0<m
    %non obstantes se tiene que tomar en cuenta los eventuales
    %intercesiones de este contorno con otros para identificar los medios
    %en contactos
    cv2     = c1*ones(length(x2),1);
    m       = cont1(im).m;
    nint    = length(cont(m,c1).int);
    k       = 1;
    xint    = [];
    for i=1:nint
        for j=2:length(cont(m,c1).int(i).xz(:,1))
            xint(k) = cont(m,c1).int(i).xz(j,1);
            k       = k+1;
        end
    end
    xint	= sort(xint);
    xint 	= [x2(1) xint x2(end)];
    xint    = unique(xint);
    n0     	= length(xint);
    xtest  	= (xint(1:(n0-1))+xint(2:(n0)))/2;
    for i=1:n0-1
        [~,ix]=min(abs(xtest(i)-x2));
        xtest0=x2(ix);
        ztest0=z2(ix);
        
        if c1==1
            ztest0=ztest0-1e-3;
        else
            ztest0=ztest0+1e-3;
        end
        m0  = inclusiontest(xtest0,ztest0,para,1);
        ix0=find(x2>=xint(i  ),1,'first');
        ix1=find(x2<=xint(i+1),1,'last');
        
        mv2(ix0:ix1,1) = m0;
    end
else
    %el nuevo contorno se construye sobre el contorno del medio m1,
    %asi que las posiciones arriba/abajo son inversadas
    c12 = [1 2];
    cv2 = c12(c12~=c1)*ones(length(x2),1);
    mv2 = m1*ones(length(x2),1);
end
x1  = cont1(im).vec.xc;
z1  = cont1(im).vec.zc;
vnx1= cont1(im).vec.vnx;
vnz1= cont1(im).vec.vnz;
cv1 = cont1(im).vec.cv; %contorno virtual (para saber si estamos arriba o abajo)
mv1 = cont1(im).vec.mv; %medio externo

if ~isempty(x1)
    if      ((x1(end)-x2(1))^2+(z1(end)-z2(1))^2)<=((x1(end)-x2(end))^2+(z1(end)-z2(end))^2) && ...
            ((x1(end)-x2(1))^2+(z1(end)-z2(1))^2)<=((x1(1  )-x2(1  ))^2+(z1(1  )-z2(1  ))^2)
        %camino esperado, nada que hacer
    elseif  ((x1(end)-x2(end))^2+(z1(end)-z2(end))^2)<=((x1(end)-x2(1))^2+(z1(end)-z2(1))^2) && ...
            ((x1(end)-x2(end))^2+(z1(end)-z2(end))^2)<=((x1(1  )-x2(1))^2+(z1(1  )-z2(1))^2)
        %x2 revuelto
        x2  =flipud(x2);
        z2  =flipud(z2);
        vnx2=flipud(vnx2);
        vnz2=flipud(vnz2);
        cv2 =flipud(cv2);
        mv2 =flipud(mv2);
    elseif  ((x1(1)-x2(1))^2+(z1(1)-z2(1))^2)<=((x1(end)-x2(1  ))^2+(z1(end)-z2(1  ))^2) && ...
            ((x1(1)-x2(1))^2+(z1(1)-z2(1))^2)<=((x1(end)-x2(end))^2+(z1(end)-z2(end))^2)
        %x1 revuelto
        x1  =flipud(x1);
        z1  =flipud(z1);
        vnx1=flipud(vnx1);
        vnz1=flipud(vnz1);
        cv1 =flipud(cv1);
        mv1 =flipud(mv1);
    end   
end

if ~(max(mv2)==0 && para.reg(m1).rho==0)
    %evita de agregar los contornos de los medios vacios en contacto con la superficie libre
    if para.reg(m1).rho==0
        indcut      =find(mv2==0);
        x2(indcut)  =[];
        z2(indcut)  =[];
        vnx2(indcut)=[];
        vnz2(indcut)=[];
        cv2(indcut) =[];
        mv2(indcut) =[];
    end
    
    
    cont1(im).vec.xc 	= [x1;  x2];
    cont1(im).vec.zc 	= [z1;  z2];
    cont1(im).vec.vnx 	= [vnx1;vnx2];
    cont1(im).vec.vnz	= [vnz1;vnz2];
    cont1(im).vec.cv	= [cv1; cv2];
    cont1(im).vec.mv    = [mv1; mv2];
    % else
    %     toto=0;
end