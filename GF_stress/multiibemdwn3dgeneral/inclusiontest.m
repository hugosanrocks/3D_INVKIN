function medio=inclusiontest(xs,zs,para,braek,dr)
%se prueba la pertenencia a los medios segun los limites de los contornos
medio   = 1; %medio de fondo
m       = para.nmed;
cont    = para.cont;
if braek==1
    dr=1e-6;
end
while m>1
    test=0;
    if      xs-min([cont(m,1).vec.xc;cont(m,2).vec.xc])<-dr || max([cont(m,1).vec.xc;cont(m,2).vec.xc])-xs<-dr || ...
            zs-min([cont(m,1).vec.zc;cont(m,2).vec.zc])<-dr || max([cont(m,1).vec.zc;cont(m,2).vec.zc])-zs<-dr
        %de plano, no pertenece a m
        test=1;
    end
    if test==0 % de cada elemento de contorno 
        zc1  = interp1(cont(m,1).vec.xc,cont(m,1).vec.zc,xs,'spline'); %mitar de arriba
        zc2  = interp1(cont(m,2).vec.xc,cont(m,2).vec.zc,xs,'spline'); %mitad de abajo
        if (zc2-zs>=-dr) && (zs-zc1>=-dr)
            if braek==1
                medio=m;
                break;
            else
                medio=[medio,m];
                if abs(zc2-zs)<dr || abs(zs-zc1)<dr
                else
                    medio(1)=[];
                    break
                end
            end
        end
    end
    m=m-1;
end

%como el medio 1 es puesto a fuerza hay que verificar que no es en realidad el medio 0
%puede ser el caso cuando el punto esta en la interfase entre la superficie y otro medio
%no hay que suprimir tampoco los puntos que tienen el principio o fin de su
%contorno en la superificie, porque en este caso la interfase es entre 1,0 y
%el medio diferente
if max(medio==1)==1
    if para.geo(1)==3
        if zs<0
            medio(medio==1)=0;
        end
    elseif para.geo(1)==2
        z1=eq_contour_m1(xs,cont(1,1),para.geo(1));
        if z1-zs>dr
            %cuando hay algun tipo de rugosidad o un elemento en la
            %superficie y zs esta lejos del z1 de la interfase del
            %medio 1
            medio(medio==1)=0;
        elseif length(medio)>1
            for i=1:length(medio)
                if medio(i)~=1
                    %si el medio diferente de 1 tiene una interfase  la superficie es de frontera igual a la del medio 1
                    %el medio 1 esta presente
                    %solo se aceptara si los 2 z a este x son iguales
                    zc1  = interp1(cont(medio(i),1).vec.xc,cont(medio(i),1).vec.zc,xs,'spline');
                    zc2  = interp1(cont(medio(i),2).vec.xc,cont(medio(i),2).vec.zc,xs,'spline');
                    if abs(zc1-zc2)>1e-10 && abs(zs-z1)<1e-10
                        medio(medio==1)=0;
                        break
                    end
                end
            end
        end
    end
end

for i=1:length(medio)
    if isfield(para,'subm')
        %cambio de medio en caso de separacion de los contornos
        if medio(i)~=0
            if length(para.subm(medio(i)).m)>1
                for im=1:length(para.subm(medio(i)).m)
                    test=0;
                    m1  =para.subm(medio(i)).m(im);
                    xmin=min(para.cont0(m1).vec.xc);
                    zmin=min(para.cont0(m1).vec.zc);
                    xmax=max(para.cont0(m1).vec.xc);
                    zmax=max(para.cont0(m1).vec.zc);
                    if ((xmin-xs)>dr) || ((xs-xmax)>dr) || ((zmin-zs)>dr) || ((zs-zmax)>dr)
                        %de plano, no pertenece a m1
                        test=1;
                    end
                    if test==0
                        % on cherche les x autour de xs, ! il peut y avoir
                        % plus de 2 groupes de solutions, en haut et en bas
                        % et on interpole z autour des 2 valeurs possibles
                        
                        n       = length(para.cont0(m1).vec.xc);
                        %####% fatigue pour bien faire
                        dr      = sqrt((para.cont0(m1).vec.xc(2)-para.cont0(m1).vec.xc(1))^2+(para.cont0(m1).vec.zc(2)-para.cont0(m1).vec.zc(1))^2);
                        dr      = 3*dr;
                        indx    = find(abs(para.cont0(m1).vec.xc-xs)<dr);
                        %partes consecutivas
                        indseg  = find(diff(indx)>1);
                        indseg(end+1)=length(indx);
                        k1      = indx(1);
                        for k = 1:length(indseg)
                            k2	= indx(indseg(k));
                            xc0	= para.cont0(m1).vec.xc(k1:k2);
                            zc0	= para.cont0(m1).vec.zc(k1:k2);
                            [xc0,ind0]  = unique(xc0);
                            zc0         = zc0(ind0);
                            zc(k)= interp1(xc0,zc0,xs,'spline','extrap');
                            if k~=length(indseg)
                                k1	= indx(indseg(k)+1);
                            end
                        end
                        
                        %                             %premier point trouve
                        %                             [~,indxm11] = min(abs(para.cont0(m1).vec.xc-xs));
                        %                         indxm11i    = mod(indxm11-1,n);
                        %                         if indxm11i==0
                        %                             indxm11i=n;
                        %                         end
                        %                         indxm11f    = mod(indxm11+1,n);
                        %                         if indxm11f==0
                        %                             indxm11f=n;
                        %                         end
                        %                         if indxm11i>indxm11f
                        %                             indxm=[(1:indxm11f),(indxm11i:n)];
                        %                         else
                        %                             indxm=indxm11i:indxm11f;
                        %                         end
                        %                         xc0         = para.cont0(m1).vec.xc(indxm);
                        %                         zc0         = para.cont0(m1).vec.zc(indxm);
                        %                         [xc0,ind0]  = unique(xc0);
                        %                         zc0         = zc0(ind0);
                        %                         zc(1)       = interp1(xc0,zc0,xs,'spline','extrap');
                        %
                        %                         %deuxieme point trouve
                        %                         %on retire les points trouves et on cherche de
                        %                         %nouveau
                        %                         xc1         = para.cont0(m1).vec.xc([1:(indxm(1)-1),indxm(end)+1:end]);
                        %                         zc1         = para.cont0(m1).vec.zc([1:(indxm(1)-1),indxm(end)+1:end]);
                        %                         n           = length(xc1);
                        %                         [~,indxm11] = min(abs(xc1-xs));
                        %                         indxm11i    = mod(indxm11-1,n);
                        %                         if indxm11i==0
                        %                             indxm11i=n;
                        %                         end
                        %                         indxm11f    = mod(indxm11+1,n);
                        %                         if indxm11f==0
                        %                             indxm11f=n;
                        %                         end
                        %                         if indxm11i>indxm11f
                        %                             indxm=[indxm11i:n,1:indxm11f];
                        %                         else
                        %                             indxm=indxm11i:indxm11f;
                        %                         end
                        %                         xc0         = xc1(indxm);
                        %                         zc0         = zc1(indxm);
                        %                         [xc0,ind0]  = unique(xc0);
                        %                         zc0         = zc0(ind0);
                        %                         zc(2)       = interp1(xc0,zc0,xs,'spline','extrap');
                        zc=sort(zc);
                        if length(zc)==1
                            if zs-zc>=-1e-5
                                medio(i)=m1;
                                test=1;
                                break;
                            end
                        elseif mod(length(zc),2)==0
                            for iz=1:2:length(zc)
                                if (zc(iz+1)-zs>=-1e-5) && (zs-zc(iz)>=-1e-5)
                                    medio(i)=m1;
                                    test=1;
                                    break;
                                end
                            end
                        else
                            %cas d un pt de debut ou fin de contour inclu
                            %ds un autre medium
                            medio(i)=m1;
                            test=1;
                            break;
                        end
                        if test==1
                            break
                        end
                    end
                    
                end
            end
        end
    end
    
    
    if isfield(para,'fus')
        %cambio de medio en caso de fusion de los contornos
        if medio(i)~=0
            m1 = find(para.fus(:,medio(i))~=0);
            for ii=1:length(m1)
                if m1(ii)>medio(i)
                    medio(i)=m1(ii);
                end
            end
        end
    end
end
toto=5;