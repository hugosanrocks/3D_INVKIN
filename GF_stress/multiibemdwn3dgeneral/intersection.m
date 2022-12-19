function para=intersection(para)

cont    = para.cont;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  1ere etape, connaitre les limites des contours  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nbpt=1e4;
% cont2=zeros(para.nmed,2);

m = 1; c=1;
if para.geo(1)==2
    %semi-espace
    cont2(m,c).xmin = para.cont(m,c).vec.xc(1);
    cont2(m,c).xmax = para.cont(m,c).vec.xc(end);
    cont2(m,c).zmin = min(para.cont(m,c).vec.zc);
    cont2(m,c).zmax = max(para.cont(m,c).vec.zc);
else
    %espace complet
    cont2(m,c).xmin=1e6;
    cont2(m,c).xmax=1e6;
    cont2(m,c).zmin=1e6;
    cont2(m,c).zmax=1e6;
end

for m=2:para.nmed
    for c=1:2
        cont2(m,c).xmin = para.cont(m,c).vec.xc(1);
        cont2(m,c).xmax = para.cont(m,c).vec.xc(end);
        cont2(m,c).zmin = min(para.cont(m,c).vec.zc);
        cont2(m,c).zmax = max(para.cont(m,c).vec.zc);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  2eme etape, trouver les points d'intersections et identifier les       %
%    segments communs                                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(cont,'int')
    cont=rmfield(cont,'int');
end
ki=ones(para.nmed,2);
ks=ones(para.nmed,2);
for m=1:para.nmed
    for c=1:2
        cont(m,c).int=[];
    end
end

%pour chaque contour
for m=1:para.nmed
    %on test avec les autres contours
    for c=1:2
        if m==1 && c==2
            break;
        end
        for m2=m+1:para.nmed
            for c2=1:2
                %critere d intersection ou de superposition
                %on vérifie d'abord si coincidence en x et z
                %!attention condition necessaire pas mais pas suffisante
                if      cont2(m2,c2).xmax-cont2(m ,c ).xmin>-1e-15 && ...
                        cont2(m ,c ).xmax-cont2(m2,c2).xmin>-1e-15 && ...
                        cont2(m2,c2).zmax-cont2(m ,c ).zmin>-1e-15 && ...
                        cont2(m ,c ).zmax-cont2(m2,c2).zmin>-1e-15
                    
                    %on rediscretise avec la meme abscisse
                    xd  = max(cont2(m,1).xmin,cont2(m2,1).xmin);
                    xf  = min(cont2(m,1).xmax,cont2(m2,1).xmax);
                    
%                     x   = linspace(xd,xf,nbpt);
                    
                    ind1=logical((para.cont(m ,c ).vec.xc>=xd).*(para.cont(m ,c ).vec.xc<=xf));
                    x1  =para.cont(m ,c ).vec.xc(ind1);
                    
                    ind2=logical((para.cont(m2,c2).vec.xc>=xd).*(para.cont(m2,c2).vec.xc<=xf));
                    x2  =para.cont(m2,c2).vec.xc(ind2);
                    
                    x   = unique([x1;x2]).';
                    
                    z1  = interp1(para.cont(m ,c ).vec.xc,para.cont(m ,c ).vec.zc,x,'linear');
                    z2  = interp1(para.cont(m2,c2).vec.xc,para.cont(m2,c2).vec.zc,x,'linear');
                    
                    %pure intersection loin des bords
                    ind0= sign(z1-z2);
                    ind2= diff(ind0);
                    indi= find(abs(ind2)>1);
                    if ~isempty(indi)
                        if indi(end)==9999
                            indi(end)=[];
                        end
                        if indi(1)==1
                            indi(1)=[];
                        end
                    end
                    indi=squeeze(indi);
                    
                    %segment commun ou intersection
                    inds=find(abs(z1-z2)<1e-15);
                    %on peut tomber sur un pt d intersection isole,
                    %il faut pouvoir les separer
                    
                    %on retire de inds les indi
                    for i=1:length(indi)
                        inds(indi(i)==inds)=[];
                    end
                    inds=squeeze(inds);
                    
                    %puis on cherche les extremites des segments communs
                    %prealablement on supprime les points d intersections
                    %qui n ont pas ete reporte dans indi
                    %certaines intersection peuvent en effet se produire
                    %quand z1=z2 pour seulement 1 x et non un segment
                    if length(inds)==1
                        indi(end+1)=inds(1);
                        inds = [];
                    elseif ~isempty(inds)
                        deb=0;
                        for i=length(inds):-1:2
                            if inds(i)~=inds(i-1)+1 && deb == 0
                                %point isole
                                %il faut donc le rajouter a indi
                                indi(end+1)=inds(i);
                                %et le supprimer de inds
                                inds(i) = [];
                                deb     = 0;
                            else
                                deb     = 1;
                            end
                        end
                        if i==2 && deb==0
                            indi(end+1)=inds(1);
                            inds(1) = [];
                        end
                        inds=squeeze(inds);
                    end
                    indi=unique(indi);
                    %on cherche les extremites des segments communs
                    j=1;
                    for i=length(inds)-1:-1:2
                        if inds(i)==inds(i-1)+1 && inds(i)==inds(i+1)-1
                            indsj(j)=i;
                            j=j+1;
                        end
                    end
                    if exist('indsj','var')
                        inds(indsj)=[];clear indsj
                        inds=squeeze(inds);
                    end
                    
                    %le z est celui du milieu d indice le plus haut
                    if m>m2
                        zok = z1;
                    else
                        zok = z2;
                    end
                    
                    %inds est maintenant representatif des segments communs
                    %seulement et indi des points d'intersections seulement
                    if ~isempty(indi)
                        cont(m,c).int(ki(m,c)).xz=[m2 c2;x(indi).' zok(indi).'];
                        ki(m,c)=ki(m,c)+1;
                        %connaissance mutuelle
                        cont(m2,c2).int(ki(m2,c2)).xz=[m c;x(indi).' zok(indi).'];
                        ki(m2,c2)=ki(m2,c2)+1;
                    end
                    if para.geo(m2)==2 && para.geo(m)==2
                        %il faut retirer les intersections au debut et
                        %a la fin qui ne sont en réalité pas des début
                        %et fin
                        inds(find(x(inds)==x(1)))     =[];
                        inds(find(x(inds)==x(end)))   =[];
                    elseif para.geo(m2)==3 && para.geo(m)==2 || para.geo(m2)==2 && para.geo(m)==3
                        %il faut retirer les intersections
                        %a la fin qui ne sont en réalité pas des début
                        %et fin
                        inds(find(x(inds)==x(end)))   =[];
                    elseif para.geo(m2)==4 && para.geo(m)==2 || para.geo(m2)==2 && para.geo(m)==4
                        %il faut retirer les intersections au debut et
                        %a la fin qui ne sont en réalité pas des début
                        %et fin
                        inds(find(x(inds)==x(1)))     =[];
                    end
                    if ~isempty(inds)
                        cont(m,c).seg(ks(m,c)).xz=[m2 c2;x(inds).' zok(inds).'];
                        ks(m,c)=ks(m,c)+1;
                        %connaissance mutuelle
                        cont(m2,c2).seg(ks(m2,c2)).xz=[m c;x(inds).' zok(inds).'];
                        ks(m2,c2)=ks(m2,c2)+1;
                    end
                    %                     figure(2);
                    %                     hold on;
                    %                     plot(x,z1,'')
                    %                     plot(x,z2,'--r')
                    %                     plot(x(inds),z1(inds),'.')
                    %                     plot(x(inds),z2(inds),'.r')
                    %                     plot(x(indi),z1(indi),'x')
                    %                     plot(x(indi),z2(indi),'xr')
                    clear inds indi
                end
            end
        end
    end
    
    %si le milieu est une plaque non bornee sur les bords, il peut y avoir
    %des intersections avec d autres milieux sur les bords du domaine avec
    %entre autre le milieu 1 (espace complet ou semi-espace) et d autres
    %plaques
%     if para.geo(m)==2
%         for m2=m+1:para.nmed
%             if para.geo(m2)==2
%                 if cont(m,c).a
%                     
%                      if c==1
%                 para.cont(m,c).vec.xc(1)    = para.cont(m,c).xa-1e-6;
%                 para.cont(m,c).vec.zc(1)    = para.cont(m,c).za+para.cont(m,2).h;
%                 para.cont(m,c).vec.r(1)     = 0;
%                 para.cont(m,c).vec.vnx(1)   = 1;
%                 para.cont(m,c).vec.vnz(1)   = 0;
%             else
%                 para.cont(m,c).vec.xc(end)    = para.cont(m,c).xa+para.cont(m,c).a+1e-6;
%                 para.cont(m,c).vec.zc(end)    = para.cont(m,c).za+para.cont(m,1).h;
%                 para.cont(m,c).vec.vnx(end)   =-1;
%                 para.cont(m,c).vec.vnz(end)   = 0;
%                      end
%             
%                      
% %                 if m==1
% %                     cont2(m2,c2).zmax-cont2(m ,c ).zmin>-1e-15 && ...
% %                         cont2(m ,c ).zmax-cont2(m2,c2).zmin>-1e-15
% %                     
% %                     cont(m,c).int(ki(m,c)).xz=[m2 c2;x(indi).' zok(indi).'];
% %                     ki(m,c)=ki(m,c)+1;
% %                     %connaissance mutuelle
% %                     cont(m2,c2).int(ki(m2,c2)).xz=[m c;x(indi).' zok(indi).'];
% %                     ki(m2,c2)=ki(m2,c2)+1;
% %                 end
%             end
%         end
%     end
end
para.cont=cont;