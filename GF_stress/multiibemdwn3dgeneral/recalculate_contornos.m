function para=recalculate_contornos(para,h)

para        = intersection(para);
cont        = para.cont;

if h ~= 0
% if para.smallscreen
    set(0, 'currentfigure', 100); cla; hold on
% else
%     figure(1);hold on;
%     set(gca,'ydir','reverse');
% end
end
c12         = [1 2];
testblack   = 0;
iim         = 0;    %#ok<NASGU> %numero du contour additionnel en cours (+1)

%juntar puntos de intercesion con puntos de superposicion
xtest10=zeros(para.nmed,2);
for m=1:para.nmed
    for c=1:2
        if isfield(cont(m,c),'seg')
            if ~isempty(cont(m,c).seg)
                ni=length(cont(m,c).int);
                ns=length(cont(m,c).seg);
                for i=1:ns
                    cont(m,c).int(ni+i).xz=cont(m,c).seg(i).xz;
                end
            end
        end
    end
    %puntos de principio y fin
    xtest10(m,1)=cont(m,c).xa;
    if para.geo(m)==2 || para.geo(m)==3 || para.geo(m)==4
        xtest10(m,2)=cont(m,c).xa+cont(m,c).a;
    elseif  para.geo(m)==1 && m>1
        xtest10(m,2)=cont(m,c).xa+2*cont(m,c).a;
    else %para.geo(1)==1 espacio completo
        xtest10(m,1)=-1e6;
        xtest10(m,2)= 1e6;
    end
end
%el ultimo medio se lleva el tereno, empezamos con el ultimo
%y se reescribe los contornos.
%A la diferencia de anteriormente, ya no se considera un contorno arriba y
%abajo por un medio, pero un solo contorno entero cerado
if para.siDesktop
    waitbarSwitch(0,h,'Check contornos');
else
    disp('Check contornos')
end
for m=para.nmed:-1:1
    icont   = 1;        % numero de contorno por el medio m (local)
    im      = m;        % numero del contorno actual
    
    cont1(im).vec.xc 	= []; %initialisation du contour
    cont1(im).vec.zc 	= [];
    cont1(im).vec.vnx 	= [];
    cont1(im).vec.vnz	= [];
    cont1(im).vec.cv	= [];
    cont1(im).vec.mv 	= [];
    cont1(im).m         = m;
    
    if (para.geo(1)==1 || para.geo(1)==3) && m==1
        break
    end
    
    for c=1:2
        ncont=[];
        
        if icont>1 || (m==1 && c==2)
            break
        end
        
        x     = cont(m,c).vec.xc;
        z     = cont(m,c).vec.zc;
        
        if testblack==1
            plot(x,z,'k');
        else
            %si le contour a des intersections on les prend en compte
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% compilation des points d intersections %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            %on commence a recherche tous les points d intersections de
            %tous les chemins possibles
            
            %nb initial d intersections qui se trouvent sur le contour en
            %question
            nint=length(cont(m,c).int);
            
            %on cherche ensuite les intersections derivees, celles qui
            %appartiennent aux contours qui ont d autres intersections avec
            %d autres contours
            
            %on commence avec les points de debut et fin du contour
            clear tmp
            if m==1 || para.geo(m)==2
                tmp(1,1:6)=[x(1),z(1),m,c,m,c];
                tmp(2,1:6)=[x(end),z(end),m,c,m,c];
            elseif m>1 && para.geo(m)==3
                tmp(1,1:6)=[x(1),z(1),m,c12(c12~=c),m,c];
                tmp(2,1:6)=[x(end),z(end),m,c,m,c];
            elseif m>1 && para.geo(m)==4
                tmp(1,1:6)=[x(1),z(1),m,c,m,c];
                tmp(2,1:6)=[x(end),z(end),m,c12(c12~=c),m,c];
            else
                tmp(1,1:6)=[x(1),z(1),m,c12(c12~=c),m,c];
                tmp(2,1:6)=[x(end),z(end),m,c12(c12~=c),m,c];
            end
            ptorigine=tmp;
            
            %verifie que les pts de depart et arrivee appartiennent au milieu
            %si il manque un des points on doit fermer le contour et relier c=1 avec c=2
            %et donc on supprime les points de depart et de fin
            if m>1 && para.geo(m)==1
                for i=2:-1:1 %on peut supprimer des termes alors on part de la fin
                    x1  = tmp(i,1);
                    z1  = tmp(i,2);
                    m01 = inclusiontest(x1,z1,para,1);
                    if m01~=m
                        tmp=[];
                        icont=2;
                        break
                        %un des pts au moins appartient a un autre milieu,
                        %on devra passer par un autre milieu et revenir au
                        %point de depart
                    end
                end
            elseif m>1 && para.geo(m)==3
                i=1;
                x1  = tmp(i,1);
                z1  = tmp(i,2);
                m01 = inclusiontest(x1,z1,para,1);
                if m01~=m
                    tmp=[];
                    icont=2;
                end
            elseif m>1 && para.geo(m)==4
                i=2;
                x1  = tmp(i,1);
                z1  = tmp(i,2);
                m01 = inclusiontest(x1,z1,para,1);
                if m01~=m
                    tmp=[];
                    icont=2;
                end
            end
            
            kk=size(tmp,1)+1;
            mlist=[];
            for k=1:nint
                %les intersections peuvent etre multiples et avec
                %plusieurs (au moins nint) contours differents
                m2 =cont(m,c).int(k).xz(1,1);
                c2 =cont(m,c).int(k).xz(1,2);
                xpb=cont(m,c).int(k).xz;
                
                zpb=xpb(2:end,2);
                xpb=xpb(2:end,1);
                
                if m2>m
                    mlist=union(mlist,m2);
                    for i=1:length(xpb)
                        tmp(kk,1:6)=[xpb(i),zpb(i),m2,c2,m,c];
                        kk=kk+1;
                    end
                    %suivre les intersections sous jacentes de ces contours
                    %avec d autres contours intermediaires
                    [tmp,kk]=follow_intersection_recursive(cont,m,m2,tmp,kk,mlist);
                end
            end
            
            %on elimine les repetitions
            [~,itmp]= unique(tmp(:,1:2),'rows');
            tmp     = tmp(itmp,:);%tri automatique
            
            %on elimine les points qui se situent exclusivement dans
            %d'autres milieux
            for i=size(tmp,1):-1:1
                x1  = tmp(i,1);
                z1  = tmp(i,2);
                dr  = 2*max(cont(tmp(i,5),tmp(i,6)).vec.r(2),cont(tmp(i,3),tmp(i,4)).vec.r(2));
                m01 = inclusiontest(x1,z1,para,0,dr);
                test=sum(m01==m);
                if test==0
                    tmp(i,:)=[];
                end
            end
            n1=size(tmp,1);
            
            
            %si c est une plaque et qu on retrouve les points d origine de
            %c==1 et c==2 alors la plaque est traversee par un element et
            %il faut donc imposer icont=2
            %             if para.geo(m)==2 && m>1 &&
            
            if n1>2  || icont == 2 %il y a plus de point que le debut et la fin
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% organisation des points d intersection %%%
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                xint    = tmp(:,1);         %abscisse des points d intersection
                mcont   = unique([tmp(:,3) tmp(:,5)]);
                xint0   = [xint;xtest10(mcont,1);xtest10(mcont,2)];   %on rajoute les debut et fin de chaque milieu en presence
                xint0   = sort(xint0);
                xtest0  = (xint0(1:(end-1))+xint0(2:end))/2;
                xtest0  = unique(xtest0);
                nk0     = length(xtest0);
                %tmp doit etre organise de maniere a refleter les sauts de contours
                %on commence avec le point de x plus petit qui verifie m et c
                tmp2        = tmp; %init
                itmp        = find((tmp(:,5)==m).*(tmp(:,6)==c));
                if ~isempty(itmp)
                    tmp2(1,:)       = tmp(itmp(1),:);
                    tmp(itmp(1),:)  = [];
                else
                    itmp=find((tmp2(:,1)==ptorigine(1,1)).*(tmp2(:,2)==ptorigine(1,2)));
                    if ~isempty(itmp)
                        tmp2(1,:)       = tmp(itmp(1),:);
                        tmp(itmp(1),:)  = [];
                    else
                        tmp2(1,:)       = tmp(1,:);
                        tmp(1,:)        = [];
                    end
                end
                tmp         =squeeze(tmp);
                m2          = tmp2(1,5);
                
                %puis on test les sauts pour que les segments n aient que
                %des points internes n appartenant qu au 2 milieux
                %consideres
                i=2;
                while ~isempty(tmp) || icont==2
                    if i==n1+(icont-1) && icont==2
                        %                         if para.geo(m)==1 && m>1
                        %le point de depart est aussi celui de la fin
                        tmp3	= tmp2(1,1:6);
                        icont   = 3;
                        %                         elseif  para.geo(m)==2 && m>1
                        %                             %le point d arrivee est le point de c==2
                        %                             x     = cont(m,2).vec.xc;
                        %                             z     = cont(m,2).vec.zc;
                        %                             tmp(1,1:6)=[x(1),z(1),m,c,m,c];
                        %                             tmp(2,1:6)=[x(end),z(end),m,c,m,c];
                        %                             tmp3	= tmp2(1,1:6);
                        %                             icont   = 3;
                        %                         end
                    else
                        %on reorganise tmp pour que le milieu d arrive soit en
                        %2eme position
                        itmp        = find((tmp(:,5)==m2).*(tmp(:,3)~=m2));
                        tmp(itmp,:) = tmp(itmp,[1 2 5 6 3 4]);
                        
                        %pour les points de debut et de fin on arrange pour que
                        %le contour d arrive soit le meme que celui de depart
                        %(utile que quand on fait le tour en entier)
                        itmp        = find((tmp(:,5)==m2).*(tmp(:,3)==m2));
                        if ~isempty(itmp)
                            itmp2       = find((tmp(itmp,6)==tmp2(i-1,6)));
                            tmp(itmp(itmp2),:) = tmp(itmp(itmp2),[1 2 5 6 3 4]);
                        end
                        
                        %on cherche tous les pts d intersection avec m2
                        tmp3        = tmp(tmp(:,3)==m2,:);
                        
                        %on arrange tmp pour que en priorite on ait les points
                        %avec le meme contour en premier
                        indtmp      = find((tmp3(:,4)==c));
                        indtmp2     = find((tmp3(:,4)~=c));
                        tmp3        = tmp3([indtmp;indtmp2],:);
                    end
                    
                    %parmi ceux la on regarde le seul (premier) contour dont les
                    %points entre les 2 pts d intersection appartiennent
                    %seulement aux 2 milieux communs
                    %on elimine de tmp3 ceux qui sont en contact avec 3
                    %milieux ou plus
                    %les points de test sont des points dont les abscisses
                    %sont situes au milieu des points d intersection
                    %comme il y a a chaque fois 2 chemins, dont 1 qui peut
                    %faire un coup a gauche un coup a droite et enfin
                    %retour a gauche, le nombre max de point a tester et 3
                    %fois la taille de xtest0
                    
                    xt1         = tmp2(i-1,1);
                    xt100       = xt1;
                    for j=1:size(tmp3,1)
                        %on teste chacun des points
                        xt2 = tmp3(j,1);
                        xt1 = xt100;
                        mt  = tmp3(j,3);%m2 milieu du contour a suivre
                        if m==mt
                            mt0 = tmp3(j,5);%milieu ext
                        else
                            mt0 = m;
                        end
                        
                        xi  = cont(mt,1).vec.xc(1);
                        xf  = cont(mt,1).vec.xc(end);
                        %                         xi  = cont(mt,1).xa;
                        %                         xf	= cont(mt,1).xa+2*cont(mt,1).a;
                        xtest=zeros(2,6*nk0);
                        ztest=zeros(2,6*nk0);
                        
                        xt10=min(xt1,xt2);
                        xt20=max(xt1,xt2);
                        xt1 =xt10;
                        xt2 =xt20;
                        if tmp3(j,4)==tmp2(i-1,6) || mt==1
                            tmp2(i,7)=1;
                            %on est sur le meme contour, soit on passe
                            %direct, soit on fait tout le tour
                            
                            %%% direct
                            ct              = tmp2(i-1,6);
                            xtest1          = xtest0(logical((xtest0>min(xt1,xt2)).*(xtest0<max(xt1,xt2))));
                            nk(1)           = length(xtest1);
                            ind             = zeros(nk(1),1);
                            for ik=1:nk(1)
                                [~,ind(ik)]=min(abs(cont(mt,ct).vec.xc-xtest1(ik)));
                            end
                            xtest(1,1:nk(1))=cont(mt,ct).vec.xc(ind);
                            ztest(1,1:nk(1))=cont(mt,ct).vec.zc(ind);
                            ind_i   = max(ind-1,1);
                            ind_f   = min(ind+1,length(cont(mt,ct).vec.zc));
                            if ~isempty(ind_i) && ~isempty(ind_f)
                                dr      = max(dr,.5*max(abs(cont(mt,ct).vec.r(ind_f)-cont(mt,ct).vec.r(ind_i))));
                            else
                                dr=1e-6;
                            end
                            if mt==1
                                ztest(2,1:nk(1))= ztest(1,1:nk(1));
                            else
                                %%% le tour
                                
                                %gauche
                                ct                  = tmp2(i-1,6);
                                xtest1              = xtest0(logical((xtest0>min(xt1,xi)).*(xtest0<max(xt1,xi))));
                                ng                  = length(xtest1);
                                ind                 = zeros(ng,1);
                                for ik=1:ng
                                    [~,ind(ik)]=min(abs(cont(mt,ct).vec.xc-xtest1(ik)));
                                end
                                xtest(2,1:ng)=cont(mt,ct).vec.xc(ind);
                                ztest(2,1:ng)=cont(mt,ct).vec.zc(ind);
                                ind_i   = max(ind-1,1);
                                ind_f   = min(ind+1,length(cont(mt,ct).vec.zc));
                                if  ~isempty(ind_i) && ~isempty(ind_f)
                                    dr      = max(dr,.5*max(abs(cont(mt,ct).vec.r(ind_f)-cont(mt,ct).vec.r(ind_i))));
                                end
                                %2eme contour
                                ct                  = c12(c12~=tmp2(i-1,6));
                                xtest1              = xtest0(logical((xtest0>min(xf,xi)).*(xtest0<max(xf,xi))));
                                nm                  = length(xtest1);
                                ind                 = zeros(nm,1);
                                for ik=1:nm
                                    [~,ind(ik)]=min(abs(cont(mt,ct).vec.xc-xtest1(ik)));
                                end
                                xtest(2,ng+1:ng+nm)=cont(mt,ct).vec.xc(ind);
                                ztest(2,ng+1:ng+nm)=cont(mt,ct).vec.zc(ind);
                                ind_i   = max(ind-1,1);
                                ind_f   = min(ind+1,length(cont(mt,ct).vec.zc));
                                if  ~isempty(ind_i) && ~isempty(ind_f)
                                    dr      = max(dr,.5*max(abs(cont(mt,ct).vec.r(ind_f)-cont(mt,ct).vec.r(ind_i))));
                                end                                %fin du tour
                                ct                  = tmp2(i-1,6);
                                xtest1              = xtest0(logical((xtest0>min(xf,xt2)).*(xtest0<max(xf,xt2))));
                                nd                  = length(xtest1);
                                ind                 = zeros(nd,1);
                                nk(2)               = ng+nm+nd;
                                for ik=1:nd
                                    [~,ind(ik)]=min(abs(cont(mt,ct).vec.xc-xtest1(ik)));
                                end
                                xtest(2,ng+nm+1:nk(2))=cont(mt,ct).vec.xc(ind);
                                ztest(2,ng+nm+1:nk(2))=cont(mt,ct).vec.zc(ind);
                                ind_i   = max(ind-1,1);
                                ind_f   = min(ind+1,length(cont(mt,ct).vec.zc));
                                if  ~isempty(ind_i) && ~isempty(ind_f)
                                    dr      = max(dr,.5*max(abs(cont(mt,ct).vec.r(ind_f)-cont(mt,ct).vec.r(ind_i))));
                                end
                            end
                            
                        else
                            tmp2(i,7)=2;
                            %on doit changer de contour (meme milieu)
                            %on doit passer soit par la gauche, soit par la
                            %droite et passer a l autre contour
                            
                            %%%% par la gauche
                            %gauche
                            if  xt1 == xt100
                                ct	= tmp2(i-1,6);
                            else
                                ct  = c12(c12~=tmp2(i-1,6));
                            end
                            xtest1              = xtest0(logical((xtest0>min(xt1,xi)).*(xtest0<max(xt1,xi))));
                            ng1                 = length(xtest1);
                            ind                 = zeros(ng1,1);
                            for ik=1:ng1
                                [~,ind(ik)]=min(abs(cont(mt,ct).vec.xc-xtest1(ik)));
                            end
                            xtest(1,1:ng1)=cont(mt,ct).vec.xc(ind);
                            ztest(1,1:ng1)=cont(mt,ct).vec.zc(ind);
                            ind_i   = max(ind-1,1);
                            ind_f   = min(ind+1,length(cont(mt,ct).vec.zc));
                            if  ~isempty(ind_i) && ~isempty(ind_f)
                                dr      = max(dr,.5*max(abs(cont(mt,ct).vec.r(ind_f)-cont(mt,ct).vec.r(ind_i))));
                            end                            %2eme contour
                            if  xt1 == xt100
                                ct  = c12(c12~=tmp2(i-1,6));
                            else
                                ct	= tmp2(i-1,6);
                            end
                            xtest1              = xtest0(logical((xtest0>min(xt2,xi)).*(xtest0<max(xt2,xi))));
                            nd1                 = length(xtest1);
                            nk(1)               = ng1+nd1;
                            ind                 = zeros(nd1,1);
                            for ik=1:nd1
                                [~,ind(ik)]=min(abs(cont(mt,ct).vec.xc-xtest1(ik)));
                            end
                            xtest(1,ng1+1:nk(1))=cont(mt,ct).vec.xc(ind);
                            ztest(1,ng1+1:nk(1))=cont(mt,ct).vec.zc(ind);
                            ind_i   = max(ind-1,1);
                            ind_f   = min(ind+1,length(cont(mt,ct).vec.zc));
                            if  ~isempty(ind_i) && ~isempty(ind_f)
                                dr      = max(dr,.5*max(abs(cont(mt,ct).vec.r(ind_f)-cont(mt,ct).vec.r(ind_i))));
                            end                            %%% par la droite
                            %droite
                            if  xt1 == xt100
                                ct	= tmp2(i-1,6);
                            else
                                ct  = c12(c12~=tmp2(i-1,6));
                            end
                            xtest1              = xtest0(logical((xtest0>min(xt1,xf)).*(xtest0<max(xt1,xf))));
                            nd2                 = length(xtest1);
                            ind                 = zeros(nd2,1);
                            for ik=1:nd2
                                [~,ind(ik)]=min(abs(cont(mt,ct).vec.xc-xtest1(ik)));
                            end
                            xtest(2,1:nd2)=cont(mt,ct).vec.xc(ind);
                            ztest(2,1:nd2)=cont(mt,ct).vec.zc(ind);
                            ind_i   = max(ind-1,1);
                            ind_f   = min(ind+1,length(cont(mt,ct).vec.zc));
                            if  ~isempty(ind_i) && ~isempty(ind_f)
                                dr      = max(dr,.5*max(abs(cont(mt,ct).vec.r(ind_f)-cont(mt,ct).vec.r(ind_i))));
                            end                            %2eme contour
                            if  xt1 == xt100
                                ct  = c12(c12~=tmp2(i-1,6));
                            else
                                ct	= tmp2(i-1,6);
                            end
                            xtest1              = xtest0(logical((xtest0>min(xt2,xf)).*(xtest0<max(xt2,xf))));
                            ng2                 = length(xtest1);
                            nk(2)               = ng2+nd2;
                            ind                 = zeros(ng2,1);
                            for ik=1:ng2
                                [~,ind(ik)]=min(abs(cont(mt,ct).vec.xc-xtest1(ik)));
                            end
                            xtest(2,nd2+1:nk(2))=cont(mt,ct).vec.xc(ind);
                            ztest(2,nd2+1:nk(2))=cont(mt,ct).vec.zc(ind);
                            ind_i   = max(ind-1,1);
                            ind_f   = min(ind+1,length(cont(mt,ct).vec.zc));
                            if  ~isempty(ind_i) && ~isempty(ind_f)
                                dr      = max(dr,.5*max(abs(cont(mt,ct).vec.r(ind_f)-cont(mt,ct).vec.r(ind_i))));
                            end
                        end
                        
                        %                         plot(xtest(1,1:nk(1)),ztest(1,1:nk(1)),'.');plot(xtest(2,1:nk(2)),ztest(2,1:nk(2)),'r.')
                        found=0;
                        for kt2=1:2 %les 2 chemins possibles
                            for ktest=1:nk(kt2) %les points internes
                                %                                 m01=contour_test(xtest(kt2,ktest),ztest(kt2,ktest),para,1e-3);
                                %                                 dr  = 2*max(cont(tmp(i,5),tmp(i,6)).vec.r(2),cont(tmp(i,3),tmp(i,4)).vec.r(2));
                                m01 = inclusiontest(xtest(kt2,ktest),ztest(kt2,ktest),para,0,dr);
                                %                                 if length(m01)==1 || length(m01)==3
                                %                                     %faux contour
                                %                                     break
                                %                                 end
                                if isempty(m01(m01==m))
                                    %chemin non direct
                                    break;
                                elseif length(m01)>2 && max(m01==0)==1 && max(m01>m)==1
                                    %chemin non direct mais avec pb
                                    break;
                                elseif ktest==nk(kt2)
                                    %chemin ok
                                    found       = 1;
                                    tmp2(i,1:6) = tmp3(j,:);
                                    tmp2(i,8)   = kt2;
                                    m2          = tmp2(i,5);
                                    
                                    %on elimine de tmp le pt choisit
                                    tmp(logical((tmp(:,1)==tmp2(i,1)).*(tmp(:,2)==tmp2(i,2))),:)=[];
                                    tmp=squeeze(tmp);
                                    
                                    %on eliminie aussi les points doubles
                                    %deja utilise (debut ou fin d un
                                    %contour croise) tant sur tmp2 que tmp
                                    indptdbl=find(tmp2(i,1)==tmp(:,1));
                                    indptdbl=sort(indptdbl);
                                    for iindpb=length(indptdbl):-1:1
                                        if abs(tmp2(i,2)-tmp(indptdbl(iindpb),2))<1e-10
                                            tmp(indptdbl(iindpb),:) = [];
                                            tmp2(end,:)=[];
                                            n1=n1-1;
                                        end
                                    end
                                    
                                    if (icont==1 || c==2 || m==1 ) && length(ncont)==0
                                        %si on arrive a la fin du contour et qu il reste des points
                                        %on eliminie le reste
                                        if tmp2(i,1:2)==ptorigine(2,1:2)
                                            tmp2(i+1:end,:)=[];
                                            tmp=[];
                                            icont=1;
                                        end
                                    end
                                    break
                                end
                            end
                            if found==1
                                break
                            end
                        end
                        if found==0 &&  j==size(tmp3,1)
                            %on a rien trouve
                            %il s agit d un nouveau contour
                            ncont(length(ncont)+1)=i;
                            tmp2(i,1:6) = tmp(1,:);
                            tmp2(i,8)   = 0;
                            m2          = tmp2(i,5);
                            tmp(1,:)    = [];
                            tmp         = squeeze(tmp);
                        elseif found==1
                            break
                        end
                    end
                    i=i+1;
                end
                ncont(ncont==n1)=[];
                ptint=tmp2;
            elseif n1==2 && icont==1
                ptint=tmp;
                ptint(1:end,7:8)=ones(size(ptint,1),2);
            else % n1==2 && icont==2
                tmp2(1,:)   = tmp(1,[1,2,5,6,3,4]);
                tmp2(2,:)   = tmp(2,:);
                tmp2(3,:)   = tmp2(1,:);
                ptint=tmp2;
                ptint(1:end,7:8)=ones(size(ptint,1),2);
            end
            
            n1 = size(ptint,1);
            if n1>0
                m1 = ptint(1,5);
                c1 = ptint(1,6);
                x1 = ptint(1,1);
                z1 = ptint(1,2);
                x10= x1;
                z10= z1;
                
                for k=2:n1
                    m2 = [ptint(k,3),ptint(k,5)];
                    c2 = [ptint(k,4),ptint(k,6)];
                    x2 = ptint(k,1);
                    z2 = ptint(k,2);
                    
                    recalculate_contornos_sub
                end
            end
        end
    end
    if para.siDesktop
        plot(cont1(im).vec.xc(1:10:end),cont1(im).vec.zc(1:10:end),'k','linewidth',3)
        %title(['wait ',num2str(para.nmed-m+1),'/',num2str(para.nmed)])
        %pause(.1)
    end
end
para.cont0  = cont1;
para.nmed0	= length(para.cont0);
%cuando hay contornos que se dividen, hay que reencontrar el medio exterior
%con la nueva numerotacion
nmed1  = length(cont1);
for m1=1:nmed1
    if ~isempty(m1)
        subm1(m1)=cont1(m1).m;
    end
end
for m=1:para.nmed
    subm(m).m=find(subm1==m);
end
para.subm1=subm1;
para.subm =subm;

% Hasta ahora se identifico los contornos resultantes de las interseciones.
% Ahora vamos a tomar en cuenta los medios que hay que fusionar y asi hacer
% desaparecer interfaces inecesarias.
% Tambien se va a indentificar los segmentos que son frontera entre 2 dos medios o de frontera libre
% De "cont1" extraemos "seg"
j   = 1;
fus = zeros(nmed1,nmed1);
if para.siDesktop
    waitbarSwitch(0,h,'Calculo segmentos de colocacion');
else
    disp('Calculo segmentos de colocacion');
end

for i = 1:nmed1
    %     waitbarSwitch(i/nmed1,h);
    m       = cont1(i).m;          %medio interior
    m1      = cont1(i).vec.mv;     %indices de los medios exterior a los puntos del contorno
    
    %Los segmentos que consideramos corresponden a una interface entre el
    %medio de indice m (el de medio actual) con uno de indice mas pequeno.
    %Los otros segmentos se encontraran con los otros contornos.
    %Asi no tenemos 2 interfases.
    
    %se prueba cada medio de indice mas pequeno
    for im = 0: m-1
        %indice de los puntos que verifican una interface entre el medio m y im
        indok   = find(m1==im);
        %verificacion de si las propiedades son diferentes
        res     = comp_mat(m,im,para);
        
        if ~isempty(indok) && res==0
            %no hay fusion, y hay puntos
            
            %se busca los indices consecutivos y por cada lista de indices
            %se crea un segmento
            indseg  = find(diff(indok)>1);
            k1      = 1;
            indseg(end+1)=length(indok);
            for k = 1:length(indseg)
                k2              = indseg(k);
                seg(j).mi       = i;            %nuevo medio de sub
                seg(j).m        = cont1(i).m;   %medio original
                seg(j).vec.xc   = cont1(i).vec.xc(indok(k1:k2));
                seg(j).vec.zc   = cont1(i).vec.zc(indok(k1:k2));
                seg(j).vec.vnx  = cont1(i).vec.vnx(indok(k1:k2));
                seg(j).vec.vnz  = cont1(i).vec.vnz(indok(k1:k2));
                seg(j).vec.cv   = cont1(i).vec.cv(indok(k1:k2));
                
                %identificacion del indice del medio exterior
                %como hubo nuevos sub-medios debido a interseciones
                %se debe de considerar los nuevos indices
                if im>0
                    if length(subm(im).m)>1
                        %este medio se dividio, hay que buscar el nuevo
                        %indice
                        nsegj=length(seg(j).vec.xc);
                        nsegj2=round(nsegj/2);
                        for isubm=1:length(subm(im).m)
                            indtest=find(cont1(subm(im).m(isubm)).vec.xc==seg(j).vec.xc(nsegj2));
                            if ~isempty(indtest)
                                if max(cont1(subm(im).m(isubm)).vec.zc(indtest)==seg(j).vec.zc(nsegj2))==1
                                    seg(j).mv = subm(im).m(isubm);
                                    break;
                                end
                            end
                        end
                    else
                        %este medio no se dividio y se queda entonces con su indice
                        seg(j).mv       = im;
                    end
                else
                    %la region externa es vacio debido a que es una
                    %frontera libre del semi-espacio
                    seg(j).mv       = im;
                end
                if para.siDesktop
                    plot( seg(j).vec.xc(1:10:end), seg(j).vec.zc(1:10:end),'r','linewidth',3)
                end
                
                k1              = k2+1;
                j               = j+1;
            end
        elseif ~isempty(indok) && res==1
            %hay fusion de 2 medios
            %hay que checar que la fusion no se haga con una subdivision de
            %otro media
            if length(subm(im).m)>1
                %como el medio puede estar varias veces en contacto con el
                %medio im por interupciones, hay que checar cada uno de los
                %segmentos
                ind2            = diff(indok);
                ind3            = find(abs(ind2)>1);
                n2              = length(ind3);
                indi            = zeros(n2+1,1);
                indf            = zeros(n2+1,1);
                indi(1)         = indok(1);
                indi(2:(n2+1))  = indok(ind3+1);
                indf(1:n2)      = indok(ind3);
                indf(n2+1)      = min(indok(end)+1,length(m1));
                nsegj2          = round((indi+indf)/2);
                for isubm=1:length(subm(im).m)
                    for ix = 1:length(nsegj2)
                        indtest=find(cont1(subm(im).m(isubm)).vec.xc==cont1(i).vec.xc(nsegj2(ix)));
                        if ~isempty(indtest)
                            if max(cont1(subm(im).m(isubm)).vec.zc(indtest)==cont1(i).vec.zc(nsegj2(ix)))==1
                                im1 = subm(im).m(isubm);
                                fus(m,im1) = 1;
                                fus(im1,m) = 1;
                            end
                        end
                    end
                end
            else
                %no hay otra subdivision
                im1=im;
                fus(m,im1) = 1;
                fus(im1,m) = 1;
            end
        end
    end
end

%heredar las fusiones entre medios
for i=1:nmed1
    if sum(fus(i,:))>1
        ind=find(fus(i,:)==1);
        n=length(ind);
        for k=1:n
            fus(ind(k),ind([1:k-1,k+1:n]))=1;
        end
    end
end

%luego se busca los segmentos contiguos con las mismas propiedades para
%juntarlos
nseg1 = j-1;
k = 1;
listm=[];
for i = 1:nseg1
    %cada segmento se puede juntar con 2 otros por cada punta
    go=0;
    if ~isempty(listm)
        if max(i==listm)==0
            %no ha sido ya llamado, es un segmento aparte
            go=1;
        end
    else
        go=1;
    end
    
    if go==1
        seg_unidos(k)        = seg(i);
        listm=[listm,i];
        
        %medios en contactos de indice min en caso de fusion
        mi(1)=seg(i).mi;
        if mi(1)~=0
            m0 = find(fus(:,mi(1))~=0);
            if ~isempty(m0)
                mi(1)=max([m0;mi(1)]);
            end
        end
        
        mi(2)=seg(i).mv;
        if mi(2)~=0
            m0 = find(fus(:,mi(2))~=0);
            if ~isempty(m0)
                mi(2)=max([m0;mi(2)]);
            end
        end
        seg_unidos(k).mi= mi(1);
        seg_unidos(k).mv= mi(2);
        
        dr = sqrt((seg(i).vec.xc(1)-seg(i).vec.xc(2))^2+(seg(i).vec.zc(1)-seg(i).vec.zc(2))^2);
        [seg_unidos(k),listm]=suit_seg(seg,seg_unidos(k),0,dr,mi,listm,fus);
        [seg_unidos(k),listm]=suit_seg(seg,seg_unidos(k),1,dr,mi,listm,fus);
        if para.siDesktop
            plot(seg_unidos(k).vec.xc(1:10:end),seg_unidos(k).vec.zc(1:10:end),'m','linewidth',3)
        end
        
        %calculo de r
        clear r
        n               = length(seg_unidos(k).vec.xc);
        r(1)            = 0;
        r(2:n)          = cumsum(sqrt(diff(seg_unidos(k).vec.xc).^2+diff(seg_unidos(k).vec.zc).^2));
        seg_unidos(k).vec.r    = r;
        
        k=k+1;
    end
end

para.fus    = fus;
if (exist('seg_unidos','var'))
    para.cont1	= seg_unidos;
    para.nmed1  = length(seg_unidos);
else
    error('No se distinguen contornos, revise velocidades de propagacion.')
end

nfus=0;
for i=1:nmed1
    nfus=nfus+sum(fus(i,i+1:nmed1));
end

nsubm=0;
for i=1:para.nmed
    nsubm=nsubm+length(subm(i).m);
end

para.nmedf=nsubm;%-nfus;
