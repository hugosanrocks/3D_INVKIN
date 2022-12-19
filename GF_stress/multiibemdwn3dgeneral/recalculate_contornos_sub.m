if ptint(k,7)==1
    %le second point appartient au meme contour que le premier point
    if  ptint(k,8)==1
        %chemin direct
        if x1>x2
            xmin=x2;
            xmax=x1;
            fp  =1;
            indpb=find((x>=x2).*(x<=x1));
            indpb=fliplr(indpb);
        else
            xmin=x1;
            xmax=x2;
            fp  =0;
        end
        cont1=nuevo_contorno(cont,m1,c1,xmin,xmax,cont1,im,fp,para);
    else
        %il faut faire le tour:
        %il faut finir un des 2 bouts du contour du milieu m1,
        %prendre le 2eme contour du meme milieu m1 et
        %revenir au premier contour pour aboutir au 2eme point d intersection
        
        %on cherche alors a determiner la direction a suivre
        %pour rester en contact avec le milieu m
        %pour cela on cherche les milieux auxquels
        %appartiennent les extremites du contours
        
        %point gauche
        xg  = cont(m1,1).xa;
        xd	= cont(m1,1).xa+2*cont(m1,1).a;
        
        if x1<x2
            %il faut aller vers la gauche
            xd1     = x1;
            xg1     = xg;
            fp1     = 1;
            xg2     = xg;
            xd2     = xd;
            fp2     = 0;
            xg3     = x2;
            xd3     = xd;
            fp3     = 1;
        else
            %il faut aller vers la droite %pas tester
            xg1     = x1;
            xd1     = xd;
            fp1     = 0;
            xg2     = xg;
            xd2     = xd;
            fp2     = 1;
            xd3     = x2;
            xg3     = xg;
            fp3     = 0;
        end
        
        %meme milieu, meme contour
        cont1=nuevo_contorno(cont,m1,c1,xg1,xd1,cont1,im,fp1,para);
        %meme milieu, 2eme contour
        cont1=nuevo_contorno(cont,m1,c12(c12~=c1),xg2,xd2,cont1,im,fp2,para);
        %meme milieu, meme contour
        cont1=nuevo_contorno(cont,m1,c1,xg3,xd3,cont1,im,fp3,para);
    end
elseif ptint(k,7)==2
    %le 2eme point appartient au 2eme contour du meme milieu du point d origine
    
    %il faut finir un des 2 bouts du contour du milieu m1,
    %prendre le 2eme contour du meme milieu m1 et aboutir au 2eme point d intersection
    
    xg=cont(m1,1).xa;
    xd=cont(m1,1).xa+2*cont(m1,1).a;

    if  ptint(k,8)==1
        %on doit passer par la gauche
        xd1     = x1;
        xg1     = xg;
        xg2     = xg;
        xd2     = x2;
        fp1     = 1;
        fp2     = 0;
    else
        %il faut aller vers la droite
        %abscisse du point droit
        xg1     = x1;
        xd1     = xd;
        xd2     = xd;
        xg2     = x2;
        fp1     = 0;
        fp2     = 1;
    end
    %meme milieu, meme contour
    cont1=nuevo_contorno(cont,m1,c1,xg1,xd1,cont1,im,fp1,para);
    %meme milieu, 2eme contour
    cont1=nuevo_contorno(cont,m1,c12(c12~=c1),xg2,xd2,cont1,im,fp2,para);
end
%     plot(cont1(im).vec.xc(1:10:end),cont1(im).vec.zc(1:10:end),'c')


%le prochain point d origine appartient au deuxieme
%milieu du point d intersection
if sum(m2==m1)>1
    tmp=find(c2==c1);
    c2(tmp(1))=[];c2=squeeze(c2);
    m2(1)=[];m2=squeeze(m2);
elseif sum(m2==m1)==1
    c2(m2==m1)=[];c2=squeeze(c2);
    m2(m2==m1)=[];m2=squeeze(m2);
else
    %we got a pb
    toto=0;
end
m1=m2;
c1=c2;
x1=x2;
z1=z2;

if sum((k+1)==ncont)>0
    %le prochain point d intersection appartient a un autre contour
    %on finit donc ce contour
    %pour trouver le point de fin, il faut regarder si on a change de
    %contour ou non
    
    ncont2              = ncont;
    ncont((k+1)==ncont) = [];%pour eviter la recursivite
    ncont               = squeeze(ncont);
    if icont==1
        if m1==m && c==2 %|| para.geo(m1)==2 && para.geo(m)==2 
            %les contours se finissent tout seul
            
            icont   =icont+1;
        elseif m1==m && c1~=c
            %on recoupe le contour oppose du meme milieu
            %il faut donc separer le contour en 2 contours distincts
            %pour finir le 1er sous nouveau contour, on revient au point d origine du
            %contour original
            %pour commencer le 2eme sous nouveau contour, on part de l'autre point
            %final
            %si il reste des points ds "ptint" apres ces 2 nouveaux sous contours
            %cela signifie qu il y a d autres sous contours (pas tester... %###)
            
            cont1   =nuevo_contorno(cont,m1,c1,x10,x1,cont1,im,1,para);

            icont   =icont+1;
        else
            %on finit le contour
            contour_rebuild_end_point;
            recalculate_contornos_sub;
        end
    else
        %les contours se finissent tout seul
    end
    if para.siDesktop
        plot(cont1(im).vec.xc(1:10:end),cont1(im).vec.zc(1:10:end),'k','linewidth',3)
    end
    %le prochain point d origine appartient au deuxieme
    %milieu du point d intersection
    ik=find(k+1==ncont2);
    if ik==length(ncont2)
        %le prochain point d origine est le dernier point de ptint
        if para.geo(m)==2 % && c==1
            x1=cont(m,2).vec.xc(end);
            z1=cont(m,2).vec.zc(end);
            m1=m;
            c1=2;
        else
            m1=ptint(n1,5);
            c1=ptint(n1,6);
            x1=ptint(n1,1);
            z1=ptint(n1,2);
        end
    else
        m1=ptint(ncont(1)-1,5);%###
        c1=ptint(ncont(1)-1,6);
        x1=ptint(ncont(1)-1,1);
        z1=ptint(ncont(1)-1,2);
    end
    ptint(k+1,7:8)=[1 1];

    ncont   = ncont2;
    iim     = iim+1;
    im      = para.nmed+iim;
    cont1(im).vec.xc 	= [];
    cont1(im).vec.zc 	= [];
    cont1(im).vec.vnx 	= [];
    cont1(im).vec.vnz	= [];
    cont1(im).vec.cv	= [];
    cont1(im).vec.mv 	= [];
    cont1(im).m         = m;
end