function para = axi_contornos(para)

%suppression de tous les contours avec x>0
for i=para.nmed1:-1:1
    if min(para.cont1(i).vec.xc)>=0
        para.nmed1=para.nmed1-1;
        para.cont1(i)=[];
    elseif max(para.cont1(i).vec.xc)>=0
        ind=logical(para.cont1(i).vec.xc>=0);
        ch =find(diff(ind)~=0);
        nch=length(ch);%numero de cambios (changements)
        if nch==1
           %on change rien
           para.cont1(i).vec.r(ind)    = [];
           ind0=true(length(para.cont1(i).vec.r),1);
           para.cont1(i).vec.r(ind0)   = para.cont1(i).vec.r(ind0) -min(para.cont1(i).vec.r(ind0));
        elseif nch==2
            %contour en boucle mais on commence  - + - ou + - +
            %il faut faire une rotation et recalculer r
            n2      = length(ind);
            dr      = zeros(1,n2);
            j      	= 2:(n2-1);
            r       = para.cont1(i).vec.r;
            dr(j)  	= (r(j+1)-r(j-1))/2;
            dr(1) 	= (r(2)-r(1))/2;
            dr(n2)	= (r(n2)-r(n2-1))/2; 
            ind0                    = -find(diff(ind)~=0,1,'first');
            para.cont1(i).vec.xc    = circshift(para.cont1(i).vec.xc,[ind0,0]);
            para.cont1(i).vec.zc	= circshift(para.cont1(i).vec.zc,[ind0,0]);
            para.cont1(i).vec.vnx   = circshift(para.cont1(i).vec.vnx,[ind0,0]);
            para.cont1(i).vec.vnz   = circshift(para.cont1(i).vec.vnz,[ind0,0]);
            para.cont1(i).vec.cv    = circshift(para.cont1(i).vec.cv,[ind0,0]);
            %             para.cont1(i).vec.r   	= cumsum((cont1.vec.xc(1:end-1)-cont1.vec.xc(2:end)).^2+(cont1.vec.zc(1:end-1)-cont1.vec.zc(2:end)).^2);
            if size(dr)~=size(para.cont1(i).vec.cv)
                dr=dr.';
            end
            dr                      = circshift(dr,[ind0,0]);
            ind                   	= circshift(ind,[ind0,0]);
            dr(ind)                 = [];
            para.cont1(i).vec.r   	= cumsum(dr);
            
        elseif nch>2
            error('perhaps several contours should be considered instead')
            %a programmer
        end
        para.cont1(i).vec.xc(ind)   = [];
        para.cont1(i).vec.zc(ind)   = [];
        para.cont1(i).vec.vnx(ind)  = [];
        para.cont1(i).vec.vnz(ind)  = [];
        para.cont1(i).vec.cv(ind)   = [];
    end
end