function ind=chercherzero(y)
%pour le calculs de courbes de dispersion
%abandonné pour cherchemin
n=length(y);
ind2 =find((y(1:n-1).*y(2:n)<=0));
ind=[];
iind=1;
for ii=1:length(ind2)
    % attention, il y a de faux 0 dus aux points singulier de tan,
    % en realité il n y a pas de passage par 0
    % on verifie que la dérivé ne change pas de signe autour des pts
    % entourant les probables zeros

    if ind2(ii)+2<=n && ind2(ii)-1>0
        s=0;
        for is=0:1
            if (y(ind2(ii)+is)-y(ind2(ii)-1+is))*(y(ind2(ii)+1+is)-y(ind2(ii)+is))<0
                s=s+1;
            end
        end
        if s<2
            ind(iind)=ind2(ii);
            iind=iind+1;
        end
    end
end