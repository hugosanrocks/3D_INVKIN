figure(154);
hold on;
coloooor={'.-','.-r'};
im=[];
for m=1:para.nmed0
    im0 = max(coord.Xim(:,m)~=0);
    if im0==1
        im = [im m];
    end
end

for m=im
    if m>1
        x1  =[];
        phi1=[];
        for c=1:2
            tmp=(coord.Xim(:,m)==c);
            ind=logical(coord.indm(m).ind.*tmp);
            drj         = coord.dr(ind);
            [x0,ind0]=sort(coord.x(ind));
            phi=abs(phi_fv(coord.phi(ind,m)));
            phi=phi(ind0).';
            if c==2
                x0  = fliplr(x0);
                phi = fliplr(phi);
            end
            x1  =[x1,x0];
            phi1=[phi1,phi];
        end
    else
        ii          = coord.indm(m).ind;
        [x1,ind1]   = sort(coord.x(ii));
        ii          = ii(ind1);
        drj         = coord.dr(ii);
        phi1        = abs(phi_fv(coord.phi(ii,m))).';
    end
    r = cumsum(drj);
    plot(r,phi1,'r.');
%     plot(r,drj,'.');
    
    %construction de fuerzas iguales
    df              = mean(phi1.*drj);
    n               = length(r);
    long            = r(n);
    %mailla 100 veces mas fina
    nf              = 100*n;
    drfin           = long/nf;
    rfin            = linspace(drfin,long,nf);
    i               = 2:(nf-1);
    drv(i)       	= (rfin(i+1)-rfin(i-1))/2;
    drv(1)        	= (rfin(2)-rfin(1))/2;
    drv(nf)         = (rfin(nf)-rfin(nf-1))/2;
    phifin          = interp1(r,phi1,rfin,'spline');
    plot(rfin,phifin,'.','markersize',2);
    j=0;
    r2=zeros(1,n);
    phi2=zeros(1,n);
    i=0;
    while i <n && j+1<nf
        i   = i+1;
        j   = j+1;
        j0  = j;
        dfc = cumsum(phifin(j0).*drv(j0));
        while dfc<df && j<nf
            j   = j+1;
            dfc = sum(phifin(j0:j).*drv(j0:j));
        end
        r2(i)   = (rfin(j0)+rfin(j))/2;
        phi2(i) = dfc/sum(drv(j0:j));
    end
    r2(i+1:end)=[];
    phi2(i+1:end)=[];
    plot(r2,phi2,'*k');
end