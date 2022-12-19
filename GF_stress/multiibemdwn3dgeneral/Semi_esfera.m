function [nbp]=Semi_esfera(n,npt)
%----------------------------------------------------------------------
%     rectificacion aproximada de una semiesfera
%----------------------------------------------------------------------

x   = zeros(npt,1);
y   = zeros(npt,1);
z   = zeros(npt,1);
vnx = zeros(npt,1);
vny = zeros(npt,1);
vnz = zeros(npt,1);
a   = zeros(npt,1);
r   = zeros(npt,1);

dtet= pi/(2*n+1);
r1  = sin(dtet/2.0);
r(1)= r1;
npn = 1;

a(1)    = pi*r1*r1;
x(1)    = 0.0;
y(1)    = 0.0;
z(1)    = 1.0;
vnx(1)  = 0.0;
vny(1)  = 0.0;
vnz(1)  =-1.0;

% atot=a(1);
nac=npn;
nfac=2.5*n;%?
for j=2:nfac
    
    if j<=n+1
        b=sin(dtet*(j-1));
    else
        b=1.0+(2*(j-n-1)+1)*r1;
    end
    
    %       npn=2*int(b/r1+0.5)
    npn=4*floor(0.5*b/r1+0.25);
    
    ntest=npn;
    
    if ntest<4
        npn=4;
    end
    
    aj=2.0*pi*b*(2.0*r1)/npn;
    rj=sqrt(aj/pi);
    
    if j<=n+1
        zj=cos(dtet*(j-1));
    else
        zj=0.0;
    end
    
    for l=1:npn
        nac=nac+1;
        %         atot=atot+aj;
        a(nac)=aj;
        r(nac)=rj;
        z(nac)=zj;
        fi=2.0*pi*(l-1)/npn;
        x(nac)=b*cos(fi);
        y(nac)=b*sin(fi);
        if j<n+1
            vnz(nac)=-zj;
            vnx(nac)=-b*cos(fi);
            vny(nac)=-b*sin(fi);
        else
            vnz(nac)=-1.0;
            vnx(nac)= 0.0;
            vny(nac)= 0.0;
        end
    end
end
nbp=nac;

