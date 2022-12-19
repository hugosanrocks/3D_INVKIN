
uw=zeros(para.nf+1,para.rec.nrec);


para0=normalizacion(para);

df  = para.fmax/(para.nf/2);     %paso en frecuencia
para0.df=df;

para0=tmppatch(para0); 

for j=0:para.nf/2
    
    fj=j*df;
    if j==0
        fj=0.1*df;
    end
    kE      = 2*pi*fj;                         %w/betE mais normalise
    ki     = kE*complex(1,-1/(2*para0.regE.qd));
    
    i=1:para0.rec.nrec;
    x=para0.rec.xri+(i-1)*para0.rec.dxr;
    z=0;
    rij=sqrt((x-para0.xs).^2+(z-para0.zs).^2);
    uw(j+1,:)=2*G22_SH(ki,rij,mu);
end

%%%%%%%%%%%%%%%
% inversion w %
%%%%%%%%%%%%%%%
para0.ninc=1;
para0.fuente=2;

utc    = inversion_w(uw,para0);

%%%%%%%%%%%%%%%%%%%%%
% dibujo resultados %
%%%%%%%%%%%%%%%%%%%%%
para0   = pos_rec(para0);
b_dib=dibujo(para0,bdessin,utc,uw,stc,sw,cont1);