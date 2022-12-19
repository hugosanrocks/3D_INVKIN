nth     = 100;
para.ninc=nth;
th      = linspace(90,.1,nth)*pi/180;
kx      = cos(th);
kx      = fliplr(2 -kx);
para.kxk=kx;
para.zs(1:para.ninc)    = 0;
para.xs(1:para.ninc)    =-1.5;
para.gam(1:para.ninc)	= 90;

nth     = 100;
para.ninc=nth;
th      = linspace(90,.1,nth)*pi/180;
kx      = cos(th);
kx      = -fliplr(2 -kx);
para.kxk=kx;
para.zs(1:para.ninc) = 0;
para.xs(1:para.ninc) = 1.5;
para.gam(1:para.ninc) = 90;

para.zs(1:para.ninc) = 1;
para.xs(1:para.ninc) = 1.5;
para.gam(1:para.ninc) = 90;



para.gam = linspace(0,90,para.ninc+1);
para.gam = para.gam(1:para.ninc);

para.gam    = linspace(-90,90,para.ninc+1);
para.gam    = para.gam(1:para.ninc);
para.kzsigno=ones(1:para.ninc,1);
para.zs(1:para.ninc) = 1;
para.xs(1:para.ninc/2+1) = 1.5;
para.xs(para.ninc/2+2:para.ninc) = -1.5;

para.ninc                       = 800;
para.kxk                         = linspace(-2,2,para.ninc);
para.zs(1:para.ninc)            = 0;
para.xs(1:para.ninc/2)          = 1.5;
para.xs((1+para.ninc/2):para.ninc)=-1.5;
para.kzsigno(1:para.ninc)       = 1;
para.gam=asin(para.kxk)*180/pi;
para.gam=real(para.gam);