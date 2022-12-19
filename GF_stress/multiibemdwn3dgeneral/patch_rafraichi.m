% para.recpos=1;
% para.rec.nrec=1;
% para.rec.nrecx=para.rec.nrec;
% para.rec.yri=0;
% para.rec.dyr=0.1;
% para.rec.zri=0;
% para.rec.dzr=0.1;

% para0=para;

para.tipo_onda=ones(1,para.ninc);

for i=1:para.nmed,para.reg(i).tipoatts=1;end
for i=1:para.nmed,para.cont(i,1).th=0;end
for i=1:para.nmed;para.cont(i,1).ruggeo=1;para.cont(i,2).ruggeo=1;end
for i=1:para.nmed;para.cont(i,1).rba=1;para.cont(i,2).rba=1;end
for i=1:para.nmed;para.cont(i,1).rh=1;para.cont(i,2).rh=1;end
for i=2:para.nmed;para.geo(i)=1;end;para.geo(1)=2;

para.cont(2,1).xa   =-1;
para.cont(2,1).za   = 0;
para.cont(1,1).xa   =-3;
para.cont(1,1).a    = 6;

for i=1:para.nmed;para.cont(i,2).xa=para.cont(i,1).xa;end
for i=1:para.nmed;para.cont(i,2).a=para.cont(i,1).a;end
for i=1:para.nmed;para.cont(i,2).za=para.cont(i,1).za;end

para.sortie.Ux	= 1;
para.sortie.Uy  = 1;
para.sortie.Uz  = 1;

para.sortie.Ut  = 1;
para.sortie.UP  = 0;
para.sortie.US  = 0;
para.sortie.UI  = 0;


para.DWNkmax=1;
para.DWNnbptkx=100;
para.DWNomei=0.0001;
