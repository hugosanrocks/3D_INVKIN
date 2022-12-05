clear all
close all
clc

nsta = 40;
ncomp = 3;

stacomp = nsta*ncomp;

for k=1:nsta
   e=sprintf('dat/obs_S%03d_C1',k);
   n=sprintf('dat/obs_S%03d_C2',k);
   v=sprintf('dat/obs_S%03d_C3',k);
   eobs=load(e);
   nobs=load(n);
   vobs=load(v);

   e2=eobs.^2;
   n2=nobs.^2;
   v2=vobs.^2;

   sume2 = sum(e2);
   sumn2 = sum(n2);
   sumv2 = sum(v2);

   ce(k) = sume2;
   cn(k) = sumn2;
   cv(k) = sumv2;
end

sumtotal = (sum(ce) + sum(cn) + sum(cv))*(0.5);
contrib = (sumtotal*2)/stacomp;
contrib = 1 / stacomp;

ce_norm = ce.^(-1)*contrib;
ce_norm = ce_norm.^(1/2);
cn_norm = cn.^(-1)*contrib;
cn_norm = cn_norm.^(1/2);
cv_norm = cv.^(-1)*contrib;
cv_norm = cv_norm.^(1/2);


for k=1:nsta
   e=sprintf('dat/obs_S%03d_C1',k);
   n=sprintf('dat/obs_S%03d_C2',k);
   v=sprintf('dat/obs_S%03d_C3',k);
   eobs=load(e);
   nobs=load(n);
   vobs=load(v);

   e2=ce_norm(k)^2*eobs.^2;
   n2=cn_norm(k)^2*nobs.^2;
   v2=cv_norm(k)^2*vobs.^2;

   sume2 = sum(e2);
   sumn2 = sum(n2);
   sumv2 = sum(v2);

   ce(k) = sume2;
   cn(k) = sumn2;
   cv(k) = sumv2;
end

sumtotal = (sum(ce) + sum(cn) + sum(cv))*(0.5);

weights = [ce' cn' cv'];

save('-ascii','weights.dat','weights')

