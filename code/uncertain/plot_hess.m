clear all
close all
clc

file=fopen('hess.full','r')
hess=fread(file,'single');
fclose(file)

n=24*12*57;
hr = reshape(hess,n,n);

pcolor(hr(1:250,1:250))
title('Hessian matrix')
set(gca,'fontsize',45)
pause


[nparam,~] = size(hess);
nparam = sqrt(nparam);

%k = 1;
%for i=1:nparam
% hess_r(:,i) = hess(k:k+nparam-1,1);
% k = k+nparam;
%end


hesscol=load('fort.26');
k = 1;
for i=1:12
 for j=1:24
  he(i,j,1:57) = hesscol(k:k+56);
  k=k+57;
 end
end
subplot(221)
p(1:12,1:24)=he(:,:,20);
maxp=max(max(p));
colormap(jet),pcolor(p),shading interp,caxis([-maxp,maxp])
subplot(222)
p(1:12,1:24)=he(:,:,18);
colormap(jet),pcolor(p),shading interp,caxis([-maxp,maxp])
subplot(223)
p(1:12,1:24)=he(:,:,22);
colormap(jet),pcolor(p),shading interp,caxis([-maxp,maxp])

for i=15:25
p(1:12,1:24)=he(:,:,i);
subplot(224)
colormap(jet),pcolor(p),shading interp,caxis([-maxp,maxp])
%pause
%plot(10,6,'bp','markersize',24,'linewidth',10)
end

