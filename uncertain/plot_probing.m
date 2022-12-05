clear all
close all
clc

src = load('source.src');
for i=1:875
 p=1;
  for j=1:12
   for k=1:24
    srcr(j,k,i) = src(i,p);
    p=p+1;
   end
  end
end
file=fopen('cube_src','w');
fwrite(file,srcr,'single');
fclose(file);

prob = load('Corr_time.out');
prob_stk = load('Corr_strike.out');
prob_dip = load('Corr_dip.out');

sampw = 25;  %n+n-1
sampwx = 7;  %n+n-1
sampwz = 7;  %n+n-1
hsampx = 4;
hsampz = 4;

winsf = 4;  %each node
winzt = 6;  %each strike row
winxt = 3;  %each dip column
sf = 288;
sf_s = 24;
sf_d = 12;
sf_t = 57;
realiz = 50;
dt = 0.25;
dx = 1.5;
dz = 1.5;

t = 0:dt:dt*(sampw)-dt;
x = 0:dx:dx*(sampwx)-dx;
z = 0:dz:dz*(sampwz)-dz;

nsamp = sampw*winsf*sf;  %total samples time corr
nsampx = sampwx*winzt*sf_d*sf_t;  %total samples strike corr
nsampz = sampwz*winxt*sf_s*sf_t;  %total samples dip corr

prob_r = reshape(prob,nsamp,realiz);
probstk_r = reshape(prob_stk,nsampx,realiz);
probdip_r = reshape(prob_dip,nsampz,realiz);

for i=1:nsamp
  prob_mean(i) = mean(prob_r(i,:));
end
display('Time done')
for i=1:nsampx
  probstk_mean(i) = mean(probstk_r(i,:));
end
display('Strike done')
for i=1:nsampz
  probdip_mean(i) = mean(probdip_r(i,:));
end
display('Dip done')

totlen = nsamp;
tot_win = totlen/sampw;
tot_winx = nsampx/sampwx;
tot_winz = nsampz/sampwz;

p=1;
for i=1:tot_win
  serie = prob_mean(p:p+sampw-1);
  seriemax(i) = max(serie);
  seriehalf(i) = seriemax(i)/2;
  val1(i)=interp1(serie(1:13),t(1:13),seriehalf(i));
  val2(i)=interp1(serie(13:end),t(13:end),seriehalf(i));
  L(i) = val2(i) - val1(i);
  sigma(i) = L(i)/(2*sqrt(2*log(2)));
  p=p+sampw;
end

p=1;
for i=1:tot_winx
  seriex = probstk_mean(p:p+sampwx-1);
  seriexmax(i) = max(seriex);
  seriexhalf(i) = seriexmax(i)/2;
  val1x(i)=interp1(seriex(1:hsampx),x(1:hsampx),seriexhalf(i));
  val2x(i)=interp1(seriex(hsampx:end),x(hsampx:end),seriexhalf(i));
  Lx(i) = val2x(i) - val1x(i);
  sigmax(i) = Lx(i)/(2*sqrt(2*log(2)));
  p=p+sampwx;
end

p=1;
for i=1:tot_winz
  seriez = probdip_mean(p:p+sampwz-1);
  seriezmax(i) = max(seriez);
  seriezhalf(i) = seriezmax(i)/2;
  val1z(i)=interp1(seriez(1:hsampz),z(1:hsampz),seriezhalf(i));
  val2z(i)=interp1(seriez(hsampz:end),z(hsampz:end),seriezhalf(i));
  Lz(i) = val2z(i) - val1z(i);
  sigmaz(i) = Lz(i)/(2*sqrt(2*log(2)));
  p=p+sampwz;
end

%time
%Arrange for cube plotting
cont = 1;
for k=1:sf_d
  for j=1:sf_s
   for i=1:winsf
    sigma_sf(j,k,i) = sigma(cont);
    cont = cont + 1;
   end
  end
end

%strike
%Arrange for cube plotting
cont = 1;
for k=1:sf_t
  for j=1:sf_d
   for i=1:winzt
    sigma_stk(i,j,k) = sigmax(cont);
    cont = cont + 1;
   end
  end
end

%dip
%Arrange for cube plotting
cont = 1;
for k=1:sf_t
  for j=1:sf_s
   for i=1:winxt
    sigma_dip(j,i,k) = sigmaz(cont);
    cont = cont + 1;
   end
  end
end



file=fopen('cube_temp','w');
fwrite(file,sigma_sf,'single');
fclose(file);

file=fopen('cube_stk','w');
fwrite(file,sigma_stk,'single');
fclose(file);

file=fopen('cube_dip','w');
fwrite(file,sigma_dip,'single');
fclose(file);


p=1;
max_y = max(max(prob_r(p:p+sampw,1:10)));
min_y = min(min(prob_r(p:p+sampw,1:10)));

xmax=t;
ymax=ones(sampw,1)*seriemax(1);
yhalf=ones(sampw,1)*seriehalf(1);

yval1 = -1:1;
xval1 = ones(3,1)*val1(1);
yval2 = -1:1;
xval2 = ones(3,1)*val2(1);


p=1;
for i=1:10
h1=  plot(t,prob_r(p:p+24,i),'b'),hold on,xlim([t(1),t(sampw)]);
end
p=1;
h2=plot(t,prob_mean(p:p+24),'*-r'),ylim([min_y,max_y]);
h3=plot(xmax,ymax,'k','linewidth',1.5);
h4=plot(xmax,yhalf,'r','linewidth',1.5);
legend([h1 h2 h3 h4],'realizations','mean','max of mean','FWHM')
plot(xval1,yval2,'r','linewidth',1.5)
plot(xval2,yval2,'r','linewidth',1.5)
set(gca,'fontsize',20)
title('First autocorrelation along time')
xlabel('Time [s]'),ylabel('Autocorrelation')

p=1;
max_y = max(max(probstk_r(p:p+6,1:10)));
min_y = min(min(probstk_r(p:p+6,1:10)));


x=0:1.5:7*1.5-1.5;
xmax=x;
ymax=ones(sampwx,1)*seriexmax(1);
yhalf=ones(sampwx,1)*seriexhalf(1);

yval1 = -1:1;
xval1 = ones(3,1)*val1x(1);
yval2 = -1:1;
xval2 = ones(3,1)*val2x(1);


p=1;
figure(2)
for i=1:10
h1=  plot(x,probstk_r(p:p+6,i),'b'),hold on,xlim([x(1),x(7)]);
end
p=1;
h2=plot(x,probstk_mean(p:p+6),'*-r'),ylim([min_y,max_y]);
h3=plot(xmax,ymax,'k','linewidth',1.5);
h4=plot(xmax,yhalf,'r','linewidth',1.5);
legend([h1 h2 h3 h4],'realizations','mean','max of mean','FWHM')
plot(xval1,yval2,'r','linewidth',1.5)
plot(xval2,yval2,'r','linewidth',1.5)
set(gca,'fontsize',20)
title('First autocorrelation along strike')
xlabel('Time [s]'),ylabel('Autocorrelation')

