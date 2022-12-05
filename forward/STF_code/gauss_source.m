close all;
clear all;

% Print a Gaussian source with an specific
% - rise time
% - origin time

% Origin time
% t0=0.25;
t0=1.0;
% Rise time
trr=0.5;
%trr=1.0;


tr=trr/4.0;
% Delta time
dt=0.001;
% Time series
t=(0:0.001:3.0);
% Computation of the gaussian
for i=1:length(t)
   gauss(i)=1/(tr*sqrt(pi))*exp(-((i-1)*dt-t0)^2/tr^2); 
end

gausst=load('slip-rate');
gauss=gausst(:,2);
dt=gausst(2,1);

% Get the spectrum
npts=length(gauss);
sgauss = dt*abs(fft(gauss,npts));
% Generate associated frequencies
f = (0:npts-1)/(dt*npts);

% Plot the Gaussian
subplot(1,2,1)
plot(gausst(:,1),gausst(:,2),'LineWidth',3,'color','r');
xlabel('Time (sec)');
%axis([0 5 -2 5])
grid on
title ('Gauss source rt = 0.05');

% Plot the spectrum
subplot(1,2,2)
%loglog(f,sgauss,'LineWidth',1,'color','r');
%semilogx(f,sgauss,'LineWidth',1,'color','r');
plot(f,sgauss,'LineWidth',3,'color','r');
xlabel('Frequency (Hz)');
axis([0 10 0 1])
grid on
title ('Spectrum');

% Print the gaussian
myfile = fopen('gauss_rt0p5.src' ,'w');
for ii=1:length(gauss)
%    fprintf(myfile, '%f %E\n',t(ii), gauss(ii));
    fprintf(myfile, '%E\n', gauss(ii));
end
fclose(myfile);
