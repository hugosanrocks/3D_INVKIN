clear all
close all
clc

nsta=17;
ncomp=6;

freq = [0.25 0.5 1.0 2.0];
order = 4;
dthe = 0.015619;
the = 0:dthe:4999*dthe;
%dtho = ;
dtho = dthe;

%Half-space file
fis = load('res_all_hspace.mat');

count = 1;
for ifile=1:17
    for j=1:3
      resi(count) = (j-1)*17 + ifile;
      count = count + 1;
    end
end

nsamp = 5000;
tseries = zeros(nsamp,1);
tseries = zeros(nsamp,1);


station = input('Which station to plot:')
subfault = input('Which subfault to plot:')
tmax = input('Maximum time to plot:');


for ista=17:1:17
  station = ista;
  filhe = sprintf('res%02d.mat',station);
  reshe = load(filhe);
  %  filho = sprintf('sol%02d.mat',ista);
  %  resho = load(filho);
  counter = (station - 1)*3 + 1;
  for cont = counter:counter%+2
   resho = fis.sol(:,:,resi(cont),:);
   resi(cont)
  end

  k = 1;
  for i=1:ncomp

   tseries = reshe.sol(:,subfault,1,i);
   tseries2 = resho(:,subfault,1,i);

   for ifil=1:4

     freqfil = freq(ifil);
     ratbhe=freqfil/(1/(2*dthe));
     ratbho=freqfil/(1/(2*dtho));
     [ahe,bhe]=butter(order,ratbhe);
     [aho,bho]=butter(order,ratbho);

     filhe = filtfilt(ahe,bhe,tseries);
     maxplot1 = max(abs(filhe));
     filho = filtfilt(aho,bho,tseries2);
     maxplot2 = max(abs(filho));
     maxamp = max([maxplot1 maxplot2]);
     figure(k)
     text = sprintf('Frequency 0 - %02.2f [Hz]',freqfil);
     subplot(2,2,ifil),plot(the,filhe),hold on
     set(gca,'FontSize',12),title(text)
     subplot(2,2,ifil),plot(the,filho,'-r','LineWidth',0.5)
     xlim([0,tmax]),ylim([-maxamp,maxamp])
     if ( ifil == 2 )
       legend('hete','hspace')
     end

   end
  figg = sprintf('green_comp_sta%02d_sub%03d_comp%02d.eps',ista,subfault,k);
  print(figg,'-depsc');
  k = k + 1;
  end
close all
end
