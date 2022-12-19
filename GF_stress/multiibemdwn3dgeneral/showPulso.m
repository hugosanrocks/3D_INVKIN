function showPulso( para )
%showPulso Mostrar la señal del pulso en la fuente
%   La señal se muestra en frecuencia y tiempo así como los valores máximos
%   calcualdos 
if para.zeropad < para.nf; para.zeropad = para.nf; end

nf      = para.nf;           disp(['nf = ',num2str(nf)])
nfN     = nf/2+1; 
df      = para.fmax/nfN;     disp(['df = ',num2str(df)])
Fq      = (0:nf/2)*df;       disp(['Fmx= ',num2str(Fq(end))])
dt = 1/(df*para.zeropad);    disp(['dt = ',num2str(dt)])
tps = (0:para.zeropad-1)*dt; disp(['tmx= ',num2str(tps(end))])

cspectre  = correction_spectre(para,nfN,df);

% crepa:
vec = zeros(1,para.zeropad);
vec(1:nfN) = cspectre(1:nfN); 
vec(end-nfN+3:end)=conj(cspectre(nfN-1:-1:2));
% escala:
  dt_nopad = 1/df/(nf-1); % dt sin zeropading
  sca = dt_nopad*para.zeropad/nf/nf;
  signal= real(sca*ifft(vec)); % inversa 

% graficar
figure(1);
set(gcf,'name','amplitud del pulso en el origen')
subplot(1,1,1)
cla
plot(Fq,real(cspectre.'),'r');hold on;
plot(Fq,imag(cspectre.'),'b');
plot(Fq,abs(cspectre.'),'k')
ax1 = gca;
ax1_pos = ax1.Position;
ax1.Box = 'off';
xlabel('frecuencia en Hertz')
ax2=axes('Position',ax1_pos,...
    'XAxisLocation','top',...
    'YAxisLocation','right',...
    'Color','none');
xlim(ax1.XLim)
nn = length(ax1.XTick);
ax2.XTickMode = 'auto';
ax2.XTickLabelMode = 'manual';
cl = cell(nn,1);
cl{1} = ' ';
for i=2:nn
    cl{i} = num2str(1/ax1.XTick(i),3);
end
ax2.XTickLabel = cl;
set(gca,'ytick',[])
title('Espectro de pulso')
xlabel('Periodo en segundos')
%hugo added this to plot separetly
figure(2),subplot(1,1,1)
cla
plot(tps,real(signal.'),'r');hold on;
plot(tps,imag(signal.'),'b');
title('Ondícula')
xlabel('tiempo en segundos')

if para.pulso.tipo==6 % butterworth
Ha = abs(cspectre);
% Hr = real(cspectre);
% Hi = imag(cspectre);
Hadb  = 20*log10(Ha/1); % Convert to dB scale
% Hrdb  = 20*log10(Hr/1); % Convert to dB scale
% Hidb  = 20*log10(Hi/1); % Convert to dB scale
figure(3333); hold on
% plot(Fq, Hrdb,'r-','LineWidth',1)
% plot(Fq, Hidb,'b-','LineWidth',1)
plot(Fq, Hadb,'k-','LineWidth',3,'DisplayName',...
  ['abs n' num2str(para.pulso.a) ' Wn' num2str(para.pulso.b)])
xlabel('Frequencia (Hertz)')
ylabel('Magnitud (db)')
grid on
end
end

