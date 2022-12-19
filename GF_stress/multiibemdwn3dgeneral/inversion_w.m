function [utc,stc]=inversion_w(uw,sw,para)
% utc : desplazamiento en tiempo corregido con la respuesta temporal de la
%       fuente
if para.zeropad < para.nf; para.zeropad = para.nf; end
nf      = para.nf;       disp(['nf = ',num2str(nf)])
nfN     = nf/2+1; 
para.df = para.fmax/(nf/2);     %paso en frecuencia
df      = para.fmax/nfN; disp(['df = ',num2str(df)])
% zeropad = para.zeropad;
% tps     = 0:(1/(df*2*(nfN+zeropad))*(2*(nfN+zeropad)/(2*(nfN+zeropad)-2))):1/df;
dt = 1/(df*para.zeropad);    disp(['dt = ',num2str(dt)])
tps = (0:para.zeropad-1)*dt; disp(['tmx= ',num2str(tps(end))])

nuw     = size(uw);
nsw     = size(sw);

% correction espectros %
cspectre    =correction_spectre(para,nfN,df);

%wrap de las dimensiones
ntot=1;
nnuw=length(nuw);
for i=2:length(nuw)
    ntot=ntot*nuw(i);
end

uw      = reshape(uw,nuw(1),ntot);
utc     = zeros(para.zeropad,ntot);

for i=1:ntot
    tmp     = uw(1:nfN,i).';
    tmp     = tmp.*cspectre;
   
   % crepa:
   vec = zeros(1,para.zeropad);
   vec(1:nfN) = tmp(1:nfN); 
   vec(end-nfN+3:end)=conj(tmp(nfN-1:-1:2));
   % escala:
  dt_nopad = 1/df/(nf-1); % dt sin zeropading
  sca = dt_nopad*para.zeropad/nf/nf;
  utc(:,i)= real(sca*ifft(vec)).*exp(para.DWNomei*tps); % inversa y frec imag
end
utc   = reshape(utc,[para.zeropad,nuw(2:nnuw)]);

ntots=1;
nnsw=length(nsw);
for i2=2:length(nsw)
    ntots=ntots*nsw(i2);
end
sw      = reshape(sw,nsw(1),ntots);
stc     = zeros(para.zeropad,ntot);

for i2=1:ntots
  tmp     = sw(1:nfN,i2).';
  tmp     = tmp.*cspectre;
  vec = zeros(1,para.zeropad);
  vec(1:nfN) = tmp(1:nfN); 
  vec(end-nfN+3:end)=conj(tmp(nfN-1:-1:2));
   % escala:
  dt_nopad = 1/df/(nf-1); % dt sin zeropading
  sca = dt_nopad*para.zeropad/nf/nf;
  stc(:,i2)= real(sca*ifft(vec)).*exp(para.DWNomei*tps); % inversa y frec imag
end
stc   = reshape(stc,[para.zeropad,nsw(2:nnsw)]);