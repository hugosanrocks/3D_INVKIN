%% variables del modelo:
%cd '/Users/marshall/Documents/DOC/SScode/AXITRA/A X I T R A - FORCE/SINGLEFORCEZIP'

RESULT=result;
dfreq = 7.81983137E-03;
omei = 2*pi/127.88;
nt = 2048;
nr = 4;

tps2=0:para.tmax/8192:para.tmax-(para.tmax/8192);
xi = [60,75,0];
x1 = [ 32.5 40.0868241 2.53798062];
d  = x1-xi;
th(1) = atan(d(2)/d(1));
thd= th(1)*180/pi;

x2 = [59.5 40.0868241 2.53798062];
d = x2-xi;
th(2) = atan(d(2)/d(1));

x3 = [61.0 40.0868241 2.53798062];
d = x3-xi;
th(3) = atan(d(2)/d(1));

x4 = [ 67.0 42.9520190 18.7873085];
d = x4-xi;
th(4) = atan(d(2)/d(1));


% leer funcion de amplitud que fue usada por source.f
aux = importdata('funcAmp.txt');
fA = aux(:,1) + aux(:,2)*1j;
nfreq = length(fA);
Fq = (0:nfreq)*dfreq; Fq(end)=[];
dt = 1/2/Fq(end);
disp(['Fmax =' num2str(Fq(end))])
tps = (0:nt)*dt; tps(end)=[];
% graficar función de amplitud:
% figure(10242); clf; hold on
% plot(Fq,real(fA),'r-')
% plot(Fq,imag(fA),'b-')
% plot(Fq,abs(fA),'k-')
% xlabel('frecuencia Hertz')
% xlim([0 50])
clear aux % no se usa.

% leer sismogramas
aux = importdata('axi.x');
resR = zeros(nt,nr);
resX= zeros(nt,nr);
auy = importdata('axi.y');
resT = zeros(nt,nr);
resY = zeros(nt,nr);
auz = importdata('axi.z');
resZ = zeros(nt,nr);
ini = 1;
fin = nt;
%sca = 1;
sca = dt/nfreq/sqrt(nt*2); % escala

for ir=1:nr
	resR(:,ir)=aux(ini:fin);%.'*sca;
	resT(:,ir)=auy(ini:fin);%.'*sca;
	resZ(:,ir)=auz(ini:fin);%.'*sca;
	fin=fin+nt;
	ini=ini+nt;

% rotar a cartesianas
resX(:,ir) = cos(th(ir))*resR(:,ir)+sin(th(ir))*resT(:,ir);
resY(:,ir) = sin(th(ir))*resR(:,ir)-cos(th(ir))*resT(:,ir);

% hacer gráficas
figure(8000+ir);clf;
set(gcf,'name',['Sismogramas receptor ' num2str(ir)])
subplot(3,1,1);hold on
plot(tps2,RESULT.utc(:,ir,1,1),'r-'); % DWN matriz global
plot(tps,resX(:,ir),'k--');
legend({'DWN matglob','AXITRA'})
ylabel('ux')
xlim([0 20])
subplot(3,1,2);hold on
plot(tps2,RESULT.utc(:,ir,1,2),'r-'); % DWN matriz global
plot(tps,resY(:,ir),'k--');
ylabel('uy')
xlim([0 20])
subplot(3,1,3);hold on
plot(tps2,RESULT.utc(:,ir,1,3),'r-'); % DWN matriz global
plot(tps,resZ(:,ir),'k--');
ylabel('uz')
xlim([0 20])
xlabel('tiempo (segundos)')
end

%cd '/Users/marshall/Documents/DOC/coco/ibem_matlab'
