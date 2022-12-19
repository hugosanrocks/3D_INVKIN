function dibujo_allSpec(para,bouton,RESULT,xind,xiind,tipolinea)
f = gcf;%figure;%(3456);
set(f,'position',[24  1 1246 726]);
uw = RESULT.uw;
sw = RESULT.sw;
df = para.fmax/(para.nf/2); %paso en frecuencia
Fq = (0:para.nf/2)*df;
nfN = para.nf/2+1;
cs    =correction_spectre(para,nfN,df);%pour supprimer la convolution temporelle de la source

if nargin > 3
r = xind;
x = [para.rec.xr(r) para.rec.yr(r) para.rec.zr(r)];
% r = x(4);
else
r = 1;
x = [para.rec.xr(r) para.rec.yr(r) para.rec.zr(r)];
end

if nargin > 4
% i = xi(4);
i = xiind;
xi = [para.xs(i) para.ys(i) para.zs(i)];
else
i = get(bouton.inc,'value');
xi = [para.xs(i) para.ys(i) para.zs(i)];
end

tx = ['x(',num2str(x(1)),...
  ',',num2str(x(2)),...
  ',',num2str(x(3)),')'...
  ' xi(',num2str(xi(1)),...
  ',',num2str(xi(2)),...
  ',',num2str(xi(3)),')'];

tmpc    =get(bouton.couleur,'string');
coul    =get(bouton.couleur,'value');
if nargin == 6
c =[tmpc{coul} tipolinea];
else
c = tmpc{coul};
end
ha = zeros(1,12);

subplot(4,3,1);  tT('u',uw(1:nfN,r,i,1),tx,c,Fq,cs); ha(1)=gca;
subplot(4,3,2);  tT('v',uw(1:nfN,r,i,2),tx,c,Fq,cs); ha(2)=gca;
subplot(4,3,3);  tT('w',uw(1:nfN,r,i,3),tx,c,Fq,cs); ha(3)=gca;
if size(sw,4) ~= 0
subplot(4,3,4);  tT('sxx',sw(1:nfN,r,i,1),tx,c,Fq,cs); ha(4)=gca;
subplot(4,3,5);  tT('sxy',sw(1:nfN,r,i,4),tx,c,Fq,cs); ha(5)=gca;
subplot(4,3,6);  tT('sxz',sw(1:nfN,r,i,5),tx,c,Fq,cs); ha(6)=gca;
% subplot(4,3,7);  title('syx'); ha(7)=gca; set(ha(7),'visible','off')
subplot(4,3,8);  tT('syy',sw(1:nfN,r,i,2),tx,c,Fq,cs); ha(8)=gca;
subplot(4,3,9);  tT('syz',sw(1:nfN,r,i,6),tx,c,Fq,cs); ha(9)=gca;
% subplot(4,3,10); title('szx'); ha(10)=gca; set(ha(10),'visible','off')
% subplot(4,3,11); title('szy'); ha(11)=gca; set(ha(11),'visible','off')
subplot(4,3,12); tT('szz',sw(1:nfN,r,i,3),tx,c,Fq,cs); ha(12)=gca;
end
linkaxes(ha,'x');

subplot(4,3,11); tT('Legend',NaN,tx,c,Fq,cs);
% h = legend; 
% if ~isempty(h)
% hS = h.String; 
% else
% hS = '';
% delete(h);
% end
legend('show'); 
LEG=legend('-DynamicLegend');
set(gca,'xtick',[]); set(gca,'ytick',[])
xlabel('');ylabel(''); set(LEG,'FontSize',15);
set(LEG,'position',[4.1573e-01 1.2121e-01 2.0546e-01 1.3774e-01]);
clear LEG
if nargin <= 4 % copiar la figura desde la pantalla principal
hg = subplot(4,3,[7 10]);
h=allchild(bouton.axe_conf_geo);
copyobj(h,hg);
% darle formato
if(para.dim>=3)
  set(gca,'zdir','reverse');
  view(3)
  xlabel('');ylabel('');zlabel('')
  grid on
  set(gca,'dataaspectratio',[1 1 1],'Projection','perspective','Box','off');
  light('Position',[2*min(get(gca,'xlim')) ...
                    2*max(get(gca,'ylim')) ...
                    5*min(get(gca,'zlim'))],'Style','infinite');
  axis tight
  set(gca, 'XColor', 'r');set(gca, 'YColor', 'g');set(gca, 'ZColor', 'b')
else
  set(gca,'ydir','reverse');
  view(2)
  xlabel('X');ylabel('Z');zlabel('');
  grid off
  set(gca,'dataaspectratio',[1 1 1],'Projection','orthographic','Box','on');
  axis image
  set(gca, 'XColor', 'k');set(gca, 'YColor', 'k');set(gca, 'ZColor', 'k')
end
end
end

function tT(tit,var,tx,c,Fq,cs)
if max(max(var)) == 0; return;end

espe = var.*cs.';
hold on
hr=plot(Fq,real(espe),c,'LineWidth',0.5);
hi=plot(Fq,imag(espe),c,'LineWidth',1.5);
set(hr,{'DisplayName'},{['fRe ',tx]});
set(hi,{'DisplayName'},{['fIm ',tx]});
xlabel('frequency')
ylabel('amplitude')
title(tit)
end

% %%
% h = allchild(gca);
% f = figure('name','uy');
% set(f,'position',[0 0 578 196]);
% copyobj(h,gca);
% clear h
% xlabel('Frecuencia [Hertz]')
% 
% %%
% h = allchild(gca);
% %%
% copyobj(h,gca);



