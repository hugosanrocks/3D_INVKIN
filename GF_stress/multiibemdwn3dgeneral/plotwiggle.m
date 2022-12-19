function plotwiggle(h,x,tr)
% PlotWig: plot section of traces in wiggle format
axes(h);
fillcolorPOS = 'k';
fillcolorNEG = 'r';
m = mean(tr);
trP = zeros(size(tr));
trN = zeros(size(tr));
for i=1:size(tr)
    if (tr(i)-m >= 0)
        trP(i)=tr(i);
        trN(i)=m;
    else
        trP(i)=m;
        trN(i)=tr(i);
    end
end
f=fill( x , trP,  fillcolorPOS);
set(f,'edgecolor','none')
f=fill( x , trN,  fillcolorNEG); 
set(f,'edgecolor','none')
% plot(h,x,tr,'w');