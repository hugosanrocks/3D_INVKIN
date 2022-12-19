function medio=contour_test(xs,zs,para,dr)
%indica los medios de cada lado de un punto que pertenece a un contorno

% dr	= 1e-6;

i=1;
medio=zeros(37,1);
for th=0:10:350 %[0 90 180 270]
    th1=th*pi/180;
    medio(i)=inclusiontest(xs+cos(th1)*dr,zs+sin(th1)*dr,para,1);
%     figure(1);plot(xs+cos(th1)*dr,zs+sin(th1)*dr,'.');
    i=i+1;
end
medio=unique(medio);

