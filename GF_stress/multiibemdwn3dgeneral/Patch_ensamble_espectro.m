function [utc,uw,stc,sw,name,b_dib,cont1]=Patch_ensamble_espectro(bdessin,b_dib)
%funcion a utilizar cuando, por alguna razon, la simulacion no se ha
%completo pero que se conta con multiplos archivos de respaldo que
%contienen cada frecuencia por separado

%Repertorio de los archivos
raiz    = 'C:\Users\MPerton\Desktop\pbsimu\';
% Lista de archivos de la simulacion con problema
list0   = ls([raiz,'*.mattmp*']);
Narch   = size(list0,1);
% Lista de archivos de la simulacion auxiliarios acabados
%se repuso simulaciones con los elementos que hacen falta
list1   = ls([raiz,'*.mat']);
Narch1  = size(list1,1);

h       = waitbarSwitch(0,'Concatenacion archivos');
ifq0    = zeros(Narch,1);

%recuperacion de la informacion de los archivos inacabados
for i = 1:Narch
    waitbarSwitch(i/Narch,h);
    %reemplaza espacio por nada
    tmparch = strrep(list0(i,1:end),' ','');
    %identificacion de la frecuencia del archivo
    indn    = strfind(tmparch,'mattmp')+6;
    ifq     = str2num(tmparch(indn:end));
    %cargar el archivo
    load([raiz,tmparch],'-mat');
    
    %inicialisacion de los campos a la primera iteracion
    if i==1
        nf      = para.nf;
        if exist('uw','var')
            [n1,n2,n3,n4]=size(uw);
            %uw  = zeros(nf/2+1,para.rec.nrec,para.ninc,ns);
            uw0=zeros(nf/2+1,n2,n3,n4);
        end
        if exist('sw','var')
            [n1,n2,n3,n4]=size(sw);
            %sw  = zeros(nf/2+1,para.rec.nrec,para.ninc,nss);
            sw0=zeros(nf/2+1,n2,n3,n4);
        end
    end
    
    if exist('uw','var')
        uw0(ifq+1,:,:,:)=uw;
    end
    if exist('sw','var')
        sw0(ifq+1,:,:,:)=sw;
    end
    
    ifq0(i)=ifq;
end

%recuperacion de la informacion de los archivos acabados 
for i = 1:Narch1
    waitbarSwitch(i/Narch1,h);
    %reemplaza espacio por nada
    tmparch = strrep(list1(i,1:end),' ','');
    %cargar el archivo
    load([raiz,tmparch],'-mat');
    
    
    if exist('uw','var')
        %identificacion de las frecuencias hechas
        for ifq=0:para.nf/2
            if uw(ifq+1,1,1,1)~=0 && uw0(ifq+1,1,1,1)==0
                uw0(ifq+1,:,:,:)=uw(ifq+1,:,:,:);
            end
        end
    end
%     if exist('sw','var')
%         for ifq=0:para.nf/2
%             if sw(ifq+1,1,1,1)~=0 && sw0(ifq+1,1,1,1)~=0
%                 sw0(ifq+1,:,:,:)=sw;
%             end
%         end
%     end
    
    ifq0(i)=ifq;
end


%copia de las variables concatenadas en su denominacion inicial
if exist('uw','var')
    uw	= uw0;
else
    uw  = zeros(nf/2+1,n2,n3,0);
end
% if exist('sw','var')
%     sw  = sw0;
% else
    sw  = zeros(nf/2+1,n2,n3,0);
% end
%nombre del archivo a guardar
name    = [raiz,tmparch(1:indn-15),'.mat'];

%recuperacion informacion del contorno (occupado en dibujo)
cont1   = para.cont1;
%inversion w
para.spct=1;
if para.spct==0
    waitbarSwitch(0,h,'Inversion de los espectros');

    [utc,stc]   = inversion_w(uw,sw,para);
    
    waitbarSwitch(0,h,'Guardando los resultados ...');
    save(name,'para','utc','uw','stc','sw','cont1');
else
    waitbarSwitch(0,h,'Guardando los resultados ...');
    save(name,'para','uw','sw','cont1');
    utc         = 0;
    stc         = 0;
end
% delete([name,'tmp*']); %prefiero en este caso que se borre a mano

%%%%%%%%%%%%%%%%%%%%%
% dibujo resultados %
%%%%%%%%%%%%%%%%%%%%%
if exist('bdessin','var')
    waitbarSwitch(0,h,'Dibujando los resultados ...');
    para.redraw     = 0;
    b_dib(1).name   = name;
    b_dib           = dibujo(para,bdessin,utc,uw,stc,sw,cont1,b_dib);
end
close(h)