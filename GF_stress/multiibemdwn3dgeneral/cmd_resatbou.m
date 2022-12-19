if (isfield(para,'cont1'))
    paratmp = para;
else
    paratmp = para;
    h = waitbarSwitch(0,'Calculo del contorno fino');
    if isfield(para,'subm') % ----------?
        paratmp = rmfield(paratmp,'subm');
        warning('deleted para.subm  correct?');
    end
    if isfield(para,'fus')
        paratmp = rmfield(paratmp,'fus');
        warning('deleted para.fus  correct?');
    end
    % se normaliza unas variables
    paratmp = normalizacion(paratmp);
    sal     = paratmp.sortie;
    nf      = paratmp.nf;
    paratmp.df = paratmp.fmax/(nf/2);     %paso en frecuencia
    df      = paratmp.df;
    %discretiza los contornos iniciales
    paratmp = malla_fina_init_2(paratmp,h);
    %calcula los verdaderos contornos, buscando los puntos de intersecion,
    %los contornos que se quedan y los que se quitan
    paratmp = recalculate_contornos(paratmp,h);
    close(h)
end

if (~isfield(paratmp,'cont1')); error('No para.cont1');end
nmed1  = length(paratmp.cont1);
xctmp = [];zctmp = []; 
    ranMed0 = 1:nmed1; ranMed=(eval(['ranMed0(',paratmp.rec.resatboundaryScaleMediaRange,')']));
for i = 1:nmed1
    auxX = paratmp.cont1(i).vec.xc(1:paratmp.rec.resatboundaryDecimate:end);
    auxZ = paratmp.cont1(i).vec.zc(1:paratmp.rec.resatboundaryDecimate:end);
    % a veces diezmar un contorno cerrado no resulta en otro concotorno
    % cerrado. Así que agregamos el punto final para que cierre.
    if (abs((paratmp.cont1(i).vec.xc(1) - paratmp.cont1(i).vec.xc(end)) < 1E-9) && ...
        abs((paratmp.cont1(i).vec.zc(1) - paratmp.cont1(i).vec.zc(end)) < 1E-9) && ...
        abs((auxX(end) - paratmp.cont1(i).vec.xc(end)) > 0))
           auxX = [auxX;paratmp.cont1(i).vec.xc(end)];
           auxZ = [auxZ;paratmp.cont1(i).vec.zc(end)];
    end
    % escala
    if (paratmp.rec.resatboundaryScale ~= 1)
       xmed = mean(auxX); zmed = mean(auxZ); 
       if (find(ranMed,i,'first')>0)
           auxX = xmed + (auxX - xmed) * paratmp.rec.resatboundaryScale;
           auxZ = zmed + (auxZ - zmed) * paratmp.rec.resatboundaryScale;
       end
    end
    % stack receivers
    xctmp = [xctmp;auxX];
    zctmp = [zctmp;auxZ];
end
nrestmp = size(xctmp,1);
clear paratmp