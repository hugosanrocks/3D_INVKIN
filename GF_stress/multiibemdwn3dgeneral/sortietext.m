function sortietext(para,utc)
%pour exporter les resultats en .txt

zerospad= para.zeropad;

para    = normalizacion(para);

df      = para.fmax/(para.nf/2);     %paso en frecuencia
nf      = para.nf;
nfN     = nf/2+1; %Nyquist
tps     = 0:(1/(df*2*(nfN+zerospad))*(2*(nfN+zerospad)/(2*(nfN+zerospad)-2))):1/df;
tps     = para.pulso.b+tps;
tps     = para.ta*tps;


[filename,pathname] = uiputfile('*.txt','Escoger el nombre y el lugar del archivo a escriber',para.nomrep);
if ~isempty(filename)
    M=tps.';
    for i=1:para.rec.nrec
        if para.fuente==2
            M=[M squeeze(utc(i,:))];
        else
            for iinc=1:para.ninc
                M=[M squeeze(utc(i,iinc,:))];
            end
        end
    end
    dlmwrite([pathname,filename],M,'delimiter','\t','precision',6)
    message=['archivo escrito'];
    hfreq=helpdlg(message,'info');
else
    errordlg('No hubo selection de archivo','Problema');
end

