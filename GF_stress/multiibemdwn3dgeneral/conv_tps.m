function txttime=conv_tps(tstart)
%funcion que convierte el tiempo CPU en s en tiempo h:min:s

tpscalcul   = toc(tstart);

timediv   =[86400 3600 60];
timedivtxt={' d, ',' h, ',' m, ',' s'};
tnj=zeros(4,1);
txttime=[];
for i=1:3
    if tpscalcul>timediv(i)
        tnj(i)=floor(tpscalcul/timediv(i));
        tpscalcul=rem(tpscalcul,timediv(i));
        txttime=[txttime,num2str(tnj(i)),timedivtxt{i}]; %#ok<AGROW>
    end
end
i=4;
tnj(i)=tpscalcul;
txttime=[txttime,num2str(tnj(i)),timedivtxt{i}];