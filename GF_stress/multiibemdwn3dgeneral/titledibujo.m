titletxt=[fieldV(ifig).name,' for '];
titletxt=[titletxt,'FP in xs:',num2str(para.xs(iinc))];
if para.dim>1
    titletxt=[titletxt,', ys:',num2str(para.ys(iinc))];
end
titletxt=[titletxt,', zs:',num2str(para.zs(iinc)),'; '];


if para.fuente==1
    if para.dim==1
        if para.pol==1
            if para.tipo_onda(iinc)==1
                titletxt=[titletxt,'OP SH'];
            else
                titletxt=[titletxt,'OP Love'];
            end
        else
            if para.tipo_onda(iinc)==1
                titletxt=[titletxt,'OP P'];
            elseif para.tipo_onda(iinc)==2
                titletxt=[titletxt,'OP SV'];
            elseif para.tipo_onda(iinc)==3
                titletxt=[titletxt,'OP Rayleigh'];
            end
        end
    else
        if para.tipo_onda(iinc)==1
            titletxt=[titletxt,'OP P'];
        elseif para.tipo_onda(iinc)==2
            titletxt=[titletxt,'OP SV'];
        elseif para.tipo_onda(iinc)==3
            titletxt=[titletxt,'OP SH'];
        elseif para.tipo_onda(iinc)==4
            titletxt=[titletxt,'OP Rayleigh'];
        elseif para.tipo_onda(iinc)==5 %&& (para.geo(1)==3 || max(para.geo(2:para.nmed)==2)==1
            titletxt=[titletxt,'OP Love'];
        end
    end
    titletxt=[titletxt,' with incidence of ',num2str(round(para.gam(iinc))),', '];
elseif para.fuente==2
    if para.dim==1 && para.pol==2
        titletxt=[titletxt,'con orientacion fx: ',num2str(para.fij(iinc,1)),' fz: ',num2str(para.fij(iinc,2)),'; '];
    elseif para.dim>=3
        titletxt=[titletxt,'con orientacion fx: ',num2str(para.fij(iinc,1)),', fy: ',num2str(para.fij(iinc,2)),', fz: ',num2str(para.fij(iinc,3)),'; '];
    end
else
    titletxt=[];
end