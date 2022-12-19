function res=comp_mat(m1,m2,para)
%res==0 indica que los medios son distintos
res=0;
if (m1==1 || m2==1) && para.geo(1)==3 
    return
end
if m1~=0 && m2~=0
    if para.reg(m1).rho==para.reg(m2).rho
        if para.reg(m2).rho==0
            %dos medios vasios en contacto
            res=1;
        end
        
        if para.reg(m1).bet==para.reg(m2).bet
            if para.reg(m1).tipoatts==para.reg(m2).tipoatts
                if para.reg(m1).qd==para.reg(m2).qd
                    if para.pol==2 || para.dim>1
                        if  para.reg(m1).alpha==para.reg(m2).alpha
                            res=1;
                        end
                    else
                        res=1;
                    end
                end
            end
        end
    end
elseif m1==0 && m2==0
    res=1;
% elseif (m1==0 && m2~=0 && para.reg(m2).rho==0) || (m1~=0 && m2==0 && para.reg(m1).rho==0)
%     res=1;
end