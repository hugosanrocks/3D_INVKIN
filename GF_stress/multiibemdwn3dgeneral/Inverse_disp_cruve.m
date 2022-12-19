function var_opt=Inverse_disp_cruve(f_in,vg_in,nlay,para)

tic
if para.pol==2
    var_guess=zeros(4*nlay-1,1);
    n0=4;
else
    var_guess=zeros(3*nlay-1,1);
    n0=3;
end

for im  = 1:nlay
    var_guess(1+(im-1)*n0)  = para.reg(1).sub(im).rho;
    var_guess(2+(im-1)*n0)  = para.reg(1).sub(im).bet;
    if para.pol==2
        var_guess(3+(im-1)*n0)  = para.reg(1).sub(im).alpha;
        if im~=nlay
            var_guess(4+(im-1)*n0)  = para.reg(1).sub(im).h;
        end
    else
        if im~=nlay
            var_guess(3+(im-1)*n0)  = para.reg(1).sub(im).h;
        end
    end
end


var_opt = fminsearchMP(@(var2opt) fonctionelle_disp_cruve(f_in,vg_in,para,var2opt),var_guess,...
    optimset('MaxIter',1000,'Display','iter','PlotFcns',@optimplotfval));
toc
var_opt