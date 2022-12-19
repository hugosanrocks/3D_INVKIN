submed     =get(bouton.submed,'value');

set(bouton.bet  ,'string',para.reg(1).sub(submed).bet);
set(bouton.alpha,'string',para.reg(1).sub(submed).alpha);
set(bouton.rho  ,'string',para.reg(1).sub(submed).rho);
set(bouton.subh ,'string',para.reg(1).sub(submed).h);
set(bouton.Q    ,'string',para.reg(1).sub(submed).qd);
if para.reg(1).sub(submed).tipoatts==0 %update field
    para.reg(1).sub(submed).tipoatts=3;
end
set(bouton.lstatts  ,'value',para.reg(1).sub(submed).tipoatts);

if para.rafraichi==0
    cmd_att_sub;
    cmd_lambda_mu_sub;
end
