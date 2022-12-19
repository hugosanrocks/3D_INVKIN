med     =get(bouton.med,'value');

alpha   = para.reg(med).alpha;
bet     = para.reg(med).bet;
rho     = para.reg(med).rho;

lambda  = rho*(alpha^2-2*bet^2);
mu      = rho*bet^2;

para.reg(med).lambda	= lambda;
para.reg(med).mu    	= mu;

set(bouton.lambda   ,'string',num2str(lambda));
set(bouton.mu       ,'string',num2str(mu));