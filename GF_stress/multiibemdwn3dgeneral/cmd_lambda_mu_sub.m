submed     =get(bouton.submed,'value');

alpha   = para.reg(1).sub(submed).alpha;
bet     = para.reg(1).sub(submed).bet;
rho     = para.reg(1).sub(submed).rho;

lambda  = rho*(alpha^2-2*bet^2);
mu      = rho*bet^2;

para.reg(1).sub(submed).lambda	= lambda;
para.reg(1).sub(submed).mu    	= mu;

set(bouton.lambda   ,'string',num2str(lambda));
set(bouton.mu       ,'string',num2str(mu));