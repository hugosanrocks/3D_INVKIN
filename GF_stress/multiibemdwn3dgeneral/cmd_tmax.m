para.tmax = (para.zeropad-1)/(para.fmax/(para.nf/2)*para.zeropad);
para.tmaxinteres = para.tmax;
set(info.tmax,'string',...
  {['tmax = ' num2str(para.tmax)];'';'T interes:'});
set(bouton.tmax,'string',para.tmaxinteres);
