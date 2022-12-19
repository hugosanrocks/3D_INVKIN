para.film.filmeRange = eval(get(bouton.bfilmeRangeBt,'string'));
set(info.filmeRangeTime,'string','');
if exist('RESULT','var')
  nf      = para.nf;
  df      = para.fmax/(nf/2);     %paso en frecuencia
%   Fq      = (0:nf/2)*df;
  nfN     = nf/2+1; %Nyquist
  zerospad= para.zeropad;
  dt = (1/(df*2*(nfN+zerospad))*...
  (2*(nfN+zerospad)/(2*(nfN+zerospad)-2)));
  tps     = 0:dt:1/df;
  tps     = para.pulso.b+tps;
  
  if para.film.filmeRange(end) > length(tps)
    para.film.filmeRange(para.film.filmeRange(:) > length(tps)) = [];
  end
  u = RESULT.utc(para.film.filmeRange,1,1,1);
  str = [num2str(tps(para.film.filmeRange(1)),2) ' : ' ...
         num2str(dt,2) '[' num2str(size(u,1)) '] : ' ...
         num2str(tps(para.film.filmeRange(end)),2)];
  set(info.filmeRangeTime,'string',str);
  clear u nf df nfN zerospad dt tps str
end