thismed = get(bouton.med,'value');
para.cont(thismed,1).piece{info.ThisPiece}.kind = get(bouton.gfThisPieceKind,'value');
% if(para.cont(thismed,1).piece{info.ThisPiece}.kind == 2) % continuidad
%   set(bouton.gfThisPieceContinuosTo,'visible','on')
% else
%   set(bouton.gfThisPieceContinuosTo,'visible','off')
%   para.cont(thismed,1).piece{info.ThisPiece}.continuosTo = 0;
% end
if(para.cont(thismed,1).piece{info.ThisPiece}.kind == 3) % pieza sólo para cerrar el polihedro
  set(bouton.gfThisPieceColor,'value',8);
  para.cont(thismed,1).piece{info.ThisPiece}.ColorIndex = 8;
end
clear thismed