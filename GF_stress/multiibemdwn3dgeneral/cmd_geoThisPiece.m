thismed = get(bouton.med,'value');
info.ThisPiece = get(bouton.gfThisPiece,'value');
if size(para.cont(thismed,1).piece,2) >= info.ThisPiece
set(bouton.geoFileSelect,'string',          para.cont(thismed,1).piece{info.ThisPiece}.fileName);
set(bouton.gfThisPieceContinuosTo,'string', para.cont(thismed,1).piece{info.ThisPiece}.continuosTo);
% if(para.cont(thismed,1).piece{info.ThisPiece}.kind == 2) % continuidad
%   set(bouton.gfThisPieceContinuosTo,'visible','on')
% else
%   set(bouton.gfThisPieceContinuosTo,'visible','off')
% end
set(bouton.gfThisPieceKind,'value',         para.cont(thismed,1).piece{info.ThisPiece}.kind);
set(bouton.gfThisPieceColor,'value',        para.cont(thismed,1).piece{info.ThisPiece}.ColorIndex);
if size(para.cont(thismed,1).piece{info.ThisPiece}.fileName,2)>1
  cla(bouton.gfPreview)
  plotSTL(bouton.gfPreview,para.cont(thismed,1).piece{info.ThisPiece}.geoFileData,...
                           para.cont(thismed,1).piece{info.ThisPiece}.ColorIndex,0.5);
end
end
clear thismed