thismed = get(bouton.med,'value');
oldNumpieces = para.cont(thismed,1).NumPieces;
para.cont(thismed,1).NumPieces = str2double(get(bouton.gfNumPieces,'string'));
set(bouton.gfThisPiece,'string',1:para.cont(thismed,1).NumPieces,'value',1);

if (para.cont(thismed,1).NumPieces > oldNumpieces)
  % para.cont(thismed,1).piece=cell(para.cont(thismed,1).NumPieces);
  for iimed = oldNumpieces+1:para.cont(thismed,1).NumPieces
    para.cont(thismed,1).piece{iimed}.fileName = '';
    para.cont(thismed,1).piece{iimed}.kind = 1;
    para.cont(thismed,1).piece{iimed}.continuosTo = 1;
    para.cont(thismed,1).piece{iimed}.ColorIndex = 1;
    para.cont(thismed,1).piece{iimed}.geoFileData = [];
  end
  clear iimed
end
clear thismed