thismed = get(bouton.med,'value');
para.cont(thismed,1).NumPieces = str2double(get(bouton.gfNumPieces,'string'));
nom = chercheSTL(thismed);
set(bouton.geoFileSelect,'string',['Piece(',num2str(info.ThisPiece),'): ',nom]);
% tomar los valores más actuales:
para.cont(thismed,1).piece{info.ThisPiece}.fileName = nom;
para.cont(thismed,1).piece{info.ThisPiece}.continuosTo = str2double(get(bouton.gfThisPieceContinuosTo,'string'));
para.cont(thismed,1).piece{info.ThisPiece}.ColorIndex = get(bouton.gfThisPieceColor,'value');
para.cont(thismed,1).piece{info.ThisPiece}.kind = get(bouton.gfThisPieceKind,'value');
if para.cont(thismed,1).piece{info.ThisPiece}.kind == 1 % free surface
  para.cont(thismed,1).piece{info.ThisPiece}.continuosTo = 0;
end
k = strfind(nom,'.txt');
if isempty(k) % Es un archivo .stl
para.cont(thismed,1).piece{info.ThisPiece}.isalist = false;
[para.cont(thismed,1).piece{info.ThisPiece}.geoFileData,flag] = ...
 previewSTL(bouton.gfPreview,para.cont(thismed,1).piece{info.ThisPiece});
else % es una lista de archivos .stl en la misma carpeta
    para.cont(thismed,1).piece{info.ThisPiece}.isalist = true;
    para.cont(thismed,1).piece{info.ThisPiece} = loadStageList(...
      para.cont(thismed,1).piece{info.ThisPiece},para.nf,nom);
    
    % El ultimo de la lista para que los receptores tengan la región
    % correcta:
    auxcont.fileName = para.cont(thismed,1).piece{info.ThisPiece}.stage(end).fileName;
    auxcont.isalist = true;
    auxcont.ColorIndex = para.cont(thismed,1).piece{info.ThisPiece}.ColorIndex;
[para.cont(thismed,1).piece{info.ThisPiece}.geoFileData,flag] = ...
 previewSTL(bouton.gfPreview,auxcont);
    clear auxcont
end

if flag == 0 % cuando no se pudo cargar
  para.cont(thismed,1).piece{info.ThisPiece}.fileName ='X';
else
  % apilar 
  para.cont(thismed,1).FV.vertices = [];
  para.cont(thismed,1).FV.faces = [];
  para.cont(thismed,1).FV.facenormals = [];
  for ip = 1:para.cont(thismed,1).NumPieces
      if size(para.cont(thismed,1).piece,2) >= ip
    if isfield(para.cont(thismed,1).piece{ip}.geoFileData,'V')
      para.cont(thismed,1).FV.vertices =    [para.cont(thismed,1).FV.vertices;    para.cont(thismed,1).piece{ip}.geoFileData.V];
      para.cont(thismed,1).FV.faces =       [para.cont(thismed,1).FV.faces;       para.cont(thismed,1).piece{ip}.geoFileData.F];
      para.cont(thismed,1).FV.facenormals = [para.cont(thismed,1).FV.facenormals; para.cont(thismed,1).piece{ip}.geoFileData.N];
    end
      end
  end
end
clear thismed nom ip k
rafraichi;
