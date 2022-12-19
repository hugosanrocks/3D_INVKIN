function [ThisPiece] = loadStageList(ThisPiece,nf,nom)
%loadStageList Load from a txt file with a list of stl files
    ThisPiece.stagelist=(0:nf/2)*0;
    data = importdata(nom);
    nn = length(data.rowheaders);
    stage = struct('fileName','','jini',0,'jfin',0);
    for idat = 1:nn
      v = data.rowheaders(idat);
      stage(idat).fileName = v{1};
      stage(idat).jini = data.data(idat,1);
      stage(idat).jfin = data.data(idat,2);
      ThisPiece. ...
        stagelist(stage(idat).jini+1:stage(idat).jfin+1) = idat;
    end
    ThisPiece.stage = stage;
end