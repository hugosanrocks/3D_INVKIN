% DX              = pi/para.DWNkmax;
% xmax            = DX*para.DWNnbptkx/2;
% DK              = para.DWNkmax/(para.DWNnbptkx*pi);
% strinfoDWN      = {['DX=',num2str(DX,2)];...
%                    ['Xmax=',num2str(xmax,2)];...
% %                    ['DK=',num2str(DK,2)];...
%                    ['L=',num2str(2*pi/DK,2)]};
if isfield(para.reg(1).sub(1),'bet')
minbeta = 100000000000000000000000;
for im = 1:para.nsubmed
minbeta = min(minbeta,para.reg(1).sub(im).bet);
end
else
  minbeta = 1;
end
% minbeta = min(para.reg(1).sub(:).bet);
para.DWNkmax = 0.9*(2*pi*para.fmax)/minbeta * 1.5;
if (para.dim >= 3)
DK = para.DWNkmax/para.DWNnbptkx;
else
DK = para.DWNkmax/para.DWNnbptkx/2;
end
strinfoDWN      = ['L=',num2str(2*pi/DK,2)];
set(info.infoDWN,'string',strinfoDWN);
