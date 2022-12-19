function [ h ] = waitbarSwitch(varargin)
% turned off waitbar
ver = false;

if ver
  if nargin == 2 
    h = waitbar(varargin{1},varargin{2});
  elseif nargin == 3
    h = waitbar(varargin{1},varargin{2},varargin{3});
  elseif nargin == 4
    h = waitbar(varargin{1},varargin{2},varargin{3},varargin{4});
  end
else
  if nargin == 2 
    h = [];
    tx = varargin{2};
  elseif nargin == 3
    h = varargin{2};
    tx = varargin{3};
  elseif nargin == 4
    h = [];
    tx = varargin{2};
  end
  disp(tx);
end
end

