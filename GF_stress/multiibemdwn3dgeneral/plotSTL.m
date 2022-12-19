function plotSTL(h,geoFileData,colorflag,EA)
% Graficar datos leidos del archivo STL
C = geoFileData.centers;
A = geoFileData.areas;
F = geoFileData.F;
V = geoFileData.V;

axes(h);
drawnow limitrate
p = patch('faces', F, 'vertices' ,V); %Cada columna es un patch distinto,
                                      %Cada renglón es vértice
%set(p, 'facec', 'b');              % Set the face color (force it)
set(p, 'facec', 'flat');            % Set the face color flat
set(p, 'EdgeColor',[0 0 0]);        % Use to see triangles, if needed.
set(p, 'FaceAlpha',0.5)               % Use for transparency
set(p, 'EdgeAlpha',EA)
normcol = 'r-';
switch colorflag % 'b','r','g','y','m','c','none'
  case 1
    set(p, 'FaceColor', 'blue');
  case 2
    set(p, 'FaceColor', 'red');
  case 3
    set(p, 'FaceColor', 'green');
  case 4
    set(p, 'FaceColor', 'yellow');
  case 5
    set(p, 'FaceColor', 'magenta');
  case 6
    set(p, 'FaceColor', 'cyan');
  case 7
    set(p, 'FaceColor', 'none');
    set(p, 'EdgeColor', [0.8 0.8 0.8]); 
    set(p, 'EdgeAlpha',EA/2)
    normcol='';
  otherwise
    set(p, 'FaceColor', 'black');
end
% set(p, 'FaceVertexCData', [0.5 1.0 0.333]);

% set(p, 'EdgeColor','none');         % Set the edge color
% light                               % add a default light
% daspect([1 1 1])                    % Setting the aspect ratio
view(3)                             % Isometric view
xlabel('X'),ylabel('Y'),zlabel('Z')
if isfield(geoFileData,'N')
N = geoFileData.N';
if (size(N,1)==3)
if ~strcmp(normcol,'')
for i=1:size(C,2)
    vn(:,1) = C(:,i);
    vn(:,2) = C(:,i)+0.5*A(i)*N(:,i);
    plot3(vn(1,:),vn(2,:),vn(3,:),normcol);
end
end
end
end
drawnow update                            %, axis manual
end
