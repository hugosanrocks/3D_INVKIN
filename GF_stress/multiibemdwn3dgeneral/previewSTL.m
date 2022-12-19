function [out,flag] = previewSTL(h,cont)
% Adaptado de CAD2MATDEMO de Don Riley (c) 2003
out = []; flag = 0; 

if cont.isalist
cd ..
cd ins
v = fullfile(pwd, cont.fileName);
if iscell(v)
nombreCompleto = v{1};
else
nombreCompleto = v;  
end
cd ..
cd multi-dwn-ibem.matlab
else
  nombreCompleto = cont.fileName;
end

disp(['loading: ' nombreCompleto]);
[F,V,N,C,T,A,name,errmsg] = rndread(nombreCompleto);
disp(['cant de puntos de colocación: ' num2str(length(C))])
if ~strcmp(errmsg,''); return; end
out.fileName = name;
out.centers = C;
out.triangles = T;
out.areas = A;
out.F = F;
out.V = V;
out.N = N';
if (h ~= 0)
axes(h); cla;
hold on
plotSTL(h,out,cont.ColorIndex,0.5);
axes(h);light;daspect([1 1 1]) 
end
flag = 1;
end

function [fout, vout, nout, centers, triangles, areas, CAD_object_name, errmsg] = rndread(filename)
% Reads CAD STL ASCII files, which most CAD programs can export.
% Used to create Matlab patches of CAD 3D data.
% Returns a vertex list and face list, for Matlab patch command.
% 
% filename = 'hook.stl';  % Example file.
%
fout = []; vout =[]; nout=[]; centers=[]; areas=[]; CAD_object_name = [];
if strcmp(filename,''); errmsg = 'no file selected jet'; disp(errmsg); return; end
[fid,errmsg]=fopen(filename, 'r'); %Open the file, assumes STL ASCII format.
if ~strcmp(errmsg,''); disp(errmsg); return; end
if fid == -1 
    error('File could not be opened, check name or path.')
end
% The first line is object name, then comes multiple facet and vertex lines.
CAD_object_name = sscanf(fgetl(fid), '%*s %s');  %CAD object name, if needed.
%                                                %Some STLs have it, some don't.   
vnum=0;       %Vertex number counter.
% report_num=0; %Report the status as we go.
% VColor = 0;
nnum=0;
%
while feof(fid) == 0                    % test for end of file, if not then do stuff
    tline = fgetl(fid);                 % reads a line of data from file.
    fword = sscanf(tline, '%s ');       % make the line a character string
% Check for color
%     if strncmpi(fword, 'c',1) == 1;    % Checking if a "C"olor line, as "C" is 1st char.
%        VColor = sscanf(tline, '%*s %f %f %f'); % & if a C, get the RGB color data of the face.
%     end                                % Keep this color, until the next color is used.
    if strncmpi(fword, 'v',1) == 1;    % Checking if a "V"ertex line, as "V" is 1st char.
       vnum = vnum + 1;                % If a V we count the # of V's
%        report_num = report_num + 1;    % Report a counter, so long files show status
%        if report_num > 249;
%            disp(sprintf('Reading vertix num: %d.',vnum));
%            report_num = 0;
%        end
       v(:,vnum) = sscanf(tline, '%*s %f %f %f'); % & if a V, get the XYZ data of it.
 %        c(:,vnum) = VColor;              % A color for each vertex, which will color the faces.
    end                                 % we "*s" skip the name "color" and get the data.                                          
    if strncmpi(fword, 'facetnormal',11) == 1;     %facet normal
        nnum = nnum + 1;
        n(:,nnum) = sscanf(tline, '%*s %*s %f %f %f');
    end
end
%   Build face list; The vertices are in order, so just number them.
%
fnum = vnum/3;      %Number of faces, vnum is number of vertices.  STL is triangles.
flist = 1:vnum;     %Face list of vertices, all in order.
F = reshape(flist, 3,fnum); %Make a "3 by fnum" matrix of face list data.

triangles = zeros(3,3,fnum); % ((x,y,z),(v1,v2,v3),1:fnum)
centers = zeros(3,fnum);
areas = zeros(fnum,1);
for i = 1:fnum
centers(:,i) = mean(v(1:3,F(1:3,i)),2);
triangles(:,1,i) = v(1:3,F(1,i));
triangles(:,2,i) = v(1:3,F(2,i));
triangles(:,3,i) = v(1:3,F(3,i));
areas(i) = HeronsArea(v(1:3,F(1:3,i)));
end
%
%   Return the faces and vertexs.
%
fout = F';  %Orients the array for direct use in patch.
vout = v';  % "
% cout = c';
if exist('n','var');nout = n;else nout=0;end
%
fclose(fid);
end

% Copyright (c) 2003, Don Riley
% All rights reserved.
%
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are
% met:
%
%     * Redistributions of source code must retain the above copyright
%       notice, this list of conditions and the following disclaimer.
%     * Redistributions in binary form must reproduce the above copyright
%       notice, this list of conditions and the following disclaimer in
%       the documentation and/or other materials provided with the distribution
%     * Neither the name of the Walla Walla University nor the names
%       of its contributors may be used to endorse or promote products derived
%       from this software without specific prior written permission.
%
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.