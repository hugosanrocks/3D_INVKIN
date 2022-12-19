function [vg,f1,ikmax]=dispersion_curve_wfix_Wathelet(para,wi,nmode)

% funcion que da las curvas de dispersion llamando el programa de Wathelet
% por un medio estratificado sobre un semi espacio
% inconvenientes: no se puede dar los wi como uno quiere
%                 consecuentemente, se resuelve mal el principio de cada modo
% ventaja       : super rapido

nomcarpeta      = pwd;
nc  = para.nsubmed;
h   = zeros(1,nc);
rho = zeros(1,nc);
vS  = zeros(1,nc);
vP  = zeros(1,nc);
for im  = 1:nc
    h(im)   = para.reg(1).sub(im).h;
    rho(im) = para.reg(1).sub(im).rho;
    vS(im)  = para.reg(1).sub(im).bet;
    vP(im)  = para.reg(1).sub(im).alpha;
end

h(nc)=0;
cd('C:\Users\MPerton\Dropbox\courbe_dispersion\DC');
dlmwrite('input_para.txt',para.nsubmed)
M=[h(1:nc);vP(1:nc);vS(1:nc);rho(1:nc);1e3*ones(2,para.nsubmed)].';
dlmwrite('input_para.txt',M,'-append','delimiter','\t')
if para.pol==1
    modetxt=['-R 0 -L ',num2str(nmode)];
else
    modetxt=['-R ',num2str(nmode)];
end
dfi=(wi(2)-wi(1))/2/pi;
command=['gpdc ',modetxt,' -min ',num2str(max(min(wi/2/pi),1e-2*dfi)),' -max ',num2str(max(wi/2/pi)), ...
    ' -s frequency -n ',num2str(length(wi)),' <input_para.txt >curve.txt'];

tic
system(command);
toc

fileID = fopen('curve.txt','r');
tline = fgetl(fileID);
L=0;
clear vR vL fL fR
fL=zeros(nmode,1e3);
fR=zeros(nmode,1e3);
vL=zeros(nmode,1e3);
vR=zeros(nmode,1e3);
ikmax=zeros(nmode,1);

while ischar(tline)
    if ~isempty(strfind(tline,'Rayleigh'))
        fgetl(fileID);
        fgetl(fileID);
        mode  = 1;j=1;L=0;
    elseif ~isempty(strfind(tline,'Love'))
        fgetl(fileID);
        fgetl(fileID);
        mode  = 1;j=1;L=1;
    elseif ~isempty(strfind(tline,'# Mode'))
        mode=str2double(tline(8:end))+1;
        if mode>1
            ikmax(mode-1)=j-1;
        end
        j=1;
    else
        if L==1
            ik=strfind(tline,' ');
            fL(mode,j)=str2double(tline(1:ik-1));
            vL(mode,j)=str2double(tline(ik+1:end));
            j=j+1;
        else
            ik=strfind(tline,' ');
            fR(mode,j)=str2double(tline(1:ik-1));
            vR(mode,j)=str2double(tline(ik+1:end));
            j=j+1;
        end
    end
    tline = fgetl(fileID);
    ikmax(mode)=j-1;
end

fclose(fileID);

if para.pol==1
    f1=fL;
    vp=1./vL;
else
    f1=fR;
    vp=1./vR;
end
figure(206);hold on;
for j=1:nmode
    plot(f1(j,1:ikmax(j)),vp(j,1:ikmax(j)),'c--')
end

vg=vp;
cd(nomcarpeta);

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %                    vg                      %
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% command=['gpdc ',modetxt,' -group -min ',num2str(max(min(wi/2/pi),0.00001)),' -max ',num2str(max(wi/2/pi)), ...
%     ' -s frequency -n ',num2str(length(wi)),' <input_para.txt >curve.txt'];
% % tic
% system(command);
% % toc
% 
% fileID = fopen('curve.txt','r');
% tline = fgetl(fileID);
% L=0;
% clear vR vL fL fR
% fL=zeros(nmode,1e3);
% fR=zeros(nmode,1e3);
% vL=zeros(nmode,1e3);
% vR=zeros(nmode,1e3);
% ikmax=zeros(nmode,1);
% 
% while ischar(tline)
%     if ~isempty(strfind(tline,'Rayleigh'))
%         fgetl(fileID);
%         fgetl(fileID);
%         mode  = 1;j=1;L=0;
%     elseif ~isempty(strfind(tline,'Love'))
%         fgetl(fileID);
%         fgetl(fileID);
%         mode  = 1;j=1;L=1;
%     elseif ~isempty(strfind(tline,'# Mode'))
%         mode=str2double(tline(8:end))+1;
%         if mode>1
%             ikmax(mode-1)=j-1;
%         end
%         j=1;
%     else
%         if L==1
%             ik=strfind(tline,' ');
%             fL(mode,j)=str2double(tline(1:ik-1));
%             vL(mode,j)=str2double(tline(ik+1:end));
%             j=j+1;
%         else
%             ik=strfind(tline,' ');
%             fR(mode,j)=str2double(tline(1:ik-1));
%             vR(mode,j)=str2double(tline(ik+1:end));
%             j=j+1;
%         end
%     end
%     tline = fgetl(fileID);
%     ikmax(mode)=j-1;
% end
% fclose(fileID);
% cd(nomcarpeta);
% 
% if para.pol==1
%     f1=fL;
%     vg=1./vL;
% else
%     f1=fR;
%     vg=1./vR;
% end
% figure(207);hold on;
% for j=1:nmode
%     plot(f1(j,1:ikmax(j)),vg(j,1:ikmax(j)),'c')
% end