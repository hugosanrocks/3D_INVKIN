nomcarpeta      = pwd;
clear h vS vP rho

para0           = para;
nc              = para0.nsubmed;
for im  = 1:nc
    h(im)   = para0.reg(1).sub(im).h;
    rho(im) = para0.reg(1).sub(im).rho;
    vS(im)  = para0.reg(1).sub(im).bet;
    vP(im)  = para0.reg(1).sub(im).alpha;
    para0.reg(1).sub(im).qd = 1000;
    para0.reg(1).sub(im).tipoatts= 0;
end

%simulation des données 
%mise en forme des données
for im  = 1:nc
    para0.reg(1).sub(im).h       = h(im);
    para0.reg(1).sub(im).rho     = rho(im);
    para0.reg(1).sub(im).bet     = vS(im);
    para0.reg(1).sub(im).alpha   = vP(im);
    para0.reg(1).sub(im).qd      = 1000;
    para0.reg(1).sub(im).tipoatts= 0;
end
%simulation des courbes de dispersion en entier
% [vg,f1,k2,w0,ikmax]=dispersion_curve_kfix_Haskel_fast_Vg(para0);
wi=linspace(0,para.fmax*2*pi,500);
dispersion_curve_wfix_Haskel(para,wi);

h(nc)=0;
cd('..\courbe_dispersion\DC');
dlmwrite('input_para.txt',para0.nsubmed)
M=[h(1:nc);vP(1:nc);vS(1:nc);rho(1:nc);1e3*ones(2,para0.nsubmed)].';
dlmwrite('input_para.txt',M,'-append','delimiter','\t')
command=['gpdc -R ',num2str(length(ikmax)),' -L ',num2str(0*length(ikmax)),' -min 0.01 -max ',num2str(para0.fmax), ...
    ' -s freq -n 500 <input_para.txt >curve.txt'];
tic
system(command);
toc

fileID = fopen('curve.txt','r');
tline = fgetl(fileID);
L=0;
clear vPR vPL fL fR
while ischar(tline)
    if ~isempty(strfind(tline,'Rayleigh'))
        tline = fgetl(fileID);
        tline = fgetl(fileID);
        mode  = 1;j=1;
    elseif ~isempty(strfind(tline,'Love'))
        tline = fgetl(fileID);
        tline = fgetl(fileID);
        mode  = 1;j=1;L=1;
    elseif ~isempty(strfind(tline,'# Mode'))
        mode=str2double(tline(8:end))+1;j=1;
    else
        if L==1
            ik=strfind(tline,' ');
            fL(mode,j)=str2double(tline(1:ik-1));
            vPL(mode,j)=str2double(tline(ik+1:end));
            j=j+1;
        else
            ik=strfind(tline,' ');
            fR(mode,j)=str2double(tline(1:ik-1));
            vPR(mode,j)=str2double(tline(ik+1:end));
            j=j+1;
        end
    end
    tline = fgetl(fileID);
end
fclose(fileID);

figure(206);hold on
if para.pol==1
    for i=1:size(fL,1)
        plot(fL(i,:),1./vPL(i,:),'-o')
    end
else
    for i=1:size(fR,1)
        plot(fR(i,:),1./vPR(i,:),'-o')
    end
end



cd(nomcarpeta);