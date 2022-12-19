% test de l inversion des courbes de dispersion

para0.pol       = 1;
para0.nsubmed   = 3;
para0.fmax      = 3;
para0.nmed      = 1;
para0.geo     	= 3;
nmodem  =2;

%données a retrouver
h(1)    = 1;
vS(1)   = 1;
vP(1)   = 2;
rho(1)  = 1;

h(2)    = 4;
vS(2)   = 1.5;
vP(2)   = 3;
rho(2)  = 1.5;

h(3)    = 30;
vS(3)   = 3;
vP(3)   = 6;
rho(3)  = 2;

h(4)    = 0;
vS(4)   = 5;
vP(4)   = 10;
rho(4)  = 2.2;

%simulation des données 
%mise en forme des données
for im  = 1:para0.nsubmed
    para0.reg(1).sub(im).h       = h(im);
    para0.reg(1).sub(im).rho     = rho(im);
    para0.reg(1).sub(im).bet     = vS(im);
    para0.reg(1).sub(im).alpha   = vP(im);
    para0.reg(1).sub(im).qd      = 1000;
    para0.reg(1).sub(im).tipoatts= 0;
end
%simulation des courbes de dispersion en entier
% [~,~,k2,w0,ikmax]=dispersion_curve_kfix_MG_fast2(para0);
tic
[vg,f1,ikmax]=dispersion_curve_kfix_Haskel_fast_Vg(para0);

toc
figure(207);hold on;
for imode=1:nmodem%length(ikmax)
    plot(f1(imode,1:ikmax(imode)),vg(imode,1:ikmax(imode)),'k')
end


%selection de quelques points
f_in    =[];
vg_in   =[];
for imode = 1:nmodem
    if imode==1
        fin=linspace(.05,2,40);
    elseif imode==2
        fin=linspace(.2,1.2,15);
    elseif imode==3
        fin=linspace(.3,2,15);
    end
    vgin=interp1(f1(imode,1:ikmax(imode)),vg(imode,1:ikmax(imode)),fin);
    f_in    =[f_in,fin];
    vg_in 	=[vg_in,vgin];
end
figure(207);hold on;
plot(f_in,vg_in,'xr')



%guess
h(1)    = 1.2;
vS(1)   = 1.1;
vP(1)   = 3;
rho(1)  = 1;


h(2)    = 4;
vS(2)   = 2.;
vP(2)   = 5;
rho(2)  = 1.5;

para1=para0;
%simulation des données de point de départ
%mise en forme des données
for im  = 1:para0.nsubmed
    para1.reg(1).sub(im).h       = h(im);
    para1.reg(1).sub(im).rho     = rho(im);
    para1.reg(1).sub(im).bet     = vS(im);
    para1.reg(1).sub(im).alpha   = vP(im);
    para1.reg(1).sub(im).qd      = 1000;
    para1.reg(1).sub(im).tipoatts= 0;
end
[vg,f1,k2,w0,ikmax]=dispersion_curve_kfix_MG_fast2(para1);
figure(207);hold on;
for imode=1:nmodem%length(ikmax)
    plot(f1(imode,1:ikmax(imode)),vg(imode,1:ikmax(imode)),'')
end
pause(.1)

indpb=find(isnan(vg_in));
vg_in(indpb)=[];
f_in(indpb)=[];
%invertion des données 
para1.fmax= max(f_in);
var_opt=Inverse_disp_cruve(f_in,vg_in,para0.nsubmed,para1);

%dessin des courbes de dispersion des données optimisées
paraopt=para1;

if paraopt.pol==2
    n0=4;
else
    n0=3;
end    
for im  = 1:paraopt.nsubmed
    paraopt.reg(1).sub(im).rho     = var_opt(1+(im-1)*n0);
    paraopt.reg(1).sub(im).bet     = var_opt(2+(im-1)*n0);
    if paraopt.pol==2
        paraopt.reg(1).sub(im).alpha   = var_opt(3+(im-1)*n0);
        if im~=paraopt.nsubmed
            paraopt.reg(1).sub(im).h       = var_opt(4+(im-1)*n0);
        end
    else
        if im~=paraopt.nsubmed
            paraopt.reg(1).sub(im).h       = var_opt(3+(im-1)*n0);
        end
    end
end
paraopt.fmax= max(f_in);

paraopt    = Vp2Cij(paraopt);
[vg,f1,ikmax]   = dispersion_curve_wfix_Haskel_4inv(paraopt,2*pi*f_in);
% [vg,f1,ikmax]   = dispersion_curve_kfix_MG_fast2(paraopt);
figure(207);hold on;
for imode=1:nmodem%length(ikmax)
    plot(f1(imode,1:ikmax(imode)),vg(imode,1:ikmax(imode)),'r')
end

%comparaison des resultats parametrique
nlay=para0.nsubmed;
disp('        target        result')
sub=para0.reg(1).sub;
for im  = 1:nlay
    disp(['rho(',num2str(im),')=    ',num2str(sub(im).rho),'        ',num2str(var_opt(1+(im-1)*n0))])
    disp(['bet(',num2str(im),')=    ',num2str(sub(im).bet),'        ',num2str(var_opt(2+(im-1)*n0))])
     if para0.pol==2
        disp(['alf(',num2str(im),')=    ',num2str(sub(im).alpha),'        ',num2str(var_opt(3+(im-1)*n0))])
         if im~=nlay
            disp(['  h(',num2str(im),')=    ',num2str(sub(im).h),'        ',num2str(var_opt(4+(im-1)*n0))])
         end
     else
         if im~=nlay
            disp(['  h(',num2str(im),')=    ',num2str(sub(im).h),'        ',num2str(var_opt(3+(im-1)*n0))])
         end
    end
end