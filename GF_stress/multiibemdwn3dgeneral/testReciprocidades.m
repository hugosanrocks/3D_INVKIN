%% receptor en 000
set(bouton.couleur,'value',2); % color negro
x  = [0 0 0 1];  % coord x, coord y, coord z, indice del receptor x
xi = [1 2 -3 1]; % coord x, coord y, coord z, indice de la fuente xi
var = RESULT;%xi12_3x000;%xi_1_2_3x000;%xi_123x000;

figure(8881);
xi(4) = 1;
dibujo_allSpec(para,bouton,var,x,xi,'-') % G_i1
%%
figure(8882);
xi(4) = 2;
dibujo_allSpec(para,bouton,var,x,xi,'-') % G_i2

figure(8883);
xi(4) = 3;
dibujo_allSpec(para,bouton,var,x,xi,'-') % G_i3

%% receptor en 123
x = [1 2 3 1]; 
xi = [0 0 0 1];
var = xi000x12_3;

figure(8891);
xi(4) = 1;
dibujo_allSpec(para,bouton,var,x,xi,'-') % G_i1

figure(8892);
xi(4) = 2;
dibujo_allSpec(para,bouton,var,x,xi,'-') % G_i2

figure(8893);
xi(4) = 3;
dibujo_allSpec(para,bouton,var,x,xi,'-') % G_i3


%% reciprocidad de desplazamientos y esfuerzos:
set(bouton.couleur,'value',3); % color rojo
x = [1 2 -3 1];
xi = [0 0 0 1];
ch = '--';
RES.uw = xi12_3x000.uw*0;
RES.sw = xi12_3x000.sw*0;
var = xi000x12_3;%xi000x_1_2_3;%xi000x_123;

figure(8881);
RES.uw(:,1,1,1) = var.uw(:,1,1,1); 
RES.uw(:,1,1,2) = var.uw(:,1,2,1);
RES.uw(:,1,1,3) = var.uw(:,1,3,1);
RES.sw(:,1,1,1) = -var.sw(:,1,1,1); %sxx
RES.sw(:,1,1,2) = -var.sw(:,1,1,2); %syy    
RES.sw(:,1,1,3) = -var.sw(:,1,1,3); %szz
RES.sw(:,1,1,4) = -var.sw(:,1,1,4); %sxy 
RES.sw(:,1,1,5) = -var.sw(:,1,1,5); %sxz
RES.sw(:,1,1,6) = -var.sw(:,1,1,6); %syz
xi(4) = 1;
dibujo_allSpec(para,bouton,RES,x,xi,ch)

figure(8882);
RES.uw(:,1,2,1) = var.uw(:,1,1,2); 
RES.uw(:,1,2,2) = var.uw(:,1,2,2);
RES.uw(:,1,2,3) = var.uw(:,1,3,2);
RES.sw(:,1,2,1) = -var.sw(:,1,2,1); %sxx
RES.sw(:,1,2,2) = -var.sw(:,1,2,2); %syy    
RES.sw(:,1,2,3) = -var.sw(:,1,2,3); %szz
RES.sw(:,1,2,4) = -var.sw(:,1,2,4); %sxy 
RES.sw(:,1,2,5) = -var.sw(:,1,2,5); %sxz
RES.sw(:,1,2,6) = -var.sw(:,1,2,6); %syz
xi(4) = 2;
dibujo_allSpec(para,bouton,RES,x,xi,ch)

figure(8883);
RES.uw(:,1,3,1) = var.uw(:,1,1,3); 
RES.uw(:,1,3,2) = var.uw(:,1,2,3);
RES.uw(:,1,3,3) = var.uw(:,1,3,3);
RES.sw(:,1,3,1) = -var.sw(:,1,3,1); %sxx
RES.sw(:,1,3,2) = -var.sw(:,1,3,2); %syy    
RES.sw(:,1,3,3) = -var.sw(:,1,3,3); %szz
RES.sw(:,1,3,4) = -var.sw(:,1,3,4); %sxy 
RES.sw(:,1,3,5) = -var.sw(:,1,3,5); %sxz
RES.sw(:,1,3,6) = -var.sw(:,1,3,6); %syz
xi(4) = 3;
dibujo_allSpec(para,bouton,RES,x,xi,ch)