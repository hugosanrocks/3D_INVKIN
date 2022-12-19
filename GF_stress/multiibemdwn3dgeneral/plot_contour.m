nmed1  = length(cont1);
ies_i = nrecx*nrecy*nrecz+1;
plusOne = 0;
BouwarpRange = eval(strrep(film.BoundaryWarpRange,'end',num2str(nmed1)));

for i = 1:nmed1
  if min(BouwarpRange-i) ~= 0; continue;end 
  if (warpBoundaries)
    ies_f = ies_i-1 + size(cont1(i).vec.xc(1:para.rec.resatboundaryDecimate:end),1);
    xctmp = cont1(i).vec.xc(1:para.rec.resatboundaryDecimate:end);
    zctmp = cont1(i).vec.zc(1:para.rec.resatboundaryDecimate:end);
    xVtmp = cont1(i).vec.vnx(1:para.rec.resatboundaryDecimate:end);
    zVtmp = cont1(i).vec.vnz(1:para.rec.resatboundaryDecimate:end);
    % a veces diezmar un contorno cerrado no resulta en otro concotorno
    % cerrado. Así que agregamos el punto final para que cierre.
    if ((abs(cont1(i).vec.xc(1) - cont1(i).vec.xc(end)) < 1E-9) && ...
        (abs(cont1(i).vec.zc(1) - cont1(i).vec.zc(end)) < 1E-9) && ...
        (abs(xctmp(end) - cont1(i).vec.xc(end)) > 0))
      xctmp = [xctmp;cont1(i).vec.xc(end)];
      zctmp = [zctmp;cont1(i).vec.zc(end)];
      xVtmp = [xVtmp;cont1(i).vec.vnx(end)];
      zVtmp = [zVtmp;cont1(i).vec.vnz(end)];
      plusOne = 1;
      % Normales positivas hacia afuera (corregir mitad de abajo)
      xVtmp(floor(end/2)+1:end)=-xVtmp(floor(end/2)+1:end);
      zVtmp(floor(end/2)+1:end)=-zVtmp(floor(end/2)+1:end);
    else
      plusOne =0;
    end
    
    if (filmStyle == 2 || filmStyle == 3)
      plot3(squeeze(xctmp)+ squeeze(u(it,ies_i:(ies_f+plusOne)))',...
        squeeze(zctmp)+ squeeze(v(it,ies_i:(ies_f+plusOne)))',...
        1+0*squeeze(xctmp),'k','linewidth',2);
    elseif filmStyle == 4
      if strcmp(nam(1:1),'U')
        if (it == 1)
          plot3(squeeze(xctmp),squeeze(zctmp),1+0*squeeze(xctmp),'Color',[1 0 0],'linewidth',2); % contorno inicial
        else
          set(hbou(1:it-1),'Color',[0.8 0.8 0.8],'LineWidth',0.1); % los contornos anteriores en color gris
        end
        hbou(it) = plot3(squeeze(xctmp)+ squeeze(u(it,ies_i:(ies_f+plusOne)))',...
          squeeze(zctmp)+ squeeze(v(it,ies_i:(ies_f+plusOne)))',...
          1+0*squeeze(xctmp),'k','linewidth',0.75); % contorno actual en color negro
      else %es un esfuerzo
        cla
        plot3(squeeze(xctmp),squeeze(zctmp),1+0*squeeze(xctmp),'Color',[1 1 1]*0.8,'linewidth',2); % contorno
        if (it == 1)
          for in=1:size(xctmp) % graficar las normales en el primer snapshot
            plot3([xctmp(in);xctmp(in)+m3max*xVtmp(in)],...
              [zctmp(in);zctmp(in)+m3max*zVtmp(in)],...
              [1;1],'Color',[0 1 0],'linewidth',2);
          end
        end
        if (strcmp(nam(1:3),'sxx') || ...
            strcmp(nam(1:3),'szz') || ...
            strcmp(nam(1:3),'sxz') || ...
            strcmp(nam(1:3),'sxy') || ...
            strcmp(nam(1:3),'syz'))
          uaux = u(it,ies_i:(ies_f+plusOne))';
        else % esfuerzo polar
          local.centroid.x = mean(squeeze(xctmp));
          local.centroid.z = mean(squeeze(zctmp));
          local.r = sqrt((squeeze(xctmp)-local.centroid.x).^2 + (squeeze(zctmp)-local.centroid.z).^2);
          local.co = xVtmp./ local.r;%(squeeze(xctmp)-local.centroid.x)./local.r;
          local.si = zVtmp./ local.r;%(squeeze(zctmp)-local.centroid.z)./local.r;
          if size(stc,4)~=3; error('not enough stress components');end
          if strcmp(nam(1:3),'stt')
            %  stt = si**2*Sxx + co**2*Szz - 2*si*co*Szx
            uaux = local.si.^2.* stc(filmeRange(it),ies_i:(ies_f+plusOne),iinc,1)'+...
                   local.co.^2.* stc(filmeRange(it),ies_i:(ies_f+plusOne),iinc,2)'-...
                   2*local.si.*local.co.* stc(filmeRange(it),ies_i:(ies_f+plusOne),iinc,3)';
%             uaux = stc(filmeRange(it),ies_i:(ies_f+plusOne),iinc,1)';
          elseif strcmp(nam(1:3),'srt')
            %  srt = si*co*(Szz-Sxx) + Szx * (co**2 - si**2)
            uaux = local.si.*local.co.* (stc(filmeRange(it),ies_i:(ies_f+plusOne),iinc,2)'-...
                                         stc(filmeRange(it),ies_i:(ies_f+plusOne),iinc,1)')+...
                   stc(filmeRange(it),ies_i:(ies_f+plusOne),iinc,3)'.* (local.co.^2 - local.si.^2);
          end
          uaux = uaux/sqrt(m3max);
          clear local 
        end
%           disp([num2str(it) '  :: ' num2str(uaux(1:3)')])
          patch([(squeeze(xctmp)+ squeeze(uaux.*xVtmp));squeeze(xctmp(end:-1:1))],...
                [(squeeze(zctmp)+ squeeze(uaux.*zVtmp));squeeze(zctmp(end:-1:1))],...
                ones(size(xctmp,1)*2,1),'FaceColor',[0.8 0.8 0.8],'FaceAlpha',0.8);
          clear uaux
      end
    end
    ies_i = ies_f +1;
  else
    plot3(cont1(i).vec.xc,cont1(i).vec.zc,1+0*cont1(i).vec.xc,'k','linewidth',2);
  end
end

if (malla)
  axis([min(xr)-para.rec.dxr ...
    max(xr)+para.rec.dxr ...
    min(zr)-para.rec.dzr ...
    max(zr)+para.rec.dxr -1 1])
  
  shading flat
  shading interp
  set(tmph,'FaceLighting','phong','FaceColor','interp','AmbientStrength',0.5)
else
  axis([min(min(cont1(:).vec.xc))-para.rec.dxr*3 ...
    max(max(cont1(:).vec.xc))+para.rec.dxr*3 ...
    min(min(cont1(:).vec.zc))-para.rec.dzr*3 ...
    max(max(cont1(:).vec.zc))+para.rec.dxr*3 -1 1])
end
set(gca,'ydir','reverse');
set(gca,'dataaspectratio',[1 1 1],'nextplot','replacechildren')
view([0 0 1])
% lightangle(30,40);
% lightangle(-30,140);
% light('Position',[0 .3 .5],'Style','infinite');
set(gcf,'Renderer','zbuffer')
set(findobj(gca,'type','surface'),...
  'FaceLighting','phong',...
  'AmbientStrength',.3,'DiffuseStrength',.8,...
  'SpecularStrength',.9,'SpecularExponent',25,...
  'BackFaceLighting','lit')
grid off;
axis off;

