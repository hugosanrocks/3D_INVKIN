%     d = dialog('Position',[150 150 320 150],'Name','¿Cerrar?','color',[0.8 .1 .6]);
%     btn1 = uicontrol('Parent',d,...
%                'Position',[85 80 150 35],...
%                'String','Close IBEM',...
%                'Callback','delete(gcf)');
%     btn2 = uicontrol('Parent',d,...
%                'Position',[85 20 150 35],...
%                'String','Back to IBEM',...
%                'Callback','delete(gcf)');
%    
%              clear d btn1 btn2
%              closereq
             %%
             % Construct a questdlg with three options
choice = questdlg('Would you like to close IBEM?', ...
	'Close?', ...
	'Yes');
closeMe = false;
% Handle response
switch choice
    case 'Yes'
        disp('I''ll bring you your check.')
        closeMe = true;
    case 'No'
        disp('back to work')
end
clear choice