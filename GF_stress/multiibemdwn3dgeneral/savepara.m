prompt={['Enter the name of config:', ...
    '                                                   .']};
name='Are you sure to save para ?';
numlines=1;
defaultanswer={'configparatmp.mat'};

options.resize='on';
options.WindowStyle='normal';


answer=inputdlg(prompt,name,numlines,defaultanswer,options);
if ~isempty(answer)
    save([para.nomrep,'/',answer{1}],'para');
end