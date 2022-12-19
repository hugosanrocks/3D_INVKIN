function [nomrep]=chercherep

str='escoger un repertorio donde guardar los resultados';
if isunix
    nomrep = uigetdir('/Users/mathieuperton/Desktop/',str);
else
    if strcmp(version,'7.4.0.287 (R2007a)')==1
        nomrep = 'C:\Users\Mathieu\Desktop\TMP_calcul';
    else
        nomrep = uigetdir('C:\Documents and Settings\MPerton\Escritorio\',str);
    end
end