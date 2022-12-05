clear all
close all
clc

%code to create several simul.info files
%required to estimate the different 
%progressive time-windows needed to 
%run a PIS

system('rm simul_*.info');

%lines that do not change in the files
message1='3.0 1.0';
message2='40 3 648';               
message3='8192 127.9844 0.015625';
message4='3';

%time samples of the limited rupture
samples=[3 5 7 9 13 19 27 35 39];
samplesint=[6 11 16 21 31 46 66 86 96];

for i=1:length(samples)

 message5=sprintf('%01d 0.25 %01d 271 351',samples(i),samplesint(i));
 file=sprintf('simul_%02d.info',i);
 fileout=fopen(file,'a')
 fprintf(fileout,'%s\n',message1);
 fprintf(fileout,'%s\n',message2);
 fprintf(fileout,'%s\n',message3);
 fprintf(fileout,'%s\n',message4);
 fprintf(fileout,'%s\n',message5);
 fclose(file);

end


