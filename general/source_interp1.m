clear all
close all
clc

sourcein=load('source.out');
[lensource,cols]=size(sourcein);
dtin=input('dt of input source file:');
dtout=input('dt of output source file:');

tmaxin=(lensource-1)*dtin;
text=sprintf('samples in input file: %03d tmaxin: %03d',lensource,tmaxin);
display(text)
t1sam=input('samples in output file:');

tin=0:dtin:(lensource-1)*dtin;
tout=0:dtout:(t1sam-1)*dtout;

for i=1:cols
  sourceout(:,i) = interp1(tin,sourcein(:,i),tout);
end

save('-ascii','source.out_interp_time','sourceout');
