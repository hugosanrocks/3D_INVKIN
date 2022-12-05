
elements = load('../out/elementsc.out');
position = load('../dat/nodesc.in');

[row,col] = size(elements);

for i=1:row
 line(1,1:2) = [position(elements(i,1),1),position(elements(i,1),2)];
 line(2,1:2) = [position(elements(i,2),1),position(elements(i,2),2)];
 line(3,1:2) = [position(elements(i,3),1),position(elements(i,3),2)];
 line(4,1:2) = [position(elements(i,1),1),position(elements(i,1),2)];
 linet=line';
 %[elements(i,1),elements(i,2),elements(i,3)]
 message=sprintf('line_%04d.bin',i);
 file=fopen(message,'w');
 fwrite(file,linet,'single');
 fclose(file);
end
