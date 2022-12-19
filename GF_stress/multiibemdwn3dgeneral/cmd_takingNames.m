prev_common_dessin; %get fieldV list of names
para.film.strmecelemlist = [];
for im = 1:size(fieldV,2)
  para.film.strmecelemlist = [para.film.strmecelemlist fieldV(im).namec(:)'];
end
clear fieldV im