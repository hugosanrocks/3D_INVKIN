function l=local_long_contorno(x,cont)

dzdx= dzdx_contorno(x,cont);
l   = sqrt((dzdx).^2+1.^2);