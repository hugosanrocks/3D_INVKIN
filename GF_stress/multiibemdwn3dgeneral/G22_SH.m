function gn22=G22_SH(ki,rij,mu)

kr  =ki*rij;
h0kr=besselh(0,2,kr);
gn22=1/(4*1i*mu).*h0kr;
