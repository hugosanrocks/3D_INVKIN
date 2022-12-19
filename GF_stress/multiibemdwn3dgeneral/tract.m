function t22=tract(g1,g3,ki,rij,vnx,vnz)

ui  =1i;
% kr  =ki*rij;
% h1kr=besselh(1,2,kr);
% dh02=kr.*h1kr;
% t22 =ui/4./rij.*dh02.*(g1.*vnx+g3.*vnz);

t22 =ui/4*ki.*besselh(1,2,ki*rij).*(g1.*vnx+g3.*vnz);
