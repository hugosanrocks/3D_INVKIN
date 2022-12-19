function ond_inc=Source_inc(theta,k,xs,zs)
%onda incidente unitaria con un angulo theta con respecto a -z

% xs,zs coordenadades del origen de la onda (fase nula)

ond_inc.kx	= k*sin(theta);
ond_inc.kz	= k*cos(theta);

ond_inc.amp = 1*exp(-1i*(ond_inc.kx.*(-xs)-ond_inc.kz.*(-zs)));