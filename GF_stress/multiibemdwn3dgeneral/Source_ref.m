function ond_ref=Source_ref(theta,k,xs,zs)
%onda incidente unitaria con un angulo theta con respecto a -z

% xs,zs coordenadades del origen de la onda (fase nula)

ond_ref.kx	= k*sin(theta);
ond_ref.kz	= k*cos(theta);

ond_ref.amp = 1*exp(-1i*(ond_ref.kx*(-xs)+ond_ref.kz*(-zs)));