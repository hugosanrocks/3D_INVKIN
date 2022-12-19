function ond_inc=OP_inc(k,kxk,xs,zs,kzsigno)
% onda incidente unitaria con un componente horizontal cual es 
% definido con respectoa k el numero de onda por un factor kxk
% asi cuando abs(akx)<1, se reencuentra theta kx=akx*k=sin(theta)*k
% y   cuando abs(akx)>1, se puede incluir la ondas evanecentes

% xs,zs coordenadades del origen de la onda (fase nula)
% importante a tomar en cuenta cuando se incluye la atenuacion

ond_inc.kx	= kxk*k;
ond_inc.kz	= kzsigno.*sqrt(k^2-ond_inc.kx.^2); %kz = + o -, la 2 soluciones son correctas
if abs(kxk)>1
    ond_inc.kz	=-ond_inc.kz.*sign(imag(ond_inc.kz));
end
    
ond_inc.amp = 1*exp(-1i*(ond_inc.kx.*(-xs)-ond_inc.kz.*(-zs)));