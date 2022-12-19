function ric0=ric(n,dt,tp)

n0=length(n);

% sigma=sqrt(tp^2/(8*log(2)));
% a=(dt*(n-1)).^2/(2*sigma^2);
% ric0=1/(sigma*sqrt(2*pi)) *(exp(-a));
% 
% a=(dt*(n0+n-1)).^2/(2*sigma^2);
% ric0=ric0+1/(sigma*sqrt(2*pi)) *(exp(-a));
% 
% a=(dt*(-n0+n-1)).^2/(2*sigma^2);
% ric0=ric0+1/(sigma*sqrt(2*pi)) *(exp(-a));

a = (pi*(dt*(n-1))/tp).^2;
ric0=(0.5-a).*exp(-a); 
a = (pi*(dt*(n0+n-1))/tp).^2;
ric0=ric0+(0.5-a).*exp(-a); 
a = (pi*(dt*(-n0+n-1))/tp).^2;
ric0=ric0+(0.5-a).*exp(-a); 
% a = (pi*(dt*(2*n0+n-1))/tp).^2;
% ric0=ric0+(0.5-a).*exp(-a); 
% a = (pi*(dt*(-2*n0+n-1))/tp).^2;
% ric0=ric0+(0.5-a).*exp(-a);