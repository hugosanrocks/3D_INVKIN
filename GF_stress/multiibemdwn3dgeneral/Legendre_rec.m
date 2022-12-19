function [Pxn,Qxn,xi]=Legendre_rec(n,x)
%Pxn:            polyn�me de Legendre de degr� n �valu� en n
%Qxn: d�riv�e du polyn�me de Legendre de degr� n �valu� en n

%evaluation des polyn�mes de Legendre par r�currence
%attention: les indices sont d�cal�s
P=zeros(n+1,1);
P(1)=1; %P0
P(2)=x; %P1
for i=2:n
    j=i-1;
    P(i+1)=1/(j+1)*((2*j+1)*x*P(i)-j*P(i-1));
end
Pxn=P(i+1);

%d�riv�e des polyn�mes de Legendre
Q=zeros(n+1,1);
Q(1)=0; %P0
Q(2)=1; %P1
for i=2:n
    j=i-1;
    Q(i+1)=1/(j+1)*((2*j+1)*(P(i)+x*Q(i))-j*Q(i-1));
end
Qxn=Q(i+1);

