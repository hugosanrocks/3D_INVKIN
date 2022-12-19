function [Pxn,Qxn,xi]=Legendre_rec(n,x)
%Pxn:            polynôme de Legendre de degré n évalué en n
%Qxn: dérivée du polynôme de Legendre de degré n évalué en n

%evaluation des polynômes de Legendre par récurrence
%attention: les indices sont décalés
P=zeros(n+1,1);
P(1)=1; %P0
P(2)=x; %P1
for i=2:n
    j=i-1;
    P(i+1)=1/(j+1)*((2*j+1)*x*P(i)-j*P(i-1));
end
Pxn=P(i+1);

%dérivée des polynômes de Legendre
Q=zeros(n+1,1);
Q(1)=0; %P0
Q(2)=1; %P1
for i=2:n
    j=i-1;
    Q(i+1)=1/(j+1)*((2*j+1)*(P(i)+x*Q(i))-j*Q(i-1));
end
Qxn=Q(i+1);

