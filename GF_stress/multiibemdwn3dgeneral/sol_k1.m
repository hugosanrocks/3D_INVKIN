function [k1]=sol_k1(C,A)
%
% fonction calculant les differences valeurs de k1
% correspondant a la solution homogene du probleme
% pour le cas k2 different de zero

a=C(1,1)*C(6,6);
b=-C(1,1)*A(4,:)-C(6,6)*A(1,:)-A(3,:).^2;
c=A(1,:).*A(4,:);

y1=(-b+sqrt(b.^2-4*a*c))/(2*a);
y2=(-b-sqrt(b.^2-4*a*c))/(2*a);


% la fonction sqrt de matlab donne
% un resultat complexe

k1(1,:)=sqrt(y1);
k1(2,:)=sqrt(y2);