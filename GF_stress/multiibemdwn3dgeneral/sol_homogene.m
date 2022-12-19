function [U1,U2]=sol_homogene(k1,C,A,n)
%
% fonction permettant d'obtenir les solutions homogenes
% du systeme considere : U1h la premiere solution
% et U2h la seconde solution

U1=zeros(2,n);
U2=zeros(2,n);

for i=1:2
    U1(i,:)=A(4,:)-C(6,6)*k1(i,:).^2;
    U2(i,:)=k1(i,:).*A(2,:); 
end
% on trouve une solution de l'equation homogene