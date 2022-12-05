clear all
close all
clc

samp = 20

x = 1:samp;
c = -2:0.02:0;
rake = 30*pi/180;
trake = ones(samp,1).*rake;
xtrake = cos(trake);
ztrake = sin(trake);


a = -pi;
b = pi;
rake_sol = (b-a).*rand(samp,1)+a;

cont =1;
while cont < length(c)+1

cc = c(cont);
rake_dis = rake_sol - trake;

cost_ang = 0;
cos2=0;
for i=1:samp
  x = cos(rake_sol(i));
  z = sin(rake_sol(i));
  vec(i,1:2) = [x, z];
  x2 = 1 / ( (z^2. / x) + x  );
  z2 = -1*z / ( z^2. + x^2);
  x3 = trake(i);
  cost_ang = cost_ang + (atan(z/x) - x3)^2;
  cos2 = cos2 + (x3-rake_sol(i))^2;
  grad1(i) = (rake_sol(i) - x3)*x2;
  grad2(i) = (rake_sol(i) - x3)*z2;
end


for i=1:samp
solx(i) = vec(i,1) - grad1(i).*cc;
solz(i) = vec(i,2) - grad2(i).*cc;
an(i) = atan(solz(i)/solx(i));
end

cost = 0;
for k=1:samp
 cost = cost + (an(k) - trake(k))^2;
end
acost(cont) = cost;

plot(acost),hold on
%if ((cont>1) && (acost(cont) > acost(cont-1)) )

%end

rake_sol = an;
cont = cont + 1;
end

x=1:samp;
%plot(x,trake),hold on,plot(x,rake_sol,'y')
%stem(x,rake_dis,'k')
%plot(x,grad1,'k')
%plot(x,grad2,'r')
%plot(x,an,'*-k')
%ylabel('Rake (rad)'),xlabel('Time samples')











