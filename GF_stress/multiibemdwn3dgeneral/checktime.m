clear all
j0=0;
j1=1:100:3000;
for j=j1
    j
    j0=j0+1;
    tic
    for i=1:30
        a=rand(j,j);
        b=rand(j,j);
        c=a.*b;
    end
    elapsedTime(j0)= toc;
end
plot(j1,sqrt(elapsedTime))
