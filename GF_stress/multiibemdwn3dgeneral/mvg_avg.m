function mvg=mvg_avg(y,l,a,rec)
%moving average 
%mvg=mvg_avg(y,l,a,rec)
%y signal
%l size of the window
%a = 1: triangle, = 0: monotone
%rec recurence

[n,n1]=size(y);
if n1<n
    y=y.';
    test=1;
else
    n=n1;
    test=0;
end

mvg=y;
l20=floor(l/2);
if a==0 %todos del mismo peso
    l21=l-l20-1;
    for j=1:rec
        for i=l20+1:n-l21
            mvg(i)=mean(y(i-l20:i+l21));
        end
        for i=1:l20
            mvg(i)=mean(y(1:2*i));
        end
        for i=(n-l21+1):n
            mvg(i)=mean(y((2*i-(n-1)):n));
        end
    end
else %peso triangular centrado
    tmp =linspace(0,l20+2,l20+2);
    tgl =[tmp(2:l20+2) tmp(l20+1:-1:2)];
    tgl =tgl/sum(tgl);
    tglg=zeros(2*l20+1,l20);
    tgld=zeros(2*l20+1,l20);
    for i=1:l20
        tmp1=linspace(0,i,i+1);
        tmp1=tmp1(2:i+1);
        tmp2=linspace(i,0,2*l20+3-i);
        tmp2=tmp2(2:(2*l20+2-i));
        tglg(:,i)=[tmp1 tmp2];
        tglg(:,i)=tglg(:,i)/sum(tglg(:,i));
        tgld(:,i)=flipud(tglg(:,i));
    end
    
    for j=1:rec
        for i=l20+1:n-l20
            mvg(i)=sum(y(i-l20:i+l20).*tgl);
        end
        for i=1:l20
            mvg(i)=sum(y(1:2*l20+1).*tglg(:,i).');
        end
        for i=n-l20+1:n
            mvg(i)=sum(y((n-2*l20):n).*tgld(:,n-i+1).');
        end
        y=mvg;
    end
end
if test==1
    mvg=mvg.';
end
