i=strfind(tmp1,'-');
for k=length(i):-1:1
    tmp1(i(k))=[];
end
if mod(length(i),2)==1
    tmp1=['-',tmp1];
end
i=strfind(tmp1,'+');
for k=length(i):-1:1
    tmp1(i(k))=[];
end

icP=strfind(tmp1,'.*cP');
isP=strfind(tmp1,'.*sP');
icS=strfind(tmp1,'.*cS');
isS=strfind(tmp1,'.*sS');
if length(icP)==2
    tmp1(icP(2):icP(2)+3)=[];
    tmp1(icP(1):icP(1)+3)=[];
    tmp1=[tmp1,'.*cP.^2'];
elseif length(isP)==2
    tmp1(isP(2):isP(2)+3)=[];
    tmp1(isP(1):isP(1)+3)=[];
    tmp1=[tmp1,'.*sP.^2'];
elseif length(icS)==2
    tmp1(icS(2):icS(2)+3)=[];
    tmp1(icS(1):icS(1)+3)=[];
    tmp1=[tmp1,'.*cS.^2'];
elseif length(isS)==2
    tmp1(isS(2):isS(2)+3)=[];
    tmp1(isS(1):isS(1)+3)=[];
    tmp1=[tmp1,'.*sS.^2'];
elseif length(isP)==1 && length(isS)==1
    i1=max(isP,isS);i2=min(isP,isS);
    tmp1(i1:i1+3)=[];
    tmp1(i2:i2+3)=[];
    tmp1=[tmp1,'.*sP.*sS'];
elseif length(isP)==1 && length(icS)==1
    i1=max(isP,icS);i2=min(isP,icS);
    tmp1(i1:i1+3)=[];
    tmp1(i2:i2+3)=[];
    tmp1=[tmp1,'.*sP.*cS'];
elseif length(icP)==1 && length(isS)==1
    i1=max(icP,isS);i2=min(icP,isS);
    tmp1(i1:i1+3)=[];
    tmp1(i2:i2+3)=[];
    tmp1=[tmp1,'.*cP.*sS'];
elseif length(icP)==1 && length(icS)==1
    i1=max(icP,icS);i2=min(icP,icS);
    tmp1(i1:i1+3)=[];
    tmp1(i2:i2+3)=[];
    tmp1=[tmp1,'.*cP.*cS'];
elseif length(icP)==1 && length(isP)==1
    i1=max(icP,isP);i2=min(icP,isP);
    tmp1(i1:i1+3)=[];
    tmp1(i2:i2+3)=[];
    tmp1=[tmp1,'.*cP.*sP'];
elseif length(icS)==1 && length(isS)==1
    i1=max(icS,isS);i2=min(icS,isS);
    tmp1(i1:i1+3)=[];
    tmp1(i2:i2+3)=[];
    tmp1=[tmp1,'.*cS.*sS'];
end


i=strfind(tmp1,'1i*');
if length(i)==2
    tmp1(i(2):i(2)+2)=[];
    tmp1(i(1):i(1)+2)=[];
    if tmp1(1)=='-'
        tmp1(1)=[];
    else
        tmp1=['-',tmp1];
    end
end
i=strfind(tmp1,'1i*');
j=strfind(tmp1,'2i*');
if length(i)==1 && length(j)==1
    tmp1(max(i,j):max(i,j)+2)=[];
    tmp1(min(i,j):min(i,j)+2)=[];
    if tmp1(1)=='-'
        tmp1=['2*',tmp1(2:end)];
    else
        tmp1=['-2*',tmp1];
    end
end
i=strfind(tmp1,'2i*');
if length(i)==2
    tmp1(i(2):i(2)+2)=[];
    tmp1(i(1):i(1)+2)=[];
    if tmp1(1)=='-'
        tmp1=['4*',tmp1(2:end)];
    else
        tmp1=['-4*',tmp1];
    end
end

i=strfind(tmp1,'2*');
for k=length(i):-1:1
    tmp1(i(k):i(k)+1)=[];
end
if ~isempty(i)
    if tmp1(1)=='-'
        tmp1=['-',num2str(2^length(i)),'*',tmp1(2:end)];
    else
        tmp1=[num2str(2^length(i)),'*',tmp1];
    end
end

i=strfind(tmp1,'4*');
for k=length(i):-1:1
    tmp1(i(k):i(k)+1)=[];
end
if ~isempty(i)
    if tmp1(1)=='-'
        tmp1=['-',num2str(4^length(i)),'*',tmp1(2:end)];
    else
        tmp1=[num2str(4^length(i)),'*',tmp1];
    end
end

%simplification
i=strfind(tmp1,'./k1P(:,ic)');
j=strfind(tmp1,'.*k1P(:,ic)');
if length(i)==1 && length(j)==1
    tmp1(i(1):i(1)+10)=[];
    j=strfind(tmp1,'.*k1P(:,ic)');
    tmp1(j(1):j(1)+10)=[];
elseif length(i)==2 && length(j)==2
    tmp1(i(2):i(2)+10)=[];
    tmp1(i(1):i(1)+10)=[];
    j=strfind(tmp1,'.*k1P(:,ic)');
    tmp1(j(2):j(2)+10)=[];
    tmp1(j(1):j(1)+10)=[];
elseif length(i)>=1 && length(j)>=1
    tmp1(i(1):i(1)+10)=[];
    j=strfind(tmp1,'.*k1P(:,ic)');
    tmp1(j(1):j(1)+10)=[];
end

i=strfind(tmp1,'./k1S(:,ic)');
j=strfind(tmp1,'.*k1S(:,ic)');
if length(i)==1 && length(j)==1
    tmp1(i(1):i(1)+10)=[];
    j=strfind(tmp1,'.*k1S(:,ic)');
    tmp1(j(1):j(1)+10)=[];
elseif length(i)==2 && length(j)==2
    tmp1(i(2):i(2)+10)=[];
    tmp1(i(1):i(1)+10)=[];
    j=strfind(tmp1,'.*k1S(:,ic)');
    tmp1(j(2):j(2)+10)=[];
    tmp1(j(1):j(1)+10)=[];
elseif length(i)>=1 && length(j)>=1
    tmp1(i(1):i(1)+10)=[];
    j=strfind(tmp1,'.*k1S(:,ic)');
    tmp1(j(1):j(1)+10)=[];
end

%regroupement et mis a la fin
i=strfind(tmp1,'./k1P(:,ic)');
j=strfind(tmp1,'./k1S(:,ic)');
if length(i)>=1 && length(j)>=1
    tmp1(i(1):i(1)+10)=[];
    j=strfind(tmp1,'./k1S(:,ic)');
    tmp1(j(1):j(1)+10)=[];
    tmp1=[tmp1,'./(k1P(:,ic).*k1S(:,ic))'];
elseif length(i)>=1
    for k=length(i):-1:1
        tmp1(i(k):i(k)+10)=[];
    end
    if length(i)>1
        tmp1=[tmp1,'./k1P(:,ic).^',num2str(length(i))];
    else
        tmp1=[tmp1,'./k1P(:,ic)'];
    end
elseif length(j)>=1
    for k=length(j):-1:1
        tmp1(j(k):j(k)+10)=[];
    end
    if length(j)>1
        tmp1=[tmp1,'./k1S(:,ic).^',num2str(length(j))];
    else
        tmp1=[tmp1,'./k1S(:,ic)'];
    end
end
    
i=strfind(tmp1,'.*k1P(:,ic)');
j=strfind(tmp1,'.*k1S(:,ic)');
if length(i)>=1 && length(j)>=1
    tmp1(i(1):i(1)+10)=[];
    j=strfind(tmp1,'.*k1S(:,ic)');
    tmp1(j(1):j(1)+10)=[];
    tmp1=[tmp1,'.*k1P(:,ic).*k1S(:,ic)'];
end


i=strfind(tmp1,'/C66');
if length(i)==1
    tmp1(i(1):i(1)+3)=[];
    tmp1=[tmp1,'/C66'];
elseif length(i)==2
    tmp1(i(2):i(2)+3)=[];
    tmp1(i(1):i(1)+3)=[];
    tmp1=[tmp1,'/C66^2'];
end
i=strfind(tmp1,'*C66');
if length(i)==1
    tmp1(i(1):i(1)+3)=[];
    if i(1)>1
        if  tmp1(i(1)-1)=='.'
            tmp1(i(1)-1)=[];
        end
    end
    tmp1=[tmp1,'*C66'];
elseif length(i)==2
    tmp1(i(2):i(2)+3)=[];
    tmp1(i(1):i(1)+3)=[];
    tmp1=[tmp1,'*C66^2'];
end
i=strfind(tmp1,'/C66*C66');
if length(i)==1
    tmp1(i:i+7)=[];
end


i=strfind(tmp1,'k2.*');
if length(i)==2
    tmp1(i(2):i(2)+3)=[];
    tmp1(i(1):i(1)+3)=[];
    tmp1=[tmp1,'.*k2.^2'];
end
i=strfind(tmp1,'k2.*');
j=strfind(tmp1,'k2*');
if length(i)==1 && length(j)==1
    if i>j
        tmp1(i:i+3)=[];
        tmp1(j:j+2)=[];
    else
        tmp1(j:j+2)=[];
        tmp1(i:i+3)=[];
    end
    tmp1=[tmp1,'.*k2.^2'];
end


i=strfind(tmp1,'ln.*');
if length(i)==2
    tmp1(i(2):i(2)+3)=[];
    tmp1(i(1):i(1)+3)=[];
    tmp1=[tmp1,'.*ln.^2'];
end



i=strfind(tmp1,'k2.^2');
j=strfind(tmp1,'ln.^2');
if length(i)>=1
    tmp1(i(1):(i(1)+4))=[];
    if i(1)>2
        if tmp1(i(1)-2:i(1)-1)=='.*'
            tmp1(i(1)-2:i(1)-1)=[];
        elseif tmp1(i(1)-1)=='*'
            tmp1(i(1)-1)=[];
        end
    end
elseif isempty(i) && length(j)>=1
    tmp1(j(1):(j(1)+4))=[];
    if j(1)>2
        if tmp1(j(1)-2:j(1)-1)=='.*'
            tmp1(j(1)-2:j(1)-1)=[];
        elseif tmp1(j(1):j(1)+1)=='.*'
            tmp1(j(1):j(1)+1)=[];
        end
    end
    tmp1=[tmp1,'.*(ln./k2).^2'];
else
    tmp1=[tmp1,'./k2.^2'];
end


if strcmp(tmp1(1:2),'.*')
    tmp1(1:2)=[];
end
if strcmp(tmp1(1:3),'-.*')
    tmp1(1:3)=[];
end
if strcmp(tmp1(1:3),'+.*')
    tmp1(1:3)=[];
end