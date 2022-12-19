pos=[...
    1 2
    1 3
    1 4
    2 3
    2 4
    3 4];
% det0        = zeros(6,6,nk2);

G11    = '2*k2.^2.*cP-ln.*cS';
G12    = '1i*k2./k1P(:,ic).*(ln.*sP+2*k1P(:,ic).*k1S(:,ic).*sS)';
G13    = '1i*k2/C66.*(cP-cS)';
G14    = '1/C66./k1P(:,ic).*(k2.^2.*sP+k1P(:,ic).*k1S(:,ic).*sS)';
G21    ='-1i*k2./k1S(:,ic).*(2*k1P(:,ic).*k1S(:,ic).*sP+ln.*sS)';
G22    ='-ln.*cP+2*k2.^2.*cS';
G23    = '1/C66./k1S(:,ic).*(k1P(:,ic).*k1S(:,ic).*sP+k2.^2.*sS)';
G24    = '1i*k2/C66.*(cP-cS)';
G31    = '2i*C66*k2.*ln.*(cP-cS)';
G32    ='-C66./k1P(:,ic).*(ln.^2.*sP+4*k2.^2.*k1P(:,ic).*k1S(:,ic).*sS)';
G33    ='-ln.*cP+2*k2.^2.*cS';
G34    = '1i*k2./k1P(:,ic).*(ln.*sP+2*k1P(:,ic).*k1S(:,ic).*sS)';
G41    ='-C66./k1S(:,ic).*(4*k2.^2.*k1P(:,ic).*k1S(:,ic).*sP+ln.^2.*sS)';
G42    = '2i*C66*k2.*ln.*(cP-cS)';
G43    ='-1i*k2./k1S(:,ic).*(2*k1P(:,ic).*k1S(:,ic).*sP+ln.*sS)';
G44    = '2*k2.^2.*cP-ln.*cS';


for jic=6%1:6
    a=pos(jic,1);
    b=pos(jic,2);
    for jic2=1:6
        c=pos(jic2,1);
        d=pos(jic2,2);
        linetxt1=['det0(',num2str(jic),',',num2str(jic2),',:)='];
        %linetxt2=['(',eval(['G',num2str(a),num2str(c)]),').*(',eval(['G',num2str(b),num2str(d)]),')-(',eval(['G',num2str(a),num2str(d)]),').*(',eval(['G',num2str(b),num2str(c)]),')'];
        
        linetxt2=eval(['G',num2str(a),num2str(c)]);
        %isoler les facteurs communs
        i0=strfind(linetxt2,'.*(');
        facac=linetxt2(1:i0+1);
        if isempty(i0)
            i0=-2;
            ie=0;
        else
            ie=1;
        end
        %isoler les termes
        i1=strfind(linetxt2,'+');
        i2=strfind(linetxt2,'-');
        it=[i1 i2];
        sign2ac =~isempty(i1)-~isempty(i2);
        term1ac=linetxt2(i0+3:(it-1));
        term2ac=linetxt2(it:end-ie);
        
        linetxt2=eval(['G',num2str(b),num2str(d)]);
        %isoler les facteurs communs
        i0=strfind(linetxt2,'.*(');
        facbd=linetxt2(1:i0+1);
        if isempty(i0)
            i0=-2;
            ie=0;
        else
            ie=1;
        end
        %isoler les termes
        i1=strfind(linetxt2,'+');
        i2=strfind(linetxt2,'-');
        it=[i1 i2];
        sign2bd =~isempty(i1)-~isempty(i2);
        term1bd=linetxt2(i0+3:(it-1));
        term2bd=linetxt2(it:end-ie);
        
        linetxt2=eval(['G',num2str(a),num2str(d)]);
        %isoler les facteurs communs
        i0=strfind(linetxt2,'.*(');
        facad=linetxt2(1:i0+1);
        if isempty(i0)
            i0=-2;
            ie=0;
        else
            ie=1;
        end
        %isoler les termes
        i1=strfind(linetxt2,'+');
        i2=strfind(linetxt2,'-');
        it=[i1 i2];
        signad =~isempty(i1)-~isempty(i2);
        term1ad=linetxt2(i0+3:(it-1));
        term2ad=linetxt2(it:end-ie);
        
        linetxt2=eval(['G',num2str(b),num2str(c)]);
        %isoler les facteurs communs
        i0=strfind(linetxt2,'.*(');
        facbc=linetxt2(1:i0+1);
        if isempty(i0)
            i0=-2;
            ie=0;
        else
            ie=1;
        end
        %isoler les termes
        i1=strfind(linetxt2,'+');
        i2=strfind(linetxt2,'-');
        it=[i1 i2];
        signbc =~isempty(i1)-~isempty(i2);
        term1bc=linetxt2(i0+3:(it-1));
        term2bc=linetxt2(it:end-ie);
        
        %linetxt2=['(',eval(['G',num2str(a),num2str(c)]),').*(',eval(['G',num2str(b),num2str(d)]),')-(',eval(['G',num2str(a),num2str(d)]),').*(',eval(['G',num2str(b),num2str(c)]),')'];
        tmp ='k2.^2.*(...';
        %!!!!!!!!!!!les sign ab  bd bc bañoihsaefldiugw4efu!!!!!!!!
        tmp1=[facac,facbd,term1ac,'.*',term1bd];simplifie_exp_det_rayleigh;tmp=[tmp,tmp1];
        tmp1=[facac,facbd,term1ac,'.*',term2bd];simplifie_exp_det_rayleigh;tmp=[tmp,' ... +',tmp1];
        tmp1=[facac,facbd,term2ac,'.*',term1bd];simplifie_exp_det_rayleigh;tmp=[tmp,' ... +',tmp1];
        tmp1=[facac,facbd,term2ac,'.*',term2bd];simplifie_exp_det_rayleigh;tmp=[tmp,' ... +',tmp1];
        tmp1=[facad,facbc,term1ad,'.*',term1bc];simplifie_exp_det_rayleigh;tmp=[tmp,' ... -',tmp1];
        tmp1=[facad,facbc,term1ad,'.*',term2bc];simplifie_exp_det_rayleigh;tmp=[tmp,' ... -',tmp1];
        tmp1=[facad,facbc,term2ad,'.*',term1bc];simplifie_exp_det_rayleigh;tmp=[tmp,' ... -',tmp1];
        tmp1=[facad,facbc,term2ad,'.*',term2bc];simplifie_exp_det_rayleigh;tmp=[tmp,' ... -',tmp1];
        tmp =[tmp,')'];
        
        i=strfind(tmp,'--');
        for k=length(i):-1:1
            tmp(i(k)+1)='+';
            tmp(i(k))=[];
        end
        
        i=strfind(tmp,'-+');
        for k=length(i):-1:1
            tmp(i(k)+1)='-';
            tmp(i(k))=[];
        end
        
        i=strfind(tmp,'+-');
        for k=length(i):-1:1
            tmp(i(k)+1)='-';
            tmp(i(k))=[];
        end
        
        i=strfind(tmp,'++');
        tmp(i)=[];
        
        i=strfind(tmp,'1.*');
        for k=length(i):-1:1
            tmp(i(k):(i(k)+2))=[];
        end
        
        i=strfind(tmp,'1i*');
        if length(i)==8
            for k=length(i):-1:1
                tmp(i(k):(i(k)+2))=[];
            end
            tmp=[tmp(1:end),'*1i'];
        end
        
        i=strfind(tmp,'/C66^2');
        if length(i)==8
            for k=8:-1:1
                tmp(i(k):i(k)+5)=[];
            end
            tmp=[tmp(1:end),'/C66^2'];
        end
        
        i=strfind(tmp,'/C66');
        if length(i)==8
            for k=8:-1:1
                tmp(i(k):i(k)+3)=[];
            end
            tmp=[tmp(1:end),'/C66'];
        end
        
        i=strfind(tmp,'/C66^2');
        if length(i)==8
            for k=8:-1:1
                tmp(i(k):i(k)+3)=[];
            end
            tmp=[tmp(1:end),'/C66^2'];
        end
        
        i=strfind(tmp,'*C66');
        j=strfind(tmp,'C66.*');
        if length(i)==8
            for k=8:-1:1
                tmp(i(k):i(k)+3)=[];
            end
            tmp=[tmp(1:end),'*C66'];
        elseif (length(i)+length(j))>=8
            for k=length(i):-1:1
                tmp(i(k):i(k)+3)=[];
            end
            j=strfind(tmp,'C66.*');
            for k=length(j):-1:1
                tmp(j(k):j(k)+4)=[];
            end
            tmp=[tmp(1:end),'*C66'];
        elseif isempty(i)
        else
            toto=0;
        end
        
        i=strfind(tmp,'*C66^2');
        if length(i)==8
            for k=8:-1:1
                tmp(i(k):i(k)+3)=[];
            end
            tmp=[tmp(1:end),'*C66^2'];
        end
        
        disp([linetxt1,tmp,';'])
    end
end