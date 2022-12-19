function [tmp,kk]=follow_intersection_recursive(cont,m,m2,tmp,kk,mlist)

if m2>=m
    for c2=1:2
        nint2=length(cont(m2,c2).int);
        for k2=1:nint2
            m3 =cont(m2,c2).int(k2).xz(1,1);
            c3 =cont(m2,c2).int(k2).xz(1,2);
            xpb=cont(m2,c2).int(k2).xz;
            zpb=xpb(2:end,2);
            xpb=xpb(2:end,1);
            if m2>=m && m3>=m
                for i=1:length(xpb)
                    tmp(kk,1:6)=[xpb(i),zpb(i),m3,c3,m2,c2];
                    kk=kk+1;
                end
                %suivre les intersections de ces contours avec d autres
                %contours intermediaires
                if ~ismember(m3,mlist)
                    mlist=union(mlist,m3);
                    [tmp,kk]=follow_intersection_recursive(cont,m,m3,tmp,kk,mlist);
                end
            end
        end
    end
end