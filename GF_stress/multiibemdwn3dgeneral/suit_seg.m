function [seg_unidos,listm]=suit_seg(seg,seg_unidos,lr,dr,mi,listm,fus)
if lr==0 %left-right
    xi=seg_unidos.vec.xc(1);
    zi=seg_unidos.vec.zc(1);
else
    xi=seg_unidos.vec.xc(end);
    zi=seg_unidos.vec.zc(end);
end
c12         = [1 2];

nseg=length(seg);
for j=1:nseg
    go=0;
    if isempty(listm)
        go=1;
    elseif max(j==listm)==0
        go=1;
    end
    if go==1
        x1j=seg(j).vec.xc(1);
        z1j=seg(j).vec.zc(1);
        x2j=seg(j).vec.xc(end);
        z2j=seg(j).vec.zc(end);
        
        mj(1)=seg(j).mi;
        if mj(1)~=0
            m0 = find(fus(:,mj(1))~=0);
            if ~isempty(m0)
                mj(1)=max([m0;mj(1)]);
            end
        end
        mj(2)=seg(j).mv;
        if mj(2)~=0
            m0 = find(fus(:,mj(2))~=0);
            if ~isempty(m0)
                mj(2)=max([m0;mj(2)]);
            end
        end
        
        if min(sort(mi)==sort(mj))==1
            if mi~=mj
                %hay que cambiar c
                seg(j).vec.cv= (seg(j).vec.cv==1)*2+(seg(j).vec.cv==2)*1;
            end
            
            drj11 = sqrt((xi-x1j)^2+(zi-z1j)^2);
            if drj11<3*dr
                if lr==0
                    seg(j).vec.xc  = flipud(seg(j).vec.xc);
                    seg(j).vec.zc  = flipud(seg(j).vec.zc);
                    seg(j).vec.vnx = flipud(seg(j).vec.vnx);
                    seg(j).vec.vnz = flipud(seg(j).vec.vnz);
                    seg(j).vec.cv  = flipud(seg(j).vec.cv);
                    seg_unidos.vec.xc 	= [seg(j).vec.xc;   seg_unidos.vec.xc   ];
                    seg_unidos.vec.zc 	= [seg(j).vec.zc;   seg_unidos.vec.zc   ];
                    seg_unidos.vec.vnx 	= [seg(j).vec.vnx;  seg_unidos.vec.vnx  ];
                    seg_unidos.vec.vnz	= [seg(j).vec.vnz;  seg_unidos.vec.vnz  ];
                    seg_unidos.vec.cv	= [seg(j).vec.cv;   seg_unidos.vec.cv   ];
                else
                    seg_unidos.vec.xc 	= [seg_unidos.vec.xc;  seg(j).vec.xc];
                    seg_unidos.vec.zc 	= [seg_unidos.vec.zc;  seg(j).vec.zc];
                    seg_unidos.vec.vnx 	= [seg_unidos.vec.vnx; seg(j).vec.vnx];
                    seg_unidos.vec.vnz	= [seg_unidos.vec.vnz; seg(j).vec.vnz];
                    seg_unidos.vec.cv	= [seg_unidos.vec.cv;  seg(j).vec.cv];
                end
                listm=[listm,j];
                [seg_unidos,listm]=suit_seg(seg,seg_unidos,lr,dr,mi,listm,fus);
            end
            drj12 = sqrt((xi-x2j)^2+(zi-z2j)^2);
            if drj12<3*dr
                if lr==1
                    seg(j).vec.xc  = flipud(seg(j).vec.xc);
                    seg(j).vec.zc  = flipud(seg(j).vec.zc);
                    seg(j).vec.vnx = flipud(seg(j).vec.vnx);
                    seg(j).vec.vnz = flipud(seg(j).vec.vnz);
                    seg(j).vec.cv  = flipud(seg(j).vec.cv);
                    seg_unidos.vec.xc 	= [seg_unidos.vec.xc;  seg(j).vec.xc];
                    seg_unidos.vec.zc 	= [seg_unidos.vec.zc;  seg(j).vec.zc];
                    seg_unidos.vec.vnx 	= [seg_unidos.vec.vnx; seg(j).vec.vnx];
                    seg_unidos.vec.vnz	= [seg_unidos.vec.vnz; seg(j).vec.vnz];
                    seg_unidos.vec.cv	= [seg_unidos.vec.cv;  seg(j).vec.cv];
                else
                    seg_unidos.vec.xc 	= [seg(j).vec.xc;   seg_unidos.vec.xc  ];
                    seg_unidos.vec.zc 	= [seg(j).vec.zc;   seg_unidos.vec.zc  ];
                    seg_unidos.vec.vnx 	= [seg(j).vec.vnx;  seg_unidos.vec.vnx ];
                    seg_unidos.vec.vnz	= [seg(j).vec.vnz;  seg_unidos.vec.vnz ];
                    seg_unidos.vec.cv	= [seg(j).vec.cv;   seg_unidos.vec.cv  ];
                end
                listm=[listm,j];
                [seg_unidos,listm]=suit_seg(seg,seg_unidos,lr,dr,mi,listm,fus);
            end
        end
    end
end