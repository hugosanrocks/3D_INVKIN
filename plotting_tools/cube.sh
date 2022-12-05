subset <cube_temp n1=12 n2=24 n3=4 >block1 
pscube <block1 n1=12 n2=24 n3=4 \
       wclip=0.26 bclip=0.31 \
       d1=1.0 d2=1.0 d3=2. size1=5 size2=10 size3=4 \
       angle=45 whls=0.666666,.5,1 bhls=0,.5,1.1 \
       d1num=5 d2num=5 d3num=2 \
       legend=1 lstyle=horibottom units="Resolution length [s]" \
       label1="Along dip [km]" label2="Along strike [km]" label3="Time [s]" \
       labelsize=45 >block1.ps
#subset <cube_temp n1=35 n2=35 n3=35 n1s=17 if3s=17 >block2 
#pscube <block2 n1=17 n2=35 n3=18 wclip=0.0 bclip=0.7 d1=5.0 d2=5.0 d3=5.0 size1=5 size2=10 size3=5 angle=45 whls=0.666666,.5,1 bhls=0,.5,1 legend=1 labelsize=45 d1num=1000 d2num=1000 d3num=1000 label2=X label3=Y >block2.ps
#~/Desktop/cwp/bin/psmerge in=block1.ps translate=0,0 in=block2.ps translate=3.53,8.53 >merge.ps 
