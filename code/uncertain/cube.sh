subset <cube_temp n1=12 n2=24 n3=45 >block1 
pscube <block1 n1=12 n2=24 n3=45 \
       wclip=0.3 bclip=0.5 \
       d1=1.5 d2=1.5 d3=0.16 size1=5 size2=10 size3=5 \
       angle=45 whls=0.666666,.5,1 bhls=0,.5,1.1 \
       d1num=5 d2num=5 d3num=2 \
       f3=0.96 f3num=1 \
       legend=1 lstyle=horibottom units="Resolution in time (s)" \
       label1="Along dip (km)" label2="Along strike (km)" label3="Time (s)" \
       labelsize=55 >block_time.ps
#subset <cube_temp n1=35 n2=35 n3=35 n1s=17 if3s=17 >block2 
#pscube <block2 n1=17 n2=35 n3=18 wclip=0.0 bclip=0.7 d1=5.0 d2=5.0 d3=5.0 size1=5 size2=10 size3=5 angle=45 whls=0.666666,.5,1 bhls=0,.5,1 legend=1 labelsize=45 d1num=1000 d2num=1000 d3num=1000 label2=X label3=Y >block2.ps
#~/Desktop/cwp/bin/psmerge in=block1.ps translate=0,0 in=block2.ps translate=3.53,8.53 >merge.ps 
