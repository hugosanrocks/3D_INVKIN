subset <cube_src n1=12 n2=24 n3=57 >block1 
pscube <block1 n1=12 n2=24 n3=57 \
       wclip=0. bclip=2.5 \
       d1=1.5 d2=1.5 d3=0.16 size1=5 size2=10 size3=3 \
       angle=45 whls=0.666666,.5,1 bhls=0,.5,1.1 \
       d1num=5 d2num=5 d3num=2 \
       legend=1 lstyle=horibottom units="Slip rate [m/s]" \
       label1="Along dip [km]" label2="Along strike [km]" label3="Time [s]" \
       labelsize=45 >block_src.ps
