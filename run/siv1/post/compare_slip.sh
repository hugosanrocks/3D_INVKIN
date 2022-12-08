
cwp_dir=$1

#epicenter file names
fileepi=epi.bin
epips=epi.ps
filenodes=nodes.bin
nodesps=nodes.ps
elems=elements.ps
nnode=500
#samples along dip and strike
ny=171
nx=351
#ny=161
#nx=341
#ny=160
#nx=340
#increment along strike and dip
dx=0.1
dy=0.1
#increment for tics
dxtic=10
dytic=5
#maximum/minimum slip value
max=1.8
min=0.0
maxw=1.8
minw=0
dcolw=0.3
dcol=0.3
#box width and height
wid=12.6
hei=6.5
nelem=982

title1='SIV1'
title2="Inversion"
title3='Test2'
title4='Test3'
titlesiz=50
outfile=../graphics/final_slip_compared.ps

#labels
label1=label1.ps
label2=label2.ps
label3=label3.ps
label4=label4.ps
label5=label5.ps
backcolor=black
textcolor=white
labelsiz=40
    $cwp_dir/pslabel t=$title1 size=$labelsiz \
    bcolor=$textcolor tcolor=black f=Helvetica-Bold > $label1
    $cwp_dir/pslabel t=$title2 size=$labelsiz \
    bcolor=$textcolor tcolor=black f=Helvetica-Bold > $label2

    $cwp_dir/pslabel t="0-6 sec rupture" size=35 \
    bcolor=$blackcolor tcolor=$textcolor f=Helvetica > $label3
    $cwp_dir/pslabel t="Prior model" size=35 \
    bcolor=$textcolor tcolor=$backcolor f=Helvetica > $label4
    $cwp_dir/pslabel t="Weight" size=35 \
    bcolor=$textcolor tcolor=$backcolor f=Helvetica > $label5



#Hypocenter
  $cwp_dir/psgraph <$fileepi \
  n=1 pairs=1 mark=7 marksize=18 \
  x1beg=0 x1end=36 x2beg=18 x2end=0 \
  f1num=-20 f2num=-20 d1num=200 d2num=200 \
  wbox=$wid hbox=$hei \
  linecolor=black >$epips

  $cwp_dir/psgraph <$fileepi \
  n=1 pairs=1 mark=7 marksize=18 \
  grid1=solid grid2=solid x1beg=0 x1end=36 x2beg=18 x2end=0 \
  f1num=-20 f2num=-20 d1num=200 d2num=200 \
  wbox=$wid hbox=$hei \
  linecolor=black >epi2.ps


#nodes
#  $cwp_dir/psgraph <$filenodes \
#  n=$nnode pairs=2 mark=6 marksize=5 \
#  lineon=1 lineoff=5000 \
#  x1beg=0 x1end=35 x2beg=17 x2end=0 \
#  f1num=-20 f2num=-20 d1num=200 d2num=200 \
#  wbox=$wid hbox=$hei \
#  linecolor=black >nodes.ps

#cp $epips $elems
#elements
#for i in $(seq 1 $nelem) ; do
#  fileelem=$(printf "line_%04i.bin" $i)
#  lineelem=$(printf "line_%04i.ps" $i)
#  $cwp_dir/psgraph <$fileelem \
#  n=4 pairs=2 linewidth=0.1 \
#  x1beg=0 x1end=35 x2beg=17 x2end=0 \
#  f1num=-20 f2num=-20 d1num=200 d2num=200 \
#  wbox=$wid hbox=$hei \
#  linecolor=white >$lineelem

#  $cwp_dir/psmerge \
#  in=$lineelem translate=0,0 \
#  in=$elems translate=0,0 >$epips
#  cp $epips $elems
#done

#$cwp_dir/psmerge \
# in=hypo.ps translate=0,0 \
# in=point.ps translate=0.0,0.0 >$epips


#1 subplot
  filein=$(printf "file_%01i.bin" 1)
  fileout=$(printf "surface%01i.ps" 1)
  fileout1=$(printf "surface_%01i.ps" 1)
  $cwp_dir/psimage <$filein \
       n1=$ny n2=$nx \
       d1=$dy d2=$dx d1num=5 d2num=10 \
       f1num=0 f2num=0 \
       whls=0.7,0.5,1 bhls=0,0.3,1.1 \
       width=$wid height=$hei  \
       lbeg=$min lend=$max ldnum=$dcol units="Cumulative slip (m)" \
       wclip=$min bclip=$max legend=1 lstyle=horibottom \
       lheight=0.4 lwidth=12 lx=9.0 ly=--0.6 \
       label1="" label2="" \
       labelsize=45 >$fileout
   $cwp_dir/psmerge \
   in=$fileout translate=0,0 \
   in=$epips translate=0.0,0.0 \
   in=$label1 translate=1.5,0.65 >$fileout1

#2 subplot
  filein=$(printf "file_%01i.bin" 2)
  fileout=$(printf "surface%01i.ps" 2)
  fileout2=$(printf "surface_%01i.ps" 2)
  $cwp_dir/psimage <$filein \
       n1=$ny n2=$nx \
       d1=$dy d2=$dx d1num=$dytic \
       d2num=10 \
       whls=0.7,0.5,1 bhls=0,0.3,1.1 \
       width=$wid height=$hei  \
       lbeg=$minw lend=$maxw ldnum=$dcolw units="Total weight" \
       wclip=$minw bclip=$maxw lstyle=vertright \
       lheight=4.5 lwidth=0.4 lx=10.8 ly=1.5 \
       labelsize=45 >$fileout
   $cwp_dir/psmerge \
   in=$fileout translate=0,0 \
   in=epi2.ps translate=0.0,0.0 \
   in=$label2 translate=1.5,0.65 >$fileout2


#for slides
   $cwp_dir/psmerge in=$fileout1 translate=0,0 \
	            in=$fileout2 translate=15,0 > $outfile

#clean all ps files
#rm *.ps line*bin
#convert to pdf
#ps2pdf -dEPSCrop $filemerge final_merge.pdf
