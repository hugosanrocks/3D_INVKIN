cwp_dir=~/Desktop/cwp/bin


#cwp_dir=$1

#number of snapshots
nsamp=8
#sample incremetn
snapleg=1

#output file
output=../graphics/snapshots.ps

#time samples of each snapshot
samples=(1 6 11 16 21 26 31 36)

#time step of inversion
timestep=0.25

#file with epicenter location
fileepi=epi.bin
#file name for epicenter plot
epips=epi.ps

#samples along dip and strike
#ny=171  #dip
#nx=351  #strike
ny=155
nx=340

#maximum value for sliprate
maxslip=1.2
#maximum value for weights
maxweig=1.2
#minimum value
minv=0.0
#increment of weight colorbar tics
dwei=1.
#increment of sliprate colorbar tics
dslip=0.25

#width and height of plot box
widt=9.0
heig=4.5

#Title
titlsiz=60
titleps=title.ps

#labels
l1="Target-SIV1"
l2="Inversion"
l3="Weight"
label1=label1.ps
label2=label2.ps
label3=label3.ps
backcolor=black
textcolor=white
labelsiz=40

#Plot hypocenter
  $cwp_dir/psgraph <$fileepi \
  n=1 pairs=1 mark=7 marksize=10 \
  x1beg=0 x1end=36 x2beg=18 x2end=0 \
  f1num=-20 f2num=-20 d1num=200 d2num=200 \
  wbox=$widt hbox=$heig \
  linecolor=black >$epips

#3 labels
#title
    $cwp_dir/pslabel t="Initial vs Final Prior model" size=$titlsiz \
    bcolor=$textcolor tcolor=$backcolor > $titleps
#labels
    $cwp_dir/pslabel t=$l1 size=$labelsiz \
    bcolor=$backcolor tcolor=$textcolor > $label1
    $cwp_dir/pslabel t=$l2 size=$labelsiz \
    bcolor=$backcolor tcolor=$textcolor > $label2
    $cwp_dir/pslabel t=$l3 size=$labelsiz \
    bcolor=$backcolor tcolor=$textcolor > $label3


#3 plots per snapshot
for i in $(seq 1 $nsamp) ; do
  samp=${samples[i-1]}
  time=$(expr $timestep*$samp | bc)
  time=$(expr $time-$timestep | bc)
  echo $time

  if [ $i -eq $snapleg ] ; then
  #write axes only for first snap shot
  #inversion snap
  filein1=$(printf "file_inv_%01i.bin" $i)
  fileout=$(printf "surface_inv%01i.ps" $i)
  fileout1=$(printf "surface_inv_%01i.ps" $i)
  $cwp_dir/psimage <$filein1 \
       n1=$ny n2=$nx \
       d1=0.1 d2=0.1 d1num=5 d2num=10 \
       lheight=0.4 lwidth=11 lx=4.5 ly=-28 \
       whls=0.666666,.5,1 bhls=0,.5,1.1 \
       width=$widt height=$heig legend=0 lstyle=horibottom \
       lbeg=$minv lend=$maxslip ldnum=$dslip units="Slip rate (m/s)" \
       wclip=$minv bclip=$maxslip \
       label1="Along dip (km)" label2="Along strike (km)" \
       labelsize=45 >$fileout
   $cwp_dir/psmerge in=$fileout \
   translate=0,0 in=$epips translate=0.0,0.0 >$fileout1
  else
  filein1=$(printf "file_inv_%01i.bin" $i)
  fileout=$(printf "surface_inv%01i.ps" $i)
  fileout1=$(printf "surface_inv_%01i.ps" $i)
  $cwp_dir/psimage <$filein1 \
       n1=$ny n2=$nx \
       d1=0.1 d2=0.1 d1num=100 d2num=100 \
       f1num=-20 f2num=-20 lheight=4.0 \
       whls=0.666666,.5,1 bhls=0,.5,1.1 \
       width=$widt height=$heig \
       lbeg=$minv lend=$maxslip ldnum=0.25 units="Weight" \
       bclip=$maxslip wclip=$minv \
       labelsize=45 >$fileout
   $cwp_dir/psmerge in=$fileout \
   translate=0,0 in=$epips translate=0.0,0.0 >$fileout1
  fi

  #prior snap
  filein2=$(printf "file_pri_%01i.bin" $i)
  fileout=$(printf "surface_pri%01i.ps" $i)
  fileout2=$(printf "surface_pri_%01i.ps" $i)
  $cwp_dir/psimage <$filein2 \
       n1=$ny n2=$nx \
       d1=0.1 d2=0.1 d1num=100 f1num=-20 \
       d2num=100 f2num=-20 \
       whls=0.666666,.5,1 bhls=0,.5,1.1 \
       bclip=$maxslip wclip=$minv \
       width=$widt height=$heig \
       labelsize=45 >$fileout
   $cwp_dir/psmerge in=$fileout \
   translate=0,0 in=$epips translate=0.0,0.0 >$fileout2

  #weight snap
  if [ $i -eq $snapleg ] ; then
  filein3=$(printf "file_wei_%01i.bin" $i)
  fileout=$(printf "surface_wei%01i.ps" $i)
  fileout3=$(printf "surface_wei_%01i.ps" $i)
  $cwp_dir/psimage <$filein3 \
       n1=$ny n2=$nx \
       legend=1 lheight=0.4 lwidth=11 lx=26.5 ly=-17.5 \
       lstyle=horibottom \
       d1=0.1 d2=0.1 \
       d1num=100 f1num=-20 d2num=100 f2num=-20 \
       whls=0.666666,.5,1 bhls=0,.5,1.1 \
       width=$widt height=$heig \
       bclip=$maxweig wclip=$minv \
       lbeg=$minv lend=$maxweig ldnum=$dwei units="Slip rate (m/s)" \
       title="t = $time sec" titlesize=50 \
       labelsize=45 >$fileout
   $cwp_dir/psmerge in=$fileout \
   translate=0,0 in=$epips translate=0.0,0.0 >$fileout3
  else
  filein3=$(printf "file_wei_%01i.bin" $i)
  fileout=$(printf "surface_wei%01i.ps" $i)
  fileout3=$(printf "surface_wei_%01i.ps" $i)
  $cwp_dir/psimage <$filein3 \
       n1=$ny n2=$nx \
       d1=0.1 d2=0.1 \
       d1num=100 f1num=-20 d2num=100 f2num=-20 \
       whls=0.666666,.5,1 bhls=0,.5,1.1 \
       width=$widt height=$heig \
       bclip=$maxweig wclip=$minv \
       lbeg=$minv lend=$maxweig ldnum=$dwei units="Weight" \
       title="t = $time sec" titlesize=50 \
       labelsize=45 >$fileout
   $cwp_dir/psmerge in=$fileout \
   translate=0,0 in=$epips translate=0.0,0.0 >$fileout3
  fi


  #merging snaps
  filemerge=$(printf "surface_merge_%01i.ps" $i)
  $cwp_dir/psmerge \
  in=$fileout3 translate=0,0 \
  in=$label3 translate=1.5,5.7 \
  in=$fileout2 translate=0,5.3 \
  in=$label2 translate=1.5,11.0 \
  in=$fileout1 translate=0,10.6 \
  in=$label1 translate=1.5,16.3 >$filemerge

done

#final merge
$cwp_dir/psmerge \
in=surface_merge_5.ps translate=0,0 \
in=surface_merge_6.ps translate=11,0.0 \
 in=surface_merge_7.ps translate=22,0.0 \
in=surface_merge_8.ps translate=33,0.0 \
in=surface_merge_1.ps translate=0.0,17 \
in=surface_merge_2.ps translate=11,17.0 \
in=surface_merge_3.ps translate=22,17.0 \
in=surface_merge_4.ps translate=33,17.0 \
in=$titleps translate=16.5,34 >$output

#clean ps files
#rm surface*.ps lab*ps tit*ps epi*ps
#convert to pdf
#ps2pdf -dEPSCrop snapshots.ps snapshots.pdf

$cwp_dir/psmerge \
 in=surface_inv_7.ps translate=0,0.0 \
 in=surface_inv_4.ps translate=0,5.0 \
 in=surface_inv_1.ps translate=0,10.0 > out.ps

