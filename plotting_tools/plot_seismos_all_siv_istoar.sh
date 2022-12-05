#!/bin/bash
nsta=40
used=40

obs1=graphics/seismo_to_plote_01.obs
obs2=graphics/seismo_to_plotn_01.obs
obs3=graphics/seismo_to_plotv_01.obs
scale1=graphics/scale1.txt
scale2=graphics/scale2.txt

# find start and end time and a reasonable step 
t_start=`grep T_START ${obs1} | awk '{print $NF}'`
t_end=`grep T_END ${obs1} | awk '{print $NF}'`
t_step=`echo $t_end $t_start | awk '{print int(($1-$2)/6)}'`

# default settings
paper="a3"
orient=""
proj=X2/6
origin="-Y1"
shift="-Xr2.2"
shifty="-Ya20"

#gmt gmtset PAPER_MEDIA $paper MEASURE_UNIT inch HEADER_FONT_SIZE 14p LABEL_FONT_SIZE 16p
gmt gmtset PS_MEDIA $paper PROJ_LENGTH_UNIT inch FONT_TITLE 20p FONT_LABEL 18p MAP_TITLE_OFFSET 0.1i

nhead=4
red=255/0/0
black=0/0/0
blue=0/0/255
green=00/200/00
paleblue=220/220/255

# set output filename
out=graphics/all_seismograms.eps
outpdf=graphics/all_seismograms.pdf

# set region
min=0
max=43
nhead=3
region=$t_start/$t_end/$min/$max

#start the map
gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -K $orient $origin > $out

#title of traces
gmt pstext -R$region -J$proj -N -M -F+f20,Times-Roman+jCT -O -K << EOF >> $out
> 15 $max 21 50p 1i l
@;$black;E-W @;;
EOF

#plot first trace
gmt psxy $obs1 -R$region -J$proj -h$nhead -W1 -O -K >> $out

nsta=20
nsta1=21
nsta2=40

#plot first column
for i in $(seq 1 $nsta) ; do
 filename=$(printf "graphics/seismo_to_plote_%02i" $i)
 basename=${filename}
 label=`basename $basename`
 #files
 obs1=${basename}.obs
 syn1=${basename}.syn
 #maximum amplitude per station
 maxamp=`grep PLOT_MAX ${obs1} | awk '{print $NF}'`
 # set region
 #region=$t_start/$t_end/$min/$max
 #Jump to put station names
 skip=2
 jump=$(expr "$i" \* "$skip")
 shift1=0.8
 shift2=0.2
 yy1=`echo $jump + $shift1 | bc`
 yy2=`echo $jump - $shift2 | bc`
 sta_i=`echo STA $i`

 if [ $i -le $used ] ; then
   #Plot the ones used for inversion
   gmt psxy $obs1 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
   gmt psxy $syn1 -R$region -J$proj -h$nhead -W0.5p,$blue,- -O -K >> $out
 else
   #Plot the ones not used for inversion
   gmt psxy $obs1 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
   gmt psxy $syn1 -R$region -J$proj -h$nhead -W0.5p,$red,- -O -K >> $out
 fi

 #To plot all
 #gmt psxy $obs1 -R$region -J$proj -h$nhead -W3 -O -K >> $out
 #gmt psxy $syn1 -R$region -J$proj -h$nhead -W0.5p,$blue -O -K >> $out

 #Station names
gmt pstext -R$region -J$proj -N -M -F+f10,Times-Roman+jLT -O -K << EOF >> $out
> 26 $yy1 7 50p 1i l
@;$black;$sta_i @;;
EOF

#maximum amplitude text
gmt pstext -R$region -J$proj -N -M -F+f10,Times-Roman+jLT -O -K << EOF >> $out
> 26 $yy2 7 50p 1i l
@;$black;$maxamp (m/s) @;;
EOF


done

#gmt pstext -R$region -J$proj -Baf -F+f40+cTC+t"Inner Title" -O -K -D0/-0.5 >> $out

#next column
gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

#TITLE
gmt pstext -R$region -J$proj -N -M -F+f20,Times-Roman+jCT -O -K << EOF >> $out
> 15 $max 21 50p 1i l
@;$black;S-N @;;
EOF

#first trace
gmt psxy $obs2 -R$region -J$proj -h$nhead -W1 -O -K >> $out

for i in $(seq 1 $nsta) ; do
filename=$(printf "graphics/seismo_to_plotn_%02i" $i)
basename=${filename}
label=`basename $basename`
#files
obs2=${basename}.obs
syn2=${basename}.syn

if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  gmt psxy $obs2 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn2 -R$region -J$proj -h$nhead -W0.5p,$blue,- -O -K >> $out
else
  #Plot the ones not used for inversion
  gmt psxy $obs2 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn2 -R$region -J$proj -h$nhead -W0.5p,$red,- -O -K >> $out
fi

#Plot all
#gmt psxy $obs2 -R$region -J$proj -h$nhead -W3 -O -K >> $out
#gmt psxy $syn2 -R$region -J$proj -h$nhead -W0.5p,$blue -O -K >> $out

done

#third column
gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out
#TITLE
gmt pstext -R$region -J$proj -N -M -F+f20,Times-Roman+jCT -O -K << EOF >> $out
> 15 $max 21 50p 1i l
@;$black;U-D @;;
EOF

#first trace
gmt psxy $obs3 -R$region -J$proj -h$nhead -W1 -O -K >> $out


for i in $(seq 1 $nsta) ; do
filename=$(printf "graphics/seismo_to_plotv_%02i" $i)
basename=${filename}
label=`basename $basename`

#files
obs3=${basename}.obs
syn3=${basename}.syn

#plot traces
if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  gmt psxy $obs3 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn3 -R$region -J$proj -h$nhead -W0.5p,$blue,- -O -K >> $out
else
  #Plot the ones not used for inversion
  gmt psxy $obs3 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn3 -R$region -J$proj -h$nhead -W0.5p,$red,- -O -K >> $out
fi

#Plot all
#gmt psxy $obs3 -R$region -J$proj -h$nhead -W3 -O -K >> $out
#gmt psxy $syn3 -R$region -J$proj -h$nhead -W0.5p,$blue -O -K >> $out

done


#repeats for 4-6 column

obs1=graphics/seismo_to_plote_21.obs
obs2=graphics/seismo_to_plotn_21.obs
obs3=graphics/seismo_to_plotv_21.obs
head=graphics/seismo_all.head


gmt psbasemap -R$region -J$proj -B${t_step}::.:"Time (s)":/S -O -K $shift >> $out

#title of traces
gmt pstext -R$region -J$proj -N -M -F+f20,Times-Roman+jCT -O -K << EOF >> $out
> 15 $max 21 50p 1i l
@;$black;E-W @;;
EOF



gmt psxy $obs1 -R$region -J$proj -h$nhead -W1 -O -K >> $out


for i in $(seq $nsta1 $nsta2) ; do
filename=$(printf "graphics/seismo_to_plote_%02i" $i)
basename=${filename}
label=`basename $basename`
#files
obs1=${basename}.obs
syn1=${basename}.syn
sta_i=`echo STA $i`

maxamp=`grep PLOT_MAX ${obs1} | awk '{print $NF}'`


#jump for names
skip=2
jump=$(expr "$i" \* "$skip")
shift1=0.8
shift2=0.2
shifty=40
yy1=`echo $jump + $shift1 | bc`
yy2=`echo $jump - $shift2 | bc`
yy1=`echo $yy1 - $shifty | bc`
yy2=`echo $yy2 - $shifty | bc`

if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  gmt psxy $obs1 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn1 -R$region -J$proj -h$nhead -W0.5p,$blue,- -O -K >> $out
else
  #Plot the ones not used for inversion
  gmt psxy $obs1 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn1 -R$region -J$proj -h$nhead -W0.5p,$red,- -O -K >> $out
fi

#To plot all
#gmt psxy $obs1 -R$region -J$proj -h$nhead -W3 -O -K >> $out
#gmt psxy $syn1 -R$region -J$proj -h$nhead -W0.5p,$blue -O -K >> $out

 #Station names
gmt pstext -R$region -J$proj -N -M -F+f10,Times-Roman+jLT -O -K << EOF >> $out
> 26 $yy1 7 50p 1i l
@;$black;$sta_i @;;
EOF

#maximum amplitude text
gmt pstext -R$region -J$proj -N -M -F+f10,Times-Roman+jLT -O -K << EOF >> $out
> 26 $yy2 7 50p 1i l
@;$black;$maxamp (m/s) @;;
EOF

done




gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

#TITLE
gmt pstext -R$region -J$proj -N -M -F+f20,Times-Roman+jCT -O -K << EOF >> $out
> 15 $max 21 50p 1i l
@;$black;S-N @;;
EOF

gmt psxy $obs2 -R$region -J$proj -h$nhead -W1 -O -K >> $out

for i in $(seq $nsta1 $nsta2) ; do
filename=$(printf "graphics/seismo_to_plotn_%02i" $i)
basename=${filename}
label=`basename $basename`
#files
obs2=${basename}.obs
syn2=${basename}.syn

if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  gmt psxy $obs2 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn2 -R$region -J$proj -h$nhead -W0.5p,$blue,- -O -K >> $out
else
  #Plot the ones not used for inversion
  gmt psxy $obs2 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn2 -R$region -J$proj -h$nhead -W0.5p,$red,- -O -K >> $out
fi

#Plot all
#gmt psxy $obs2 -R$region -J$proj -h$nhead -W3 -O -K >> $out
#gmt psxy $syn2 -R$region -J$proj -h$nhead -W0.5p,$blue -O -K >> $out

done



gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out
#TITLE
gmt pstext -R$region -J$proj -N -M -F+f20,Times-Roman+jCT -O -K << EOF >> $out
> 15 $max 21 50p 1i l
@;$black;U-D @;;
EOF
gmt psxy $obs3 -R$region -J$proj -h$nhead -W1 -O -K >> $out

for i in $(seq $nsta1 $nsta2) ; do
filename=$(printf "graphics/seismo_to_plotv_%02i" $i)
basename=${filename}
label=`basename $basename`

#files
obs3=${basename}.obs
syn3=${basename}.syn

if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  gmt psxy $obs3 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn3 -R$region -J$proj -h$nhead -W0.5p,$blue,- -O -K >> $out

else
  #Plot the ones not used for inversion
  gmt psxy $obs3 -R$region -J$proj -h$nhead -W1,$black -O -K >> $out
  gmt psxy $syn3 -R$region -J$proj -h$nhead -W0.5p,$red,- -O -K >> $out
fi

#Plot all
#gmt psxy $obs3 -R$region -J$proj -h$nhead -W3 -O -K >> $out
#gmt psxy $syn3 -R$region -J$proj -h$nhead -W0.5p,$blue -O -K >> $out

done


#scale
#gmt psxy $scale1 -R$region -J$proj -W3,$black -O -K -V>> $out
#gmt psxy $scale2 -R$region -J$proj -W0.5p,$blue -O -K -V>> $out

#gmt pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 15 42 13 0 0 LT 1 10 l 
#@;$black;OBS @;;
#END
#gmt pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 15 43.5 13 0 0 LT 1 10 l 
#@;$black;SYN @;;
#END

#END THE PLOT
gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O >> $out


#gmt ps2pdf $out $outpdf
#echo $out $outpdf


