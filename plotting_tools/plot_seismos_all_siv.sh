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
paper="a4+"
orient="-P"
proj=X2/6
origin="-Y16"
shift="-X2.2"
shifty="-Y20"

gmtset PAPER_MEDIA $paper MEASURE_UNIT inch HEADER_FONT_SIZE 14p LABEL_FONT_SIZE 16p
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
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -K $orient $origin -Y1 > $out

#title of traces
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 $max 21 0 0 LT 1 10 l 
@;$black;E-W @;;
END

#plot first trace
psxy $obs1 -R$region -J$proj -H$nhead -W3 -O -K >> $out

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

 if [ $i -le $used ] ; then
   #Plot the ones used for inversion
   psxy $obs1 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
   psxy $syn1 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out
 else
   #Plot the ones not used for inversion
   psxy $obs1 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
   psxy $syn1 -R$region -J$proj -H$nhead -W0.5p,$red -O -K >> $out
 fi

 #To plot all
 #psxy $obs1 -R$region -J$proj -H$nhead -W3 -O -K >> $out
 #psxy $syn1 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out

 #Station names
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 26 $yy1 7 0 0 LT 1 10 l 
@;$black;STA $i @;;
END
#maximum amplitude text
#pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 26 $yy2 7 0 0 LT 1 10 l 
#@;$black;$maxamp (m/s) @;;
#END

pstext -R$region -J$proj -N -O -K <<  END >> $out
26 $yy2 7 0 0 LT $maxamp (m/s)
END

done



#next column
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 $max 21 0 0 LT 1 10 l 
@;$black;S-N @;;
END

#first trace
psxy $obs2 -R$region -J$proj -H$nhead -W3 -O -K >> $out

for i in $(seq 1 $nsta) ; do
filename=$(printf "graphics/seismo_to_plotn_%02i" $i)
basename=${filename}
label=`basename $basename`
#files
obs2=${basename}.obs
syn2=${basename}.syn

if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  psxy $obs2 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out
else
  #Plot the ones not used for inversion
  psxy $obs2 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$red -O -K >> $out
fi

#Plot all
#psxy $obs2 -R$region -J$proj -H$nhead -W3 -O -K >> $out
#psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out

done

#third column
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out
#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 $max 21 0 0 LT 1 10 l 
@;$black;U-D @;;
END
#first trace
psxy $obs3 -R$region -J$proj -H$nhead -W3 -O -K >> $out


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
  psxy $obs3 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out
else
  #Plot the ones not used for inversion
  psxy $obs3 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$red -O -K >> $out
fi

#Plot all
#psxy $obs3 -R$region -J$proj -H$nhead -W3 -O -K >> $out
#psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out

done


#repeats for 4-6 column

obs1=graphics/seismo_to_plote_21.obs
obs2=graphics/seismo_to_plotn_21.obs
obs3=graphics/seismo_to_plotv_21.obs
head=graphics/seismo_all.head


psbasemap -R$region -J$proj -B${t_step}::.:"Time (s)":/S -O -K $shift >> $out
#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 $max 21 0 0 LT 1 10 l 
@;$black;E-W @;;
END

psxy $obs1 -R$region -J$proj -H$nhead -W3 -O -K >> $out


for i in $(seq $nsta1 $nsta2) ; do
filename=$(printf "graphics/seismo_to_plote_%02i" $i)
basename=${filename}
label=`basename $basename`
#files
obs1=${basename}.obs
syn1=${basename}.syn

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
  psxy $obs1 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn1 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out
else
  #Plot the ones not used for inversion
  psxy $obs1 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn1 -R$region -J$proj -H$nhead -W0.5p,$red -O -K >> $out
fi

#To plot all
#psxy $obs1 -R$region -J$proj -H$nhead -W3 -O -K >> $out
#psxy $syn1 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out

#STATION
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 26 $yy1 7 0 0 LT 1 10 l 
@;$black;STA $i @;;
END
#NORMALIZATION
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 26 $yy2 7 0 0 LT 1 10 l 
@;$black;$maxamp (m/s) @;;
END

done




psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 $max 21 0 0 LT 1 10 l 
@;$black;S-N @;;
END

psxy $obs2 -R$region -J$proj -H$nhead -W3 -O -K >> $out

for i in $(seq $nsta1 $nsta2) ; do
filename=$(printf "graphics/seismo_to_plotn_%02i" $i)
basename=${filename}
label=`basename $basename`
#files
obs2=${basename}.obs
syn2=${basename}.syn

if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  psxy $obs2 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out
else
  #Plot the ones not used for inversion
  psxy $obs2 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$red -O -K >> $out
fi

#Plot all
#psxy $obs2 -R$region -J$proj -H$nhead -W3 -O -K >> $out
#psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out

done



psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out
#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 $max 21 0 0 LT 1 10 l 
@;$black;U-D @;;
END
psxy $obs3 -R$region -J$proj -H$nhead -W3 -O -K >> $out

for i in $(seq $nsta1 $nsta2) ; do
filename=$(printf "graphics/seismo_to_plotv_%02i" $i)
basename=${filename}
label=`basename $basename`

#files
obs3=${basename}.obs
syn3=${basename}.syn

if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  psxy $obs3 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out

else
  #Plot the ones not used for inversion
  psxy $obs3 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$red -O -K >> $out
fi

#Plot all
#psxy $obs3 -R$region -J$proj -H$nhead -W3 -O -K >> $out
#psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out

done


#scale
psxy $scale1 -R$region -J$proj -W3,$black -O -K -V>> $out
psxy $scale2 -R$region -J$proj -W0.5p,$blue -O -K -V>> $out

pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 15 42 13 0 0 LT 1 10 l 
@;$black;OBS @;;
END
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 15 43.5 13 0 0 LT 1 10 l 
@;$black;SYN @;;
END

#END THE PLOT
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O >> $out


#ps2pdf $out $outpdf
#echo $out $outpdf


