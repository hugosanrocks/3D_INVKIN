#!/bin/bash
nsta=1 #30 all traces
used=0

obs1=graphics/seismo_to_plote_01.obs
obs2=graphics/seismo_to_plote_23.obs
obs3=graphics/seismo_to_plotv_01.obs
head=graphics/seismo_all.head
esc1=graphics/escaobs_1.xy
esc2=graphics/escasyn_1.xy
esc3=graphics/escasyn2_1.xy
esc4=graphics/escasyn3.xy
esc5=graphics/escasyn4.xy

# find start and end time and a reasonable step 
t_start=`grep T_START ${head} | awk '{print $NF}'`
t_end=`grep T_END ${head} | awk '{print $NF}'`
t_step=`echo $t_end $t_start | awk '{print int(($1-$2)/3)}'`
t_step=20

#all stations distance and used
array=(    8  9  10   17    3    5    4    1    2   15   16   22   26   21   23   24   27   25 30   7     11   13   12   14    6   19   18   20   28   29 )
nom=(KMMH14 KMMH16 KMM003 TMC KMMH03 KMMH09 KMMH06 KMMH01 KMMH02 MYZ001 TKD FKOH10 MYZH04 OITH08 SAGH04 FKOH03 MYZH08 SAGH02 KGSH08 KMMH14 KMM012 KMM014 KMM013 KMM018 KMMH11 KMMH13 KMMH12 KMMH15 NGSH06 KGSH04)


# default settings
paper="a4+"
orient="-P"
proj=X3/1.5
origin="-Y16"
shift="-X3.3"

gmtset PAPER_MEDIA $paper MEASURE_UNIT inch HEADER_FONT_SIZE 12p LABEL_FONT_SIZE 12p ANNOT_FONT_SIZE_PRIMARY 12p
nhead=4
red=255/0/0
black=0/0/0
blue=0/0/255
green=00/200/00
paleblue=220/220/255
aqua=153/50/204

# set output filename
out=graphics/all_seismograms.eps
outpdf=graphics/all_seismograms.pdf

# set region
#max=`grep PLOT_MAX ${obs1} | awk '{print $NF}'`
min=1
max=3.9 #66 all traces
region=$t_start/$t_end/$min/$max

psbasemap -R$region -J$proj -B${t_step}::.${label}:"Time (s)":/S -K $orient $origin -Y1 > $out

#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 16 3.7 18 0 0 LT 1 10 l 
@;$black;Inverted @;;
END

pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 3.3 15 0 0 LT 1 10 l 
@;$black;E-W @;;
END


psxy $obs1 -R$region -J$proj -H$nhead -W3 -O -K >> $out


for i in $(seq 1 $nsta) ; do
filename=$(printf "graphics/seismo_to_plote_%02i" $i)
filename2=$(printf "graphics/seismo_to_plote2_%02i" $i)
basename=${filename}
basename2=${filename2}
label=`basename $basename`
#files
obs1=${basename}.obs
syn1=${basename}.syn
synn1=${basename2}.syn

maxamp=`grep PLOT_MAX ${obs1} | awk '{print $NF}'`


# set region
#max=`grep PLOT_MAX ${obs1} | awk '{print $NF}'`
region=$t_start/$t_end/$min/$max
if [ $i -le 1 ] ; then
skip=1.8
else
skip=1.6
fi
#jump=$(expr "$i" \* "$skip")
jump=$(expr $i*$skip | bc)
shift1=0.70
shift2=0.4
yy1=`echo $jump + $shift1 | bc`
yy2=`echo $jump - $shift2 | bc`


if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  psxy $obs1 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
#  psxy $syn1 -R$region -J$proj -H$nhead -W0.5p,$red,- -O -K >> $out
#  psxy $synn1 -R$region -J$proj -H$nhead -W0.5p,$aqua,- -O -K >> $out

else
  #Plot the ones not used for inversion
  psxy $obs1 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn1 -R$region -J$proj -H$nhead -W0.5p,$red,- -O -K >> $out
  psxy $synn1 -R$region -J$proj -H$nhead -W0.5p,$aqua,- -O -K >> $out
fi


#STATION
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 45 $yy1 15 0 0 LT 1 10 l 
@;$black; @;;
END
#NORMALIZATION
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 55 $yy2 15 0 0 LT 1 10 l 
@;$black;$maxamp @;;
END
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 $yy2 15 0 0 LT 1 10 l 
@;$black;${nom[$i-1]} @;;
END

done




psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 16 3.7 18 0 0 LT 1 10 l 
@;$black;Predicted @;;
END


#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 3.3 15 0 0 LT 1 10 l
@;$black;E-W @;;
END



psxy $obs2 -R$region -J$proj -H$nhead -W3 -O -K >> $out

nsta=23
#for i in $(seq 23 $nsta) ; do
filename=$(printf "graphics/seismo_to_plote_23")
filename2=$(printf "graphics/seismo_to_plote2_23")
basename=${filename}
basename2=${filename2}
label=`basename $basename`
#files
obs2=${basename}.obs
syn2=${basename}.syn
synn2=${basename2}.syn

#################################
# do plot for STA/LTA
#################################
nhead=4

# set region
region=$t_start/$t_end/$min/$max

#if [ $i -le $used ] ; then
  #Plot the ones used for inversion
#  psxy $obs2 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
#  psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$blue,- -O -K >> $out
#  psxy $synn2 -R$region -J$proj -H$nhead -W0.5p,$green,- -O -K >> $out
#else
  #Plot the ones not used for inversion
  psxy $obs2 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$red,- -O -K >> $out
  psxy $synn2 -R$region -J$proj -H$nhead -W0.5p,$aqua,- -O -K >> $out
#fi

#STATION
i=23
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 45 $yy1 15 0 0 LT 1 10 l 
@;$black; @;;
END
#NORMALIZATION
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 55 $yy2 15 0 0 LT 1 10 l 
@;$black;0.02 @;;
END
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 $yy2 15 0 0 LT 1 10 l 
@;$black;${nom[$i-1]} @;;
END


#done


psxy $esc1 -R$region -J$proj -W3,$black -O -K -V >> $out
psxy $esc2 -R$region -J$proj -W0.5p,$red,- -O -K -V >> $out
psxy $esc3 -R$region -J$proj -W0.5p,$aqua,- -O -K -V >> $out
#psxy $esc4 -R$region -J$proj -W0.5p,$green,- -O -K -V >> $out
#psxy $esc5 -R$region -J$proj -W0.5p,$aqua,- -O -K -V >> $out


pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 55 3.3 12 0 0 LT 1 10 l 
@;$black;Obs @;;
END
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 55 2.9 12 0 0 LT 1 10 l 
@;$black;PIS2 @;;
END
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 55 2.5 12 0 0 LT 1 10 l 
@;$red;SIS2 @;;
END



psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O >> $out


ps2pdf $out $outpdf
echo $out $outpdf



