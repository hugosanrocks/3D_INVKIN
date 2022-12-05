#!/bin/bash
nsta=30 #30 all traces
used=19

obs1=graphics/seismo_to_plote_01.obs
obs2=graphics/seismo_to_plotn_01.obs
obs3=graphics/seismo_to_plotv_01.obs
head=graphics/seismo_all.head
esc1=graphics/escaobs.xy
esc2=graphics/escasyn.xy
esc3=graphics/escasyn2.xy
esc4=graphics/escasyn3.xy
esc5=graphics/escasyn4.xy
line=graphics/linexy.xy


# find start and end time and a reasonable step 
t_start=`grep T_START ${head} | awk '{print $NF}'`
t_end=`grep T_END ${head} | awk '{print $NF}'`
t_step=`echo $t_end $t_start | awk '{print int(($1-$2)/6)}'`


#all statation in distance order
#array=(8 7 10 17 3 9 11 5 4 1 13 12 2 15 14 6 16 22 26 19 18 21 20 23 28 24 27 25 29 30)
#all stations distance and used
array=(    8  9  10   17    3    5    4    1    2   15   16   22   26   21   23   24   27   25 30   7     11   13   12   14    6   19   18   20   28   29 )
nom=(KMMH16 KMM003 KMM05 TMC KMMH03 KMMH09 KMMH06 KMMH01 KMMH02 MYZ001 TKD FKOH10 MYZH04 OITH08 SAGH04 FKOH03 MYZH08 SAGH02 KGSH08 KMMH14 KMM012 KMM014 KMM013 KMM018 KMMH11 KMMH13 KMMH12 KMMH15 NGSH06 KGSH04)

#array=(1 3 5 7 eleven)
#echo ${array[0]}

#for i in "${array4[@]}" 
#do 

# find number of windows
#nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`

# default settings
paper="a4+"
orient="-P"
proj=X2/7.5
origin="-Y16"
shift="-X2.3"

gmtset PAPER_MEDIA $paper MEASURE_UNIT inch HEADER_FONT_SIZE 14p LABEL_FONT_SIZE 16p
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
max=66 #66 all traces
region=$t_start/$t_end/$min/$max

psbasemap -R$region -J$proj -B${t_step}::.${label}:"Time (s)":/S -K $orient $origin -Y1 > $out
#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 65 21 0 0 LT 1 10 l 
@;$black;E-W @;;
END


psxy $line -R$region -J$proj -W0.5p,$aqua,*- -O -K >> $out


psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 65 21 0 0 LT 1 10 l 
@;$black;N-S @;;
END

#TITLE
pstext -R$region -J$proj -N -X-1i -Y0.3i -m -O -K <<  END >> $out
> 0 66 24 0 0 LT 1 10 l 
@;$black; 2016 Kumamoto misfit @;;
END


psxy $obs2 -R$region -J$proj -H$nhead -W3 -O -X1i -Y-0.3i -K >> $out

for i in $(seq 1 $nsta) ; do
filename=$(printf "graphics/seismo_to_plotn_%02i" $i)
filename2=$(printf "graphics/seismo_to_plotn2_%02i" $i)
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
#max=`grep PLOT_MAX ${obs2} | awk '{print $NF}'`
region=$t_start/$t_end/$min/$max

if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  psxy $obs2 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$blue,- -O -K >> $out
  psxy $synn2 -R$region -J$proj -H$nhead -W0.5p,$green,- -O -K >> $out
else
  #Plot the ones not used for inversion
  psxy $obs2 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$blue,- -O -K >> $out
#  psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$red,- -O -K >> $out
#  psxy $synn2 -R$region -J$proj -H$nhead -W0.5p,$aqua,- -O -K >> $out
  psxy $synn2 -R$region -J$proj -H$nhead -W0.5p,$green,- -O -K >> $out

fi

#Plot all
#psxy $obs2 -R$region -J$proj -H$nhead -W3 -O -K >> $out
#psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out

done

psxy $line -R$region -J$proj -W0.5p,$aqua,*- -O -K >> $out


#################################
# plot envelopes with windows
#################################

psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out
#TITLE
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 0 65 21 0 0 LT 1 10 l 
@;$black;U-D @;;
END
psxy $obs3 -R$region -J$proj -H$nhead -W3 -O -K >> $out

nhead=4
for i in $(seq 1 $nsta) ; do
filename=$(printf "graphics/seismo_to_plotv_%02i" $i)
basename=${filename}
filename2=$(printf "graphics/seismo_to_plotv2_%02i" $i)
basename2=${filename2}
label=`basename $basename`

#files
obs3=${basename}.obs
syn3=${basename}.syn
synn3=${basename2}.syn

# set region
#max=`grep PLOT_MAX ${obs3} | awk '{print $NF}'`
region=$t_start/$t_end/$min/$max
skip=2
jump=$(expr "$i" \* "$skip")
shift1=0.8
shift2=0.2
yy1=`echo $jump + $shift1 | bc`
yy2=`echo $jump - $shift2 | bc`


if [ $i -le $used ] ; then
  #Plot the ones used for inversion
  psxy $obs3 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$blue,- -O -K >> $out
  psxy $synn3 -R$region -J$proj -H$nhead -W0.5p,$green,- -O -K >> $out
else
  #Plot the ones not used for inversion
  psxy $obs3 -R$region -J$proj -H$nhead -W3,$black -O -K >> $out
  psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$blue,- -O -K >> $out
#  psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$red,- -O -K >> $out
#  psxy $synn3 -R$region -J$proj -H$nhead -W0.5p,$aqua,- -O -K >> $out
  psxy $synn3 -R$region -J$proj -H$nhead -W0.5p,$green,- -O -K >> $out

fi



#STATION
#pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 55 $yy1 9 0 0 LT 1 10 l 
#@;$black;${nom[$i-1]} @;;
#END
#NORMALIZATION
#pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 55 $yy2 9 0 0 LT 1 10 l 
#@;$black;$maxamp (m/s) @;;
#END



#Plot all
#psxy $obs3 -R$region -J$proj -H$nhead -W3 -O -K >> $out
#psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$blue -O -K >> $out

done



psxy $esc1 -R$region -J$proj -W3,$black -O -K -V >> $out
psxy $esc2 -R$region -J$proj -W0.5p,$blue,- -O -K -V >> $out
psxy $esc3 -R$region -J$proj -W0.5p,$green,- -O -K -V >> $out
#psxy $esc4 -R$region -J$proj -W0.5p,$green,- -O -K -V >> $out
#psxy $esc5 -R$region -J$proj -W0.5p,$aqua,- -O -K -V >> $out


#pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 31 67 15 0 0 LT 1 10 l 
#@;$black;Observed @;;
#END
#pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 31 67 15 0 0 LT 1 10 l 
#@;$black;Standard Inv. @;;
#END
#pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 31 67 15 0 0 LT 1 10 l 
#@;$black;Progressive Inv. @;;
#END
#pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 20.5 64.5 11 0 0 LT 1 10 l 
#@;$black;SIS2 @;;
#END
#pstext -R$region -J$proj -N -m -O -K <<  END >> $out
#> 45.5 64.5 11 0 0 LT 1 10 l 
#@;$black;PIS2 @;;
#END
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 40 67.0 11 0 0 LT 1 10 l 
@;$black;Observed @;;
END
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 40 65.0 11 0 0 LT 1 10 l 
@;$black;PIS-KUMA @;;
END
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 40 63.0 11 0 0 LT 1 10 l 
@;$black;Asanao & @;;
END
pstext -R$region -J$proj -N -m -O -K <<  END >> $out
> 40 62.0 11 0 0 LT 1 10 l 
@;$black;Iwata 2016 @;;
END


psxy $line -R$region -J$proj -W0.5p,$aqua,*- -O -K >> $out



psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O >> $out
#psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O >> $out


ps2pdf $out $outpdf
echo $out $outpdf

#===================================================


