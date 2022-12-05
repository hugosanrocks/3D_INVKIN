#!/bin/bash

#usage: windows.sh  NAME  Station  Component
#./windows.sh S002_C1 2 1

#directory of plotting tools
PLOT_DIR=~/svn/3D_INVKIN/3D_INVKIN/branches/predict_prior/plotting_tools
 
basename=$1
label=`basename $basename`
obs1=../out/obs_${basename}.a
syn1=../out/syn_${basename}.ascii
obshead=../out/obs_${basename}.head
windows=../dat/windows.info
station=$2
comp=$3
echo 'plotting station' $station 'component' $comp
escobs=escala.obs
escsyn=escala.syn

max=`grep MAX ${obshead} | awk '{print $NF}'`
maxamp=$(printf "%02.3f" $max)
max=$maxamp

t_start=0
t_end=18
t_step=`echo $t_end $t_start | awk '{print int(($1-$2)/18)}'`
echo $max

# default settings
paper="a4+"
orient="-P"
proj=X5/1.5
origin="-Y8"
shift="-Y-2.0"

gmtset PAPER_MEDIA $paper MEASURE_UNIT inch HEADER_FONT_SIZE 20p LABEL_FONT_SIZE 15p ANNOT_FONT_SIZE_PRIMARY 15p
nhead=0
red=255/0/0
black=0/0/0
blue=0/0/255
green=00/200/00
paleblue=220/220/255

nline=22
ncols=3
exec 5< $PLOT_DIR/jetneg.ascii
for i in `seq 1 ${nline}`
do
read -a arr -u 5
  for j in `seq 0 ${ncols}`
  do
  color[$i,$j]=$(echo ${arr[$j]})
  c1=${color[$j,0]}
  c2=${color[$j,1]}
  c3=${color[$j,2]}
  #c1=255
  #c2=$c1
  #c3=$c1
  c1=255 #23
  c2=255 #255
  c3=255 #255
  colors[$i]=$(printf "%01i/%01i/%01i" $c1 $c2 $c3)
  colors2[$i]=255/255/255
  done
done

#################################
# plot seismograms with windows
#################################

# set output filename
out=${basename}.eps
outpdf=${basename}.pdf


region=$t_start/$t_end/-$max/$max

psbasemap -R$region -J$proj -B${t_step}::."":"Time (sec)":/S -K $orient $origin > $out
nwin=11
#11
nsta=1
#2-21 windows

#for i in $(seq 1 $nsta) ; do
i=$station
skip=2
for j in $(seq 3 $nwin) ; do
 k=$j
 k=$(expr "$k" \* "$skip")
 time2=$(awk -v i=$i -v j=$j 'NR==i{print $j}' $windows)
 samp=$(expr $j-1| bc)
 time1=$(awk -v i=$i -v j=$samp 'NR==i{print $j}' $windows)
 if [ $j -le 11 ] ; then
a=2
echo $time1 $time2
 echo $time1 $time2 | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $1,-1*m,$1,m,$2,m,$2,-1*m}' m=$max | psxy -R$region -J$proj -m -W1,black -G${colors[$k]} -O -K >> $out 
 else 
a=1
 echo $time1 $time2 | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $1,-1*m,$1,m,$2,m,$2,-1*m}' m=$max | psxy -R$region -J$proj -m -W1,black -G${colors2[$k]} -O -K >> $out 
fi
 finish=$(awk -v i=$i -v j=1 'NR==i{print $j}' $windows)
done
#done



psxy $obs1 -R$region -J$proj -W1.0 -O -K >> $out
psxy $syn1 -R$region -J$proj -W1.5p,$red,- -O -K >> $out

#psxy $escobs -R$region -J$proj -W3 -O -K >> $out
#psxy $escsyn -R$region -J$proj -W1.0p,$red,- -O -K >> $out


if [ $comp -eq 1 ] ; then
pstext -R$region -J$proj -Gblack -Ni -O -K <<  END >> $out
$t_start $max 18 0 0 LT E-W
END
elif [ $comp -eq 2 ] ; then
pstext -R$region -J$proj -Gblack -Ni -O -K <<  END >> $out
$t_start $max 18 0 0 LT N-S
END
elif [ $comp -eq 3 ] ; then
pstext -R$region -J$proj -Gblack -Ni -O -K <<  END >> $out
$t_start $max 18 0 0 LT Z
END
fi

region=$t_start/$t_end/0/1
#pstext -R$region -J$proj -Gblack -N -O -K <<  END >> $out
#32 0.86 12 0 0 LT Syn
#32 0.74 12 0 0 LT Obs
#28 1.00 14 0 0 LT Station i$station
#34 0.20 14 0 0 RT max=$max (m/s)
#END

pstext -R$region -J$proj -Gblack -N -O -K <<  END >> $out
28 1.00 14 0 0 LT Station i$station
34 0.20 14 0 0 RT max=$max (m/s)
END

psbasemap -R$region -J$proj -B${t_step}:"":/S -O >> $out

ps2pdf $out $outpdf


