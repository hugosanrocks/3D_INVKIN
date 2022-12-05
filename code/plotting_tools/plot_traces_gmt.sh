#!/bin/bash

#name of file to plot
basename=$1
station=$2
label=`basename $basename`

#3 components for each station
obs1=${basename}_C1.obs
syn1=${basename}_C1.syn
obs2=${basename}_C2.obs
syn2=${basename}_C2.syn
obs3=${basename}_C3.obs
syn3=${basename}_C3.syn

#legends
escal1=graphics/escala1.dat
escal2=graphics/escala2.dat

#file with current data window information used
win=${basename}_C1.win

# find start and end time and a reasonable step 
t_start=`grep T_START ${obs1} | awk '{print $NF}'`
t_end=`grep T_END ${obs1} | awk '{print $NF}'`
t_step=`echo $t_end $t_start | awk '{print int(($1-$2)/5)}'`
t_legend=75

# If you want to fix these values
#t_start=
#t_end=
#t_step=


# find number of windows
nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`

# default settings
paper="a4+"
orient="-P"
proj=X8/1.5
origin="-Y6"
shift="-Y-2.5"

gmtset PAPER_MEDIA $paper MEASURE_UNIT inch HEADER_FONT_SIZE 20p LABEL_FONT_SIZE 18p HEADER_OFFSET 0.1i
#line inside headers
nhead=4

#some colors
red=255/0/0
black=0/0/0
blue=0/0/255
green=00/200/00
paleblue=220/220/255

#################################
# plot seismograms with windows
#################################

# set output filename
out=${basename}.seis.eps
outpdf=${basename}.seis.pdf

# set region
max=`grep PLOT_MAX ${obs1} | awk '{print $NF}'`
region=$t_start/$t_end/-$max/$max

# do plot

#basemap
psbasemap -R$region -J$proj -B${t_step}::."Station #$station":"Time (s)":/S -K $orient $origin > $out

#colored window
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | psxy -R$region -J$proj -m -W1 -G$paleblue -O -K >> $out

#seismograms and legend
psxy $obs1 -R$region -J$proj -H$nhead -W2 -O -K >> $out
psxy $syn1 -R$region -J$proj -H$nhead -W0.5p,$red,- -O -K >> $out

#plot legend
psxy $escal2 -R$region -J$proj -W0.5p,$red,- -O -K >> $out
psxy $escal1 -R$region -J$proj -W2 -O -K >> $out

#some descriptive text
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $max 18 0 0 LT E-W
END

maxstring=$(printf "%03.3f" $max)

region=$t_start/$t_end/0/1
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_end 1 18 0 0 RT max=$maxstring (m/s)
$t_legend 0.82 16 0 0 LT Obs
$t_legend 0.68 16 0 0 LT Syn
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out

#################################
# do plot for next component
#################################
nhead=4

# set region
max=`grep PLOT_MAX ${obs2} | awk '{print $NF}'`
region=$t_start/$t_end/-$max/$max

#basemap
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

#colored window
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | psxy -R$region -J$proj -m -W1 -G$paleblue -O -K >> $out

#seismograms
psxy $obs2 -R$region -J$proj -H$nhead -W2 -O -K >> $out
psxy $syn2 -R$region -J$proj -H$nhead -W0.5p,$red,- -O -K >> $out
#psxy $syn2 -R$region -J$proj -H$nhead -W2/$red -O -K >> $out

maxstring=$(printf "%03.3f" $max)

#some descriptive text
region=$t_start/$t_end/0/1
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start 1 18 0 0 LT S-N
END
region=$t_start/$t_end/0/1
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_end 1 18 0 0 RT Max=$maxstring (m/s)
END
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out

#################################
# plot last component
#################################

nhead=4

# set region
max=`grep PLOT_MAX ${obs3} | awk '{print $NF}'`
region=$t_start/$t_end/-$max/$max

#basemap
psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

#colored window
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | psxy -R$region -J$proj -m -W1 -G$paleblue -O -K >> $out

#seismogrmas
psxy $obs3 -R$region -J$proj -H$nhead -W2 -O -K >> $out
psxy $syn3 -R$region -J$proj -H$nhead -W0.5p,$red,- -O -K >> $out

maxstring=$(printf "%03.3f" $max)

#some descriptive text
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_start $max 18 0 0 LT U-D
END
region=$t_start/$t_end/0/1
pstext -R$region -J$proj -N -O -K <<  END >> $out
$t_end 1 18 0 0 RT Max=$maxstring (m/s)
END

#End of plot
psbasemap -R$region -J$proj -O -B >> $out


#convert to pdf
#ps2pdf $out $outpdf

