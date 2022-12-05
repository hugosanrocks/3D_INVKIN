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
escal=graphics/escala1.dat
escal2=graphics/escala2.dat

#file with current data window information used
win=${basename}_C1.win

# find start and end time and a reasonable step 
t_start=`grep T_START ${obs1} | awk '{print $NF}'`
t_end=`grep T_END ${obs1} | awk '{print $NF}'`
t_step=`echo $t_end $t_start | awk '{print int(($1-$2)/15)}'`
t_legend=30

# If you want to fix these values
#t_start=
#t_end=
#t_step=


# find number of windows
nwin=`grep NUM_WIN ${win} | awk '{print $NF}'`

# default settings
paper="a3"
orient="-P"
proj=X8/1.5
origin="-Y6"
t_max=24
shift="-Y-2.5"

gmt gmtset PS_MEDIA $paper PROJ_LENGTH_UNIT inch FONT_TITLE 20p FONT_LABEL 18p MAP_TITLE_OFFSET 0.1i
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
gmt psbasemap -R$region -J$proj -B${t_step}::."Station #$station":"Time (s)":/S -K $orient $origin > $out

#colored window
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | gmt psxy -R$region -J$proj -W1 -G$paleblue -O -K >> $out

#seismograms and legend
gmt psxy $obs1 -R$region -J$proj -h$nhead -W2 -O -K >> $out
gmt psxy $syn1 -R$region -J$proj -h$nhead -W0.5p,$blue,- -O -K >> $out

#plot legend
gmt psxy $escal -R$region -J$proj -W0.5p,$blue,- -O -K  >> $out
gmt psxy $escal2 -R$region -J$proj -W2 -O -K  >> $out

#some descriptive text
gmt pstext -R$region -J$proj -F+f18p,Helvetica -N -O -K <<  END >> $out
$t_start $max E-W
END

maxstring=`echo $max`

region=$t_start/$t_end/0/1
gmt pstext -R$region -J$proj -F+f16p,Helvetica -N -O -K <<  END >> $out
$t_max 1 max=$maxstring (m/s)
$t_legend 0.82 Syn
$t_legend 0.68 Obs
END
gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out

#################################
# do plot for next component
#################################
nhead=4

# set region
max=`grep PLOT_MAX ${obs2} | awk '{print $NF}'`
region=$t_start/$t_end/-$max/$max

#basemap
gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

#colored window
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | gmt psxy -R$region -J$proj -W1 -G$paleblue -O -K >> $out

#seismograms
gmt psxy $obs2 -R$region -J$proj -h$nhead -W2 -O -K >> $out
gmt psxy $syn2 -R$region -J$proj -h$nhead -W0.5p,$blue,- -O -K >> $out
#gmt psxy $syn2 -R$region -J$proj -h$nhead -W2/$blue -O -K >> $out

maxstring=`echo $max`

#some descriptive text
region=$t_start/$t_end/0/1
gmt pstext -R$region -J$proj -F+f18p,Helvetica -N -O -K <<  END >> $out
$t_start 1 S-N
END
region=$t_start/$t_end/0/1
gmt pstext -R$region -J$proj -F+f16p,Helvetica -N -O -K <<  END >> $out
$t_max 1 max=$maxstring (m/s)
END
gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K >> $out

#################################
# plot last component
#################################

nhead=4

# set region
max=`grep PLOT_MAX ${obs3} | awk '{print $NF}'`
region=$t_start/$t_end/-$max/$max

#basemap
gmt psbasemap -R$region -J$proj -B${t_step}:"Time (s)":/S -O -K $shift >> $out

#colored window
tail -$nwin $win | awk '{printf "%f %f\n%f %f\n%f %f\n%f %f\n>\n", $2,-1*m,$2,m,$3,m,$3,-1*m}' m=$max | gmt psxy -R$region -J$proj -W1 -G$paleblue -O -K >> $out

#seismogrmas
gmt psxy $obs3 -R$region -J$proj -h$nhead -W2 -O -K >> $out
gmt psxy $syn3 -R$region -J$proj -h$nhead -W0.5p,$blue,- -O -K >> $out

maxstring=`echo $max`

#some descriptive text
gmt pstext -R$region -J$proj -F+f18p,Helvetica -N -O -K <<  END >> $out
$t_start $max U-D
END
region=$t_start/$t_end/0/1
gmt pstext -R$region -J$proj -F+f16p,Helvetica -N -O -K <<  END >> $out
$t_max 1 max=$maxstring (m/s)
END

#End of plot
gmt psbasemap -R$region -J$proj -O -B >> $out


#convert to pdf
#ps2pdf $out $outpdf

