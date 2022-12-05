#!/bin/bash

# NAME OF THE FIGURE
# set output filename
out=snaps_siv1


# default settings
paper="a4+"
orient="-P"
proj=X5/4.15
proj2=x3/2.5
origin="-Y1i"
shifty1="-Y1i"
shifty2="-Y1.5i"
backy="-Y-3.5i"
shiftx="-X2.5i"
region="-R0.5/36.5/0.5/18.5"
projection="-JX1.8i/0.9i"
ticks="-Ba5:"":/a4:"":SW"
ticksx="-Ba5:"":/a50:"":SW"
ticksy="-Ba50:"":/a4:"":SW"
#label="-Ba50:"Along strike (km)":/a4:"Along dip (km)":SW"


gmtset PAPER_MEDIA $paper MEASURE_UNIT inch HEADER_FONT_SIZE 20p LABEL_FONT_SIZE 15p ANNOT_FONT_SIZE_PRIMARY = 12p
gmtset ANNOT_OFFSET_PRIMARY    = 0.07i
red=255/0/0
black=0/0/0
blue=0/0/255
green=00/200/00
paleblue=220/220/255

#################################
# plot seismograms with windows
#################################

outpdf=${out}.pdf
ps=${out}.eps

rm $ps $outpdf

#T=`minmax -T9/1 table_9.ascii`
makecpt -Cjet -T-4e-8/1.6/0.1 > topo.cpt
#pscontour -R -J table_2.ascii -B2f1WSne -Wthin -Ctopo.cpt -Lthinnest,- -G1i/0 -X-3.25i -Y-3.65i -O -K \
#	-U"Example 12 in Cookbook" >> $ps

#Snapshots of slip-rate evolution
#1st column
pscontour $region $projection $origin table_1.ascii $ticks -Ctopo.cpt -I -K >> $ps
pscontour $region $projection tableinv_1.ascii $ticksy -Ctopo.cpt -I $shifty1 -O -K >> $ps

gmtset ANNOT_FONT_SIZE_PRIMARY = 12p

pscontour $region $projection table_2.ascii $ticks -Ctopo.cpt -I $shifty2 -O -K >> $ps
pscontour $region $projection tableinv_2.ascii $ticksy -Ctopo.cpt -I $shifty1 -O -K >> $ps


#2 column
pscontour $region $projection table_3.ascii $ticks -Ctopo.cpt -I $backy $shiftx -O -K >> $ps
pscontour $region $projection tableinv_3.ascii $ticksy -Ctopo.cpt -I $shifty1 -O -K >> $ps

pscontour $region $projection table_4.ascii $ticks -Ctopo.cpt -I $shifty2 -O -K >> $ps
pscontour $region $projection tableinv_4.ascii $ticksy -Ctopo.cpt -I $shifty1 -O -K >> $ps


#3 column
pscontour $region $projection table_5.ascii $ticks -Ctopo.cpt -I $backy $shiftx -O -K >> $ps
pscontour $region $projection tableinv_5.ascii $ticksy -Ctopo.cpt -I $shifty1 -O -K >> $ps

pscontour $region $projection table_6.ascii $ticks -Ctopo.cpt -I $shifty2 -O -K >> $ps
pscontour $region $projection tableinv_6.ascii $ticksy -Ctopo.cpt -I $shifty1 -O -K >> $ps

gmtset LABEL_OFFSET = -0.1i
gmtset HEADER_OFFSET = -0.1i
##gmtset ANNOT_OFFSET_SECONDARY  = 0.07i

#4 column
pscontour $region $projection table_7.ascii $ticks -Ctopo.cpt -I $backy $shiftx -O -K >> $ps
pscontour $region $projection tableinv_7.ascii -Ba50:"":/a40:"Along dip [km]":SE -Ctopo.cpt -I -O -K >> $ps
gmtset LABEL_OFFSET = 0.16i
gmtset HEADER_OFFSET = 0.16i
pscontour $region $projection tableinv_7.ascii -Ba50:"Along strike [km]":/a40:"":SE -Ctopo.cpt -I -O -K >> $ps


pscontour $region $projection tableinv_7.ascii $ticksy -Ctopo.cpt -I $shifty1 -O -K >> $ps

pscontour $region $projection table_8.ascii $ticks -Ctopo.cpt -I $shifty2 -O -K >> $ps
pscontour $region $projection tableinv_8.ascii $ticksy -Ctopo.cpt -I $shifty1 -O -K >> $ps

reg2="-R0/1/0/1"
proj2="-JX10i/4i"
#legend for rupture time
#legend for solutions
pstext $reg2 $proj2 -X-8i -Y-3.4i -N -O -K <<  END >> $ps
0. 0.53 18 0 0 LT t=0 [s]
0. 1.15 18 0 0 LT t=1.28 [s]
0.25 0.53 18 0 0 LT t=2.56 [s]
0.25 1.15 18 0 0 LT t=3.84 [s]
0.5 0.53 18 0 0 LT t=5.12 [s]
0.5 1.15 18 0 0 LT t=6.40 [s]
0.75 0.53 18 0 0 LT t=7.68 [s]
0.75 1.15 18 0 0 LT t=8.96 [s]

0.00 0.475 14 0 0 LT SIV1
0.00 0.23 14 0 0 LT PIS2
0.25 0.475 14 0 0 LT SIV1
0.25 0.23 14 0 0 LT PIS2
0.50 0.475 14 0 0 LT SIV1
0.50 0.23 14 0 0 LT PIS2
0.75 0.475 14 0 0 LT SIV1
0.75 0.23 14 0 0 LT PIS2
0.00 0.85 14 0 0 LT PSI2
0.00 1.095 14 0 0 LT SIV1
0.25 0.85 14 0 0 LT PSI2
0.25 1.095 14 0 0 LT SIV1
0.50 0.85 14 0 0 LT PSI2
0.50 1.095 14 0 0 LT SIV1
0.75 0.85 14 0 0 LT PSI2
0.75 1.095 14 0 0 LT SIV1
END

gmtset ANNOT_FONT_SIZE_PRIMARY = 12p 
gmtset LABEL_OFFSET = 0.05i

psscale -Ctopo.cpt -D4.6/-0.45/2.5/0.2h -B0.2::."":"Slip rate [m/s]":/S -O -V >> $ps


#echo "0.5 8.5 4.0 0.0 0.0 0.0 0.0" | awk '{printf "%f %f %f %f %f %f %f>\n", $1,$2,$3,$4,$5,$6,$7}' | psvelo -Se0.2i/0.1i/18 -R0.5/73/0.5/37 -JX7.2i/3.6i -G255/0/0 -X-3.5i -Y-2i -A0.05/0.1/0.15 -O -K -V >> $out

#echo "0.5 5 4.0 0.0 0.0 0.0 0.0" | awk '{printf "%f %f %f %f %f %f %f>\n", $1,$2,$3,$4,$5,$6,$7}' | psvelo -Se0.2i/0.1i/18 -R0.5/73/0.5/37 -JX7.2i/3.6i -G255/0/0 -A0.05/0.1/0.15 -O -V >> $out
#echo "3.16 8 30 0 1 BC Delaunay Triangulation" | \
#	pstext -R0/8/0/11 -Jx1i -O -X-3.25i >> $ps



ps2pdf $ps $outpdf


