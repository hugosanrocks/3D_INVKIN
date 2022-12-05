#!/usr/bin/perl

`rm *.ps`;
$plate_boundaries="plate_boundaries_peter_2002";
$scale="-98/17.25/12";
#Region for all stations
$region="-R-101/-97/14/19.5";
$region2="-R-117/-86/14/32";

$projection="-Jm$scale";
$projection2="-JM20";
$psfile="map.ps";
$ticks="-Ba0.5f0.5/a0.5f0.5SW";
$ticks2="-Ba60f60/a60f60";
$etopo2="temp2.grd";
$etopo2_int="temp2.grad";
$relief="relief_sea.cpt";
$lon_min="-99.5";
$lon_max="-96.5";
$lat_min="15";
$lat_max="17";
$topo="ETOPO1_Bed_c_gmt4.grd";
$edge_fault="edge_fault.out";

$portrait     = " -P";
$misc_begin   = "-K -V -X1 -Y1 $portrait";
$misc_middle  = "-K -O -V $portrait";
$misc_end     = "-O -V $portrait";
$res          = "h";
$earth        = "-G50/200/50 -S100/100/240";

$asano="~/MEGA/INV3DKIN/run/jap/latlon.dat";
$gpssta="~/svn/3D_INVKIN/3D_INVKIN/run/pue/plot/stations.sta";
$gpssta2="~/MEGA/INV3DKIN/run/jap/prepare/stations.sta2";
$gpssta2="~/MEGA/INV3DKIN/run/jap/prepare/stations.stause";
$coorsta="~/svn/3D_INVKIN/3D_INVKIN/run/pue/plot/stations.xy";
$coorsta2="~/MEGA/INV3DKIN/run/jap/prepare/stations.coor2";
$coorstaf="~/MEGA/INV3DKIN/run/jap/prepare/stations.cooruse";


$epi="pue.epi";
$after="sym.xyz";
$after2="~/MEGA/INV3DKIN/run/jap/prepare/aftershocks.coor";
$elat="16.254";
$elon="-98.531 16.254";

$portrait     = "-P ";
$misc_begin   = "-K -V -X1 -Y1 $portrait";
$misc_middle  = "-K -O -V $portrait";
$misc_end     = "-O -V $portrait";
$res          = "h";
$earth        = "-G50/200/50 -S100/100/240";

`gmtset PAPER_MEDIA a0`;
`gmtset BASEMAP_TYPE plain`;
`gmtset GRID_PEN_PRIMARY 0.25p,. `;
`gmtset ANNOT_FONT_SIZE_PRIMARY = 70p`;
`gmtset LABEL_FONT_SIZE = 70p`;
`gmtset PAGE_ORIENTATION landscape`;


$t_min="1.9";
$t_max="8";
$y_min="0";
$y_max="45";
$reg1="$t_min/$t_max/$y_min/$y_max";
$proj1="X9/-11.7";


`pscoast $region $projection -Na -Dh -X7.5 -A1000 -K -Ba2.5f2.5:%::."":WSne -P -V > $psfile`;


`grdcut $topo -Gtemp.grd $region -V`;
`grd2xyz -V temp.grd>topo.xyz `;
`surface $region -V -I0.5m topo.xyz -H2 -T0.25 -C0.1 -Gtemp2.grd`;
`grdsample temp.grd -Gtemp2.grd -I0.5m -V`;
`grdgradient temp2.grd -A300 -Gtemp2.grad -Ne0.6 -V`;
`grdimage -R -Jm temp2.grd -C$relief -Q -Itemp2.grad -O -K -V>>$psfile`;

#`pscoast $region $projection $ticks -W0.5p -P -K -O -V >> $psfile`;

#`psxy $projection $region -M -O -K -V -W2p,red $plate_boundaries >>$psfile`;


#Edge of the fault
#1
`psxy $projection $region -M -O -K -V -W10p,red edge_fault.out >> $psfile`;
#2
#`psxy $projection $region -M -O -K -V -W5p,blue edge_fault2.out >> $psfile`;

#all stations
`psxy $coorsta $region $projection -Si2.5 -Gwhite -W15,black -O -K -V >> $psfile`;
#`psxy $coorsta2 $region $projection -Sc0.4 -Gwhite -W15,black -O -K -V >> $psfile`;
`pstext $gpssta $region $projection -O -K -Dj0.0/1 -Gblack -V >> $psfile`;

`psxy $epi $region $projection -Sa1.5 -Gwhite -W1,black -O -K -V >> $psfile`;

#Plotting basemap and coast and country lines
`pscoast $region2 $projection2 -X26 $ticks2 -P -N1 -Ggreen -Wblue -K -O -V -Dh -Q >> $psfile`;

`pscoast $region2 $projection2 -P -W1 -Gblack -Sgray -N1 -O -K -V >> $psfile`;
`psxy $projection2 $region2 -M -O -K -V -W30p,yellow square.xy >> $psfile`;
`psxy $projection2 $region2 -M -O -V -W2p,red $plate_boundaries >>$psfile`;


#`rm ../fig11/kuma*ps`;
`ps2eps $psfile`;

