#rm fort.*

#modify paths
BIN_DIR=~/svn/3D_INVKIN/3D_INVKIN/branches/predict_prior/bin
PLOT_DIR=~/svn/3D_INVKIN/3D_INVKIN/branches/predict_prior/plotting_tools
CWP_DIR=~/Desktop/cwp/bin

#remove output files and graphics
#rm out/*
#rm post/*bin post/*out post/*ps post/*pdf graphics/*

#INPUTS:
#inversion options
sis_pis=1      #1=SIS         , 2=PIS
option=2       #1=1D inversion, 2=2D inversion

#Variables to post process results and plot
node2see=422   #node to plot time history
finaltime=9.0  #final time for intergration (seconds)
dt=0.25        #time sampling of output sliprate file (seconds)

#Grid option
optgrid=1      #1=regular, 2=iregular

#prepare files according to inversion options
if [ $sis_pis -eq 1 ] ; then
 cp dat/syn.infobup dat/syn.info
 cp initial_model_sis.src dat/vitesse.out 
# cp dat/model_target_triang.dat dat/vitesse.out
 cp dat/fwioption.info_sis dat/fwioption.info
  if [ $option -eq 1 ] ; then
   cp dat/focal_1d.info dat/focal.info
  else
   cp dat/focal_2d.info dat/focal.info
  fi
else
 cp dat/syn.info_pis dat/syn.info
 cp dat/fwioption.info_pis dat/fwioption.info
  if [ $option -eq 1 ] ; then
   cp dat/focal_1d.info dat/focal.info
   cp initial_model_pis.src dat/modelpri1d.dat
  else
   cp dat/focal_2d.info dat/focal.info
   cp initial_model_pis.src dat/modelpri.dat
  fi
fi

#Run inversion
$BIN_DIR/INV3DKIN

#Plot syn vs obs single traces
$PLOT_DIR/plot_single_traces.sh

#Plot syn vs obs all traces together
octave $PLOT_DIR/makeplot_siv.m
$PLOT_DIR/plot_seismograms_siv.sh

#Plot slip-rate snapshots
cd post/
octave out2sliprate.m <<!
$option
$node2see
!

file="opt_grid.dat"
cat <<EOM >$file
$optgrid
$finaltime
$dt
EOM


matlab -nodisplay -nosplash < plot_snaps.m
./plot_snaps.sh $CWP_DIR
cd ../

#Plot comparison of final slip
cd post/
matlab -nodisplay -nosplash < sliprate_integration.m
./compare_slip.sh $CWP_DIR


