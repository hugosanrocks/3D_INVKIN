
#modify paths if necessary inside "paths.sh"

source paths.sh
HOST=$(echo $HOSTNAME) 
echo $HOST

#======================================#
#prepare files according to options    #
#======================================#

#info for synthetic seismograms
 cp dat/syn.infobup dat/syn.info

#source sliprate model
 #zero inital model
 #cp initial_model_sis.src dat/vitesse.out

 #target SIV1 model
 cp dat/model_target.dat_bup dat/vitesse.out

 #with unstructured mesh
 #target SIV1 model
 #cp dat/model_target_triang.dat dat/vitesse.out

m#SIV1 model but shifted maximum slip to the left end of the fault
 #cp source_left.src dat/vitesse.out


#other options related to the length of forward modeling
 cp dat/fwioption.info_sis dat/fwioption.info

#focal mechanism is known
 cp dat/focal_1d.info dat/focal.info

#information about arrays
 cp dat/simul.infobup dat/simul.info

#Forward modeling
$BIN_DIR/FORWARD

#Plot syn vs obs single traces
#hugo-PC
#$PLOT_DIR/plot_single_traces.sh
#ist-oar Cluster
#$PLOT_DIR/plot_single_traces_istoar.sh


#Plot syn vs obs all traces together
octave $PLOT_DIR/makeplot_siv.m
#hugo-PC
#$PLOT_DIR/plot_seismograms_siv.sh
#ist-oar Cluster
#$PLOT_DIR/plot_seismograms_siv_istoar.sh


#INPUTS:
#inversion options
sis_pis=1      #1=SIS         , 2=PIS
option=1       #1=1D inversion, 2=2D inversion

#Variables to post process results and plot
node2see=422   #node to plot time history
finaltime=9.0  #final time for intergration (seconds)
dt=0.25        #time sampling of output sliprate file (seconds)

#Grid option
optgrid=1      #1=regular, 2=iregular


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


octave  plot_snaps.m
./plot_snaps.sh $CWP_DIR
cd ../

#Plot comparison of final slip
cd post/
#matlab -nodisplay -nosplash < sliprate_integration.m
#./compare_slip.sh $CWP_DIR


