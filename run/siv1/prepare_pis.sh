
BIN_DIR=~/svn/3D_INVKIN/3D_INVKIN/branches/predict_prior/bin
PLOT_DIR=~/svn/3D_INVKIN/3D_INVKIN/branches/predict_prior/plotting_tools

#clean graphics directory
rm graphics/*

#prepare simul.info files for every limited source time-space zone
octave pis_simul.m

#set correct initial files to estimate correctly
#copy model used to estimate time windows
cp dat/model_target.dat dat/vitesse.out
#copy correct synthetic information
cp dat/syn.infobup dat/syn.info

#number of windows
nstep=9
#time samples at each window
samples=(3 5 7 9 13 19 27 35 39)
#time sampling
timestep=0.25

  k=1
  for i in $(seq 1 $nstep) ; do
    samp=${samples[i-1]}
    time=$(expr $timestep*$samp | bc)
    time=$(expr $time-$timestep | bc)
    j=$(expr $i*$k | bc)
    cp $(printf "simul_%02i.info dat/simul.info" $j)
    $BIN_DIR/FORWARD
    echo $time $j
    octave windows_station.m <<!
$time
!
  done

#save window limits in the correct directory and format
octave write_windows.m

#Plot seismograms and time limits for each receiver location
cd graphics/
$PLOT_DIR/plot_windows.sh
cd ../

#Set the right synwin.info to perform a PIS exercise
cp dat/syn.info_pis dat/syn.info

#clean directory of tmp files
rm simul_*.info window_*.** out/*
