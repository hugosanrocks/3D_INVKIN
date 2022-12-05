
#directory of plotting tools
PLOT_DIR=/data/cycle/sanchezh/3D_INVKIN/branches/predict_prior/plotting_tools

#stations
nsta=$1
#components
ncomp=3

echo $(printf "  Script to plot traces")
echo $(printf "  GMT is required!")

echo $(printf "  Total receivers to plot: " $nsta)


#clean graphics directory
rm graphics/syn_*ps

#prepare files and plot them
  for i in $(seq 1 $nsta) ; do
    for j in $(seq 1 $ncomp) ; do
    #prepare file with headers
    cat $(printf "out/syn_S%03i_C%01i.head  out/syn_S%03i_C%01i.ascii" $i $j $i $j) > $(printf "graphics/syn_S%03i_C%01i.syn" $i $j)
    cat $(printf "out/obs_S%03i_C%01i.head  out/obs_S%03i_C%01i.a" $i $j $i $j) > $(printf "graphics/syn_S%03i_C%01i.obs" $i $j)
    #move file with time window information
    cp  $(printf "out/syn_S%03i.win  graphics/syn_S%03i_C1.win" $i $i)
    done
octave $PLOT_DIR/escala_trace.m <<! 
$i
!
    $PLOT_DIR/plot_traces_gmt.sh $(printf "graphics/syn_S%03i %03i" $i $i)
    echo $(printf "  Plotting station: %02i" $i)
  done

cd graphics/
rm *win *syn *obs

