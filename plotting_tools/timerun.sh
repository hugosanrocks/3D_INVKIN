nsta=30
#15
ncomp=3
dt=0.04
fac=0.0
wdCount=$(wc -w <<<$winget)
rm graphics/syn_*ps
rm input.graphics
#first line of input file
echo $(printf "%02i" $nsta) >> input.graphics

  for i in $(seq 1 $nsta) ; do
    for j in $(seq 1 $ncomp) ; do
    cat $(printf "out/syn_S%03i_C%01i.head  out/syn_S%03i_C%01i.ascii" $i $j $i $j) > $(printf "graphics/syn_S%03i_C%01i.syn" $i $j)
    cat $(printf "out/obs_S%03i_C%01i.head  out/obs_S%03i_C%01i.a" $i $j $i $j) > $(printf "graphics/syn_S%03i_C%01i.obs" $i $j)
    cp  $(printf "out/syn_S%03i.win  graphics/syn_S%03i_C1.win" $i $i)
    cp  $(printf "out/syn_S%03i.win.qual  graphics/syn_S%03i_C1.win.qual" $i $i)
    done
    ~/svn/3D_INVKIN/3D_INVKIN/branches/slip_interp/plotting_tools/plot_seismos_gmt.sh $(printf "graphics/syn_S%03i" $i)
    echo $(printf "Plotting station: %02i" $i)
  done

cd graphics/
rm *win *win.* *syn *obs *head

# NPTS =          875
# PLOT_MAX =    0.075
# T_START =      0.00
# T_END =       34.96

