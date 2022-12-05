
#directory of plotting tools
#PLOT_DIR=~/svn/3D_INVKIN/3D_INVKIN/branches/predict_prior/plotting_tools

nsta=40
ncomp=3

#octave makeplot_siv.m
rm graphics/all_s*eps

  for i in $(seq 1 $nsta) ; do
    #OBSERVATIONS
    cat $(printf "seismo_to_plote_%02i.obs" $i) > $(printf "graphics/seismo_to_plote_%02i.obs" $i)
    cat $(printf "seismo_to_plotn_%02i.obs" $i) > $(printf "graphics/seismo_to_plotn_%02i.obs" $i)
    cat $(printf "seismo_to_plotv_%02i.obs" $i) > $(printf "graphics/seismo_to_plotv_%02i.obs" $i)
    #SYNTHETICS
    cat $(printf "seismo_to_plote_%02i.syn" $i) > $(printf "graphics/seismo_to_plote_%02i.syn" $i)
    cat $(printf "seismo_to_plotn_%02i.syn" $i) > $(printf "graphics/seismo_to_plotn_%02i.syn" $i)
    cat $(printf "seismo_to_plotv_%02i.syn" $i) > $(printf "graphics/seismo_to_plotv_%02i.syn" $i)
  done

  rm  $(printf "seismo_to_plot*")
  $PLOT_DIR/plot_seismos_all_siv_istoar.sh
  rm graphics/seismo_to_plot*
  ps2pdf graphics/all_seismograms.eps
  pdfcrop all_seismograms.pdf
  mv all_seismograms-crop.pdf graphics/all_seismograms.pdf

