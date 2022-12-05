nsta=30
#15
ncomp=3

#mv seismo_to*obs seismo_to*syn graphics/
#octave makeplot.m
rm graphics/seismo_to*eps

  for i in $(seq 1 $nsta) ; do
#    cat $(printf "seismo_all.head seismo_to_plote_%02i.obs" $i) > $(printf "graphics/seismo_to_plote_%02i.obs" $i)
#    cat $(printf "seismo_all.head seismo_to_plotn_%02i.obs" $i) > $(printf "graphics/seismo_to_plotn_%02i.obs" $i)
#    cat $(printf "seismo_all.head seismo_to_plotv_%02i.obs" $i) > $(printf "graphics/seismo_to_plotv_%02i.obs" $i)
    mv $(printf "seismo_to_plote_%02i.obs graphics/seismo_to_plote_%02i.obs" $i $i)
    mv $(printf "seismo_to_plotn_%02i.obs graphics/seismo_to_plotn_%02i.obs" $i $i)
    mv $(printf "seismo_to_plotv_%02i.obs graphics/seismo_to_plotv_%02i.obs" $i $i)

    mv $(printf "seismo_to_plote_%02i.syn graphics/seismo_to_plote_%02i.syn" $i $i)
    mv $(printf "seismo_to_plotn_%02i.syn graphics/seismo_to_plotn_%02i.syn" $i $i)
    mv $(printf "seismo_to_plotv_%02i.syn graphics/seismo_to_plotv_%02i.syn" $i $i)

    mv $(printf "seismo_to_plote2_%02i.syn graphics/seismo_to_plote2_%02i.syn" $i $i)
    mv $(printf "seismo_to_plotn2_%02i.syn graphics/seismo_to_plotn2_%02i.syn" $i $i)
    mv $(printf "seismo_to_plotv2_%02i.syn graphics/seismo_to_plotv2_%02i.syn" $i $i)


#    cat $(printf "seismo_all.head seismo_to_plote_%02i.syn" $i) > $(printf "graphics/seismo_to_plote_%02i.syn" $i)
#    cat $(printf "seismo_all.head seismo_to_plotn_%02i.syn" $i) > $(printf "graphics/seismo_to_plotn_%02i.syn" $i)
#    cat $(printf "seismo_all.head seismo_to_plotv_%02i.syn" $i) > $(printf "graphics/seismo_to_plotv_%02i.syn" $i)
  done

#    rm  $(printf "seismo_to_plot*")
#    rm graphics/seismo_to_plot*

#Only two stations
#    ~/svn/3D_INVKIN/3D_INVKIN/branches/slip_interp/plotting_tools/plot_1sta.sh

# All stations
    ~/svn/3D_INVKIN/3D_INVKIN/branches/predict_prior/plotting_tools/plot_seismos_all.sh

    rm graphics/seismo_to_plot*

