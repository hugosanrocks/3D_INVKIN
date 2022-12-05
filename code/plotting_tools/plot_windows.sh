
#directory of plotting tools
PLOT_DIR=/data/cycle/sanchezh/3D_INVKIN/branches/predict_prior/plotting_tools
cwp_dir=/soft/cwpsu

sta=$1
ncomp=3
nsta=2

for j in $(seq 1 $nsta) ; do
  sta=$j
  for i in $(seq 1 $ncomp) ; do
   $PLOT_DIR/windows.sh $(printf "S%03i_C%01i" $sta $i) $sta $i
  name1=$(printf "S%03i_C%01i.eps" $sta 1)
  name2=$(printf "S%03i_C%01i.eps" $sta 2)
  name3=$(printf "S%03i_C%01i.eps" $sta 3)
  out=$(printf "windows_S%03i.eps" $sta)
  outpdf=$(printf "S%03i.pdf" $sta)
  done
  $cwp_dir/bin/psmerge \
   in=$name1 translate=0,5 \
   in=$name2 translate=0,2.5 \
   in=$name3 translate=0,0 > $out
   ps2pdf -dEPSCrop $out $outpdf
done
rm *C1.eps *C2.eps *C3.eps *pdf

