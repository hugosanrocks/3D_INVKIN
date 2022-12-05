nstep=19
samples=(3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 33 35 37 39)
timestep=0.25

  for i in $(seq 1 $nstep) ; do
    samp=${samples[i-1]}
    time=$(expr $timestep*$samp | bc)
    time=$(expr $time-$timestep | bc)
    cp $(printf "s%02i dat/simul.info" $i)
    ./for2.sh
    echo $time
    cat $time file.in
    octave windows_sta.m <<!
$time
!
  done
