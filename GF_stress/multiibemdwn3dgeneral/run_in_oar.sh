#!/bin/bash
#OAR -n greenfun
#OAR -l /nodes=1/core=12,walltime=40:30:00
#OAR -O result.out
#OAR -E error.out
#OAR --project iste-equ-cycle
#OAR --notify mail:hugo.geofisica@gmail.com

####OAR -p network_address='ist-calcul1.ujf-grenoble.fr'

source /soft/env.bash
module load MATLAB/R2016a

matlab -nodisplay -nosplash  < run_stress_invedkin.m

