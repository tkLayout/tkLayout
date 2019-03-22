#!/bin/bash

# Following line set up full LCG view: nice but hide what is strictly necessary.
# source /cvmfs/sft.cern.ch/lcg/views/LCG_95/$ARCH-opt/setup.sh 

# Following 2 lines (1 is enough) would be nicer to set up ROOT env.
# source $RELEASE_LCG/ROOT/6.16.00/$ARCH-opt/ROOT-env.sh 
# export LD_LIBRARY_PATH=$VIEW_LCG/$ARCH-opt/lib64/:$VIEW_LCG/$ARCH-opt/lib/:$LD_LIBRARY_PATH

# Unfortunately, all solutions above break a lot of things on slc6 (fonts, emacs...), as described in 
# https://root-forum.cern.ch/t/root-production-release-v6-14-00-is-out/29360/23 .

# Only source strictly necessary ROOT dependencies, if ever you find a better solution for slc6, please PR it!
source $RELEASE_LCG/vdt/0.4.2/$ARCH-opt/vdt-env.sh
source $RELEASE_LCG/Davix/0.7.1/$ARCH-opt/Davix-env.sh
source $RELEASE_LCG/png/1.6.17/$ARCH-opt/png-env.sh
source $RELEASE_LCG/tbb/2019_U1/$ARCH-opt/tbb-env.sh
source $RELEASE_LCG/sqlite/3210000/$ARCH-opt/sqlite-env.sh
