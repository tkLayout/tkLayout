#!/bin/bash

ARCH=x86_64-centos8-gcc10
RELEASE=/cvmfs/sft.cern.ch/lcg/releases
RELEASE_LCG=$RELEASE/LCG_99
VIEW_LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_99


# COMPILER
source $RELEASE/gcc/10.1.0/x86_64-centos8/setup.sh

# ROOT
source $RELEASE_LCG/ROOT/v6.22.06/$ARCH-dbg/bin/thisroot.sh
export LD_LIBRARY_PATH=$VIEW_LCG/$ARCH-opt/lib64/:$VIEW_LCG/$ARCH-opt/lib/:$LD_LIBRARY_PATH

# BOOST
export BOOST_INCLUDE=$RELEASE_LCG/Boost/1.73.0/$ARCH-opt/include
export BOOST_LIB=$RELEASE_LCG/Boost/1.73.0/$ARCH-opt/lib
export BOOST_SUFFIX=
export LD_LIBRARY_PATH=$BOOST_LIB:$LD_LIBRARY_PATH

# DOXYGEN
export DOXYGEN_PATH=$RELEASE_LCG/doxygen/1.8.18/$ARCH-opt/bin/
export PATH=${DOXYGEN_PATH}:${PATH}

# GRAPHVIZ
source $RELEASE_LCG/graphviz/2.40.1/$ARCH-opt/graphviz-env.sh

# UPDATE PATH
export PATH=`pwd`/bin:$PATH

