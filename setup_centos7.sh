#!/bin/bash

ARCH=x86_64-centos7-gcc8
RELEASE=/cvmfs/sft.cern.ch/lcg/releases
RELEASE_LCG=$RELEASE/LCG_95
VIEW_LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_95


# COMPILER
source $RELEASE/gcc/8.2.0/x86_64-centos7/setup.sh

# ROOT
source $RELEASE_LCG/ROOT/6.16.00/$ARCH-dbg/bin/thisroot.sh
export LD_LIBRARY_PATH=$VIEW_LCG/$ARCH-opt/lib64/:$VIEW_LCG/$ARCH-opt/lib/:$LD_LIBRARY_PATH

# BOOST
export BOOST_INCLUDE=$RELEASE_LCG/Boost/1.69.0/$ARCH-opt/include
export BOOST_LIB=$RELEASE_LCG/Boost/1.69.0/$ARCH-opt/lib
export BOOST_SUFFIX=
export LD_LIBRARY_PATH=$BOOST_LIB:$LD_LIBRARY_PATH

# DOXYGEN
export DOXYGEN_PATH=$RELEASE_LCG/doxygen/1.8.11/$ARCH-opt/bin/
export PATH=${DOXYGEN_PATH}:${PATH}

# GRAPHVIZ
export $RELEASE_LCG/graphviz/2.28.0/$ARCH-opt/graphviz-env.sh

# UPDATE PATH
export PATH=`pwd`/bin:$PATH

