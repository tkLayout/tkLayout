#!/bin/bash

ARCH=x86_64-centos7-gcc8
RELEASE=/cvmfs/sft.cern.ch/lcg/releases
RELEASE_LCG=$RELEASE/LCG_95
VIEW_LCG=/cvmfs/sft.cern.ch/lcg/views/LCG_95


# CMAKE
export PATH=$CONTRIB/CMake/3.13.4/Linux-x86_64/bin/:$PATH

# COMPILER
source $RELEASE/gcc/8.2.0/x86_64-centos7/setup.sh
export CC=`which gcc`
export CXX=`which g++`

# ROOT
source $RELEASE_LCG/ROOT/6.16.00/$ARCH-dbg/bin/thisroot.sh
export LD_LIBRARY_PATH=$VIEW_LCG/$ARCH-opt/lib64/:$VIEW_LCG/$ARCH-opt/lib/:$LD_LIBRARY_PATH

# BOOST
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$RELEASE_LCG/Boost/1.69.0/$ARCH-opt/include:$RELEASE_LCG/Boost/1.69.0/$ARCH-opt/lib

# DOXYGEN
export DOXYGEN_PATH=$RELEASE_LCG/doxygen/1.8.11/$ARCH-opt/bin/
export PATH=${DOXYGEN_PATH}:${PATH}

# GRAPHVIZ
source $RELEASE_LCG/graphviz/2.28.0/$ARCH-opt/graphviz-env.sh

# UPDATE PATH
export PATH=`pwd`/build/bin:$PATH
