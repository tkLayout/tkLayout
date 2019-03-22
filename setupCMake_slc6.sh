#!/bin/bash

ARCH=x86_64-slc6-gcc8
CONTRIB=/cvmfs/sft.cern.ch/lcg/contrib
RELEASE_LCG=/cvmfs/sft.cern.ch/lcg/releases/LCG_95
TK_DIRECTORY=$(dirname $BASH_SOURCE)


# CMAKE
export PATH=$CONTRIB/CMake/3.13.4/Linux-x86_64/bin/:$PATH

# COMPILER
source $CONTRIB/gcc/8.2.0/$ARCH-opt/setup.sh
export CC=`which gcc`
export CXX=`which g++`

# ROOT
source $RELEASE_LCG/ROOT/6.16.00/$ARCH-dbg/bin/thisroot.sh
source $TK_DIRECTORY/ROOT-env.sh

# BOOST
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$RELEASE_LCG/Boost/1.69.0/$ARCH-opt/include

# DOXYGEN
export DOXYGEN_PATH=$RELEASE_LCG/doxygen/1.8.11/$ARCH-opt/bin/
export PATH=${DOXYGEN_PATH}:${PATH}

