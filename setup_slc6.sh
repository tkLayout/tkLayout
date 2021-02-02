#!/bin/bash

ARCH=x86_64-slc6-gcc8
CONTRIB=/cvmfs/sft.cern.ch/lcg/contrib
RELEASE_LCG=/cvmfs/sft.cern.ch/lcg/releases/LCG_95
TK_DIRECTORY=$(dirname $BASH_SOURCE)


# COMPILER
source $CONTRIB/gcc/8.2.0/$ARCH-opt/setup.sh

# ROOT
source $RELEASE_LCG/ROOT/6.16.00/$ARCH-dbg/bin/thisroot.sh
source $TK_DIRECTORY/ROOT-env.sh

# BOOST
export BOOST_INCLUDE=$RELEASE_LCG/Boost/1.69.0/$ARCH-opt/include
export BOOST_LIB=$RELEASE_LCG/Boost/1.69.0/$ARCH-opt/lib
export BOOST_SUFFIX=
export LD_LIBRARY_PATH=$BOOST_LIB:$LD_LIBRARY_PATH

# DOXYGEN
export DOXYGEN_PATH=$RELEASE_LCG/doxygen/1.8.11/$ARCH-opt/bin/
export PATH=${DOXYGEN_PATH}:${PATH}

# UPDATE PATH
export PATH=`pwd`/bin:$PATH
