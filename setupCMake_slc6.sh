#!/bin/bash

# EXPLICIT, MINIMAL DPEENDENCIES FOR SLC6
source setup_slc6.sh


# CMAKE
export PATH=$CONTRIB/CMake/3.13.4/Linux-x86_64/bin/:$PATH

# BOOST
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$RELEASE_LCG/Boost/1.69.0/$ARCH-opt/include:$RELEASE_LCG/Boost/1.69.0/$ARCH-opt/lib

# UPDATE PATH
export PATH=`pwd`/build/bin:$PATH
