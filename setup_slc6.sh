# COMPILER
source /cvmfs/sft.cern.ch/lcg/external/gcc/4.7.2/x86_64-slc6-gcc47-opt/setup.sh

# ROOT
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc47-dbg/root/bin/thisroot.sh

# BOOST
export BOOST_INCLUDE=/cvmfs/sft.cern.ch/lcg/external/Boost/1.50.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_50
export BOOST_LIB=/cvmfs/sft.cern.ch/lcg/external/Boost/1.50.0_python2.7/x86_64-slc6-gcc47-opt/lib
export BOOST_SUFFIX=-gcc47-mt-1_50
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/external/Boost/1.50.0_python2.7/x86_64-slc6-gcc47-opt/lib:$LD_LIBRARY_PATH

# DOXYGEN
export DOXYGEN_PATH=/cvmfs/sft.cern.ch/lcg/external/doxygen/1.8.2/x86_64-slc6-gcc47-opt/bin
if [ -d $DOXYGEN_PATH ] ; then
  export PATH=${DOXYGEN_PATH}:${PATH}
fi

# UPDATE PATH
export PATH=`pwd`/bin:$PATH
