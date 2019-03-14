#source /cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc8-opt/setup.sh 


# COMPILER
source /cvmfs/sft.cern.ch/lcg/contrib/gcc/8.2.0/x86_64-slc6-gcc8-opt/setup.sh
#source /cvmfs/sft.cern.ch/lcg/contrib/gcc/8binutils/x86_64-slc6-gcc8-opt/setup.sh
export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc8-opt/lib64/:/cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc8-opt/lib/:$LD_LIBRARY_PATH




# ROOT
#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/5.34.09/x86_64-slc6-gcc47-dbg/root/bin/thisroot.sh
#source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.03.02/x86_64-slc6-gcc48-dbg/root/bin/thisroot.sh




#source /cvmfs/sft.cern.ch/lcg/releases/LCG_95/ROOT/6.16.00/x86_64-slc6-gcc62-opt/bin/thisroot.sh
#source /cvmfs/sft.cern.ch/lcg/releases/LCG_95/ROOT/6.16.00/x86_64-slc6-gcc62-opt/bin/thisroot.sh
source /cvmfs/sft.cern.ch/lcg/releases/LCG_95/ROOT/6.16.00/x86_64-slc6-gcc8-opt/bin/thisroot.sh





# BOOST
#export BOOST_INCLUDE=/cvmfs/sft.cern.ch/lcg/external/Boost/1.50.0_python2.7/x86_64-slc6-gcc47-opt/include/boost-1_50
#export BOOST_INCLUDE=/cvmfs/sft.cern.ch/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/include/boost-1_53

#export BOOST_LIB=/cvmfs/sft.cern.ch/lcg/external/Boost/1.50.0_python2.7/x86_64-slc6-gcc47-opt/lib
#export BOOST_SUFFIX=-gcc47-mt-1_50
#export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/external/Boost/1.50.0_python2.7/x86_64-slc6-gcc47-opt/lib:$LD_LIBRARY_PATH
#export BOOST_LIB=/cvmfs/sft.cern.ch/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/lib
#export BOOST_SUFFIX=-gcc48-mt-1_53
#export LD_LIBRARY_PATH=/cvmfs/sft.cern.ch/lcg/external/Boost/1.53.0_python2.7/x86_64-slc6-gcc48-opt/lib:$LD_LIBRARY_PATH






#export BOOST_INCLUDE=/cvmfs/sft.cern.ch/lcg/releases/LCG_95/Boost/1.69.0/x86_64-slc6-gcc62-opt/include
#export BOOST_LIB=/cvmfs/sft.cern.ch/lcg/releases/LCG_95/Boost/1.69.0/x86_64-slc6-gcc62-opt/lib
#export BOOST_SUFFIX=
#export LD_LIBRARY_PATH=$BOOST_LIB:/cvmfs/sft.cern.ch/lcg/views/LCG_rootext20180517/x86_64-slc6-gcc62-opt/lib64:/cvmfs/sft.cern.ch/lcg/views/LCG_rootext20180517/x86_64-slc6-gcc62-opt/lib:$LD_LIBRARY_PATH


export BOOST_INCLUDE=/cvmfs/sft.cern.ch/lcg/releases/LCG_95/Boost/1.69.0/x86_64-slc6-gcc8-opt/include
export BOOST_LIB=/cvmfs/sft.cern.ch/lcg/releases/LCG_95/Boost/1.69.0/x86_64-slc6-gcc8-opt/lib
export BOOST_SUFFIX=
export LD_LIBRARY_PATH=$BOOST_LIB:$LD_LIBRARY_PATH







#---then Graphviz
#if [ -f $thisdir/bin/dot ]; then
#   GVBINDIR=$(dirname $(dirname $(readlink $thisdir/bin/dot)))/lib/graphviz; export GVBINDIR
#fi

#if [ -f $thisdir/etc/fonts/fonts.conf ]; then
#    export FONTCONFIG_PATH=$thisdir/etc/fonts
#fi

export FONTCONFIG_PATH=/cvmfs/sft.cern.ch/lcg/views/LCG_95/x86_64-slc6-gcc8-opt/etc/fonts/:$FONTCONFIG_PATH




# DOXYGEN
export DOXYGEN_PATH=/cvmfs/sft.cern.ch/lcg/releases/LCG_95/doxygen/1.8.11/x86_64-slc6-gcc8-opt/bin/
if [ -d $DOXYGEN_PATH ] ; then
  export PATH=${DOXYGEN_PATH}:${PATH}
fi

# UPDATE PATH
export PATH=`pwd`/bin:$PATH
