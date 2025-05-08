
env=LCG_105a
env_version=x86_64-el9-gcc13-opt
lcg_dir=/cvmfs/sft.cern.ch/lcg/views/$env/$env_version
export BOOST_LIB=$lcg_dir/lib
source $lcg_dir/setup.sh
