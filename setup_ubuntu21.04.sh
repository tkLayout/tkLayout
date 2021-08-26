#!/bin/bash

if ! dpkg -s git cmake libboost1.74-all-dev doxygen graphviz  >/dev/null 2>&1 ; then
  echo Some packages appear to be missing
  echo please run sudo apt install git cmake libboost1.74-all-dev doxygen graphviz
fi

ROOT_SCRIPT=`root-config --bindir`/thisroot.sh

if [ ! -f "$ROOT_SCRIPT" ] ; then
  echo root appears to be missing. Please install it manually
fi

source $ROOT_SCRIPT 

