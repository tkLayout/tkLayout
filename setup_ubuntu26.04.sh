#!/bin/bash

if ! dpkg -s git cmake libboost-filesystem-dev libboost-program-options-dev libboost-json-dev >/dev/null 2>&1 ; then
  echo Some packages appear to be missing
  echo please run:
  echo sudo apt install git cmake libboost-filesystem-dev libboost-program-options-dev libboost-json-dev 
fi

ROOT_SCRIPT=`root-config --bindir`/thisroot.sh

if [ ! -f "$ROOT_SCRIPT" ] ; then
  echo root appears to be missing. Please install it manually
fi

source $ROOT_SCRIPT 

