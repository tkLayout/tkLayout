#!/bin/bash

source ~/.tkgeometryrc 
curdir=`dirname $0`
sourceDir=$curdir/../config
diffCommand=`which colordiff || which diff || echo NOTFOUND`
[ $diffCommand == "NOTFOUND" ] && {
  echo "I cannot find diff not colordiff. I quit."
  exit 0
}

[ -d $sourceDir ] || {
  echo "ERROR: Source file config directory ($sourceDir) not found"
  exit 0
}

materialFile=Materials.cpp
referenceFile=$TKG_STANDARDDIRECTORY/config/Materials.cfg
sourceFile=$sourceDir/Materials.cfg
tempMaterialFile=tempMaterials.cfg

gcc -I modules -I parts -C -E $materialFile > $tempMaterialFile && {
  realText=`cat $tempMaterialFile | egrep -v '^#'`
  echo "$realText" > $tempMaterialFile
  echo $diffCommand $referenceFile $tempMaterialFile
  $diffCommand $referenceFile $tempMaterialFile
  echo -n "Do you want to export these changes to the active configuration directory? (y/n) [y]: "
  read answer
  if [ "$answer" == "y" ] || [ "$answer" == "" ] ; then 
    echo New configuration file will be copied to $referenceFile
    cp $tempMaterialFile $referenceFile 
  fi
  echo
  echo $diffCommand $sourceFile $tempMaterialFile
  $diffCommand $sourceFile $tempMaterialFile
  echo    "Do you want to export these changes to the source code directory?"
  echo -n "This will override the active configuration file at the next install (y/n) [n]: "
  read answer
  if [ "$answer" == "y" ] ; then 
    echo New configuration file will be copied to $sourceFile
    cp $tempMaterialFile $sourceFile
  fi 
}
rm -f $tempMaterialFile
