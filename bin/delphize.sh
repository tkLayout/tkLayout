#!/bin/bash
source ~/tklayout.env 
source $TKGEOMETRYRC
source $TKG_BINDIRECTORY/setup.sh

cd $TKG_LAYOUTDIRECTORY
fileList=`ls */linptres_*_MS*.root 2> /dev/null`

iFile=0;
for aFile in $fileList; do
  # echo file is $aFile
  ((iFile++))
  fileName[$iFile]=$aFile
done
nFiles=${#fileName[@]}
if (( nFiles == 0 )); then
  echo ERROR: I could not find any pt resolution root file in $TKG_LAYOUTDIRECTORY
  exit -1
fi

echo
echo "Available resolution files:"
for ((iFile=1; iFile<=$nFiles; iFile++)); do
  echo "File $iFile = ${fileName[iFile]}"
done
echo

selected=""
while [ "$selected" == "" ]; do
  echo -n "Select the desired resolution: "
  read resNum
  selected="${fileName[resNum]}"
  if [ "$selected" == "" ]; then
    echo "ERROR: invalid selection"
  fi
done

# then choose which to delphize
# then start
# root -l ~/

root -l "$TKG_BINDIRECTORY/delphize.cpp(\"${selected}\")"

