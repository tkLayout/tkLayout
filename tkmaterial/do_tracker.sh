#!/bin/bash

if [ "$2" == "" ]; then
    tks=180
else
    tks="$2"
fi

dir=$(dirname "$1")
tkn=$(basename "$1")

echo Directory $dir 
echo Tracker chosen $tkn
echo N of bins $tks

# Destination directories
cmssw_dir="$dir/CMSSW"
mkdir -p "$cmssw_dir"
summaries_dir="$dir/summaries"
mkdir -p "$summaries_dir"
matsum_dir="$dir/matsum"
mkdir -p "$matsum_dir"
store_dir="$dir/store"
mkdir -p "$store_dir"

mkdir -p matsum
mkdir -p xml 

../tkLayout "$1.cfg" "$1_Types.cfg"
./tkmaterial -um "$1.cfg" "$1_Types.cfg" "$1_Materials.cfg" -t $tks -h -x "$tkn"

# Moving the entire directories named after the tracker
mv xml/"$tkn"/* "$cmssw_dir"
rmdir xml/"$tkn"/
mv summaries/"$tkn"/* "$summaries_dir"
rmdir summaries/"$tkn"

# Moving the single output files
mv matsum/"$tkn".* "$matsum_dir"
mv store/"$tkn".root "$store_dir"

