#!/bin/bash

myDir=`dirname $0`
TKG_STYLEDIRECTORY=$myDir/style
TKG_MAIN=$myDir/bin/tklayout
TKG_SETUP_BIN=$myDir/bin/setup.bin

if [ ! -d $TKG_STYLEDIRECTORY ] ; then
  echo I cannot find the installation style directory $TKG_STYLEDIRECTORY
  echo Are you running the install script from the source directory?
  exit -1
fi

if [ ! -f $TKG_MAIN ] ; then
  echo I cannot find the main program binary $TKG_MAIN
  echo Did you run \'make\'?
  exit -1
fi

if [ ! -f $TKG_SETUP_BIN ]; then
  echo I cannot find the setyp program binary $TKG_SETUP_BIN
  echo Did you run \'make\'?
  exit -1
fi

echo "Where should I install the binary?"
echo -n "(Choose a directory in you PATH) [ ~/bin ]: "
read TKG_BIN_TARGET
if [ "${TKG_BIN_TARGET}" == "" ] ; then TKG_BIN_TARGET="$HOME/bin" ; fi

if [ ! -d $TKG_BIN_TARGET ] ; then
  echo I cannot find the target directory $(TKG_BIN_TARGET)
  exit -1
fi

if $TKG_SETUP_BIN ; then
  source ~/.tkgeometryrc
  cp -rf  $TKG_STYLEDIRECTORY $TKG_LAYOUTDIRECTORY \
    && echo Style directory created/updated \
    || echo Failed copying the style directory $TKG_STYLEDIRECTORY to $TKG_LAYOUTDIRECTORY
  cp -f $TKG_MAIN $TKG_BIN_TARGET \
    && echo Main program installed/updated \
    || echo Failed copying the style directory $TKG_MAIN to $TKG_BIN_TARGET
else
  echo The setup program could not read/create the configuration file ~/.tkgeometryrc properly
fi
