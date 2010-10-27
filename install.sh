#!/bin/bash

myDir=`dirname $0`
TKG_MAIN=$myDir/bin/tklayout
TKG_MATSHOW=$myDir/bin/materialShow
TKG_SETUP_BIN=$myDir/bin/setup.bin

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
echo -n "(Choose a directory in your PATH) [ $HOME/bin ]: "
read TKG_BIN_TARGET
if [ "${TKG_BIN_TARGET}" == "" ] ; then TKG_BIN_TARGET="$HOME/bin" ; fi

if [ ! -d $TKG_BIN_TARGET ] ; then
    echo I cannot find the target directory $(TKG_BIN_TARGET)
    exit -1
fi

# Create/read the configuration file
if $TKG_SETUP_BIN ; then
    # Get the variables from the program
    eval "`$TKG_SETUP_BIN --dirNames`"
    # Get the rest of the configuration from the config file itself
    source $TKG_CONFIGFILE
 
    # Create the empty destination directories
    mkdir -p $TKG_DESTINATION_ROOT
    mkdir -p $TKG_DESTINATION_GRAPH
    mkdir -p $TKG_DESTINATION_SUMMARY

    # Create the destination directories that should be filled
    mkdir -p $TKG_DESTINATION_MATTAB
    mkdir -p $TKG_DESTINATION_XML
    mkdir -p $TKG_DESTINATION_STYLE
    # Copy the installation directories
    cp -rf $myDir/$TKG_SOURCE_MATTAB/* $TKG_DESTINATION_MATTAB \
	&& echo Material config directory created/updated \
	|| echo Failed copying the style directory $TKG_SOURCE_MATTAB to $TKG_DESTINATION_MATTAB
    cp -rf $myDir/$TKG_SOURCE_XML/* $TKG_DESTINATION_XML \
	&& echo Xml directory created/updated \
	|| echo Failed copying the style directory $TKG_SOURCE_XML to $TKG_DESTINATION_XML
    cp -rf $myDir/$TKG_SOURCE_STYLE/* $TKG_DESTINATION_STYLE \
	&& echo Style directory created/updated \
	|| echo Failed copying the style directory $TKG_SOURCE_STYLE to $TKG_DESTINATION_STYLE
    cp -f $TKG_MAIN $TKG_BIN_TARGET \
	&& echo Main program installed/updated \
	|| echo Failed copying the style directory $TKG_MAIN to $TKG_BIN_TARGET
    cp -f $TKG_MATSHOW $TKG_BIN_TARGET \
	&& echo Material helper program installed/updated \
	|| echo Failed copying the style directory $TKG_MATSHOW to $TKG_BIN_TARGET
    if ! $TKG_SETUP_BIN --checkDir ; then
	echo Problem during installation
	exit -1
    else
	echo Installation successful
    fi
else
    echo The setup program could not read/create the configuration file $HOME/.tkgeometryrc properly
fi
