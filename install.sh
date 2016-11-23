#!/bin/bash

myDir=`dirname $0`
TKG_MAIN=$myDir/bin/tklayout
#TKG_MATSHOW=$myDir/bin/materialShow
#TKG_TUNE=$myDir/bin/tunePtParam
TKG_SETUP_BIN=$myDir/bin/setup

if [ ! -f $TKG_MAIN ] ; then
    echo I cannot find the main program binary $TKG_MAIN
    echo Did you run \'make\'?
    exit -1
fi

if [ ! -f $TKG_SETUP_BIN ]; then
    echo I cannot find the setup program binary $TKG_SETUP_BIN
    echo Did you run \'make\'?
    exit -1
fi

# Create/read the configuration file
if $TKG_SETUP_BIN ; then
    # Get the variables from the program
    eval "`$TKG_SETUP_BIN --dirNames`"
    # Get the rest of the configuration from the config file itself
    source $TKG_CONFIGFILE

    mkdir -p $TKG_BINDIRECTORY

    if [ ! -d $TKG_BINDIRECTORY ] ; then
        echo I cannot find the target directory $(TKG_BINDIRECTORY)
        exit -1
    fi

    cp -f $TKG_MAIN $TKG_BINDIRECTORY \
	&& echo Main program installed/updated \
	|| echo Failed copying the main program $TKG_MAIN to $TKG_BINDIRECTORY

    # Copying the directory with .css and all that jazz
    cp -R $TKG_SOURCE_STYLE $TKG_LAYOUTDIRECTORY
          

    if ! $TKG_SETUP_BIN --checkDir ; then
	echo Problem during installation
	exit -1
    else
	echo Installation successful
    fi
else
    echo The setup program could not read/create the configuration file $HOME/.tkgeometryrc properly
fi

cd $myDir
if [ "$TKG_STANDARDDIRECTORY" != `pwd` ]
then echo -e "ERROR: your config file $TKG_CONFIGFILE should contain the following line:\nTKG_STANDARDDIRECTORY=\"`pwd`\""
fi

