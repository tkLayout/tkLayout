#!/bin/bash

myDir=`dirname $0`
TKG_MAIN=$myDir/bin/tklayout
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

# Set current folder as standard directory for the setup script
export TKLAYOUTDIRECTORY=`realpath -s $myDir`

# Create/read the configuration file
if $TKG_SETUP_BIN; then
    # Get the variables from the program
    eval "`$TKG_SETUP_BIN --dirNames`"

    # Get the rest of the configuration from the config file itself
    source $TKG_CONFIGFILE

    mkdir -p $TKG_BINDIRECTORY
    if [ ! -d $TKG_BINDIRECTORY ] ; then
        echo I cannot find the target directory $TKG_BINDIRECTORY
        exit -1
    fi

    echo "Installing main program -> $TKG_BINDIRECTORY"
    cp -f $TKG_MAIN $TKG_BINDIRECTORY \
        && echo "+ Main program installed/updated" \
        || echo "- Failed copying the main program $TKG_MAIN to $TKG_BINDIRECTORY"

    # Copying the directory with .css and all that jazz
    mkdir -p $TKG_LAYOUTDIRECTORY
    if [ ! -d $TKG_LAYOUTDIRECTORY ] ; then
        echo I cannot find the target directory $TKG_LAYOUTDIRECTORY
        exit -1
    fi
    cp -R $TKG_SOURCE_STYLE $TKG_LAYOUTDIRECTORY

    # Copying the xml and the config files to the standard directory
    mkdir -p $TKG_STANDARDDIRECTORY
    if [ ! -d $TKG_STANDARDDIRECTORY ] ; then
        echo I cannot find the target directory $TKG_STANDARDDIRECTORY
        exit -1
    fi
    cp -R --parents $myDir/xml $TKG_STANDARDDIRECTORY/
    cp -R --parents $myDir/config $TKG_STANDARDDIRECTORY/

    if ! $TKG_SETUP_BIN --checkDir ; then
        echo "- Problem during installation"
        exit -1
    else
        echo "+ Installation successful"
    fi
else
    echo "- The setup program could not read/create the configuration file $HOME/.tkgeometryrc properly"
fi

