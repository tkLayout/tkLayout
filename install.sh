#!/bin/bash

myDir=`dirname $0`
TKG_MAIN=$myDir/bin/tklayout
#TKG_MATSHOW=$myDir/bin/materialShow
#TKG_TUNE=$myDir/bin/tunePtParam
TKG_SETUP_BIN=$myDir/bin/setup.bin
SVNBIN=`which svn`

if [ ! -f $SVNBIN ] ; then
    echo I cannot find svn binary in your path
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

    SVNURL=`$SVNBIN info | grep URL | cut -d' ' -f2-`
    echo $SVNURL | grep -q -e '^http' || {
      echo ERROR: the current directory is not under SVN revision
      exit -1
    }
    dirlist="$TKG_SOURCE_MATTAB $TKG_SOURCE_GEOMETRIES $TKG_SOURCE_XML $TKG_SOURCE_STYLE"
    for myDir in $dirlist; do
       if [ -d $TKG_STANDARDDIRECTORY/$myDir ] ; then
          otherBase=`$SVNBIN info $TKG_STANDARDDIRECTORY/$myDir | grep URL | cut -d' ' -f2-`
          if [ "$otherBase" != "$SVNURL/$myDir" ]; then
            echo ERROR: directory $TKG_STANDARDDIRECTORY/$myDir is not under the same version control as `pwd`
            exit -1
          fi
          echo -n Updating $myDir...
          $SVNBIN up $TKG_STANDARDDIRECTORY/$myDir
       else
          echo -n Checking out $myDir...
          mkdir -p $TKG_STANDARDDIRECTORY/$myDir
          $SVNBIN checkout $SVNURL/$myDir $TKG_STANDARDDIRECTORY/$myDir
       fi
    done

    if ! $TKG_SETUP_BIN --checkDir ; then
	echo Problem during installation
	exit -1
    else
	echo Installation successful
    fi
else
    echo The setup program could not read/create the configuration file $HOME/.tkgeometryrc properly
fi
