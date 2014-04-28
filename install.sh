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

    SVNURL=`$SVNBIN info | egrep ^URL: | cut -d' ' -f2-`
    echo $SVNURL | grep -q -e '^http' || {
      echo ERROR: the current directory is not under SVN revision
      exit -1
    }
    dirlist="$TKG_SOURCE_MATTAB $TKG_SOURCE_GEOMETRIES $TKG_SOURCE_XML $TKG_SOURCE_STYLE deployStyle"
    for myDir in $dirlist; do
       if [ "$myDir" == "deployStyle" ]; then
          myDir=$TKG_SOURCE_STYLE
          targetDirectory=$TKG_LAYOUTDIRECTORY/$myDir
       else
          targetDirectory=$TKG_STANDARDDIRECTORY/$myDir
       fi
       if [ -d $targetDirectory ] ; then
          echo -n Updating $targetDirectory...
          if [ -d $targetDirectory/.svn ]; then
            otherBase=`$SVNBIN info $targetDirectory | egrep ^URL: | cut -d' ' -f2-`
            if [ "$otherBase" != "$SVNURL/$myDir" ]; then
              echo ""
              echo Directory $targetDirectory is under a different version control than `pwd`. 
              echo Overwrite the old version control with the new one?
              select yn in "Yes" "No"; do
                case $yn in
                  Yes ) rm -rf $targetDirectory; $SVNBIN checkout $SVNURL/$myDir $targetDirectory; break;;
                  * ) echo Installation aborted by user; exit;;
                esac
              done
            fi
            $SVNBIN up $targetDirectory
          else
            echo ""
            echo Directory $targetDirectory is not under a version control. 
            echo Overwrite it with the latest version from svn?
            select yn in "Yes" "No"; do
              case $yn in
                Yes ) rm -rf $targetDirectory; $SVNBIN checkout $SVNURL/$myDir $targetDirectory; break;;
                * ) echo Installation aborted by user; exit;;
              esac
            done
          fi
       else
          echo -n Checking out $targetDirectory...
          mkdir -p $targetDirectory
          $SVNBIN checkout $SVNURL/$myDir $targetDirectory
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
