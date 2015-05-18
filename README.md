# Getting the code

   git clone https://github.com/tkLayout/tkLayout.git
   cd tkLayout

# Before the compilation:

You need a working version of root active (root and root-config should be in your path) and a few
libraries are required:
  * boost_filesystem
  * boost_regex
  * the root library set

If you are on lxplus you just have to run a bash shell and then

     source setup_slc6.sh


# Compilation
  Simply type:

      make

This will build the needed programs and put them in the ./bin directory.


# Install
  If the make command runs properly you can install the program with the script

     make install

## First-time install
If this is the first time that you install tkLayout, a few questions will be asked and a configuration
fle will be created in $HOME/.tkgeometry:

1. You will need to provide the destination directory for the www output of the program (proabably
  something like /afs/cern.ch/user/y/yourname/www/layouts) this directory needs to be writable by you
  and readable by the web server.
     The style files will be copied during the installation to that directory (for example
  /afs/cern.ch/user/y/yourname/www/layouts/style ). If the style files get changed during the development
  a new ./install.sh will be needed to get this propagated to the output.
     Alternatively, you can choose to make a symbolic link between the source directory style directory and
  the layout directory. This avoids the need of repeating the ./install.sh at every update of the style files
  but in order to do this the source directory should also be within the reach of the web server.
1. The set of momentum values to be used for the simulation
1. [...] t.b.d.

# Update
To get the latest development version (usually UNSTABLE) you simply type

   git fetch
   git chechout master
   make
   make install
  
