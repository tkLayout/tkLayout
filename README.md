# Getting the repository

```
git clone https://github.com/tkLayout/tkLayout.git
cd tkLayout
```

The latest official version is at `tkLayout/main`.


# Environment

**tkLayout** requires:

```
ROOT version 6.30.04 (root and root-config should be in user's path)
BOOST version 1.82.0
gcc version 13.1.0
```

If you are on AlmaLinux-9 and your environment can see the CVMFS, you can run a bash shell and then from inside the repository home directory:

```
source setup_el9.sh
```

**tkLayout** is also compatible with Ubuntu, though it is not the main supported version.

# Compilation

Compilation using make:

```
make -j $(nproc)
make install
make doc # Optional, generate Doxygen-based documentation in doc directory
```

This will build the binaries and place them in the `tkLayout/bin` directory. If this is your first time installing **tkLayout**, please note that you must follow the instructions in [**First time install**](#first-time-install) when running `make install`.


### First-time install

If this is the first time installing **tkLayout**, you will need to provide:

1. The directory where you want to store the **tkLayout** executable (e.g. `/eos/user/y/yourname/bin`).
2. The destination directory for the web output of the program (e.g. `/eos/user/y/yourname/www/layouts`). This directory needs to be writable by you and readable by the web server. The style files will be copied during the installation to that directory. If the style files get changed during the development, a new `./install.sh` will be needed to get this propagated to the output. Alternatively, you can choose to make a symbolic link between the `tkLayout/style` directory and the `www/layouts` directory. This avoids the need of repeating the `./install.sh` at every update of the style files (but in order to do this, the `tkLayout/style` directory should also be within the reach of the web server).
3. The path to your **tkLayout** source directory (e.g. `/eos/user/y/yourname/tkLayout`).
4. The set of momentum values to be used for the simulation.

**tkLayout** writes the installation configuration to `~/.tkgeometryrc`. In the case of a misinput, you can directly modify the file and rerun `make install`.

# Command line options

Once **tklayout** is built and installed, you can of course see the command line options with:

```
tklayout --help
```

# Have fun!
