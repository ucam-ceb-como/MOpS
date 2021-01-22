# <img align="top" src="https://user-images.githubusercontent.com/32139876/92795750-8767bf00-f3a8-11ea-921a-17fac1435171.png" alt="como_group" width="150"/> MOpS Particle Simulator

## Installation:

The current installation of the **MOpS particle simulator** uses [cmake](https://cmake.org/) to compile the program. The installation procedure below refers to installation on Linux systems, but can easily be modified for other operating systems.

### Requirements:

- [Boost C++ Libraries](https://www.boost.org/) v1.47 or later.
- [CMake](https://cmake.org/) v2.6 or later.
- [PkgConfig](https://www.freedesktop.org/wiki/Software/pkg-config/)
- [OpenBabel](http://openbabel.org/wiki/Main_Page) v2.3.1 or later (only for PAH and carbonaceous particle simulations)

### Download the source code

Download the latest version of the code from the [MOpS Github Repository](https://github.com/ucam-ceb-como/MOpS).

Alternatively use [git](https://git-scm.com/) to set up your download:

	git clone https://github.com/ucam-ceb-como/MOpS
	
### Set git variables

Set the following variables (substitute for your name/email):

    git config --global user.name "John Doe"
    git config --global user.email "jdoe@cam.ac.uk"

### Build MOpS

#### Build OpenBabel

The PAH part of the carbonaceous particle model requires [OpenBabel](http://openbabel.org/wiki/Main_Page) to be built as a statically linked library. OpenBabel has been added as a submodule of the MOpS Particle Simulator. 

Before installing OpenBabel make sure you have [PkgConfig](https://www.freedesktop.org/wiki/Software/pkg-config/) on your system. This is used with cmake/modules/FindOpenBabel2.cmake for MOpS to use the [OpenBabel API](http://openbabel.org/dev-api/cmake_project.shtml).
 
The following commands build OpenBabel 2.3.1 as a submodule of MOpS and should be run **before** building MOpS:

    git submodule init
    git submodule update
    cd openbabel-2.3.1-Git
    bash build.sh

You should now set up the OpenBabel environmental variables in your bashrc. A reminder of these variables is shown by the build.sh script. Remember to source the updated bashrc before continuing.

#### Build MOpS (DEBUG)

The following commands should be called from the MOpS Particle Simulator main folder to complete the DEBUG build of the MOpS Particle Simulator:

    mkdir debug
    cd debug
    cmake -D CMAKE_BUILD_TYPE=Debug ..
    make -j 4   #for 4 threads of compilation (more if you can!)
    ctest -j 4  #for 4 threads of testing

#### Build MOpS (RELEASE)

The following commands should be called from the MOpS Particle Simulator main folder to complete the RELEASE build of the MOpS Particle Simulator:

    mkdir release
    cd release
    cmake -D CMAKE_BUILD_TYPE=Release ..
    make -j 4   #for 4 threads of compilation (more if you can!)
    ctest -j 4  #for 4 threads of testing

### Run a case

Examples of different simulations are found in the folder /test/mopsc/

Use the MOpS Particle Simulator help to see the different options for running the program. 

	mops-app --help
	
	Generic options:
	  -h [ --help ]         print usage message
	  -v [ --version ]      print version number
	  -w [ --new ]          use new network interface

	Input file options:
	  -r [ --mops ] arg (=mops.inx)         path to main input file
	  -c [ --chem ] arg (=chem.inp)         path to chemical mechanism
	  -t [ --therm ] arg (=therm.dat)       path to thermochemical data
	  -n [ --trans ] arg (=tran.dat)        path to transport data
	  -g [ --gasphase ] arg (=gasphase.inp) path to gas phase profile
	  -s [ --sweep ] arg (=sweep.xml)       path to particle mechanism
	  -q [ --sensi ] arg (=sensi.xml)       path to sensitivity analysis
	  --schem arg (=surfchem.inp)           path to surface chemical mechanism
	  --stherm arg (=surftherm.dat)         path to surface thermochemical data

	Solver options:
	  -e [ --rand ] arg (=456) adjust random seed value
	  --surf                   turn-on surface chemistry
	  --opsplit                use (simple) opsplit solver
	  --strang                 use strang solver
	  --predcor                use predcor solver
	  --flamepp                use flamepp solver

	Output options:
	  -p [ --postproc ]     postprocess files
	  -o [ --only ]         postprocess files only (don't solve)
	  --diag arg (=0)       set diagnostics level (0-4)
	  --ensemble            write full ensembles to binary files
	  --ppah                write full PAHPP data
	  --ppri                write full primary particle data
	  --jumps               write stochastic jumps data
	  --wdotA4              postprocess based on the molar rate of production 
							                  by 	chemical reaction of the inception species


######################################################################################\
[Computational Modelling Group](https://como.ceb.cam.ac.uk), [University of Cambridge](https://www.cam.ac.uk/). \
10 September 2020 \
######################################################################################
