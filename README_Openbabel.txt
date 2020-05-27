The following lines are needed in .bashrc for Mops and Openbabel to work together. These are defined in Openbabel documentation.

export BABEL_LIB="/home/USERID/Codes/Openbabel/lib/"
export BABEL_LIBDIR="/home/USERID/Codes/Openbabel/lib"
export BABEL_SOURCE_DIR="/home/USERID/Codes/Openbabel/"
export BABEL_DATADIR="/home/USERID/Codes/Openbabel/share/openbabel/2.3.1/"
export LD_LIBRARY_PATH=$BABEL_LIB:$LD_LIBRARY_PATH
#export PKG_CONFIG_PATH=$PKG_CONFIG_PATH:"/home/USERID/Codes/Openbabel/install/lib/pkgconfig/"
export LIBRARY_PATH=$BABEL_LIB:$LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:"/home/USERID/Codes/Openbabel/lib/"
export LIBRARY_PATH=$LIBRARY_PATH:"/home/USERID/Codes/Openbabel/lib/"

