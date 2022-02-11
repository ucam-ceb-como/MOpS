#
cd openbabel-2.3.1
mkdir build
cd build
cmake ../ -DCMAKE_INSTALL_PREFIX=./ -DBUILD_GUI=OFF -DBUILD_SHARED=OFF
make -j 
make install
tailMessage="printf '
Be sure to add the following to your ~/.bashrc:\n
\n
# OpenBabel Install:\n
export BABEL_LIB=\"$(pwd)/lib/\"\n
export BABEL_LIBDIR=\"$(pwd)/lib\"\n
export BABEL_SOURCE_DIR=\"$(pwd)/\"\n
export BABEL_DATADIR=\"$(pwd)/share/openbabel/2.3.1/\"\n
\n
export PKG_CONFIG_PATH=$(pwd)/lib/pkgconfig/:\$PKG_CONFIG_PATH\n
export LD_LIBRARY_PATH=\$BABEL_LIB:\$LD_LIBRARY_PATH\n
export LIBRARY_PATH=\$BABEL_LIB:\$LIBRARY_PATH\n
'"
eval $tailMessage
