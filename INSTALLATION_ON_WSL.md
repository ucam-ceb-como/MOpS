# <img align="top" src="https://user-images.githubusercontent.com/32139876/92795750-8767bf00-f3a8-11ea-921a-17fac1435171.png" alt="como_group" width="150"/> MOpS Particle Simulator

This document details steps needed to install and run the MOpS Particle Simulator on a Windows Supported Linux (WSL) environment. The installation requires admin (sudo) privileges.

## Install WSL (Ubuntu)
The main procedure can be found on https://docs.microsoft.com/en-us/windows/wsl/install. On an elevated Powershell run

    wsl --install

Setup an user account and password and then allow for the distribution to download and run. The following commands will update the required repositories and should be input on the Linux terminal.

## General Ubuntu updates
    sudo apt update && sudo apt upgrade

## Essential commands and compilers
    sudo apt-get install build-essential

## CMake
Necessary to produce executables for both OpenBabel and MOpS.

    sudo apt install cmake

## wget
Useful to get previous releases

    sudo apt-get update
    sudo apt-get install wget
    wget --version

## Boost
Although Openbabel can compile on modern versions of Boost, MOpS has compatibility issues. 
The current version (tested February 2022) recommended is Boost 1.55.0. 
Download Boost 1.55.0 from the main website https://www.boost.org/users/history/version_1_55_0.html.
Alternatively it can be downloaded with the commands:

    wget -O boost_1_55_0.tar.gz https://sourceforge.net/projects/boost/files/boost/1.55.0/boost_1_55_0.tar.gz/download
    tar xzvf boost_1_55_0.tar.gz
    rm boost_1_55_0.tar.gz

Unzip and install boost:
    
    cd boost_1_55_0
    sudo apt-get install build-essential g++ python-dev autotools-dev libicu-dev libbz2-dev libboost-all-dev
    ./bootstrap.sh
    ./b2
    sudo ./b2 install
    cat /usr/local/include/boost/version.hpp | grep "BOOST_LIB_VERSION"

## gcc4
The Linux distribution probably already contains a version of gcc. However, OpenBabel does not appear to compile with any distribution of gcc>4. For this reason gcc4 is recommended.
The repositories for gcc4 are not easily obtained. <!---The following method can become obsolete at some point.-->

Create a backup for the file /etc/apt/sources.list:

    sudo cp /etc/apt/sources.list /etc/apt/sources.list.bak
    
Add xenial repositories to /etc/apt/sources.list where gcc4 is available.

    echo deb http://dk.archive.ubuntu.com/ubuntu/ xenial main | sudo tee -a /etc/apt/sources.list
    echo deb http://dk.archive.ubuntu.com/ubuntu/ xenial universe | sudo tee -a /etc/apt/sources.list

Update and install gcc4

    sudo apt-get update
    sudo apt update && sudo apt install g++-4.9

Possibly remove xenial repositories from /etc/apt/sources.list at this point.

Assign higher priorities for gcc4. This will probably make gcc4 the preferred compiler for every executable. It is recommended that the update-alternatives parameters can be changed later.

    sudo update-alternatives --install /usr/bin/g++ g++ /usr/bin/g++-4.9 100
    sudo update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-4.9 100
    sudo update-alternatives --install /usr/bin/cpp cpp-bin /usr/bin/cpp-4.9 100
    sudo update-alternatives --set g++ /usr/bin/g++-4.9
    sudo update-alternatives --set gcc /usr/bin/gcc-4.9
    sudo update-alternatives --set cpp-bin /usr/bin/cpp-4.9

## pkg-config
Necessary for the linking between MOpS and Openbabel

    sudo apt-get update -y
    sudo apt-get install -y pkg-config

## Installation of the MOpS Particle Simulator
From this point all the necessary libraries and packages should be installed. Follow the [Installation](./INSTALLATION.md) procedure.

## Possible errors
ZLIB: OpenBabel uses zlib to convert to different file formats. However, the zlib is also used by the gcc linker and can cause issues installing OpenBabel or linking it to MOpS. The typical error will point say "undefined reference to `crc32'". This can be fixed commenting out the line #find_package(ZLIB) in the CMakeLists.txt file contained in the openbabel-2.3.1 folder. The compatibility exporting to those files will be lost but the program will be able to compile and run succesfully.

Boost version: Many distributions of Linux come with Boost libraries installed that are not easily removed. One way to fix this is to modify CMakeLists.txt contained in the MOpS folder by forcing MOpS to look for the exact Boost distribution: find_package(Boost 1.55.0 EXACT REQUIRED COMPONENTS regex serialization filesystem system program_options)

## Debugging using Visual Studio (VS) Code on Windows
The recommended debugger for the MOpS Particle Simulator is Visual Studio Code. It can be downloaded from https://code.visualstudio.com/download.
Some of the steps listed here can be found on https://docs.microsoft.com/en-us/windows/wsl/tutorials/wsl-vscode.
Once downloaded and installed on a Windows machine, VS Code can be set up with the WSL.
Install the Remote Development pack on VS Code: https://marketplace.visualstudio.com/items?itemName=ms-vscode-remote.vscode-remote-extensionpack.
Once installed proceed to set it up on the Linux distribution.

### Update Linux

    sudo apt-get update
    sudo apt-get install wget ca-certificates

### Open VS Code on a Linux folder.
Navigate to the MOpS main folder and open VS Code with:

    code .

### Debugging with GDB
GDB is one of the recommended debuggers and can be installed with the following commands on either the Linux or VS Code terminals

    sudo apt update
    sudo apt install gdb

Add the C++ extension for VS Code https://code.visualstudio.com/docs/languages/cpp and start debugging C++ code. For example, see https://code.visualstudio.com/docs/cpp/config-wsl.

Try debugging with the MOpS particle Simulator. Assuming that all the programs are installed and that MOpS is compiled in Debug mode:
From the main menu, choose Debug: Start debugging, and then choose gdb.
This will create a launch.json file where the debugger can be called. An example of such a file is shown below:
Pay attention to the "program", "args" and "cwd" sections where the MOpS binary file is selected followed by the program arguments and working directory.

<pre>
#### launch.json file example
{
    // Use IntelliSense to learn about possible attributes.
    // Hover to view descriptions of existing attributes.
    // For more information, visit: https://go.microsoft.com/fwlink/?linkid=830387
    "version": "0.2.0",
    "configurations": [
        {
            "name": "(gdb) Launch",
            "type": "cppdbg",
            "request": "launch",
            "program": "${workspaceFolder}/debug/bin/mops-app",
            "args": ["-p", "--flamepp", "-g", "gp_init.inp", "-r", "mops_fixedpop.inx", "--wdotA4", "--ppah", "--ensemble", "-e", "165"],
            "stopAtEntry": false,
            "cwd": "${workspaceFolder}/../../Simulations/DEBUG",
            "environment": [],
            "externalConsole": false,
            "MIMode": "gdb",
            "setupCommands": [
                {
                    "description": "Enable pretty-printing for gdb",
                    "text": "-enable-pretty-printing",
                    "ignoreFailures": true
                },
                {
                    "description": "Set Disassembly Flavor to Intel",
                    "text": "-gdb-set disassembly-flavor intel",
                    "ignoreFailures": true
                }
            ]
        }
    ]
}
#### end of launch.json file example
</pre>

######################################################################################\
[Computational Modelling Group](https://como.ceb.cam.ac.uk), [University of Cambridge](https://www.cam.ac.uk/). \
February 2022 \
######################################################################################
