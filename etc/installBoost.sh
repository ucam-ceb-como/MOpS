#!/bin/bash
# Install Boost in ~.
# Written by Laurence McGlashan lrm29@cam.ac.uk

usage() 
{
    while [ "$#" -ge 1 ]; do echo "$1"; shift; done
    cat<<USAGE
usage: ${0##*/} [OPTION]
  
options:
  -root		   Install as root
  -local	   Install locally
  -help        This usage

USAGE
    exit 1
}

root=0
local=0

while [ "$#" -gt 0 ]
do
   case "$1" in
   -h | -help)
      usage
      ;;
   -root)
      root=1
      shift 1
      ;;
   -local)
      local=1
      shift 1
      ;;
   --)
      shift
      break
      ;;
   -*)
      usage "invalid option '$1'"
      ;;
   *)
      break
      ;;
   esac
done

cd ~
mkdir boost
cd boost

wget http://sourceforge.net/projects/boost/files/boost/1.47.0/boost_1_47_0.tar.gz/download -O boost_1_47_0.tar.gz

tar -zxvf boost_1_47_0.tar.gz

cd boost_1_47_0

if [ $local -eq 1 ]; then
    ./bootstrap.sh --prefix=../
    ./b2 install
fi

if [ $root -eq 1 ]; then
    ./bootstrap.sh
    sudo ./b2 install
fi
