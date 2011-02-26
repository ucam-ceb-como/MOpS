#!/bin/bash
# Install script.
# lrm29@cam.ac.uk

usage() 
{
    while [ "$#" -ge 1 ]; do echo "$1"; shift; done
    cat<<USAGE
usage: ${0##*/} [OPTION]
  
options:
  -test        Run Make and test running
  -help        This usage

USAGE
    exit 1
}

test=0

# parse options
while [ "$#" -gt 0 ]
do
   case "$1" in
   -h | -help)
      usage
      ;;
   -test)
      test=1
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

echo "Configuring local git settings."
cp -f .gitconfig .git/config

if [ $test -eq "1" ]; then
    echo "Make all libraries in release mode."
    ./etc/makeAllLibs.sh
    ./etc/makeDocs.sh
    echo "Run camflow in release mode."
    ./etc/runCamflow.sh
    if [ $? -gt 0 ] ; then
        echo "Running camflow failed."
        exit
    else
        echo "Install Successful."
    fi
fi
