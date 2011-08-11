#!/bin/bash
# Install script.
# lrm29@cam.ac.uk

usage() 
{
    while [ "$#" -ge 1 ]; do echo "$1"; shift; done
    cat<<USAGE
usage: ${0##*/} [OPTION]
  
options:
  -help        This usage

USAGE
    exit 1
}

# parse options
while [ "$#" -gt 0 ]
do
   case "$1" in
   -h | -help)
      usage
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

echo "Make all Documentation."
./etc/makeDocs.sh

echo "Install Successful."
