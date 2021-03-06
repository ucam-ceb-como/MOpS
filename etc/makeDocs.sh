#!/bin/bash
# Make all the documentation in doc/.
# Written by Laurence McGlashan lrm29@cam.ac.uk

usage() 
{
    while [ "$#" -ge 1 ]; do echo "$1"; shift; done
    cat<<USAGE
usage: ${0##*/} [OPTION]
  
options:
  -cleanall    Clean all documentation in doc/
  -help        This usage

USAGE
    exit 1
}

cleanall=0

# parse options
while [ "$#" -gt 0 ]
do
   case "$1" in
   -h | -help)
      usage
      ;;
   -cleanall)
      cleanall=1
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

#if [ $0 != "./etc/makeDocs.sh" ]; then 
#    echo "Error: Must run this script from mops-c-Git/"
#    exit 1
#fi

folders=(utils camxml sprogc geometry camflow sweepc mopsc brush)

if [ $cleanall -eq 1 ]; then
    echo "Cleaning all documentation."
    for element in ${folders[@]}
    do
        rm -rf doc/$element
        echo "$element cleaned."
    done
    rm -rf doc/mopssuite
    rm -rf doc/index.html
    echo "Documentation cleaned successfully."
    exit
fi

echo "Create the documentation:"
echo ${folders[@]}
for element in ${folders[@]}
do
    folderLocation=`find src -name $element`
    doxygen $folderLocation/Doxyfile
    if [ $? -gt 0 ] ; then
        echo "Doxygen on $element failed."
        exit
    else
        echo "Doxygen on $element successful."
    fi
done

doxygen doc/Doxyfile
if [ $? -gt 0 ] ; then
    echo "Doxygen on docs/Doxyfile failed."
    exit
else
    echo "Doxygen on docs/Doxyfile successful."
fi

rm -f doc/documentationIndex.html
cp -f doc/mopssuite/html/index.html doc/documentationIndex.html
cp -f doc/mopssuite/html/doxygen.css doc/doxygen.css
cp -f doc/mopssuite/html/tabs.css doc/tabs.css
cp -f doc/mopssuite/html/*.png doc/

echo "Documentation successfully created and stored in doc/"
