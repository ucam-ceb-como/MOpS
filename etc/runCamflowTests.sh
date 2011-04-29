#!/bin/bash
# Run camflow as a test.
usage() 
{
    while [ "$#" -ge 1 ]; do echo "$1"; shift; done
    cat<<USAGE
usage: ${0##*/} [OPTION]
  
options:
  -cleanall    Clean all
  -profile     Compile in profile mode
  -debug       Compile in debug mode
  -solver      Name of the solver to compile
  -help        This usage

USAGE
    exit 1
}

debug=0
profile=0
cleanall=0
exeDir=release
solver=camflow

# parse options
while [ "$#" -gt 0 ]
do
   case "$1" in
   -h | -help)
      usage
      ;;
   -debug)
      debug=1
      shift 1   
      ;;
   -profile)
      profile=1
      shift 1   
      ;;
   -cleanall)
      cleanall=1
      shift 1
      ;;
   -solver)
      [ "$#" -ge 2 ] || usage "'-solver' option requires an argument"
      solver=$2
      shift 2
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

if [ $0 != "./etc/runCamflowTests.sh" ]; then 
    echo "Error: Must run this script from mops-c-Git/"
    exit 1
fi

if [ $debug -eq "1" ]; then
    exeDir=debug
fi

if [ $profile -eq "1" ]; then
    exeDir=profile
fi

if [ $cleanall -eq "1" ]; then
    echo "Cleaning the solver."
    make --directory=applications/solvers/$solver -s --file="Makefile" clean
    rm -f bin/debug/$solver
    rm -f bin/release/$solver
    echo "Solver cleaned."
    exit
fi

make --directory=applications/solvers/$solver -s --file="Makefile" clean all debug=$debug profile=$profile
if [ $? -gt 0 ] ; then
    echo make $solver solver failed
    exit
else
    echo make $solver solver success
fi

cd test/camflow/hydrogenFlamelet
time ../../../bin/$exeDir/$solver

test=`diff profile.dat originalProfile.dat`
if [ -z "$test" ] ; then
    echo "Test success."
else
    echo "Test returns $test."
    echo "Test has failed."
    exit
fi
