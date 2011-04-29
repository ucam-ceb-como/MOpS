#!/bin/bash
# Run some mops tests.
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
exeDir=debug
solver=brush

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

if [ $0 != "./etc/runBrushTests.sh" ]; then 
    echo "Error: Must run this script from top level of project"
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

cd "test/brush"
executable="../../bin/$exeDir/$solver"

# Run the tests and store the return value that should indicate success or failure
# Success is indicated by 0, so initialise all the values to 0
regress1=1
regress2=1
regress3=1
regress4=1

# It is convenient to put the shortest tests first

testregress1="./regress1.sh $executable"
echo "$testregress1"
time $testregress1 > /dev/null
regress1=$?
if [ $regress1 -ne 0 ]; then
    echo "**regress1 failed**"
    echo "continuing with other tests"
fi

testregress2="./regress2.sh $executable"
echo "$testregress2"
time $testregress2 > /dev/null
regress2=$?
if [ $regress2 -ne 0 ]; then
    echo "**regress2 failed**"
    echo "continuing with other tests"
fi

testregress3="./regress3.sh $executable"
echo "$testregress3"
time $testregress3 > /dev/null
regress3=$?
if [ $regress3 -ne 0 ]; then
    echo "**regress3 failed**"
    echo "continuing with other tests"
fi

testregress4="./regress4.sh $executable"
echo "$testregress4"
time $testregress4 > /dev/null
regress4=$?
# No more tests, when adding extra tests add a warning message as for the preceding tests.

if [ $regress1 -eq 0 ]; then
    echo "regress1 passed"
else
    echo "**REGRESS1 FAILED**"
fi

if [ $regress2 -eq 0 ]; then
    echo "regress2 passed"
else
    echo "**REGRESS2 FAILED**"
fi

if [ $regress3 -eq 0 ]; then
    echo "regress3 passed"
else
    echo "**REGRESS3 FAILED**"
fi

if [ $regress4 -eq 0 ]; then
    echo "regress4 passed"
else
    echo "**REGRESS4 FAILED**"
fi
