#include "swp_kmc_typedef.h"
#include "swp_kmc_structure_comp.h"
#include <list>
#include <cmath>

using namespace Sweep::KMC_ARS;
using namespace std;

Carbon::Carbon(): C1(NULL),
    C2(NULL),
    C3(NULL),
    edge(true),
    bridge(false),
    A(0),
    bondAngle1(0),
    bondAngle2(0),
    coords(0,0)
{
}
Carbon::Carbon(const Carbon& copy){

}
Carbon::~Carbon(){
}