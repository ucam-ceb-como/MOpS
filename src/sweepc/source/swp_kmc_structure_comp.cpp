#include "swp_kmc_typedef.h"
#include "swp_kmc_structure_comp.h"
#include <list>
#include <cmath>

using namespace Sweep::KMC_ARS;
using namespace std;

Carbon::Carbon(): C1(NULLC),
    C2(NULLC),
    C3(NULLC),
    edge(true),
    bridge(false),
    A(0),
    bondAngle1(0),
    bondAngle2(0),
    coords(0,0)
{
}

Carbon::~Carbon(){
}
