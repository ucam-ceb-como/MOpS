#include "swp_kmc_typedef.h"
#include "swp_kmc_structure_comp.h"

using namespace Sweep::KMC_ARS;

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

/*!
 * Returns the next carbon atom after current, coming from the previous
 *
 * @param previous    Pointer to previous carbon atom
 * @param current     Pointer to current carbon atom
 * @return            Pointer to next carbon atom
 */
Cpointer Carbon::MoveCPointer(Cpointer &previous, Cpointer &current) {
    if (current == NULL) throw std::runtime_error("No current carbon atom supplied.");
    // Check for bridge
    if(!(current->bridge)) {
    // If not bridge continue to next C
        return current->C2;
    } else {
        if(previous == current->C1) {
            // If coming from main PAH, move iterator to bridged PAH
            return current->C3;
        } else if(previous == current->C3) {
            // If coming from bridged PAH, move to next on main PAH
            return current->C2;
        }
    }
    return current->C2;
}
