/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    This file defines the PAH particle sub-model used to provide some PAH
    structural information about soot particles.  This is a specialised model
    in sweep which requires certain elements to be defined in the mechanism file.
*/

#ifndef SWEEP_PAHMODEL_H
#define SWEEP_PAHMODEL_H

#include "swpmechanism.h"

namespace Sweep
{
class PAHModel
{
public:
    static int iED; // Index of PAH edges in particle value vector.
    static int iZZ; // Index of PAH zig-zags in particle value vector.
    static int iR5; // Index of PAH R5 rings in a zig-zag in particle value vector.
    static int iAC; // Index of PAH armchairs in particle value vector.
    static int iHL; // Index of PAH holes in particle value vector.

    // Indices of required species in mechanism species list for the calculation of
    // radical site fractions.
    static int C2H2, O2, OH, H, H2, H2O;
public:
    static int Initialise(const Mechanism &mech)
    {
        // This function initialises the PAH model variables for the
        // given mechanism.  Returns negative on error.

        // Get indices of PAH structure values in mechanism value list.
        iED = mech.GetValueIndex("Edges");
        iZZ = mech.GetValueIndex("ZigZags");
        iR5 = mech.GetValueIndex("ZigR5s");
        iAC = mech.GetValueIndex("Armchairs");
        iHL = mech.GetValueIndex("Holes");

        // Check that all the indices are valid, return an error if not.
        if ((iED<0) || (iZZ<0) || (iR5<0) || (iAC<0) || (iHL<0)) {
            return -1;
        }

        // Get indices of species required by model.
        SpeciesList &sp = mech.GetSpeciesList();
        C2H2  = sp.GetIndex("C2H2");
        O2    = sp.GetIndex("O2");
        OH    = sp.GetIndex("OH");
        H     = sp.GetIndex("H");
        H2    = sp.GetIndex("H2");
        H2O   = sp.GetIndex("H2O");

        return 0;
    }

    static real EdgeRadicalFraction(const vector<real> &chem, const real T, const real P)
    {
        real r1f, r1b, r2f, r2b, r3f, r4f, r5f, rdenom, RT;
        RT = RCAL * T;

        // Calculate the forward and back reaction rates.
        r1f = 4.2e+13 * exp(-13.0/RT)				  * chem[H];
        r1b = 3.9e+12 * exp(-11.0/RT)				  * chem[H2];
        r2f = 1.0e+10 * exp(-1.43/RT) * pow(T,0.734)  * chem[OH];
        r2b = 3.68e+8 * exp(-17.1/RT) * pow(T,1.139)  * chem[H2O];
        r3f = 2.0e+13								  * chem[H];
        r4f = 1.1e+07 * exp( -12.96/RT) * pow(T,1.71) * chem[C2H2];
        r5f = 2.2e+12 * exp( -7.5/RT)				  * chem[O2];
        rdenom = r1b+r2b+r3f+r4f+r5f;

        if (rdenom > 0.0) {
            return (r1f+r2f) / rdenom;
        } else {
            return 0.0;
        }
    }

    static real ArmchairRadicalFraction(const vector<real> &chem, const real T, const real P)
    {
        real r1f, r1b, r2f, r2b, r3f, r4f, r5f, rdenom, RT;
        RT = RCAL * T;

        // Calculate the forward and back reaction rates.
        r1f = 4.2e+13 * exp(-13.0/RT)				 * chem[H];
        r1b = 3.9e+12 * exp(-11.0/RT)				 * chem[H2];
        r2f = 1.0e+10 * exp(-1.43/RT) * pow(T,0.734) * chem[OH];
        r2b = 3.68e+8 * exp(-17.1/RT) * pow(T,1.139) * chem[H2O];
        r3f = 2.0e+13								 * chem[H];
        r4f = 8.0e+07 * exp( -3.8/RT) * pow(T,1.56)  * chem[C2H2];
        r5f = 2.2e+12 * exp( -7.5/RT)				 * chem[O2];
        rdenom = r1b+r2b+r3f+r4f+r5f;

        if (rdenom > 0.0) {
            return (r1f+r2f) / rdenom;
        } else {
            return 0.0;
        }
    }

    static real ZigZagRadicalFraction(const vector<real> &chem, const real T, const real P)
    {
        real r1f, r1b, r2f, r2b, r3f, r4f, r5f, rdenom, RT;
        RT = RCAL * T;

        // Calculate the forward and back reaction rates.
        r1f = 4.2e+13 * exp(-13.0/RT)				 * chem[H];
        r1b = 3.9e+12 * exp(-11.0/RT)				 * chem[H2];
        r2f = 1.0e+10 * exp(-1.43/RT) * pow(T,0.734) * chem[OH];
        r2b = 3.68e+8 * exp(-17.1/RT) * pow(T,1.139) * chem[H2O];
        r3f = 2.0e+13								 * chem[H];
        r4f = 6.8e11  * exp(-22.0/RT)                * chem[C2H2];
        r5f = 2.2e+12 * exp(-7.5/RT)				 * chem[O2];
        rdenom = r1b+r2b+r3f+r4f+r5f;

        if (rdenom > 0.0) {
            return (r1f+r2f) / rdenom;
        } else {
            return 0.0;
        }
    }
};

};

#endif