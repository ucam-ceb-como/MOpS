/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    This class defines the ABF Hydrogen Abstraction - C2H2 Addtion (HACA) model
    for soot particles as discussed by Appel et al. (2000).  This model includes
    routines for calculating the fraction of radical sites on soot particle surfaces
    using the steady-state assumption and a routine for for returning the number of
    active sites using the alpha correlation given by Appel et al. (2000).  These
    routines are only valid for the surface growth model described in that paper.
*/

#ifndef SWEEP_ABFMECH_H
#define SWEEP_ABFMECH_H

#include "swpsystem.h"
#include "swpparams.h"
#include "gpcspecieslist.h"
#include <cmath>
#include <vector>

namespace Sweep
{
namespace ABF
{
class ABFMech
{
public:
    static int A4, C2H2, O2, OH, CO, H, H2, H2O, Alpha;
public:
    /* Initialises the HACA mechanism by saving the indices of
       species required by the functions. */
    inline static void InitHACA(SpeciesList &sp)
    {
        A4    = sp.GetIndex("A4");
        C2H2  = sp.GetIndex("C2H2");
        O2    = sp.GetIndex("O2");
        OH    = sp.GetIndex("OH");
        CO    = sp.GetIndex("CO");
        H     = sp.GetIndex("H");
        H2    = sp.GetIndex("H2");
        H2O   = sp.GetIndex("H2O");
        Alpha = sp.GetIndex("Alpha");
    };

    /* Returns the fraction of active sites which are radical sites for the
       HACA mechanism. */
    inline static real RadicalSiteFraction(const real t, const std::vector<real> &chem, 
                                           const real T, const real P)
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
    };

    /* Calculates the number of active radical sites per unit surface area of
       particles given the chemical conditions and particle sums. */
    static real HACASites(const real t, const System &sys, 
                          const vector<real> &chem, const real T, 
                          const real P, const vector<real> &sums)
    {
        if (Alpha >= 0) {
            // 2.3e15 is the surface site concentration given by Frenklach.  This
            // must be scaled be the fraction of those sites which are active (alpha)
            // and the fraction of those sites which are radicalised (hydrogen abstracted).
            return 2.3e15 * chem[Alpha] * RadicalSiteFraction(t, chem, T, P);
        } else {
            return 0.0;
        }
    };
};
};
};

#endif