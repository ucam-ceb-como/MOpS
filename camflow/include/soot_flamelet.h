/*!
 * \file   soot_flamelet.h
 * \author Robert I A Patterson
 *
 * \brief Interface for flamelet calculations including soot.
 *
 *  Copyright (C) 2009 Robert I A Patterson.
 *

 Licence:
    This file is part of "camflow".

    brush is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Prof Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
 */
#ifndef SOOT_FLAMELET_H
#define SOOT_FLAMELET_H

#include "reset_chemistry.h"
#include "reactor1d.h"

#include "interface.h"

/*!
 *@brief    Wrapper to combine chemistry and particle solvers for sooty flamelets
 *
 * Chemistry is calculated by Camflow, with Brush (and Sweep) used for the particles.
 *
 * To run a chemistry only calculation, without soot related input files, it should
 * be possible to modify the constructor by removing all of its contents.  The Brush
 * related code in the run method would also need to be removed.  The preferred
 * way of doing such a calculation is to use an instance of Interface, instead of
 * an instance of SootFlamelet.
 */
class SootFlamelet {
public:
    //! Read settings from camflow input files
    SootFlamelet();

    //! Copy the dynamically allocated particle system
    //SootFlamelet(const SootFlamelet &rhs);

    //! Release heap memory
    virtual ~SootFlamelet();

    //! Copy the dynamically allocated particle system
    //SootFlamelet& operator=(const SootFlamelet &rhs);

    //! Advance flamelet time
    void run(const real t_stop, const vector<real>& data_times,
                                const std::vector<real>& mix_frac_diff,
                                const std::vector<real>& grad_mix_frac,
                                const std::vector<real>& lapl_mix_frac,
                                const std::vector<real>& grad_T);

protected:
    //! Build a brush object that will set the gas phase mixture details on the 1d reactor
    Brush::ResetChemistry buildResetChemistry(const vector<real>& data_times,
                                              const std::vector<real>& mix_frac_diff,
                                              const std::vector<real>& grad_mix_frac,
                                              const std::vector<real>& lapl_mix_frac,
                                              const std::vector<real>& grad_T);

    //! Build a brush style geometry object
    Geometry::Geometry1d buildBrushGeometry();

private:
    //! Flamelet chemistry
    Camflow::Interface mChemistry;

    //! Particle system
    Brush::Reactor1d *mParticles;

    //! Has calculation been started
    bool mCalcStarted;

};

#endif //SOOT_FLAMELET_H
