/*
  Author(s):      Matthew Celnik (msc37)
  Project:        none

  File purpose:
    Header file to include in any program that uses sweep.  This file includes
    all header files in the sweep code.
*/

#ifndef SWEEP_H
#define SWEEP_H

#include "swpparams.h"
#include "swpsolver.h"
#include "swpmechanism.h"
#include "swpensemble.h"
#include "swpcomponent.h"
#include "swpparticlechanger.h"

/* Particle types. */
#include "swpparticle1d.h"

/* Classes which define different systems that can be solved using sweep. */
#include "swpsystem.h"
#include "swphomogas.h"
#include "swppremixflame.h"

/* Classes which define different particle processes. */
#include "swpprocess.h"
#include "swpinception.h"
#include "swpcoagulation.h"
#include "swpsurfacereaction.h"
#include "swpactivesitesreaction.h"
#include "swpcondensation.h"

/* Input and output code. */
#include "swpxmlio.h"

/* Hardcoded particle and chemical models. */
#include "swpabf.h"

#endif