/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    A class which reads a mechanism from an XML file into a Mechanism
    object.
*/

#ifndef SWEEP_XMLIO_H
#define SWEEP_XMLIO_H

#include <string>
#include "swpmechanism.h"
#include "camxml.h"

namespace Sweep
{
class XMLIO
{
private:
    enum RXNTYPE {Surf,ActiveSites,Cond};
    enum ACTSURFMODEL {Const,ABF,Profile};
public:
    XMLIO(void);
    ~XMLIO(void);
public:
    /* Reads the given xml file into the given mechanism object. */
    static int ReadMechanism(const string &filename, Mechanism &mech);
private:
    /* Reads a version 1 sweep mechanism.  This is also the default is no
       version is specified. */
    static int readVersion1(CamXML::Document &xml, Mechanism &mech);
};
};

#endif