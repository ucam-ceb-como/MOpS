#include "gpc_thermo.h"
#include "gpc_params.h"

using namespace Sprog;
using namespace Sprog::Thermo;


// CONSTRUCTORS AND DESTRUCTORS.

ThermoInterface::ThermoInterface(void)
{
}

ThermoInterface::~ThermoInterface(void)
{
}


// INTERNAL ENERGY.

real ThermoInterface::CalcBulkU(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkU(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkU(Sprog::real T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &U) const
{
    return CalcBulkU(T, &x[0], x.size(), U);
}

real ThermoInterface::CalcBulkU(Sprog::real T, 
                                const Sprog::real *const x, 
                                unsigned int n) const
{
    return CalcBulkU(T, x, n, m_tmpvec);
}


// ENTHALPY.

real ThermoInterface::CalcBulkH(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkH(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkH(Sprog::real T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &H) const
{
    return CalcBulkH(T, &x[0], x.size(), H);
}

real ThermoInterface::CalcBulkH(Sprog::real T, 
                                const Sprog::real *const x, 
                                unsigned int n) const
{
    return CalcBulkH(T, x, n, m_tmpvec);
}


// ENTROPY.

real ThermoInterface::CalcBulkS(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkS(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkS(Sprog::real T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &S) const
{
    return CalcBulkS(T, &x[0], x.size(), S);
}

real ThermoInterface::CalcBulkS(Sprog::real T, 
                                const Sprog::real *const x, 
                                unsigned int n) const
{
    return CalcBulkS(T, x, n, m_tmpvec);
}


// GIBBS FREE ENERGY.

real ThermoInterface::CalcBulkG(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkG(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkG(Sprog::real T, 
                                const Sprog::fvector &x, 
                                Sprog::fvector &G) const
{
    return CalcBulkG(T, &x[0], x.size(), G);
}

real ThermoInterface::CalcBulkG(Sprog::real T, 
                                const Sprog::real *const x, 
                                unsigned int n) const
{
    return CalcBulkG(T, x, n, m_tmpvec);
}


// CONSTANT PRESSURE HEAT CAPACITY.

real ThermoInterface::CalcBulkCp(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkCp(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkCp(Sprog::real T, 
                                 const Sprog::fvector &x, 
                                 Sprog::fvector &Cp) const
{
    return CalcBulkCp(T, &x[0], x.size(), Cp);
}

real ThermoInterface::CalcBulkCp(Sprog::real T, 
                                 const Sprog::real *const x, 
                                 unsigned int n) const
{
    return CalcBulkCp(T, x, n, m_tmpvec);
}


// CONSTANT VOLUME HEAT CAPACITY.

real ThermoInterface::CalcBulkCv(Sprog::real T, const Sprog::fvector &x) const
{
    return CalcBulkCv(T, &x[0], x.size());
}

real ThermoInterface::CalcBulkCv(Sprog::real T, 
                                 const Sprog::fvector &x, 
                                 Sprog::fvector &Cv) const
{
    return CalcBulkCv(T, &x[0], x.size(), Cv);
}

real ThermoInterface::CalcBulkCv(Sprog::real T, 
                                 const Sprog::real *const x, 
                                 unsigned int n) const
{
    return CalcBulkCv(T, x, n, m_tmpvec);
}
