/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a default particle used by sweep.  Particles
    are assumed spherical and are only described by their volume.
*/

#ifndef SWEEP_DEFAULTPARTICLE_H
#define SWEEP_DEFAULTPARTICLE_H

#include <vector>
#include "swpparams.h"
#include "swpcomponent.h"

using namespace std;

namespace Sweep
{
class DefaultParticle
{
public:
    static const unsigned int NCACHE = 11; // Number of properties returned by particle.
    /* Indices in property vectors. */
    static const int iAS      = 0;
    static const int iD       = 1;
    static const int iD2      = 2;
    static const int iD_1     = 3;
    static const int iD_2     = 4;
    static const int iM_1_2   = 5;
    static const int iD2M_1_2 = 6;
    static const int iV       = 7;
    static const int iM       = 8;
    static const int iS       = 9;
    static const int iMOM1    = 10;
protected:
    vector<real> m_comp;              // Particle composition.
    vector<real> m_values;            // other particle properties.
    real m_createt;                   // Time at which particle was created.
    vector<Component*> *m_components; // Components used to define particle.
    vector<real> m_cache;             // Precalculated cache of particle properties.
private:
    real m_time; // Last time particle was updated.  Required for LPDA.
public: // Constructors and destructors.
    DefaultParticle(void);
    DefaultParticle(vector<Component*> &components, const unsigned int nvals);
    ~DefaultParticle(void);
    void Initialise(vector<Component*> &components, const unsigned int nvals);
public: // Override routines from BaseParticle.
    void GetProperties(vector<real> &p) const;
    bool IsValid() const;
    DefaultParticle &CreateCopy();
public: // Property get/set routines.
    real GetProperty(const unsigned int i) const;
    /* Returns a reference to the particle composition. */
    const vector<real> &Composition() const;
    /* Sets the particle composition. */
    void SetComposition(const vector<real> &comp); 
    /* Adjusts the particle with the given changes in composition and values. */
    virtual void Adjust(const vector<real> &dcomp, const vector<real> &dvalues);
    /* Same as Adjust but applies it n times. */
    virtual void Adjust(const vector<real> &dcomp, const vector<real> &dvalues, const unsigned int n);
    /* Returns a reference to the particle values. */
    const vector<real> &Values() const;
    /* Sets the particle values. */
    void SetValues(const vector<real> &values);
    /* Time at which particle was created (s). */
    real CreateTime() const;
    /* Sets the particle create time (s). */
    void SetCreateTime(real t);
    /* Returns the last time at which a particle was updated. */
    real Time(void) const;
    /* Sets the last time a particle was updated. */
    void SetTime(const real t);
public: // Particle actions and interactions.
    virtual DefaultParticle &operator=(const DefaultParticle &sp);
    virtual DefaultParticle &operator+=(const DefaultParticle &sp);
    virtual const DefaultParticle operator+(const DefaultParticle &sp) const;
    virtual DefaultParticle &Coagulate(const DefaultParticle &sp);
public: // Physical properties required by Sweep processes.
    real Mass() const; // Particle mass (g).
    real Volume() const; // Particle volume (cm3).
    real SurfaceArea() const; // Particle surface area (cm2).
    real SphereSurface() const; // Equivalent sphere surface area (cm2).
    real CollisionDiameter() const; // Collision diameter (cm).
    real SphereDiameter() const; // Equivalent sphere diameter (cm).
    real CollDiamSquared() const; // Collision diameter squared (cm2).
    real ActiveSurfaceArea() const; // Active surface area (cm2).
    real InvCollDiam() const; // Inverse collision diameter (cm-1).
    real InvCollDiamSquared() const; // Inverse squared collision diameter (cm-2).
    real InvSqrtMass() const; // Inverse of square root of mass (g-1/2).
    real CollDiamSqrdInvSqrtMass() const; // Collision diameter squared times the inverse square root of mass.
public: // Physical properties calculated from property vector (from Ensemble).
    static real Mass(const vector<real> &p); // Particle mass (g).
    static real Volume(const vector<real> &p); // Particle volume (cm3).
    static real SurfaceArea(const vector<real> &p); // Particle surface area (cm2).
    static real SphereSurface(const vector<real> &p); // Equivalent sphere surface area (cm2).
    static real CollisionDiameter(const vector<real> &p); // Collision diameter (cm).
    static real SphereDiameter(const vector<real> &p); // Equivalent sphere diameter (cm).
    static real CollDiamSquared(const vector<real> &p); // Collision diameter squared (cm2).
    static real ActiveSurfaceArea(const vector<real> &p); // Active surface area (cm2).
    static real InvCollDiam(const vector<real> &p); // Inverse collision diameter (cm-1).
    static real InvCollDiamSquared(const vector<real> &p); // Inverse squared collision diameter (cm-2).
    static real InvSqrtMass(const vector<real> &p); // Inverse of square root of mass (g-1/2).
    static real CollDiamSqrdInvSqrtMass(const vector<real> &p); // Collision diameter squared times the inverse square root of mass.
public:
    virtual void CalcCache(void); // Calculates the cache of particle properties.
};

/* Physical properties of the particle required by Sweep processes. */

inline real DefaultParticle::Mass() const {return m_cache[iM];}
inline real DefaultParticle::Volume() const {return m_cache[iV];}
inline real DefaultParticle::SurfaceArea() const {return m_cache[iS];}
inline real DefaultParticle::SphereSurface() const {return PI * pow(6.0 * m_cache[iV] / PI, TWO_THIRDS);}
inline real DefaultParticle::CollisionDiameter() const {return m_cache[iD];}
inline real DefaultParticle::SphereDiameter() const {return pow(6.0 * m_cache[iV] / PI, ONE_THIRD);}
inline real DefaultParticle::CollDiamSquared() const {return m_cache[iD2];}
inline real DefaultParticle::ActiveSurfaceArea() const {return m_cache[iAS];};
inline real DefaultParticle::InvCollDiam() const {return m_cache[iD_1];}
inline real DefaultParticle::InvCollDiamSquared() const {return m_cache[iD_2];}
inline real DefaultParticle::InvSqrtMass() const {return m_cache[iM_1_2];}
inline real DefaultParticle::CollDiamSqrdInvSqrtMass() const {return m_cache[iD2M_1_2];}

/* Physical properties calculated from an ensemble property vector. */

inline real DefaultParticle::Mass(const vector<real> &p) {if (iM < (int)p.size()) return p[iM]; else return 0.0;}
inline real DefaultParticle::Volume(const vector<real> &p) {if (iV < (int)p.size()) return p[iV]; else return 0.0;}
inline real DefaultParticle::SurfaceArea(const vector<real> &p) {if (iS < (int)p.size()) return p[iS]; else return 0.0;}
inline real DefaultParticle::SphereSurface(const vector<real> &p) {return PI * pow(6.0 * p[iV] / PI, TWO_THIRDS);}
inline real DefaultParticle::CollisionDiameter(const vector<real> &p) {if (iD < (int)p.size()) return p[iD]; else return 0.0;}
inline real DefaultParticle::SphereDiameter(const vector<real> &p) {return pow(6.0 * p[iV] / PI, ONE_THIRD);}
inline real DefaultParticle::CollDiamSquared(const vector<real> &p) {if (iD2<(int)p.size()) return p[iD2]; else return 0.0;}
inline real DefaultParticle::ActiveSurfaceArea(const vector<real> &p) {if (iAS<(int)p.size()) return p[iAS]; else return 0.0;};
inline real DefaultParticle::InvCollDiam(const vector<real> &p) {if (iD_1<(int)p.size()) return p[iD_1]; else return 0.0;}
inline real DefaultParticle::InvCollDiamSquared(const vector<real> &p) {if (iD_2<(int)p.size()) return p[iD_2]; else return 0.0;}
inline real DefaultParticle::InvSqrtMass(const vector<real> &p) {if (iM_1_2<(int)p.size()) return p[iM_1_2]; else return 0.0;}
inline real DefaultParticle::CollDiamSqrdInvSqrtMass(const vector<real> &p) {if (iD2M_1_2<(int)p.size()) return p[iD2M_1_2]; else return 0.0;}

};

#endif