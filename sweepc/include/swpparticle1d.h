/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    Definition of a default particle used by sweep.  Particles
    are assumed spherical and are only described by their volume.

    This class serves as the base class for all other stochastic particle types
    defined in sweep.
*/

#ifndef SWEEP_PARTICLE_H
#define SWEEP_PARTICLE_H

#include <vector>
#include "swpparams.h"
#include "swpcomponent.h"
#include "swpparticlecache.h"

using namespace std;

namespace Sweep
{
class Particle
{
public:
	// Constructors.
    Particle(void);   // Default Constructor.
    Particle(vector<Component*> &components, const unsigned int nvals);  // Parameterised constructor.

	// Destructor.
    ~Particle(void);

	// Initialise Particle - the only function used in the parameterised constructor.
	void Initialise(vector<Component*> &components, const unsigned int nvals);

	// Override routines from BaseParticle.
    void GetProperties(vector<real> &p) const;
    bool IsValid() const;
    Particle &CreateCopy();

	// Property get/set routines.
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

	// Particle actions and interactions.
    virtual Particle &operator=(const Particle &sp);
    virtual Particle &operator+=(const Particle &sp);
    virtual const Particle operator+(const Particle &sp) const;
    virtual Particle &Coagulate(const Particle &sp);

	// Physical properties required by Sweep processes.
    real Mass() const;                    // Particle mass (g).
    real Volume() const;                  // Particle volume (cm3).
    real SurfaceArea() const;             // Particle surface area (cm2).
    real SphereSurface() const;           // Equivalent sphere surface area (cm2).
    real CollisionDiameter() const;       // Collision diameter (cm).
    real SphereDiameter() const;          // Equivalent sphere diameter (cm).
    real CollDiamSquared() const;         // Collision diameter squared (cm2).
    real ActiveSurfaceArea() const;       // Active surface area (cm2).
    real InvCollDiam() const;             // Inverse collision diameter (cm-1).
    real InvCollDiamSquared() const;      // Inverse squared collision diameter (cm-2).
    real InvSqrtMass() const;             // Inverse of square root of mass (g-1/2).
    real CollDiamSqrdInvSqrtMass() const; // Collision diameter squared times the inverse square root of mass.

	// Physical properties calculated from property vector (from Ensemble).
    static real Mass(const vector<real> &p);                    // Particle mass (g).
    static real Volume(const vector<real> &p);                  // Particle volume (cm3).
    static real SurfaceArea(const vector<real> &p);             // Particle surface area (cm2).
    static real SphereSurface(const vector<real> &p);           // Equivalent sphere surface area (cm2).
    static real CollisionDiameter(const vector<real> &p);       // Collision diameter (cm).
    static real SphereDiameter(const vector<real> &p);          // Equivalent sphere diameter (cm).
    static real CollDiamSquared(const vector<real> &p);         // Collision diameter squared (cm2).
    static real ActiveSurfaceArea(const vector<real> &p);       // Active surface area (cm2).
    static real InvCollDiam(const vector<real> &p);             // Inverse collision diameter (cm-1).
    static real InvCollDiamSquared(const vector<real> &p);      // Inverse squared collision diameter (cm-2).
    static real InvSqrtMass(const vector<real> &p);             // Inverse of square root of mass (g-1/2).
    static real CollDiamSqrdInvSqrtMass(const vector<real> &p); // Collision diameter squared times the inverse square root of mass.

	virtual void CalcCache(void); // Calculates the cache of particle properties.

protected:
    vector<real> m_comp;              // Particle composition.
    vector<real> m_values;            // other particle properties.
    real m_createt;                   // Time at which particle was created.
    vector<Component*> *m_components; // Components used to define particle.
    ParticleCache m_cache;            // Precalculated cache of particle properties.

private:
    real m_time; // Last time particle was updated.  Required for LPDA.

	
};



/* Physical properties of the particle required by Sweep processes. */

inline real Particle::Mass() const {return m_cache[m_cache.iM];}
inline real Particle::Volume() const {return m_cache[m_cache.iV];}
inline real Particle::SurfaceArea() const {return m_cache[m_cache.iS];}
inline real Particle::SphereSurface() const {return PI * pow(6.0 * m_cache[m_cache.iV] / PI, TWO_THIRDS);}
inline real Particle::CollisionDiameter() const {return m_cache[m_cache.iD];}
inline real Particle::SphereDiameter() const {return pow(6.0 * m_cache[m_cache.iV] / PI, ONE_THIRD);}
inline real Particle::CollDiamSquared() const {return m_cache[m_cache.iD2];}
inline real Particle::ActiveSurfaceArea() const {return m_cache[m_cache.iAS];};
inline real Particle::InvCollDiam() const {return m_cache[m_cache.iD_1];}
inline real Particle::InvCollDiamSquared() const {return m_cache[m_cache.iD_2];}
inline real Particle::InvSqrtMass() const {return m_cache[m_cache.iM_1_2];}
inline real Particle::CollDiamSqrdInvSqrtMass() const {return m_cache[m_cache.iD2M_1_2];}

/* Physical properties calculated from an ensemble property vector. */

inline real Particle::Mass(const vector<real> &p) {if (m_cache.iM < (int)p.size()) return p[m_cache.iM]; else return 0.0;}
inline real Particle::Volume(const vector<real> &p) {if (m_cache.iV < (int)p.size()) return p[m_cache.iV]; else return 0.0;}
inline real Particle::SurfaceArea(const vector<real> &p) {if (m_cache.iS < (int)p.size()) return p[m_cache.iS]; else return 0.0;}
inline real Particle::SphereSurface(const vector<real> &p) {return PI * pow(6.0 * p[m_cache.iV] / PI, TWO_THIRDS);}
inline real Particle::CollisionDiameter(const vector<real> &p) {if (m_cache.iD < (int)p.size()) return p[m_cache.iD]; else return 0.0;}
inline real Particle::SphereDiameter(const vector<real> &p) {return pow(6.0 * p[m_cache.iV] / PI, ONE_THIRD);}
inline real Particle::CollDiamSquared(const vector<real> &p) {if (m_cache.iD2<(int)p.size()) return p[m_cache.iD2]; else return 0.0;}
inline real Particle::ActiveSurfaceArea(const vector<real> &p) {if (m_cache.iAS<(int)p.size()) return p[m_cache.iAS]; else return 0.0;};
inline real Particle::InvCollDiam(const vector<real> &p) {if (m_cache.iD_1<(int)p.size()) return p[m_cache.iD_1]; else return 0.0;}
inline real Particle::InvCollDiamSquared(const vector<real> &p) {if (m_cache.iD_2<(int)p.size()) return p[m_cache.iD_2]; else return 0.0;}
inline real Particle::InvSqrtMass(const vector<real> &p) {if (m_cache.iM_1_2<(int)p.size()) return p[m_cache.iM_1_2]; else return 0.0;}
inline real Particle::CollDiamSqrdInvSqrtMass(const vector<real> &p) {if (m_cache.iD2M_1_2<(int)p.size()) return p[m_cache.iD2M_1_2]; else return 0.0;}


};

#endif