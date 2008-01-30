#include "gpc_reaction_set.h"
#include "gpc_reaction.h"
#include "gpc_mech.h"
#include "gpc_stoich.h"
#include <math.h>

using namespace Sprog;
using namespace Sprog::Kinetics;
using namespace std;

// CONSTRUCTORS AND DESTRUCTORS.

// Default constructor.
ReactionSet::ReactionSet()
{
    m_mech = NULL;
}

// Copy constructor.
ReactionSet::ReactionSet(const Sprog::Kinetics::ReactionSet &rxn)
{
    m_mech = NULL;
    *this = rxn;
}

// Stream-reading constructor.
ReactionSet::ReactionSet(std::istream &in)
{
    Deserialize(in);
}

// Destructor.
ReactionSet::~ReactionSet()
{
    releaseMemory();
}


// OPERATOR OVERLOADING.

// Assignment operator.
ReactionSet &ReactionSet::operator=(const ReactionSet &rxns)
{
    // Check for self assignment!
    if (this != &rxns) {
        // Clear current memory.
        releaseMemory();

        // Copy the reaction list.  Use the Clone() member function
        // to ensure reactions of the correct type are added.
        RxnPtrVector::const_iterator i;
        for (i=rxns.m_rxns.begin(); i!=rxns.m_rxns.end(); i++) {
            m_rxns.push_back((*i)->Clone());
        }

        // Build reversible reaction map.  Loop over incoming map to
        // get the reaction indices, but remember to use the pointers
        // to the new reactions!
        unsigned int j;
        RxnMap::const_iterator jrxn;
        for (jrxn=rxns.m_rev_rxns.begin(); jrxn!=rxns.m_rev_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_rev_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build forward Landau Teller reaction map.
        for (jrxn=rxns.m_lt_rxns.begin(); jrxn!=rxns.m_lt_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_lt_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build reverse Landau Teller reaction map.
        for (jrxn=rxns.m_revlt_rxns.begin(); jrxn!=rxns.m_revlt_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_revlt_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build third-body reaction map.
        for (jrxn=rxns.m_tb_rxns.begin(); jrxn!=rxns.m_tb_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_tb_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build fall-off reaction map.
        for (jrxn=rxns.m_fo_rxns.begin(); jrxn!=rxns.m_fo_rxns.end(); jrxn++) {
            j = (*jrxn).first;
            m_fo_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }
    }

    return *this;
}

// Compound assignment operator:  Adds the contents of one reaction set
// to this one.
ReactionSet &ReactionSet::operator+=(const ReactionSet &rxns)
{
    // It is currently easier to not allow self-compounding here.
    if (this != &rxns) {
        // Save the current number of reactions, we'll need it later.
        int n = m_rxns.size();

        // Loop over all incoming reactions and copy them to the vector.
        RxnPtrVector::const_iterator i;
        for (i=rxns.m_rxns.begin(); i!= rxns.m_rxns.end(); i++) {
            m_rxns.push_back((*i)->Clone());
        }

        // Build reversible reaction map.  Loop over incoming map to
        // get the reaction indices, but remember to use the pointers
        // to the new reactions!
        unsigned int j;
        RxnMap::const_iterator jrxn;
        for (jrxn=rxns.m_rev_rxns.begin(); jrxn!=rxns.m_rev_rxns.end(); jrxn++) {
            j = n + (*jrxn).first; // Note n!
            m_rev_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build forward Landau Teller reaction map.
        for (jrxn=rxns.m_lt_rxns.begin(); jrxn!=rxns.m_lt_rxns.end(); jrxn++) {
            j = n + (*jrxn).first;
            m_lt_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build reverse Landau Teller reaction map.
        for (jrxn=rxns.m_revlt_rxns.begin(); jrxn!=rxns.m_revlt_rxns.end(); jrxn++) {
            j = n + (*jrxn).first;
            m_revlt_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build third-body reaction map.
        for (jrxn=rxns.m_tb_rxns.begin(); jrxn!=rxns.m_tb_rxns.end(); jrxn++) {
            j = n + (*jrxn).first;
            m_tb_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }

        // Build fall-off reaction map.
        for (jrxn=rxns.m_fo_rxns.begin(); jrxn!=rxns.m_fo_rxns.end(); jrxn++) {
            j = n + (*jrxn).first;
            m_fo_rxns.insert(RxnMap::value_type(j, m_rxns[j]));
        }
    }

    return *this;
}

// Addition operator:  Adds the contents of two reaction sets together.
const ReactionSet ReactionSet::operator+(const ReactionSet &rxns) const
{
    ReactionSet rs(*this);
    rs += rxns;
    return rs;
}

// Subscripting operator:  Provides a different way to access a particular
// reaction by index in the list.
Reaction *const ReactionSet::operator[](unsigned int i)
{
    return m_rxns.at(i);
}

// Subscripting operator:  Provides a different way to access a particular
// reaction by index in the list.
const Reaction *const ReactionSet::operator[](unsigned int i) const
{
    return m_rxns.at(i);
}


// REACTIONS.

// Returns the number of reactions in the set.
unsigned int ReactionSet::Count(void) const
{
    return m_rxns.size();
}

// Returns the vector of all reactions.
const RxnPtrVector &ReactionSet::Reactions() const
{
    return m_rxns;
}

// Returns a pointer to the ith reaction.  Returns NULL if i is invalid.
const Reaction *const ReactionSet::Reactions(unsigned int i) const
{
    if (i < m_rxns.size()) {
        return m_rxns[i];
    } else {
        return NULL;
    }
}

// Adds a reaction to the set.
Reaction *const ReactionSet::AddReaction(const Sprog::Kinetics::Reaction &rxn)
{
    // Clone the reaction and add it to the vector.
    Reaction *pr = rxn.Clone();
    m_rxns.push_back(pr);


    // Check for reverse parameters.
    if (rxn.RevArrhenius() != NULL) {
        m_rev_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    // Check for forward LT parameters.
    if (rxn.LTCoeffs() != NULL) {
        m_lt_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    // Check for reverse LT parameters.
    if (rxn.RevLTCoeffs() != NULL) {
        m_revlt_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    // Check for third-bodies.
    if (rxn.ThirdBodies().size() > 0) {
        m_tb_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    // Check for fall-off parameters.
    if (rxn.FallOffType() != None) {
        m_fo_rxns.insert(RxnMap::value_type(m_rxns.size()-1, pr));
    }

    return pr;
}


// TIDYING UP.

// Clears all reactions from the set.
void ReactionSet::Clear()
{
    releaseMemory();
}


// MOLAR PRODUCTION RATES RATES.

// Calculates the molar production rates of all species given the rate
// of progress of each reaction.
void ReactionSet::GetMolarProdRates(const fvector &rop,
                                    fvector &wdot) const 
{
    int k;
    const RxnStoichMap *mu;
    RxnStoichMap::const_iterator i;

    // Assign sufficient memory for output.
    wdot.resize(m_mech->Species().size());

    // Loop over all species in mechanism.
    for (k=0; k<m_mech->Species().size(); k++) {
        // Reset prod. rate of this species to zero.
        wdot[k] = 0.0;

        // Get the map of reaction stoichiometry for this species.
        mu = &m_mech->GetStoichXRef(k);

        // Sum up the production rate of this species from all reactions.
        //   first  = index of reaction.
        //   second = stoichiometry of species in reaction.
        for (i=mu->begin(); i!= mu->end(); i++) {
            wdot[k] += (*i).second * rop[(*i).first];
        }
    }
}

// REACTION RATES OF PROGRESS.

// Calculates the rate of progress of each reaction given the
// precalculated rate constants.
void ReactionSet::GetRatesOfProgress(real density,
                                     const real *const x,
                                     unsigned int n,
                                     const fvector &kforward,
                                     const fvector &kreverse,
                                     fvector &rop) const
{
    int i, j, k;
    real *rev = new real[m_rxns.size()];

    // Resize output vector to sufficient length.
    rop.resize(m_rxns.size());

    // Loop over all reactions.
    for (i=0; i<m_rxns.size(); i++) {
        // Use rop to store forward rates of production, 
        // and rev to store reverse rates.
        rop[i] = kforward[i];
        rev[i] = kreverse[i];

        // Integer reactants.
        for (k=0; k<m_rxns[i]->ReactantCount(); k++) {
            // As the stoichiometry is integer, it is more computationally efficient
            // to multiply the values together than to use the pow() function.
            for (j=0; j<m_rxns[i]->Reactant(k).Mu(); j++) {
                rop[i] *= density * x[m_rxns[i]->Reactant(k).Index()];
            }
        }

        // Integer products.
        for (k=0; k<m_rxns[i]->ProductCount(); k++) {
            // As the stoichiometry is integer, it is more computationally efficient
            // to multiply the values together than to use the pow() function.
            for (j=0; j<m_rxns[i]->Product(k).Mu(); j++) {
                rev[i] *= density * x[m_rxns[i]->Product(k).Index()];
            }
        }

        // Real reactants.
        for (k=0; k<m_rxns[i]->FReactantCount(); k++) {
            // Now the stoichiometry is non-integer, we must use the pow() function.
            rop[i] *= pow(density * x[m_rxns[i]->FReactant(k).Index()], 
                          m_rxns[i]->FReactant(k).Mu()); 

        }

        // Real products.
        for (k=0; k<m_rxns[i]->FProductCount(); k++) {
            // Now the stoichiometry is non-integer, we must use the pow() function.
            rev[i] *= pow(density * x[m_rxns[i]->FProduct(k).Index()], 
                          m_rxns[i]->FProduct(k).Mu()); 

        }

        // Calculate the net rates of production.
        rop[i] -= rev[i];
    }

    delete [] rev;
}

// Returns the rates of progress of all reactions given the mixture
// object.
void ReactionSet::GetRatesOfProgress(const Sprog::Thermo::GasPhase &mix, 
                                     const fvector &kforward, 
                                     const fvector &kreverse, 
                                     fvector &rop) const
{
    GetRatesOfProgress(mix.Density(), &(mix.MoleFractions()[0]), 
                       m_mech->Species().size(),
                       kforward, kreverse, rop);
}


// RATE CONSTANTS.

// Calculates the forward and reverse rate constants of all reactions
// given the mixture temperature, density and species mole fractions.
void ReactionSet::GetRateConstants(real T,
                                   real density,
                                   const real *const x,
                                   unsigned int n,
                                   const fvector &Gs,
                                   fvector &kforward,
                                   fvector &kreverse) const 
{
    real lnT, invRT, Patm_RT, T_1_3, T_2_3, lowk, pr, logpr;
    RxnPtrVector::const_iterator i;
    RxnMap::const_iterator im;
    int j, k;
    real *tbconcs;
    const Reaction *rxn;

    // Allocate temporary memory.
    tbconcs = new real[m_rxns.size()];

    // Resize output arrays to sufficient length.
    kforward.resize(m_rxns.size());
    kreverse.resize(m_rxns.size());

    // Precalculate some temperature parameters.
    lnT = log(T);
    switch (m_mech->Units()) {
        case SI :
            invRT = 1.0 / (R * T);
            Patm_RT = 101325.0 * invRT;
            break;
        case CGS :
            invRT = 1.0 / (R_CGS * T);
            Patm_RT = 1013250.0 * invRT;
            break;
    }
    if (m_lt_rxns.size() > 0) {
        T_1_3 = 1.0 / pow(T, ONE_THIRD);
        T_2_3 = T_1_3 * T_1_3;
    }

    // Calculate classic Arrhenius forward rate expression.
    for (i=m_rxns.begin(),j=0; i!=m_rxns.end(); i++,j++) {
        kforward[j] = (*i)->Arrhenius().A * 
                      exp(((*i)->Arrhenius().n * lnT) - 
                          ((*i)->Arrhenius().E * invRT));
    }

    // Landau-Teller rate expressions.
    for (im=m_lt_rxns.begin(); im!=m_lt_rxns.end(); im++) {
        j = (*im).first;
        kforward[j] *= exp(((*im).second->LTCoeffs()->B / T_1_3) + 
                           ((*im).second->LTCoeffs()->C / T_2_3));
    }

    // Reverse rate constants.
    for (i=m_rxns.begin(),j=0; i!=m_rxns.end(); i++,j++) {
        if ((*i)->RevArrhenius() != NULL) {
            // This reaction has explicit reverse rate parameters.
            kreverse[j] = (*i)->RevArrhenius()->A * 
                          exp(((*i)->RevArrhenius()->n * lnT) - 
                              ((*i)->RevArrhenius()->E * invRT));
        } else if ((*i)->IsReversible()) {
            // Must find the reverse rate constants by equilibrium.
            kreverse[j] = 0.0;

            // Calculate the Gibbs free energy change for reaction i and sum up
            // the stoichiometric coefficients.
            for (k=0; k<(*i)->ReactantCount(); k++) {
                // Integer Reactants.
                kreverse[j] += (*i)->Reactant(k).Mu() * Gs[(*i)->Reactant(k).Index()];
            }
            for (k=0; k<(*i)->ProductCount(); k++) {
                // Integer Products.
                kreverse[j] -= (*i)->Product(k).Mu() * Gs[(*i)->Product(k).Index()];
            }
            for (k=0; k<(*i)->FReactantCount(); k++) {
                // Real Reactants.
                kreverse[j] += (*i)->FReactant(k).Mu() * Gs[(*i)->FReactant(k).Index()];
            }
            for (k=0; k<(*i)->FProductCount(); k++) {
                // Real Products.
                kreverse[j] -= (*i)->Product(k).Mu() * Gs[(*i)->Product(k).Index()];
            }

            // Calculate the reverse rate constant.
            kreverse[j]  = exp(min(kreverse[j], log(1.0e250)));
            kreverse[j] *= pow(Patm_RT, (*i)->TotalStoich());
            kreverse[j]  = kforward[j] / max(kreverse[j], 1.0e-250);
        }
    }

    // Explicit reverse Landau-Teller parameters.
    for (im=m_revlt_rxns.begin(); im!= m_revlt_rxns.end(); im++) {
        j = (*im).first;
        kreverse[j] *= exp(((*im).second->RevLTCoeffs()->B / T_1_3) + 
                           ((*im).second->RevLTCoeffs()->C / T_2_3));
    }

    // Third body concentrations for each reaction.
    for (i=m_rxns.begin(),j=0; i!=m_rxns.end(); i++,j++) {
        if ((*i)->UseThirdBody()) {
            // Calculate enhanced third body concentration using the enhancement
            // factors defined for this reaction.
            tbconcs[j] = 0.0;
            for (k=0; k<(*i)->ThirdBodyCount(); k++) {
                tbconcs[j] += ((*i)->ThirdBody(k).Mu() - 1.0) * 
                              x[(*i)->ThirdBody(k).Index()];
            }
            tbconcs[j] = density * (1.0 + tbconcs[j]);
        } else {
            // This reaction has no third body requirement.
            tbconcs[j] = 1.0;
        }
    }

    // Pressure dependent fall-off reactions.
    for (im=m_fo_rxns.begin(); im!=m_fo_rxns.end(); im++) {
        j   = (*im).first;  // Reaction index of imth fall-off reaction.
        rxn = (*im).second;

        // Calculate low pressure limit.
        lowk = rxn->FallOffParams().LowP_Limit.A * 
               exp((rxn->FallOffParams().LowP_Limit.n * lnT) - 
                    (rxn->FallOffParams().LowP_Limit.E * invRT));

        // Calculate reduced pressure.
        if (rxn->FallOffParams().ThirdBody >= 0) {
            // A particular species is to be used as the third body.
            pr = lowk * density * x[rxn->FallOffParams().ThirdBody] / kforward[j];
        } else {
            // Use all species as third bodies.
            pr = lowk * tbconcs[j] / kforward[j];
            tbconcs[j] = 1.0;
        }

        // Calculate rate constants based on equation form.
        logpr = log10(pr);
        pr    = pr / (1.0 + pr);
        switch (rxn->FallOffType()) {
            case Troe3: // 3-parameter Troe form.
                pr *= rxn->FTROE3(T, logpr);
                kforward[j] *= pr;
                kreverse[j] *= pr;
                break;
            case Troe4: // 4-parameter Troe form.
                pr *= rxn->FTROE4(T, logpr);
                kforward[j] *= pr;
                kreverse[j] *= pr;
                break;
            case SRI: // SRI form.
                pr *= rxn->FSRI(T, logpr);
                kforward[j] *= pr;
                kreverse[j] *= pr;
                break;
            case Custom: // A custom function is defined to calculate the fall-off form.
                rxn->FallOffFn()(*rxn, lowk, tbconcs[j], T, kforward[j], kreverse[j]);
                break;
        }
    }

    // Apply third-body concentrations to rate constants.
    for (im=m_tb_rxns.begin(); im!=m_tb_rxns.end(); im++) {
        j = (*im).first;
        kforward[j] *= tbconcs[j];
        kreverse[j] *= tbconcs[j];
    }

    // Deallocate temporary memory.
    delete [] tbconcs;
}

// Calculates the forward and reverse rate constants 
// of all reactions given a mixture object.
void ReactionSet::GetRateConstants(const Sprog::Thermo::GasPhase &mix, 
                                   const std::vector<real> &Gs,
                                   std::vector<real> &kforward, 
                                   std::vector<real> &kreverse) const
{
    GetRateConstants(mix.Temperature(), mix.Density(), &(mix.MoleFractions()[0]), 
                     m_mech->Species().size(), Gs, kforward, kreverse);
}


// PARENT MECHANISM.

// Returns a pointer to the parent mechanism.
const Sprog::Mechanism *const ReactionSet::Mechanism() const
{
    return m_mech;
}

// Sets the parent mechanism.
void ReactionSet::SetMechanism(Sprog::Mechanism &mech)
{
    m_mech = &mech;

    // Set mechanism on all reactions as well.
    for (RxnPtrVector::iterator i=m_rxns.begin(); i!=m_rxns.end(); i++) {
        (*i)->SetMechanism(mech);
    }
}


// READ/WRITE/COPY FUNCTIONS.

// Writes the reaction set to a binary data stream.
void ReactionSet::Serialize(std::ostream &out) const
{
    if (out.good()) {
        // Write the serialize version to the stream.
        const unsigned int version = 0;
        out.write((char*)&version, sizeof(version));
        
        // Write the number of reactions to the stream.
        unsigned int n = m_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write the Reaction objects to the stream.
        for (RxnPtrVector::const_iterator i=m_rxns.begin(); i!=m_rxns.end(); i++) {
            (*i)->Serialize(out);
        }

        // Write number of reactions with explicit reverse parameters.
        n = m_rev_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of reactions with explicit reverse parameters.
        for (RxnMap::const_iterator i=m_rev_rxns.begin(); i!=m_rev_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

        // Write number of third body reactions.
        n = m_tb_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of third body reactions.
        for (RxnMap::const_iterator i=m_tb_rxns.begin(); i!=m_tb_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

        // Write number of fall-off reactions.
        n = m_fo_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of fall-off reactions.
        for (RxnMap::const_iterator i=m_fo_rxns.begin(); i!=m_fo_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

        // Write number of reactions with Landau Teller parameters.
        n = m_lt_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of Landau Teller reactions.
        for (RxnMap::const_iterator i=m_lt_rxns.begin(); i!=m_lt_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

        // Write number of reactions with explicit reverse LT parameters.
        n = m_revlt_rxns.size();
        out.write((char*)&n, sizeof(n));

        // Write map of reverse Landau Teller reactions.
        for (RxnMap::const_iterator i=m_revlt_rxns.begin(); i!=m_revlt_rxns.end(); i++) {
            // Write the reaction index in the main vector.  The pointer doesn't
            // need to be written.
            unsigned int ix = (*i).first;
            out.write((char*)&ix, sizeof(ix));
        }

    } else {
        throw invalid_argument("Output stream not ready (Sprog, ReactionSet::Serialize).");
    }
}

// Reads the reaction set data from a binary data stream.
void ReactionSet::Deserialize(std::istream &in)
{
    // Clear the current reaction et.  We do this before checking
    // the stream condition to avoid confusion in the calling code.
    // Even if the possible exception is handled incorrectly, the
    // set will still be empty.
    releaseMemory();

    if (in.good()) {
        // Read the serialized mechanism version.
        unsigned int version = 0;
        in.read(reinterpret_cast<char*>(&version), sizeof(version));

        unsigned int n = 0;

        switch (version) {
            case 0:
                // Read the number of reactions.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read the reactions.
                try {
                    for (unsigned int i=0; i<n; i++) {
                        // Create a reaction object using the stream-reading
                        // constructor.
                        Kinetics::Reaction *rxn = NULL;
                        rxn = new Kinetics::Reaction(in);
                        rxn->SetMechanism(*m_mech);

                        // Add the reaction to the vector.
                        m_rxns.push_back(rxn);
                    }
                } catch (exception &e) {
                    // Clear reaction set memory before throwing error to
                    // higher level.
                    releaseMemory();
                    throw e;
                }

                // Read number of reactions with explicit reverse parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Read map of reactions with explicit reverse parameters.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_rev_rxns[ix] = m_rxns[ix];
                }

                // Write number of third body reactions.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Write map of third body reactions.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_tb_rxns[ix] = m_rxns[ix];
                }

                // Write number of fall-off reactions.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Write map of fall-off reactions.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_fo_rxns[ix] = m_rxns[ix];
                }

                // Write number of reactions with Landau Teller parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Write map of Landau Teller reactions.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_lt_rxns[ix] = m_rxns[ix];
                }

                // Write number of reactions with explicit reverse LT parameters.
                in.read(reinterpret_cast<char*>(&n), sizeof(n));

                // Write map of reverse Landau Teller reactions.
                for (unsigned int i=0; i<n; i++) {
                    // Read reaction index from the stream.
                    unsigned int ix = 0;
                    in.read(reinterpret_cast<char*>(&ix), sizeof(ix));

                    // Add the reaction to map.
                    m_revlt_rxns[ix] = m_rxns[ix];
                }

                break;
            default:
                throw runtime_error("Serialized version number is unsupported (Sprog, ReactionSet::Deserialize).");
        }
    } else {
        throw invalid_argument("Input stream not ready (Sprog, ReactionSet::Deserialize).");
    }
}


// MEMORY MANAGEMENT.

// Clears all memory used by the reaction set.
void ReactionSet::releaseMemory()
{
    // Wipe the cross-referencing maps.
    m_rev_rxns.clear();
    m_lt_rxns.clear();
    m_revlt_rxns.clear();
    m_tb_rxns.clear();
    m_fo_rxns.clear();

    // Delete the reactions.
    RxnPtrVector::iterator i;
    for (i=m_rxns.begin(); i!=m_rxns.end(); i++) {
        delete *i; // Remember to delete memory!
    }
    m_rxns.clear();
}