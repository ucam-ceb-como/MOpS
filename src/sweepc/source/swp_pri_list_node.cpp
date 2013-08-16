#include "swp_pri_list_node.h"
#include <memory.h>

namespace Sweep {
namespace AggModels {

// Default constructor (private)
PrimaryListNode::PrimaryListNode(void):
    m_pmodel(NULL),
    m_comp() {}

/*!
 * Constructor with particle model
 *
 * @param model    Particle model
 */
PrimaryListNode::PrimaryListNode(const Sweep::ParticleModel &model):
    m_pmodel(&model),
    m_comp(model.ComponentCount(), 0.0)
{}

/*!
 * Constructor with composition vector
 *
 * @param model    Particle model
 * @param comp     Vector of component values to initialise with
 */
PrimaryListNode::PrimaryListNode(
    const Sweep::ParticleModel &model,
    const fvector &comp):
        m_pmodel(&model),
        m_comp(comp) {
            if (comp.size() > model.ComponentCount())
                throw std::runtime_error("Component vector too large.");
}

/*!
 * Copy constructor
 *
 * @param copy    Object to copy
 */
PrimaryListNode::PrimaryListNode(const PrimaryListNode &copy):
        m_pmodel(copy.m_pmodel),
        m_comp(copy.m_comp) {}


/*!
 * Destructor
 */
PrimaryListNode::~PrimaryListNode(void) {
    m_pmodel = NULL;
}

/*!
 * Equals operator definition
 *
 * @param rhs    Node to equate with
 * @return       Copied node
 */
PrimaryListNode &PrimaryListNode::operator=(const PrimaryListNode &rhs) {
    if (rhs.m_pmodel == m_pmodel)
        memcpy(&m_comp[0], &rhs.m_comp[0], sizeof(double)*m_comp.size());
    else {
        m_pmodel = rhs.m_pmodel;
        m_comp.assign(rhs.m_comp.begin(), rhs.m_comp.end());
    }
    return *this;
}

/*!
 * cout << operator definition
 *
 * @param os     Output stream
 * @param node   Node to output
 * @return       Output stream
 */
std::ostream& operator<<(
    std::ostream &os,
    const Sweep::AggModels::PrimaryListNode &node) {
        os << "[Primary List Node], ";
        for (size_t i(0); i != node.m_comp.size(); ++i) {
            os << node.m_pmodel->Components(i)->Name() << "=" << node.m_comp[i] << ", ";
        }
        os << "d=" << node.Diameter();
        os << ", add=" << static_cast<void const *>(&node) << std::endl;
        return os;
}

/*!
 * Return the value of the ith component
 *
 * @param i   Index to get value of
 * @return    Value at index
 */
double PrimaryListNode::Component(unsigned int i) const {
    return m_comp[i];
}

/*!
 * @return    Volume of node
 */
double PrimaryListNode::Volume() const {
    double vol(0.0);
    for (size_t i (0); i != m_comp.size(); ++i) {
        vol += m_comp[i] * m_pmodel->Components(i)->MolWt()
                / (Sweep::NA * m_pmodel->Components(i)->Density());
    }
    return vol;
}

/*!
 * @return    Surface area of node
 */
double PrimaryListNode::SurfaceArea() const {
    double d = Diameter();
    return SurfaceArea(d);
}

/*!
 * Get surface area of a sphere, given a diameter
 *
 * @param diam    Diameter
 * @return        Surface area
 */
double PrimaryListNode::SurfaceArea(double diam) {
    return Sweep::PI * diam * diam;
}

/*!
 * @return    Mass of node (kg)
 */
double PrimaryListNode::Mass() const {
    double mass(0.0);
    for (size_t i (0); i != m_comp.size(); ++i) {
        mass += m_comp[i] * m_pmodel->Components(i)->MolWt() / Sweep::NA;
    }
    return mass;
}

/*!
 * @return    Diameter of node (m)
 */
double PrimaryListNode::Diameter() const {
    return Diameter(Volume());
}

/*!
 * Get diameter of a sphere, given a volume
 *
 * @param vol    Volume
 * @return       Diameter
 */
double PrimaryListNode::Diameter(double vol) {
    return pow(6.0 * vol / Sweep::PI, Sweep::ONE_THIRD);
}

/*!
 * Adjust the particle state space after a surface reaction
 *
 * @param dcomp    Component changes
 * @param n        Number of times to do the adjustment
 * @return         Number of times the adjustment was actually done
 */
unsigned int PrimaryListNode::Adjust(
    const fvector &dcomp,
    unsigned int n) {
    for (size_t i(0); i != m_comp.size(); ++i) {
        // TODO: protect against negative dcomps
        m_comp[i] += dcomp[i] * (double) n;
    }
    return n;
}

/*!
 * @param rhs    Add the rhs's state space to this node
 */
void PrimaryListNode::Merge(PrimaryListNode &rhs) {
    // Add the component counters to this primary
    for (size_t i(0); i != m_comp.size(); ++i) {
        m_comp[i] += rhs.m_comp[i];
    }
}


} //AggModels
} //Sweep

