#include "swp_pri_list_node.h"
#include <memory.h>

namespace Sweep {
namespace AggModels {

// Default constructor
PrimaryListNode::PrimaryListNode(void):
    mPModel(NULL),
    mComp() {}

PrimaryListNode::PrimaryListNode(const Sweep::ParticleModel &model):
    mPModel(&model),
    mComp(model.ComponentCount(), 0.0) 
{}

PrimaryListNode::PrimaryListNode(
    const Sweep::ParticleModel &model,
    const fvector &comp):
        mPModel(&model),
        mComp(comp) {
            if (comp.size() > model.ComponentCount()) 
                throw std::runtime_error("Component vector too large.");
}

PrimaryListNode::PrimaryListNode(const PrimaryListNode &copy):
        mPModel(copy.mPModel),
        mComp(copy.mComp) {}


PrimaryListNode::~PrimaryListNode(void) {}

PrimaryListNode &PrimaryListNode::operator=(const PrimaryListNode &rhs) {
	if (rhs.mPModel == mPModel)
		memcpy(&mComp[0], &rhs.mComp[0], sizeof(double)*mComp.size());
	else {
		mPModel = rhs.mPModel;
		mComp.assign(rhs.mComp.begin(), rhs.mComp.end());
	}
	return *this;
}

std::ostream& operator<<(
    std::ostream &os,
    const Sweep::AggModels::PrimaryListNode &node) {
        os << "[Primary List Node], ";
        for (size_t i(0); i != node.mComp.size(); ++i) {
            os << node.mPModel->Components(i)->Name() << "=" << node.mComp[i] << ", ";
        }
        os << "d=" << node.Diameter() << std::endl;
        return os;
}

double PrimaryListNode::Volume() const {
    double vol(0.0);
    for (size_t i (0); i != mComp.size(); ++i) {
        vol += mComp[i] * mPModel->Components(i)->MolWt() 
                / (Sweep::NA * mPModel->Components(i)->Density());
    }
    return vol;
}

double PrimaryListNode::SurfaceArea() const {
    double d = Diameter();
    return Sweep::PI * d * d;
}

double PrimaryListNode::Mass() const {
    double mass(0.0);
    for (size_t i (0); i != mComp.size(); ++i) {
        mass += mComp[i] * mPModel->Components(i)->MolWt() / Sweep::NA;
    }
    return mass;
}

double PrimaryListNode::Diameter() const {
    return pow(6.0 * Volume() / Sweep::PI, Sweep::ONE_THIRD);
}

unsigned int PrimaryListNode::Adjust(
    const fvector &dcomp,
    unsigned int n) {
    for (size_t i(0); i != mComp.size(); ++i) {
        // TODO: protect against negative dcomps
        mComp[i] += dcomp[i] * (double) n;
    }
    return n;
}



} //AggModels
} //Sweep

