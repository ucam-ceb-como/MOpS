/*
  Author(s):      Matthew Celnik (msc37)
  Project:        sweep (population balance solver)

  File purpose:
    This file defines classes which evalulate different functional
    forms.  These are useful if, for example, you want to implement
    different simple models for parameters at run time, instead
    of compile time.
*/

#ifndef SWEEP_MATHS_FUNCTIONAL_H
#define SWEEP_MATHS_FUNCTIONAL_H

#include "swp_params.h"
#include <iostream>

namespace Sweep
{
namespace Maths
{
// Forward declare functional forms.
class Linear;

// Enumeration of functions forms.
enum FuncType {
    LinearFunc = 0
};

// The Functional class defines the interface for any
// function which has one independent variable.
class Functional
{
public:
    // INTERFACE.
    // Sets the value of the ith parameter.
    virtual void SetParam(unsigned int i, real val) = 0;
    // Returns the function evaluation.
    virtual real Evaluate(real x) const = 0;
    // Outputs the functional to a binary stream.
    virtual void Serialize(std::ostream &out) const = 0;
    // Read the functional from a binary stream.
    virtual void Deserialize(std::istream &in) = 0;
    // Returns the functional type ID.
    virtual FuncType Type(void) const = 0;
    // FACTORY METHODS.
    // Writes a linear functional to a binary stream and also writes
    // the functional type so that it can be successfully
    // recreated.
    static void Write(std::ostream &out, const Functional &fun);
    // Reads a functional from a binary stream, also reads
    // the functional type from the stream and creates an
    // object of that type.
    static Functional *const Read(std::istream &in);
};

// The Linear class defines a linear functional form.
class Linear : public Functional
{
public:
    // Default constructor.
    Linear(void) : m_m(0.0), m_c(0.0) {}
    
    // Copy constructor.
    Linear(const Linear &copy) {*this = copy;}
    
    // Destructor.
    ~Linear(void) {}

    // Operators.
    Linear &operator=(const Linear &rhs) {
        if (this != &rhs) {
            m_m = rhs.m_m;
            m_c = rhs.m_c;
        }
    }

    // PARAMETERS.

    // Returns the gradient.
    inline real Gradient(void) const {return m_m;}

    // Sets the gradient.
    inline void SetGradient(real m) {m_m = m;}

    // Returns the constant term.
    inline real Constant(void) const {return m_c;}

    // Sets the constant term.
    inline void SetConstant(real c) {m_c = c;}

    // FUNCTIONAL INTERFACE.

    // Sets the value of the ith parameter.
    inline void SetParam(unsigned int i, real val) {
        if (i==0)
            m_m = val;
        else
            m_c = val;
    }

    // Returns the function evaluation.
    inline real Evaluate(real x) const {return (m_m*x)+m_c;}

    // Outputs the functional to a binary stream.
    inline void Serialize(std::ostream &out) const {
        if (out.good()) {
            double val = (double)m_m;
            out.write((char*)&val, sizeof(val));
            val = (double)m_c;
            out.write((char*)&val, sizeof(val));
        }
    }

    // Read the functional from a binary stream.
    inline void Deserialize(std::istream &in) {
        if (in.good()) {
            double val = 0.0;
            in.read(reinterpret_cast<char*>(&val), sizeof(val));
            m_m = (real)val;
            in.read(reinterpret_cast<char*>(&val), sizeof(val));
            m_c = (real)val;
        }
    }

    // Returns the functional type ID.
    inline FuncType Type(void) const {return LinearFunc;}

private:
    // Gradient and constant parameters.
    real m_m, m_c;
};
};
};

// Definitions of static Functional functions.

inline void Sweep::Maths::Functional::Write(std::ostream &out, const Functional &fun)
{
    Sweep::Maths::FuncType t = fun.Type();
    out.write((char*)&t, sizeof(t));
    fun.Serialize(out);
}

inline Sweep::Maths::Functional *const Sweep::Maths::Functional::Read(std::istream &in)
{
    Sweep::Maths::FuncType t; 
    Sweep::Maths::Functional *fun = NULL;
    in.read(reinterpret_cast<char*>(&t), sizeof(t));
    switch(t) {
        case Sweep::Maths::LinearFunc:
        default:
            fun = new Sweep::Maths::Linear();
    }
    fun->Deserialize(in);
    return fun;
}

#endif
