#ifndef _LEWIS_H
#define _LEWIS_H

#include "array.h"

namespace Camflow
{

class LewisNumber
{

    const int mCord_;
    const Sprog::Mechanism *const mech_;
    Array2D Le;
    int lewisType_;

    void loadSettings(const std::string& inputFileName);
    void readFromFile(const std::string& fixedLewisFile);
    void calculateLe();

public:

    //! Enumerator for Unity or Fixed Lewis number.
    enum
    {
        UNITY,
        FIXEDFROMFILE,
        CALCULATED
    };

    LewisNumber
    (
        const std::string& inputFileName,
        const Sprog::Mechanism *const mech,
        const int& mCord,
        const int& nSpc
    );

    ~LewisNumber(){}

    inline const int& type() const {return lewisType_;}

    doublereal operator()(const int& Z, const int& species) const;
    doublereal& calcLewis(const int& Z, const int& species);

};

} // End Namespace Camflow

#endif  /* _LEWIS_H */
