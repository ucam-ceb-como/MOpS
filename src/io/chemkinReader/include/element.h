/*
 * element.h
 *
 *  Created on: Jun 22, 2011
 *      Author: lrm29
 */

#ifndef ELEMENT_H_
#define ELEMENT_H_

#include <string>
#include <iostream>

namespace IO
{

    class Element
    {

            std::string name_;
            double atomicWeight_;

        public:

            Element
            (
                const std::string name,
                const double atomicWeight
            );

            ~Element(){}

            const std::string& getName() const;

            const double& getAtomicWeight() const;

            friend std::ostream& operator<<(std::ostream& output, const Element& element);

    };


} // namespace IO

#endif /* ELEMENT_H_ */
