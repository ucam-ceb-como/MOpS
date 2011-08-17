/*
 * reaction.h
 *
 *  Created on: Jun 25, 2011
 *      Author: lrm29
 */

#ifndef REACTION_H_
#define REACTION_H_

#include <string>
#include <iostream>
#include <map>
#include <vector>

namespace IO
{

    struct Arrhenius
    {
        //! Forward and reverse Arrhenius parameters.
        double A; // Pre-exponential factor.
        double n; // Temperature exponent.
        double E; // Activation energy.
    };

    class Reaction
    {

            //! Is the reaction reversible or not?
            bool flagReversible_;

            //! Are the Arrhenius parameters for the revers reaction explicitally given.
            bool flagHasREV_;

            //! Is this reaction a duplicate?
            bool flagDuplicate_;

            //! reactant & product stoichiometry.
            std::multimap<std::string, double> reactants_, products_;

            //! Total stoichiometry changes.
            //double dstoich_, dreac_, dprod_;

            Arrhenius forwardArrhenius_, reverseArrhenius_;

            // Third bodies.
            //! Set to true if this reaction requires third bodies.
            bool flagThirdBody_;
            bool flagLOW_, flagTROE_, flagSRI_;
            //! Set if (+M) or e.g. (+H2O) is found.
            bool flagPressureDependent_;
            //! Reaction third bodies and their coefficients.
            std::multimap<std::string, double> thirdBodies_;

            std::string fallOffBody_;

            std::vector<double> LOW_, TROE_, SRI_;

        public:

            Reaction();

            ~Reaction(){}

            void setReversible(const bool flag);
            const bool& isReversible() const;
            const bool& hasREV() const;

            void setArrhenius(double A, double n, double E, bool reverse=false);
            const Arrhenius& getArrhenius(bool reverse=false) const;

            void setReactants(std::multimap<std::string, double> reactants);
            const std::multimap<std::string, double>& getReactants() const;

            void setProducts(std::multimap<std::string, double> products);
            const std::multimap<std::string, double>& getProducts() const;

            void setThirdBodies(const std::multimap<std::string, double>& thirdBodies);
            const std::multimap<std::string, double>& getThirdBodies() const;
            void checkForThirdBody(std::multimap<std::string, double>& species);
            bool hasThirdBody() const {return flagThirdBody_;}

            void setFallOffBody(const std::string& fallOffBody);
            const std::string& getFallOffBody() const;

            void setLOW(const std::vector<double>& LOW);
            const std::vector<double>& getLOW() const;
            const bool& hasLOW() const;

            void setTROE(const std::vector<double>& TROE);
            const std::vector<double>& getTROE() const;
            const bool& hasTROE() const;

            void setSRI(const std::vector<double>& SRI);
            const std::vector<double>& getSRI() const;
            const bool& hasSRI() const;

            void setPressureDependent();
            const bool& isPressureDependent() const;

            void setDuplicate();
            const bool& hasDuplicate() const;

            friend std::ostream& operator<<(std::ostream& output, const Reaction& reaction);

    };


} // namespace IO

#endif /* REACTION_H_ */
