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
#include <string>
using namespace std;

namespace IO
{

    struct Arrhenius
    {
        //! Forward and reverse Arrhenius parameters.
        double A; // Pre-exponential factor.
        double n; // Temperature exponent.
        double E; // Activation energy.
    };

	struct Cov
	{
		double Eta; // Eta.
		double Miu; // Miu.
		double Epsilon; // e.
		std::string spName; // species name associated with the coverage param. 
		 
		template<class Archive>
		void serialize(Archive & ar, const unsigned int /* file_version */)
		{
          ar & spName & Eta & Miu & Epsilon;
		}
	};
	struct Ford
	{
		double F_k; // Forward reaction order.
		std::string spName; // species name associated with the coverage param. 
		
		template<class Archive>
		void serialize(Archive & ar, const unsigned int /* file_version */)
		{
          ar & spName & F_k;
		}
	}; 
	
	typedef std::vector<Cov> CovVector;
	typedef std::vector<Ford> FordVector; 
	
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
			bool flagFORD_, flagSTICK_, flagCOV_, flagMWON_;
            //! Set if (+M) or e.g. (+H2O) is found.
            bool flagPressureDependent_;
            //! Reaction third bodies and their coefficients.
            std::multimap<std::string, double> thirdBodies_;

            std::string fallOffBody_;
			
			Cov coverageParam_; 
			Ford fordParam_; 

			//! Set to true if this surface reaction has flag
			
            std::vector<double> LOW_, TROE_, SRI_; // No need MWOFF_ or MWON_, STICK_ since both using ARRHENIUS PARAMS (NO ADDITIONAL PARAMS). 
			CovVector COV_; 
			FordVector FORD_; 
			
		
			
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

	    void setFORD(double F, std::string &name);
            const std::vector<Ford>& getFORD() const;
	    //Ford *const getFORDElement(unsigned int a) const;
            const bool& hasFORD() const;
	    unsigned int fordParamCount() const;  
			
	    void setCOV(double e, double m, double eps, std::string &name);
            const std::vector<Cov>& getCOV() const;
	    //Cov *getCOVElement(unsigned int a) const;
            const bool& hasCOV() const;
	    unsigned int coverageParamCount() const;  
			
	    void setSTICK();
            const bool& hasSTICK() const;
			
	    void setMWON();
            const bool& hasMWON() const;
			
            void setPressureDependent();
            const bool& isPressureDependent() const;

            void setDuplicate();
            const bool& hasDuplicate() const;

            friend std::ostream& operator<<(std::ostream& output, const Reaction& reaction);

    };


} // namespace IO

#endif /* REACTION_H_ */
