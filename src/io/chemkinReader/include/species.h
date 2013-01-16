/*
 * species.h
 *
 *  Created on: Jun 22, 2011
 *      Author: lrm29
 */

#ifndef SPECIES_H_
#define SPECIES_H_

#include "boost/regex.hpp"
#include <string>
#include <iostream>

#include "transport.h"
#include "thermo.h"
#include "element.h"
#include "phase.h" 

namespace IO
{

    class Species
    {

            std::string name_;
	    int siteOccupancy_; // Added by mm864 
            double molecularWeight_;	
	    std::string phaseName_;     
            Transport transport_;
            Thermo thermo_;

            std::map<std::string,double> speciesComposition_;

        public:

            explicit Species
            (
                const std::string name
            );

            ~Species();

            std::string name() const 
            {return name_;}

	    std::string phasename() const
            {return phaseName_;}
			
	    void setSiteOccupancy(const int OccupancyNo);  
		
	    void setPhaseName(const std::string &phName);

	    const int& getSiteOccupancy() const; 

            Transport& transport();
            const Transport& transport() const;

            Thermo& thermo();
            const Thermo& thermo() const;
           

			
            //void checkElementsInSpecies(const std::vector<IO::Element>& elements);

            friend std::ostream& operator<<(std::ostream& output, const Species& element);

    };


} // namespace IO

#endif /* SPECIES_H_ */
