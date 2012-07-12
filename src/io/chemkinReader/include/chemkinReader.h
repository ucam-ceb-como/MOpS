/*
 * chemkinReader.h
 *
 *  Created on: Jun 22, 2011
 *      Author: lrm29
 *
 *  \todo The element and species parsers currently cannot handle:
 *  SPEC H2
 *  spec O2
 *
 *  ELEM H
 *  elem H
 *
 */

#ifndef CHEMKINREADER_H_
#define CHEMKINREADER_H_

#include "boost/regex.hpp"
#include <string>
#include <vector>
#include <iostream>
#include "phase.h"
#include "element.h"
#include "species.h"
#include "reaction.h"

namespace IO
{

    class ChemkinReader
    {
			/*
			* For surface (and bulk - if included) input file, the element has been defined in chem.inp file 
			* so we only need to consider species 
			* reaction format is also the same as in gas phase
			*/ 
            static const boost::regex elementListRegex;
            static const boost::regex elementSingleRegex;
			static const boost::regex mainSurfaceListRegex;
            static const boost::regex speciesListRegex;
			static const boost::regex surfaceSpeciesListRegex;
            static const boost::regex speciesSingleRegex;
			static const boost::regex surfaceSpeciesSingleRegex;
			static const boost::regex surfaceSpeciesSingleSubRegex; 
            static const boost::regex reactionListRegex;
            static const boost::regex unitsRegex;
	    static const boost::regex surfaceUnitsRegex;

            const std::string chemfile_;
			const std::string chemSurfFile_;
            const std::string thermfile_;
			const std::string thermSurfFile_;
            const std::string transfile_;

            const std::string chemfilestring_;
			const std::string chemSurfFilestring_; 
            std::vector<Element> elements_;
			std::vector<Phase> phase_; 
            std::vector<Species> species_;
            std::vector<Reaction> reactions_;
            std::string globalUnits_;
	    std::string surfUnits_; // added by mm864	

            bool checkChemFile();
			bool checkChemSurfFile();
            void readElements();
			void readPhase();
            void readSpecies();
			void readSurfaceSpecies(); 
            void readReactions();
			void readSurfReactions();
            void readGlobalUnits();
	    void readSurfaceUnits();

        public:

			
			 ChemkinReader
            (
                const std::string chemfile,
                const std::string thermfile,
                const std::string transfile = "NOT READ"
            );
			
			
            ChemkinReader
            (
                const std::string chemfile,
				const std::string chemSurfFile, 
                const std::string thermfile,
				const std::string thermSurfFile,
                const std::string transfile = "NOT READ"
            );

            ~ChemkinReader(){}

            void check();

            void read();

            const std::vector<Element>& elements() const
            {
                return elements_;
            }

            std::vector<Element>& setElements()
            {
                return elements_;
            }

            const std::vector<Species>& species() const
            {
                return species_;
            }

            std::vector<Species>& setSpecies()
            {
                return species_;
            }

			const std::vector<Phase>& phase() const
            {
                return phase_;
            }

            std::vector<Phase>& setPhase()
            {
                return phase_;
            }
			
            const std::vector<Reaction>& reactions() const
            {
                return reactions_;
            }

            std::vector<Reaction>& setReactions()
            {
                return reactions_;
            }

    };

} // namespace IO

#endif /* CHEMKINREADER_H_ */
