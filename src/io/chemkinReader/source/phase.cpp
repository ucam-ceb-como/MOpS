/*
 * phase.cpp
 *
 *  Created on: Jun 15, 2012
 *      Author: mm864
 */
 
#include "phase.h"
#include <iomanip>
#include <stdexcept>
#include "stringFunctions.h"

using namespace std;
using namespace boost;

IO::Phase::Phase
(
    const string phaseName,
	const string phaseID,
    const double siteDensity
)
:
    phaseName_(phaseName),
    phaseID_(phaseID),
    siteDensity_(siteDensity), 
	species_map_()
{}

IO::Phase::~Phase(){}

void IO::Phase::setSpecies(std::map<std::string, int> species_map) {
    species_map_ = species_map;
}

const std::map<std::string, int>& IO::Phase::getSpecies() const {
    return species_map_;
}
	
const std::string& IO::Phase::getPhaseName() const
{
return this-> phaseName_;
}

const std::string& IO::Phase::getPhaseID() const
{
return this-> phaseID_;
}


const double& IO::Phase::getSiteDensity() const
{
return this-> siteDensity_; 
}   

namespace IO

{
	ostream& operator<<(ostream& output, const Phase &phase)  

	{
	output << "    Phase Data:\n"
                << "    (\n"
                << "        Phase  : " << phase.phaseName_ << "\n"
				<< "        Phase ID : " << phase.phaseID_ << "\n"
				<< "        Site Density : " << phase.siteDensity_ << "\n"
                << "        Species : {";
	
		std::map<std::string, int>::const_iterator iter, final_iter;
        final_iter = phase.species_map_.end();
        --final_iter;
        for (iter = phase.species_map_.begin(); iter != phase.species_map_.end(); ++iter) {
            output << iter->first << ":" << iter->second;
            if (iter != final_iter) {
                output << ", ";
            }
        }
	
        output << "}" << "\n";
        return output;
    }

}