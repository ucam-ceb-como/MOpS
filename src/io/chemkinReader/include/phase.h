/*
 *phase.h
 *
 * Created on June, 15 2012
 * Author: M. Martin (mm864)
 *
 */

#ifndef PHASE_H_ 
#define PHASE_H_


#include "boost/regex.hpp"
#include <string>
#include <iostream>
#include <map>
#include <vector>

namespace IO
{

  class Phase 
  {
    std::string phaseName_;
    std::string phaseID_;
    double siteDensity_;

	public:

	Phase(){}

    Phase
      (
       const std::string phaseName,
       const std::string phaseID,
       const double siteDensity
       );
    ~Phase();
	/*
    void setPhaseName(const std::string &name); 
	void setphaseID(const std::string &id);
	void setSiteDensity(const double den); 
    */
	
    const std::string& getPhaseName() const;
    const std::string& getPhaseID() const;
    const double& getSiteDensity() const;   
	
	// Set and return the species in the phase
    void setSpecies(std::map<std::string, int> species_map);
	const std::map<std::string, int>& getSpecies() const;
	
    friend std::ostream& operator<<(std::ostream& output, const Phase &phase);   

	private:
	std::map<std::string, int> species_map_;
	
	
  };

  } //namespace IO
#endif /* PHASE_H_ */
