/*
 * thermo.h
 *
 *  Created on: Jun 24, 2011
 *      Author: gigadot
 */

#ifndef THERMO_H_
#define THERMO_H_

#include <string>
#include <iostream>
#include <vector>
#include <map>

namespace IO {

    class Thermo {

    public:

        explicit Thermo(const std::string species_name);

        ~Thermo(){}

        const std::string& getSpeciesName() const;

        void setNote(std::string note);
        const std::string& getNote() const;

        void setPhase(std::string phase);
        const std::string& getPhase() const;

        void setTLow(double T_low);
        const double& getTLow() const;

        void setTCommon(double T_common);
        const double& getTCommon() const;

        void setTHigh(double T_high);
        const double& getTHigh() const;

        void setElements(std::map<std::string, int> elements_map);
        const std::map<std::string, int>& getElements() const;

        void setUpperTemperatureCoefficients(
                double ah1, double ah2, double ah3,
                double ah4, double ah5, double ah6,
                double ah7);
        const std::vector<double>& getUpperTemperatureCoefficients() const;

        void setLowerTemperatureCoefficients(
                double al1, double al2, double al3,
                double al4, double al5, double al6,
                double al7);
        const std::vector<double>& getLowerTemperatureCoefficients() const;

        friend std::ostream& operator<<(std::ostream& output, const Thermo& thermo);

    private:
        std::string species_name_;
        std::string note_;
        std::string phase_;
        double T_low_;
        double T_common_;
        double T_high_;
        std::vector<double> al_;
        std::vector<double> ah_;
        std::map<std::string, int> elements_map_;

    };

} // namespace IO

#endif /* THERMO_H_ */
