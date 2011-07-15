/*
 * thermo.cpp
 *
 *  Created on: Jun 24, 2011
 *      Author: gigadot
 */

#include "thermo.h"

IO::Thermo::Thermo(const std::string species_name)
:
species_name_(species_name),
note_("EMPTY"),
phase_("EMPTY"),
T_low_(-1),
T_common_(-1),
T_high_(-1),
al_(7, -1),
ah_(7, -1),
elements_map_() {
}

const std::string& IO::Thermo::getSpeciesName() const {
    return species_name_;
}

void IO::Thermo::setNote(std::string note) {
    note_ = note;
}

const std::string& IO::Thermo::getNote() const {
    return note_;
}

void IO::Thermo::setPhase(std::string phase) {
    phase_ = phase;
}

const std::string& IO::Thermo::getPhase() const {
    return phase_;
}

void IO::Thermo::setTLow(double T_low) {
    T_low_ = T_low;
}

const double& IO::Thermo::getTLow() const {
    return T_low_;
}

void IO::Thermo::setTCommon(double T_common) {
    T_common_ = T_common;
}

const double& IO::Thermo::getTCommon() const {
    return T_common_;
}

void IO::Thermo::setTHigh(double T_high) {
    T_high_ = T_high;
}

const double& IO::Thermo::getTHigh() const {
    return T_high_;
}

void IO::Thermo::setElements(std::map<std::string, int> elements_map) {
    elements_map_ = elements_map;
}

const std::map<std::string, int>& IO::Thermo::getElements() const {
    return elements_map_;
}

void IO::Thermo::setUpperTemperatureCoefficients(
        double ah1, double ah2, double ah3,
        double ah4, double ah5, double ah6,
        double ah7) {
    ah_.clear();
    ah_.push_back(ah1);
    ah_.push_back(ah2);
    ah_.push_back(ah3);
    ah_.push_back(ah4);
    ah_.push_back(ah5);
    ah_.push_back(ah6);
    ah_.push_back(ah7);
}

void IO::Thermo::setLowerTemperatureCoefficients(
        double al1, double al2, double al3,
        double al4, double al5,
        double al6, double al7) {
    al_.clear();
    al_.push_back(al1);
    al_.push_back(al2);
    al_.push_back(al3);
    al_.push_back(al4);
    al_.push_back(al5);
    al_.push_back(al6);
    al_.push_back(al7);
}

const std::vector<double>& IO::Thermo::getUpperTemperatureCoefficients() const{
    return ah_;
}

const std::vector<double>& IO::Thermo::getLowerTemperatureCoefficients() const{
    return al_;
}

namespace IO {

    std::ostream& operator<<(std::ostream& output, const Thermo& thermo) {
        output << "    Thermo Data:\n"
                << "    (\n"
                << "        Species  : " << thermo.species_name_ << "\n"
                << "        Elements : {";


        std::map<std::string, int>::const_iterator iter, final_iter;
        final_iter = thermo.elements_map_.end();
        --final_iter;
        for (iter = thermo.elements_map_.begin(); iter != thermo.elements_map_.end(); ++iter) {
            output << iter->first << ":" << iter->second;
            if (iter != final_iter) {
                output << ", ";
            }
        }
        output << "}" << "\n"
                << "        Note     : " << thermo.note_ << "\n"
                << "        Phase    : " << thermo.phase_ << "\n"
                << "        T_low    : " << thermo.T_low_ << "\n"
                << "        T_common : " << thermo.T_common_ << "\n"
                << "        T_high   : " << thermo.T_high_ << "\n"
                << "        AHs      : {" << thermo.ah_[0] << ", " << thermo.ah_[1] << ", " << thermo.ah_[2] << ", " << thermo.ah_[3] << ", " << thermo.ah_[4] << ", " << thermo.ah_[5] << ", " << thermo.ah_[6] << "}" << "\n"
                << "        ALs      : {" << thermo.al_[0] << ", " << thermo.al_[1] << ", " << thermo.al_[2] << ", " << thermo.al_[3] << ", " << thermo.al_[4] << ", " << thermo.al_[5] << ", " << thermo.al_[6] << "}" << "\n"
                << "    )";
        return output;
    }
}
