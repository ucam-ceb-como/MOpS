/*
 * transport.h
 *
 *  Created on: Jun 24, 2011
 *      Author: Laurence R. McGlashan
 */

#ifndef TRANSPORT_H_
#define TRANSPORT_H_

#include <iostream>

#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

namespace IO
{

    class Transport
    {
        public:

            Transport();

            ~Transport(){}

            void setMoleculeIndex(const int moleculeIndex);
            const int& getMoleculeIndex() const;

            void setCollisionDiameter(const double collisionDiameter);
            const double& getCollisionDiameter() const;

            void setDipoleMoment(const double dipoleMoment);
            const double& getDipoleMoment() const;

            void setPolarizability(const double polarizability);
            const double& getPolarizability() const;

            void setPotentialWellDepth(const double potentialWellDepth);
            const double& getPotentialWellDepth() const;

            void setRotRelaxationNumber(const double rotRelaxationNumber);
            const double& getRotRelaxationNumber() const;

            void setReducedDipoleMoment();
            void setReducedDipoleMoment(const double reducedDipoleMoment_);
            const double& getReducedDipoleMoment();

            friend std::ostream& operator<<(std::ostream& output, const Transport& transport);

        private:

            template<class Archive>
            void serialize(Archive & ar, const unsigned int /* file_version */)
            {
                ar & moleculeIndex_ & potentialWellDepth_ & collisionDiameter_
                   & dipoleMoment_ & polarizability_ & rotRelaxationNumber_
                   & reducedDipoleMoment_;
            }

            friend class boost::serialization::access;

            int moleculeIndex_;
            //! L-J potential well depth, (e/kB) *kB
            double potentialWellDepth_;
            //! L-J collision dia- convert to Angstroms.
            double collisionDiameter_;
            //! Dipole converting to C-m.
            double dipoleMoment_;
            double polarizability_;
            //! Rotational relaxation number, Zrot at 298K.
            double rotRelaxationNumber_;
            double reducedDipoleMoment_;

    };

} // namespace IO

#endif /* TRANSPORT_H_ */
