/* 
 * File:   cam_profile.h
 * Author: vinod
 *
 * Copyright (C) 2008 Vinod M Janardhanan.
 *
 * File purpose:
 *  This class contains the implementation of initial profiles
 * Licence:
 *  This file is part of "Camflow".
 *
 *  Camflow is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
 *
 * Contact:
 *  Dr Markus Kraft
 *  Dept of Chemical Engineering
 *  University of Cambridge
 *  New Museum Site
 *  Pembroke Street
 *  Cambridge
 *  CB2 3RA
 *  UK
 *
 *  Email   :   mk306@cam.ac.uk
 *  Website :   http://como.cheng.cam.ac.uk
 *
 *Created on January 18, 2009, 10:43 AM
 */

#ifndef _CAM_PROFILE_H
#define	_CAM_PROFILE_H
#include "cam_params.h"
#include "cam_converter.h"
#include "cam_conc.h"
#include "gpc.h"
#include "array.h"
#include "cam_geometry.h"
#include "cam_boundary.h"
#include <string>
#include <map>
#include "linear_interpolator.hpp"

using namespace Sprog;

namespace Camflow{
    class CamProfile : public CamConc {

        std::vector<double> u_pos, u_temp, u_species_pos, u_frac;
        std::vector<std::string> u_species;
        std::map<std::string, double> list_prdt, list_intmd;
        std::vector<double> m_prdt, m_intmd;
        double mWidth, mCenter;
        int fracType;
        bool flag_loadFracs, flag_loadTemp;
        Array2D start;
        CamGeometry& geom;


    public:

        CamProfile(CamGeometry& cg);

        virtual ~CamProfile();

        /*
         *set the geometry object
         */
        void setGeometryObj(CamGeometry& cg);
        /*
         *set mixing length and mixing width
         */
        void setMixingCenter(double len);
        void setMixingWidth(double len);
        /*
         *populte the product and intermediate mass fracs
         */
        void populateProducts(Mechanism& mech);
        void populateIntermdts(Mechanism& mech);
        /*
         *set the start profile
         */
        void setStartProfile(CamBoundary& cb, Mechanism& mech);
        /*
         *set start profile given 2 inlets
         */
        void setStartprofile(CamBoundary &left, CamBoundary &right, Mechanism &mech);
        /*
         *set the gaussian for intermediates
         */
        void setGaussian(Mechanism &mech);
        /*
         *return the start profile
         */
        Array2D& getStartProfile();

        //set guess for products
        void setProductSpecies(std::map<std::string, double> spec);
        //set guess for intermediates
        void setIntermediateSpecies(std::map<std::string, double> spec);

        //return the initial species guess
        void getmassFracs(std::map<std::string,double>& spec, Mechanism& mech, std::vector<double>& frac);

        //set the user defined temp profile
        void setUserTemp(double pos, double temp);

        void setUserFrac(double pos, double temp, std::string species);

        const bool flagLoadFracs() const
        {
            return flag_loadFracs;
        };

        const bool flagLoadTemp() const
        {
            return flag_loadTemp;
        };

        //return the user defined temperature based on the interpolation
        double getUserDefTemp(const double& pos);
        double getUserDefFracs(const double& pos, const std::string species);
        //set temperature based on a gaussian
        void setGaussTempProfile(std::vector<double>& vTemp);

        //return the position vector
        std::vector<double>& getPosition();


    };




}



#endif	/* _CAM_PROFILE_H */

