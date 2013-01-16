
/*
 * File:   cam_geometry.cpp
 * Author: vinod
 *
 * Created on January 17, 2009, 3:29 PM
 *  File purpose:
 *  This class contains the implementation of geometry details
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
 */

#include "cam_geometry.h"

using namespace std;
using namespace Camflow;
using namespace Strings;
using namespace Sprog;

void CamGeometry::setGridFile(string name){
    gridFile = name;

}

CamGeometry::CamGeometry()
:
length(0),
nCell(0),
a_slope(0.8),
curve(0.8),
minRange(0.01),
mPrune(-0.001)
{

}

/*void CamGeometry::discretize(){
    vector<double> grid;
    ifstream inf;
    inf.open(gridFile.c_str(),ios::in);
    if(inf.good()){
        std::string position;

        while(!inf.eof()){
            getline(inf,position);
            if(! isEmpty(position)){
                grid.push_back(cdble(position));
            }
        }
    }

    inf.close();
    int len = grid.size()-1;
    for(int i=0; i<len; i++){
        dz.push_back(grid[i+1]-grid[i]);
    }
    nCell = dz.size();
    //axial position info cell centers
    double pos =0;
    for(int i=0; i<nCell; i++){
        pos += dz[i];
        axPos.push_back((pos - dz[i]/2.0));
    }
    length = grid[len];
}*/

void CamGeometry::setGeometry(const vector<double>& dz_){
    dz = dz_;
    nCell=dz.size();
    //axial position info cell centers
    double pos =0;
    for(int i=0; i<nCell; i++){
        pos += dz[i];
        axPos.push_back((pos - dz[i]/2.0));
    }
}

void CamGeometry::setLength(double len){
    this->length = len;
}


double CamGeometry::getLenth() const{
    if(length == 0)
        throw CamError("Length not specified\n");
    else
        return this->length;
}

int CamGeometry::getnCells() const{
    return this->nCell;
}

void CamGeometry::setDia(double d){
    dia = d;
}

std::string CamGeometry::getGridFileName() const{
    return gridFile;
}

double CamGeometry::getDia() const{
    return dia;
}

double CamGeometry::getArea() const{
    return (PI*dia*dia/4.0);
}

double CamGeometry::getSurfArea() {
    //if(length==0)discretize();
    return (PI*dia*length);
} 

double CamGeometry::getSurfAres_l() const{
    return (PI*dia);
}

const vector<double>& CamGeometry::getGeometry() const {
    return dz;
}

const vector<double>& CamGeometry::getAxpos() const {
    return axPos;
}

void CamGeometry::addZeroWidthCells()
{
    vector<double> temp = dz;
    dz.clear();
    dz.push_back(1e-08);
    int len = temp.size();
    for (int i = 0; i < len; i++) {
        dz.push_back(temp[i]);
    }
    dz.push_back(1.e-08);

    //axpos reset
    temp = axPos;
    len = temp.size();
    axPos.clear();
    axPos.push_back(0.0);
    for(int i=0; i< len; i++){
        axPos.push_back(temp[i]);
    }
    axPos.push_back(length);

    //reset the number of cells
    nCell = dz.size();
}

static double mcPrec(){
    double x = 1.0;
    while(1.0+x != 1.0) x*= 0.5;
    return sqrt(x);

}

/*
 *grid refinement
 */
void CamGeometry::refine(double* y, const int nVar, const int nSpec, int ptrT){

    std::cout << "Refining Grid" << endl;
   /* std::cout << "Grid Size = " << axPos.size()
              << " ,nCell = " << nCell
              << " ,nVar = " << nVar
              << " ,nSpec = " << nSpec
              << endl;
    for(unsigned int i =0; i<axPos.size(); ++i)
    {
        cout << axPos[i] << endl;
    }*/

    vector<double> v,slope;
    double vmin, vmax, smin, smax;
    double maxAbsV, maxAbsS;
    double thresh;
    thresh = mcPrec();

    v.resize(nCell,0);
    slope.resize(nCell-1,0);
    /*
     *species
     */
    for(int l=0; l<nSpec; l++){


        for(int i=0; i<nCell; i++){
            v[i] = y[i*nVar+l];
        }
        for(int i=0; i<nCell-1;i++){
            slope[i] = (v[i+1]-v[i])/(axPos[i+1]-axPos[i]);
        }
        vmin = *min_element(v.begin(),v.end());
        vmax = *max_element(v.begin(),v.end());
        smin = *min_element(slope.begin(),slope.end());
        smax = *max_element(slope.begin(),slope.end());

        /*
         *max absolute of values and slope
         */
        maxAbsV = max(fabs(vmin),fabs(vmax));
        maxAbsS = max(fabs(smin),fabs(smax));

        /*
         *refine based on this component only if its range is
         *greater than a fraction of min range
         */
        if( (vmax - vmin) > maxAbsV*minRange ){

            //maximum allowable difference in value between
            //adjucent points

            double dmax = a_slope*(vmax-vmin)+thresh;
            for(int i=0; i<nCell-1; ++i){
                double r = fabs(v[i+1]-v[i])/dmax;
                if(r>1.0){
                    z_loc[i] = 1;
                    z_keep[i] = 1;
                }
                if( r>= mPrune){
                    z_keep[i] =1;
                    z_keep[i+1] =1;
                }else{
                    z_keep[i]=0;
                    z_keep[i]=-1;
                }
            }
        }

        /*
         *refine based on slope of the component
         */
        if( (smax-smin) > minRange*maxAbsS ){

            double dmax = curve*(smax-smin);

            for(int i=0; i<nCell-2; ++i){
                double r = fabs(slope[i+1]-slope[i])/(dmax + thresh/dz[i]);
                if( r>1.0){
                    z_loc[i] =1;
                    z_loc[i+1] =1;
                }
                if(r>= mPrune){
                    z_keep[i+1] =1;
                }else{
                    if(z_keep[i+1]==0){
                        z_keep[i+1]=-1;
                    }
                }
            }
        }

    }

    int nCells = 0;
    vector<double> temp = axPos;
    axPos.clear();
    for(int j=0; j<nCell-1; ++j){
        nCells++;
        axPos.push_back(temp[j]);
        if(z_loc.find(j) != z_loc.end()) {
            nCells++;
            axPos.push_back(0.5*(temp[j]+temp[j+1]));

        }
    }
    nCell = nCells;

    /*std::ofstream gridOut("gridNew.inp");

    for(int j=0; j<nCell-1; ++j){
        gridOut << axPos[j] << endl;
    }*/

}


