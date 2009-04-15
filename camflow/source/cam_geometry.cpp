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
#include <vector>
#include <fstream>
#include <iostream>
#include "cam_geometry.h"
#include "comostrings.h"
using namespace std;
using namespace Camflow;
using namespace Strings;
void CamGeometry::setGridFile(string name){
    gridFile = name;

}


void CamGeometry::descretize(){
    vector<doublereal> grid;
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
    doublereal pos =0;
    for(int i=0; i<nCell; i++){
        pos += dz[i];
        axPos.push_back((pos - dz[i]/2.0));
        
    }
}

void CamGeometry::setLength(doublereal len){
    this->length = len;
}


doublereal CamGeometry::getLenth() const{
    if(length == 0)
        throw CamError("descretisation not invoked\n");
    else
        return this->length;
}

int CamGeometry::getnCells() const{
    if(nCell == 0)
        throw CamError("descretisation not invoked\n");
    else
        return this->nCell;
}

void CamGeometry::setDia(doublereal d){
    dia = d;   
}

doublereal CamGeometry::getDia() const{
    return dia;
}

doublereal CamGeometry::getArea() const{
    return (pi*dia*dia/4.0);
}

doublereal CamGeometry::getSurfArea() {
    if(length==0)descretize();
    return (pi*dia*length);
}

doublereal CamGeometry::getSurfAres_l() const{
    return (pi*dia);
}

vector<doublereal>& CamGeometry::getGeometry(){
    return dz;
}

vector<doublereal>& CamGeometry::getAxpos(){
    return axPos;
}

void CamGeometry::addZeroWidthCells(){
    vector<doublereal> temp = dz;
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
