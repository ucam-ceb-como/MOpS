/*!
 * \file:   cam_read.cpp
 * \author: vinod
 *
 * Created on January 17, 2009, 5:38 PM
 *
 * \brief This class contains the implementation IO functions
 *
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

#include "cam_read.h"

using namespace Camflow;
using namespace Strings;

void CamRead::readInput(const std::string fileName,
                            CamControl& cc,
                            CamGeometry& cg,
                            CamConverter& convert,
                            CamAdmin& ca,
                            CamBoundary& cb,
                            CamProfile& cp,
                            CamConfiguration& config, CamSoot &cSoot){

    CamXML::Document doc;
    const CamXML::Element* root;

    ca.setInputFile(fileName);

    if(doc.Load(fileName) == 0){
        root = doc.Root();
        readGrid(cg,*root);
        readGeometry(cg,config,convert,*root);
        readProcessConditions(config, convert,ca,*root);
        readBoundary(ca,cb,convert,*root);
        readControl(cc,*root);
        readInitialGuess(ca,cp,convert,*root);
        readReport(ca,*root);
    }
}

//this function reads in the information concerning the model geometry
void CamRead::readGeometry(CamGeometry& cg,CamConfiguration& config,
                                CamConverter& convert,const CamXML::Element& node){
    CamXML::Element *reactorNode, *subnode;
    std::vector<CamXML::Element*> subnodes;
    //vector<CamXML::Element*>::iterator p;
    const CamXML::Attribute *attr;
    std::string attrValue;
    double factor = 1.0;

    reactorNode = node.GetFirstChild("reactor");
    if(reactorNode != NULL){
        attr    =   reactorNode->GetAttribute("model");
        if(attr != NULL){
            attrValue   =   attr->GetValue();
            if(!convertToCaps(attrValue).compare("PREMIX")) config.setConfiguration(config.PREMIX);
            if(!convertToCaps(attrValue).compare("PLUG")) config.setConfiguration(config.PLUG);
            if(!convertToCaps(attrValue).compare("COUNTERFLOW")) config.setConfiguration(config.COUNTERFLOW);
            if(!convertToCaps(attrValue).compare("STAGFLOW"))config.setConfiguration(config.STAGFLOW);
            if(!convertToCaps(attrValue).compare("BATCH_CV"))config.setConfiguration(config.BATCH_CV);
            if(!convertToCaps(attrValue).compare("FLAMELET"))config.setConfiguration(config.FLAMELET);
            if(!convertToCaps(attrValue).compare("UNSTEADYFLAMELET"))config.setConfiguration(config.UNSTEADYFLAMELET);
        }else{
            throw CamError("Model remains undefined\n");
        }


        subnode = reactorNode->GetFirstChild("diameter");
        if(subnode != NULL){
            attrValue = subnode->GetAttributeValue("unit");
            factor = convert.getConvertionFactor(attrValue);
            cg.setDia(cdble(subnode->Data())*factor);
        }else{
            cg.setDia(0.0);
        }

        subnode = reactorNode->GetFirstChild("length");
        if(subnode != NULL){
            attrValue = subnode->GetAttributeValue("unit");
            factor = convert.getConvertionFactor(attrValue);
            cg.setLength(cdble(subnode->Data())*factor);
        }else{
            cg.setLength(0.0);
        }

        // Read the grid in from the file, usually grid.inp.
        std::vector<double> grid;
        std::vector<double> dz;
        std::ifstream inf;
        inf.open(cg.getGridFileName().c_str(), std::ios::in);
        if (inf.good())
        {
            std::string position;

            while(!inf.eof())
            {
                getline(inf,position);

                if(! isEmpty(position))
                {
                    grid.push_back(cdble(position)*factor);
                }
            }
        }
        inf.close();

        // Calculate cell widths.
        int len = grid.size() - 1;
        cg.setLength(grid[len]);
        for (int i=0; i<len; i++)
        {
            dz.push_back(grid[i+1] - grid[i]);
        }
        // Call setGeometry in CamGeometry to assign length, axial positions
        // and cell widths.
        cg.setGeometry(dz);

        if
        (
            config.getConfiguration() == config.FLAMELET
         || config.getConfiguration() == config.UNSTEADYFLAMELET
         || config.getConfiguration() == config.STAGFLOW
         || config.getConfiguration() == config.COUNTERFLOW
        )
        {
            cg.addZeroWidthCells();
        }

    }else{
        throw CamError("Reactor definition missing\n");
    }

}

//this function reads in the process conditions
void CamRead::readProcessConditions
(
    CamConfiguration& config,
    CamConverter& convert,
    CamAdmin& ca,
    const CamXML::Element& node
)
{

    CamXML::Element *subnode, *opNode;
    const CamXML::Attribute *atr;
    std::string atrVal;

    opNode = node.GetFirstChild("op_condition");
    if(opNode == NULL){
        throw CamError("op_condition not defined\n");
    }else{
    //read temperature
        subnode = opNode->GetFirstChild("temperature");
        if(subnode != NULL)
            ca.setEnergyModel(subnode->Data());
        else
            throw CamError("temperature calculation not defined\n");

    //read the wall temperature
        subnode = opNode->GetFirstChild("twall");
        if(subnode!= NULL){
            std::string val = subnode->GetAttributeValue("unit");
            double factor = convert.getConvertionFactor(val);
            ca.setWallTemp(cdble(subnode->Data())+factor);
        }

    //read pressure
        subnode = opNode->GetFirstChild("pressure");
        if(subnode !=NULL){
            std::string unit = subnode->GetAttribute("unit")->GetValue();
            double fact = convert.getConvertionFactor(unit);
            double pre = cdble(subnode->Data())*fact;
            ca.setPressure(pre);
        }else{
            throw CamError("operating pressure not defined\n");
        }

        //read the ignition step for temperature
        subnode = opNode->GetFirstChild("step_ignite");
        if(subnode != NULL){
            ca.setIgnitionStep(cdble(subnode->Data()));
        }else{
            ca.setIgnitionStep(0.0);
        }

        subnode = opNode->GetFirstChild("radiation");

        if(subnode == NULL){
            throw CamError("op_condition::radiation is undefined\n");
        }else{
            atr = subnode->GetAttribute("activate");
            if(atr != NULL){
                atrVal = atr->GetValue();
                if(!convertToCaps(atrVal).compare("ON"))
                {
                    ca.setRadiationModel(true);
                }
                else
                {
                    ca.setRadiationModel(false);
                }
            }
        }

        if (config.getConfiguration() == config.FLAMELET)
        {
            subnode = opNode->GetFirstChild("flameletEquation");
            if(subnode == NULL){
                throw CamError("op_condition::flameletEquation is undefined\n");
            }else{
                std::string type = convertToCaps(subnode->Data());
                ca.setFlameletEquationType(type);
            }
        }

    }

}

//this function reads the boundary conditions
void CamRead::readBoundary(CamAdmin& ca,
                                CamBoundary& cb,
                                CamConverter& convert,
                                const CamXML::Element& node){

    CamXML::Element *subnode, *inletNode;

    inletNode = node.GetFirstChild("inlet");
    if(inletNode != NULL){
        subnode = inletNode->GetFirstChild("fuel");
        if(subnode != NULL){
            try{
                readNozzle(cb,convert,*subnode);
            }catch(CamError &ce){
                throw ce;
            }

            ca.setLeftBoundary(cb);

        }else{
            throw CamError("fuel definition missing in inlet element\n");
        }
        subnode = inletNode->GetFirstChild("oxidizer");
        if(subnode != NULL){
            try{
                readNozzle(cb,convert,*subnode);
            }catch(CamError &ce){
                throw ce;
            }
            ca.setRightBoundary(cb);
        }
    }else{
        throw CamError("inlet definition missing\n");
    }


}

//actual implementation of boundary conditions
void CamRead::readNozzle(CamBoundary& cb,
                                CamConverter convert,
                                const CamXML::Element& node){

    CamXML::Element *subnode;
    const CamXML::Attribute *atr;
    std::map<std::string, double> fracs;
    std::string atrVal;


    subnode = node.GetFirstChild("velocity");
    if(subnode != NULL){
        atr = subnode->GetAttribute("unit");
        atrVal = atr->GetValue();
        double fact = convert.getConvertionFactor(atrVal);
        double vel = cdble(subnode->Data())*fact;
        cb.setVelocity(vel);
    }else{
        cb.setVelocity(0.0);
    }

    subnode = node.GetFirstChild("temperature");
    if(subnode != NULL){
        atr = subnode->GetAttribute("unit");
        atrVal = atr->GetValue();
        double fact = convert.getConvertionFactor(atrVal);
        double temp = cdble(subnode->Data()) + fact;
        cb.setTemperature(temp);
    }else{
        throw CamError("Temperature at the inlet must be specified\n");
    }


    subnode = node.GetFirstChild("flowrate");
    if(subnode != NULL){
        atr = subnode->GetAttribute("unit");
        atrVal = atr->GetValue();
        double fact = convert.getConvertionFactor(atrVal);
        double flow = cdble(subnode->Data())*fact;
        cb.setFlowRate(flow);
    }else{
        cb.setFlowRate(0.0);
    }
    std::string member = "species";
    subnode = node.GetFirstChild("massfrac");
    if(subnode != NULL){
        cb.setFracType(cb.MASS);
        readFrac(member,fracs,*subnode);
        cb.setSpecies(fracs);
    }
    subnode = node.GetFirstChild("molefrac");
    if(subnode != NULL){
        cb.setFracType(cb.MOLE);
        readFrac(member,fracs,*subnode);
        cb.setSpecies(fracs);
    }
}

//function to read the mole or mass fractions for any element
void CamRead::readFrac(std::string& member, std::map<std::string,double>& fracs, const CamXML::Element& subnode){

    std::vector<CamXML::Element*> subnodes;
    std::vector<CamXML::Element*>::const_iterator p;
    const CamXML::Attribute *atr;
    std::string atrVal;
    std::string frac,finalSpecies;
    double sumfrac = 0.0;
    subnode.GetChildren(member,subnodes);
    for(p=subnodes.begin(); p<subnodes.end(); ++p){
        atr = (*p)->GetAttribute("name");
        if(atr!=NULL){
            atrVal = atr->GetValue();
            frac = (*p)->Data();
            if(!frac.compare("*")){
                finalSpecies = atrVal;
            }else{
                sumfrac += cdble(frac);
                fracs.insert(make_pair(atrVal,cdble(frac)));
                //Check to see that it's reading input file values correctly -- prints expected values. 
                //std::cout << "Species fraction:  " << frac  << std::endl;
            }
        }
    }
    if(! isEmpty(finalSpecies))
        fracs.insert(make_pair(finalSpecies,1-sumfrac));
}

//function to read solver control parameters
void CamRead::readControl(CamControl& cc, const CamXML::Element& node){
    CamXML::Element *solverNode, *subnode, *tols;
    const CamXML::Attribute *atr;
    std::string atrVal;
    solverNode = node.GetFirstChild("solver");
    if(solverNode != NULL){
        atr = solverNode->GetAttribute("mode");
        if(atr != NULL){
            atrVal = atr->GetValue();
            if(!convertToCaps(atrVal).compare("COUPLED")){
                cc.setSolutionMode(cc.COUPLED);
            }else if(!convertToCaps(atrVal).compare("SEGREGATED")){
                cc.setSolutionMode(cc.SEGREGATED);
            }else{
                throw CamError(" integration mode not specified in solver element\n");
            }
        }else{
            throw CamError(" integration mode not specified in solver element\n");
        }
        atr = solverNode->GetAttribute("solver");
        if(atr != NULL){
            atrVal = atr->GetValue();
            if(!convertToCaps(atrVal).compare("CVODE"))
                cc.setSolver(cc.CVODE);
            if(!convertToCaps(atrVal).compare("RADAU"))
                cc.setSolver(cc.RADAU);
            if(!convertToCaps(atrVal).compare("KINSOL"))
                cc.setSolver(cc.KINSOL);
            if(!convertToCaps(atrVal).compare("IDA"))
                cc.setSolver(cc.IDA);
            if(!convertToCaps(atrVal).compare("NEWTON"))
                cc.setSolver(cc.NEWTON);
        }else{
            throw CamError(" solver need to be specified\n");
        }

        atr = solverNode->GetAttribute("residual");
        if(atr != NULL){

            atrVal = atr->GetValue();
            if(!convertToCaps(atrVal).compare("ON"))
                cc.setResidualMonitor(true);
            else
                cc.setResidualMonitor(false);
        }

        subnode = solverNode->GetFirstChild("iterations");
        if(subnode != NULL)
            cc.setNumIterations(int(cdble(subnode->Data())));

        subnode = solverNode->GetFirstChild("iniStep");
        if(subnode != NULL)
            cc.setIniStep(cdble(subnode->Data()));
        subnode = solverNode->GetFirstChild("minStep");
        if(subnode!=NULL)
            cc.setMinStep(cdble(subnode->Data()));
        subnode = solverNode->GetFirstChild("maxStep");
        if(subnode!=NULL)
            cc.setMaxStep(cdble(subnode->Data()));
        subnode = solverNode->GetFirstChild("maxTime");
        if(subnode!=NULL)
            cc.setMaxTime(cdble(subnode->Data()));

        //read tolerences
        subnode = solverNode->GetFirstChild("tols");
        if(subnode != NULL){
            tols = subnode->GetFirstChild("resTol");
            if(tols != NULL){
                cc.setResTol(cdble(tols->Data()));
            }
            double atol, rtol;
            tols = subnode->GetFirstChild("species");
            if(tols != NULL) {
                readTol(*tols,atol,rtol);
                cc.setSpeciesAbsTol(atol);
                cc.setSpeciesRelTol(rtol);
            }

            tols = subnode->GetFirstChild("temperature");
            if(tols != NULL){
                readTol(*tols,atol,rtol);
                cc.setTempAbsTol(atol);;
                cc.setTempRelTol(rtol);
            }
            tols = subnode->GetFirstChild("flow");
            if(tols != NULL){
                readTol(*tols,atol,rtol);
                cc.setFlowAbsTol(atol);
                cc.setFlowRelTol(rtol);
            }

        }

    }else{
        throw CamError("solver definition missing\n");
    }

}
//funtion to read tolerence
void CamRead::readTol(const CamXML::Element& node, double& atol, double& rtol){

    CamXML::Element *subnode;
    atol = 0;
    rtol = 0;
    subnode = node.GetFirstChild("aTol");
    if(subnode!= NULL)
        atol = cdble(subnode->Data());
    subnode = node.GetFirstChild("rTol");
    if(subnode!= NULL)
        rtol = cdble(subnode->Data());
}

//function to read initial guess
void CamRead::readInitialGuess
(
    CamAdmin& ca,
    CamProfile& cp,
    CamConverter& convert,
    const CamXML::Element& node
)
{

    CamXML::Element *initialize, *subsubnode;
    std::vector<CamXML::Element*> subsubnodes;
    std::vector<CamXML::Element*>::const_iterator p,q,r;
    initialize = node.GetFirstChild("initialize");
    if(initialize != NULL)
    {
        //restart
        subsubnode = initialize->GetFirstChild("restart");
        if(subsubnode != NULL)
        {
            const CamXML::Attribute *file;
            file = subsubnode->GetAttribute("file");
            ca.setRestartType(subsubnode->Data());
            ca.setRestartFile(file->GetValue());
        }
        else
        {
            ca.setRestartType("NONE");
        }

        //mixing center
        subsubnode = initialize->GetFirstChild("mCenter");
        if(subsubnode != NULL)
        {
            const CamXML::Attribute *length;
            length = subsubnode->GetAttribute("unit");
            double convertL = convert.getConvertionFactor(length->GetValue());
            cp.setMixingCenter(cdble(subsubnode->Data())*convertL);
        }
        //mixing width
        subsubnode = initialize->GetFirstChild("mWidth");
        if(subsubnode != NULL)
        {
            const CamXML::Attribute *length;
            length = subsubnode->GetAttribute("unit");
            double convertL = convert.getConvertionFactor(length->GetValue());
            cp.setMixingWidth(cdble(subsubnode->Data())*convertL);
        }

        //temperature profile
        subsubnode = initialize->GetFirstChild("Tprofile");
        if(subsubnode != NULL)
        {
            const CamXML::Attribute *length, *temp;
            length = subsubnode->GetAttribute("unit_L");
            temp = subsubnode->GetAttribute("unit_T");
            double convertL = convert.getConvertionFactor(length->GetValue());
            double convertT = convert.getConvertionFactor(temp->GetValue());

            subsubnode->GetChildren("position",subsubnodes);

            for(p=subsubnodes.begin(); p<subsubnodes.end(); ++p){
                double pos = cdble((*p)->GetAttributeValue("x"))*convertL;
                double temp = cdble((*p)->Data())+convertT;
                cp.setUserTemp(pos,temp);
                std::cout << pos << " " << temp << std::endl;
            }
        }

        //intermediate and product species
        std::string mem1 = "product";
        std::string mem2 = "intrmdt";
        std::map<std::string,double> fracs;
        subsubnode = initialize->GetFirstChild("massfrac");
        if(subsubnode!=NULL)
        {
            cp.setFracType(cp.MASS);
            readFrac(mem1,fracs,*subsubnode);
            cp.setProductSpecies(fracs);
            readFrac(mem2,fracs,*subsubnode);
            cp.setIntermediateSpecies(fracs);
        }

        //mass frac profiles
        std::vector<CamXML::Element*> massfracs;
        initialize->GetChildren("massfracs",massfracs);
        subsubnode = initialize->GetFirstChild("massfracs");
        if(subsubnode != NULL)
        {
            for (q=massfracs.begin(); q<massfracs.end(); ++q)
            {
                if ((*q) != NULL)
                {
                    const CamXML::Attribute *species;
                    species = (*q)->GetAttribute("species");

                    (*q)->GetChildren("position",subsubnodes);

                    for (r=subsubnodes.begin(); r<subsubnodes.end(); ++r)
                    {
                        double pos = cdble((*r)->GetAttributeValue("x"));
                        double temp = cdble((*r)->Data());
                        cp.setUserFrac(pos,temp,species->GetValue());
                    }
                }
            }
        }

    }

}

void CamRead::readReport(CamAdmin& ca, const CamXML::Element& node){

    CamXML::Element *subnode;
    std::string atrVal;
    subnode = node.GetFirstChild("report");
    if(subnode!=NULL){
        atrVal = subnode->GetAttributeValue("species");
        if(!convertToCaps(atrVal).compare("MOLE")){
            ca.setSpeciesOut(ca.MOLE);
        }else{
            ca.setSpeciesOut(ca.MASS);
        }
    }else{
        throw CamError("Report information missing\n");
    }

}

void CamRead::readGrid(CamGeometry& cg, const CamXML::Element& node){
    CamXML::Element *subnode;
    subnode = node.GetFirstChild("grid");
    if(subnode!=NULL){
        cg.setGridFile(subnode->Data());
    }
}


