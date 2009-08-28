
#include "cam_residual.h"


#include "cam_configuration.h"


#include "cam_control.h"

#include "interface.h"
#include "flamelet.h"
using namespace Camflow;

Mechanism Interface::mech;
/*
 *this interface is for external programs to call any of the models
 *within camflow. The interface constructor does read the input files.
 *all the files required for the stand alone program is also required
 *here. Once the solution is complete, the interface object can be
 *queried to get the information on dependent variables
 */
Interface::Interface(){
    
    string fChem("chem.inp");
    string fThermo("therm.dat");
    string fTrans("tran.dat");
    string fCamFlow("camflow.xml");

    try{
        cm.readInput(fCamFlow,cc,cg,convert,ca,cb,cp,config,cSoot);
    }catch(CamError &ce){
        cout << ce.errorMessge;
    }
    

    //read mechanism, thermo and trasnport data
    IO::MechanismParser::ReadChemkin(fChem,mech,fThermo,fTrans);

    //get the number of species
    nSpecies = mech.SpeciesCount();
    //create the mixture 
    Thermo::Mixture mix(mech.Species());
    const SpeciesPtrVector *spv = mix.Species();
    /*
     *populate the species names
     */
    for(int l=0; l<nSpecies; l++){
        speciesNames.push_back((*spv)[l]->Name());
    }

    flmlt = NULL;
}
/*
 *return the number of species
 */
const int Interface::getNumberOfSpecies() const {
    return nSpecies;
}
/*
 *return the argument vector with the species names
 */
void Interface::getSpeciesNames(vector<string>& names){
    names.clear();
    names = speciesNames;
}
/*
 *this function is the interface for calling the flamelet code
 *for the implementation of interactive falmelets. The function
 *takes the scalar dissipation rate as input. The integration time
 *cam be passed as an optional argument. If its not passed, the
 *integration will be carried upto the time specified in the input
 *file/default value. If steady state is obtained before the specified
 *integration time, the program return with the converged solution.
 */
void Interface::flamelet(doublereal sdr, doublereal intTime, bool continuation){
    
    if(intTime!=0)cc.setMaxTime(intTime);
    if(flmlt == NULL ) flmlt = new FlameLet();
    flmlt->setExternalScalarDissipationRate(sdr);
    try{
        if(!continuation){
            flmlt->solve(cc,ca,cg,cp,mech,true);
        }else{
            flmlt->restart(cc);
        }
        /*
         *store the results for lookup
         */
        flmlt->getDensityVector(rhoVector);
        flmlt->getSpeciesMassFracs(spMassFracs);
        flmlt->getTemperatureVector(TVector);
        flmlt->getIndepedantVar(indVar);
        stMixtureFrac = flmlt->stoichiometricMixtureFraction();

    }catch(CamError &ce){
        cout << ce.errorMessge;
    }
}
/*
 *return the stoichiometric mixture fraction
 */
doublereal Interface::getStMixtureFrac(){
    return stMixtureFrac;
}
/*
 *
 *return the density
 */
doublereal Interface::getDensity(doublereal axpos){

    doublereal dens = getVariableAt(axpos,rhoVector);
    return dens;
}
/*
 *return the mass fractions
 */
doublereal Interface::getMassFrac(int spIndex, doublereal axpos){
    vector<doublereal> mf;
    int len = indVar.size();
    for(int i=0; i<len; i++){
        mf.push_back(spMassFracs(i,spIndex));
    }

    doublereal massfrac = getVariableAt(axpos,mf);
    return massfrac;
}
/*
 *return the temperature
 */
doublereal Interface::getTemperature(doublereal axpos){
    
    doublereal temp = getVariableAt(axpos,TVector);
    return temp;
}
/*
 *return the viscosity
 */
doublereal Interface::getViscosity(doublereal axpos){
    
    doublereal temp = getVariableAt(axpos,muVector);
    return temp;
}
/*
 *private function to do the interpolation of the solution variables
 */
doublereal Interface::getVariableAt(const doublereal& pos,
                                        vector<doublereal>& var){

    doublereal vu, vl, xu, xl, temp=0.0;
    int len = indVar.size();



    for(int i=0; i<len; i++){
        if(pos == indVar[i]) {
            temp= var[i];
            break;
        }else if( i>0 && (pos > indVar[i-1]) && (pos < indVar[i]) ){
            vu = var[i];
            xu = indVar[i];
            //cout << "location " << i << " pos " << pos << endl;
            vl = var[i-1];
            xl = indVar[i-1];
            //cout << tu << "  " << xu << endl;
            //cout << tl << "  " << xl << endl;
            doublereal slope = (vu-vl)/(xu-xl);
            doublereal intersect = vu- (slope*xu);
            //cout << "slope " << slope << endl;
            //cout << "intersect " << intersect << endl;
            temp= slope*pos + intersect;
            break;
        }
    }

    return temp;

}


