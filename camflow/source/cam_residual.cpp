
#include <cmath>
#include <map>
#include <vector>

#include "array.h"

#include "cam_params.h"
#include "cam_residual.h"
#include "cam_math.h"
using namespace Camflow;

void CamResidual::speciesResidual(const doublereal& time, doublereal* y, doublereal* f){


    /*
     *preperation for array offsets
     */
    int ptrSP, ptrSW;
    int ptrCE, ptrCP, ptrCW;

    /*
     *prepare flux terms
     */
    doublereal convection, diffusion, source;
    map<int, vector<doublereal> >::iterator p;
    

    for(int i= loopBegin; i< loopEnd; i++ ){
        ptrCE = nVar*(i+1) + ptrC;
        ptrCP = nVar * i + ptrC;
        ptrCW = nVar *(i-1) + ptrC;
        
        
        for(int l=0; l<nSpc; l++){
            ptrSP = nVar*i + l;
            ptrSW = nVar*(i-1) +l;
            /*
             *convection term
             */
            convection = (-y[ptrCP]/m_rho[i])*dydx(y[ptrSP],y[ptrSW],dz[i]);
            /*
             *diffusion term
             */            
            diffusion = dydx(s_jk(i,l),s_jk(i+1,l),dz[i])/m_rho[i];
            /*
             *reaction source
             */
            source = s_Wdot(i,l)*(*spv)[l]->MolWt()/m_rho[i];            

            f[ptrSP] = convection + diffusion + source ;
           
        }

        
    }
   
  
}

void CamResidual::energyResidual(const doublereal& time, doublereal* y, doublereal* f){


    if(admin->getEnergyModel() == admin->ADIABATIC){
        /*
         *array offsets
         */
        int ptrTP, ptrTW, ptrTE;
        doublereal sumHSource, flxk, sumFlx;;
        for(int i=loopBegin; i<loopEnd; i++){
            /*
             *evaluate the summation terms
             */
            sumHSource = 0.0;
            sumFlx = 0.0;
            for(int l=0; l<nSpc; l++){
                sumHSource += s_Wdot(i,l)*s_H(i,l);
                flxk = 0.5*(s_jk(i-1,l)+s_jk(i,l));
                sumFlx += flxk*s_cp(i,l)/(*spv)[l]->MolWt();
            }
            /*
             *prepate offsets
             */
            ptrTP = i*nVar + ptrT;
            ptrTE = (i+1)*nVar + ptrT;
            ptrTW = (i-1)*nVar + ptrT;

            /*
             *convection
             */
            doublereal convection =  -m_u[i]*dydx(y[ptrTP],y[ptrTW],dz[i]);
            /*
             *conduction
             */
            doublereal conduction = dydx(m_q[i+1],m_q[i],dz[i]);
            /*
             *flux transport
             */
            sumFlx *= dydx(y[ptrTP],y[ptrTW],dz[i]);

            /*
             *residual
             */
            f[ptrTP] = convection + (conduction - sumHSource - sumFlx)/(m_rho[i]*m_cp[i]);
            

        }
    }else{
        int ptrTP,ptrTW;
        for(int i=loopBegin; i< loopEnd; i++){
            ptrTP = nVar*i + ptrT;
            ptrTW = nVar*(i-1) + ptrT;
            f[ptrTP] = 0;
            //f[ptrTP] = y[ptrTW]-y[ptrTP];

        }

    }


}


void CamResidual::massFlowResidual(const doublereal& time, doublereal* y, doublereal* f){

    /*
     *prepare array offsets
     */
    int ptrCP, ptrCW, ptrCE;
    for(int i=loopBegin; i<loopEnd; i++){
        ptrCP = nVar*i + ptrC;
        ptrCE = nVar*(i+1) + ptrC;
        ptrCW = nVar*(i-1) + ptrC;

        f[ptrCP] = -(y[ptrCP]/m_rho[i])*(y[ptrCP]-y[ptrCW])/dz[i];
        //f[ptrCP] = y[ptrCW]-y[ptrCP];
        
    }
}


void CamResidual::saveMixtureProp(doublereal* y, bool thermo){
    //clear all the vectors and maps before saving
    //s_Diff.clear();
    //s_H.clear();
    //s_mf.clear();
    m_T.clear();
    m_cp.clear();
    m_k.clear();
    m_rho.clear();
    m_u.clear();
    

    s_H.resize(cellEnd,nSpc);
    s_cp.resize(cellEnd,nSpc);
    s_mf.resize(cellEnd,nSpc);
    s_Diff.resize(cellEnd,nSpc);
    s_Wdot.resize(cellEnd, nSpc);
    //first and the last cell are imaginary cells. This is done
    //in order to be able to calulate the inlet species composition
    //as they are hardly kept constant
    vector<doublereal> mf; //mas fraction
    vector<doublereal> temp, htemp, cptemp;  //diffusion coefficient
    doublereal temperature, dens;

    

    for(int i=cellBegin; i< cellEnd;i++){
        mf.clear();
        for(int l=0; l<nSpc; l++){
            mf.push_back(y[i*nVar+l]);
        }
        temperature = y[i*nVar+ptrT];
        //store the temperature
        m_T.push_back(temperature);
        camMixture->SetMassFracs(mf);
        camMixture->SetTemperature(temperature);
        dens = opPre*camMixture->getAvgMolWt()/(R*temperature);
        camMixture->SetMassDensity(dens);
        camMech->Reactions().GetMolarProdRates(*camMixture,wdot);
                                
        m_rho.push_back(dens);                
        m_u.push_back(y[i*nVar+ptrC]/dens);
        //store the diffusion coefficient
        temp = camMixture->getMixtureDiffusionCoeff(opPre);
        //the following properties are needed only when the
        //energy equation is solved. They may not be needed
        //for jacobian avaluation as well.
        htemp.resize(nSpc,0.0);
        cptemp.resize(nSpc,0.0);
        if(thermo){
            //store the molar enthalpy (J/mol)
            htemp = camMixture->getMolarEnthalpy();
            //molar specific heats
            cptemp = camMixture->getMolarSpecificHeat();
            //store the the thermal conductivity (J/m-s-K)
            m_k.push_back(camMixture->getThermalConductivity(opPre));
            //store the specific heat capacity (K/kg K)
            m_cp.push_back(camMixture->getSpecificHeatCapacity());
        }

        for(int l=0; l<nSpc; l++){
            s_Wdot(i,l)=wdot[l];
            s_Diff(i,l)=temp[l];
            s_mf(i,l) = mf[l];
            s_H(i,l) = htemp[l];
            s_cp(i,l) = cptemp[l];
        }


    }
}



//
//calculate the species diffusion fluxes in kg/m2s
void CamResidual::updateDiffusionFluxes(){

    s_jk.resize(cellEnd+1,nSpc);

    
    vector<doublereal> flx;
    doublereal avgD, grad, avgRho;
    doublereal delta;


    //preperation for flux correction
    doublereal jCorr;

    for(int i=loopBegin; i<loopEnd; i++){
    
        avgRho = (m_rho[i]+m_rho[i-1])/2.0;
        delta = (dz[i]+dz[i-1])/2.0;
        flx.resize(nSpc,0.0);

        jCorr = 0.0;
        
        for(int l=0; l<nSpc; l++){            
            grad = dydx(s_mf(i,l),s_mf(i-1,l),delta);            
            avgD = (s_Diff(i-1,l)+s_Diff(i,l))/2.0;
            flx[l] = -avgD*avgRho*grad;
            //preparing for flux correction
            jCorr += flx[l];
            
        }
        
        //cantera copy
        for(int l=0; l<nSpc; l++){
            flx[l] -= s_mf(i,l)*jCorr;
            s_jk(i,l) = flx[l];
        }

    }

}

void CamResidual::updateThermo(){
    /*
     *the conduction fluxes are stored for the
     *left face of each cell.
     */
    m_q.resize(cellEnd+1,0.0);
    doublereal delta;
    for (int i = loopBegin; i < loopEnd; i++) {
        delta = (dz[i]+dz[i-1])/2.0;
        doublereal value = 0.5*(m_k[i]+m_k[i-1])*(m_T[i]-m_T[i-1])/delta;
        m_q[i] = value;
        
    }
    
}

//return the number of species
const int& CamResidual::getNSpecies() const{
    return this->nSpc;
}
//return the number of variables
const int& CamResidual::getNVar() const{
    return this->nVar;
}
//return the total number of equations
const int& CamResidual::getNEqn() const{
    return this->nEqn;
}


