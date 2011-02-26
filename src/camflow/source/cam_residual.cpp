
#include <cmath>
#include <map>
#include <vector>

#include "array.h"

#include "cam_params.h"
#include "cam_residual.h"
#include "cam_math.h"
#include "cam_read.h"
using namespace Camflow;


CamResidual::CamResidual()
:
admin(NULL),
reporter(NULL),
reacGeom(NULL),
sootMom(NULL),
control(NULL)
{}

CamResidual::~CamResidual()
{
    //if (camMixture != NULL) delete camMixture;
   // if (camMech != NULL) delete camMech;
   // if (control != NULL) delete control;
   // if (sootMom != NULL) delete sootMom;
   // if (reacGeom != NULL) delete reacGeom;
   // if (reporter != NULL) delete reporter;

}


void CamResidual::speciesResidual(const doublereal& time, doublereal* y, doublereal* f){


    /*
     *prepare flux terms
     */
    doublereal convection, diffusion, source;


    for(int i= iMesh_s; i< iMesh_e; i++ ){


        for(int l=0; l<nSpc; l++){
            /*
             *convection term
             */
            //convection = -m_u[i]*dydx(s_mf(i,l),s_mf(i-1,l),dz[i]);
            convection = -m_u[i]*dydx(i,s_mf(i+1,l),s_mf(i,l),s_mf(i-1,l),dz[i]);
            /*
             *diffusion term
             */
            diffusion = dydx(s_jk(i,l),s_jk(i+1,l),dz[i])/m_rho[i];

            source = s_Wdot(i,l)   *(*spv)[l]->MolWt()/m_rho[i];
            f[i*nSpc+l] = convection + diffusion + source ;

        }


    }


}



void CamResidual::energyResidual(const doublereal& time, doublereal* y, doublereal* f, bool flxTrans){


    if(admin->getEnergyModel() == admin->ADIABATIC){
        doublereal sumHSource, flxk, sumFlx;;
        for(int i=iMesh_s; i<iMesh_e; i++){
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
             *convection
             */
            //doublereal convection =  -m_u[i]*dydx(m_T[i],m_T[i-1],dz[i]);
            doublereal convection =  -m_u[i]*dydx(i,m_T[i+1],m_T[i],m_T[i-1],dz[i]);
            /*
             *conduction
             */
            doublereal conduction = dydx(m_q[i+1],m_q[i],dz[i]);
            /*
             *flux transport
             */
            sumFlx *= dydx(m_T[i],m_T[i-1],dz[i]);

            /*
             *residual
             */
            if(flxTrans){
                f[i] = convection + (conduction - sumHSource - sumFlx)/(m_rho[i]*m_cp[i]);
            }else{
                f[i] = convection + (conduction - sumHSource )/(m_rho[i]*m_cp[i]);
            }


        }
    }else{

        for(int i=iMesh_s; i< iMesh_e; i++){
            f[i] = 0;
            //f[ptrTP] = y[ptrTW]-y[ptrTP];

        }

    }


}


void CamResidual::massFlowResidual(const doublereal& time, doublereal* y, doublereal* f){


    for(int i=iMesh_s; i<iMesh_e; i++){

        f[i] = -m_u[i]*(m_flow[i]-m_flow[i-1])/dz[i];
        //f[ptrCP] = y[ptrCW]-y[ptrCP];

    }
}

doublereal CamResidual::getResidual() const {

    return 0;

}

void CamResidual::saveMixtureProp(const doublereal time,
                        doublereal* y, bool thermo, bool mom){
    //clear all the vectors and maps before saving
    m_T.clear();
    m_cp.clear();
    m_k.clear();
    m_rho.clear();
    m_mu.clear();
    avgMolWt.clear();
    //m_flow.clear();

    s_H.resize(cellEnd,nSpc);
    s_cp.resize(cellEnd,nSpc);
    s_mf.resize(cellEnd,nSpc);
    s_Diff.resize(cellEnd,nSpc);
    s_Wdot.resize(cellEnd, nSpc);
    //first and the last cell are imaginary cells. This is done
    //in order to be able to calulate the inlet species composition
    //as they are hardly kept constant
    std::vector<doublereal> mf; //mas fraction
    std::vector<doublereal> temp, htemp, cptemp;  //diffusion coefficient
    doublereal temperature, dens, mwt;
    /*
     *check if particle sources are present
     */
    int npSource = std::max( s_ParticleBegin.size(), s_ParticleEnd.size());
    doublereal slope, intersect;
    for(int i=cellBegin; i< cellEnd;i++){
        mf.clear();
        for(int l=0; l<nSpc; l++){
            mf.push_back(y[i*nVar+l]);
        }
        temperature = y[i*nVar+ptrT];
        if(temperature < 0 || temperature > 3500){
            std::cout << "Invalid temperature " << temperature << std::endl;
            throw CamError("Invalid temperature");
        }

        //store the temperature
        m_T.push_back(temperature);
        camMixture->SetMassFracs(mf);
        camMixture->SetTemperature(temperature);
        mwt = camMixture->getAvgMolWt();
        dens = opPre*mwt/(R*temperature);
        avgMolWt.push_back(mwt);
        camMixture->SetMassDensity(dens);
        camMech->Reactions().GetMolarProdRates(*camMixture,wdot);

        m_rho.push_back(dens);
        //store the diffusion coefficient
        temp = camMixture->getMixtureDiffusionCoeff(opPre);
        if(mom) m_mu.push_back(camMixture->getViscosity());
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
            //store the specific heat capacity (J/kg K)
            m_cp.push_back(camMixture->getSpecificHeatCapacity());
        }


        for(int l=0; l<nSpc; l++){
            /*
             *if there is a particle process source present
             *add that into the species source terms to
             *account for the particle processes
             */
            doublereal pSource =0;
            if(npSource != 0){
                 /*
                 *the particle process source terms are assumed to be
                 *piece-wise linear. At a given time, the source term
                 *is evaluated by interpolating between the initial and
                 *final value
                 */

                slope = ( s_ParticleBegin(i,l) - s_ParticleEnd(i,l) )/control->getMaxTime();
                intersect = s_ParticleBegin(i,l);
                pSource = slope * time + intersect;
            }
            s_Wdot(i,l)=wdot[l]+pSource;
            s_Diff(i,l)=temp[l];
            s_mf(i,l) = mf[l];
            if(mf[l] > 1.1) {
                std::cout << "Species " << "  " << mf[l] << std::endl;
                std::cout << "Temperature " << temperature << std::endl;
                std::cout << "Pressure " << opPre << std::endl;
                throw CamError("invalid mass frac\n");
            }
            s_H(i,l) = htemp[l];
            s_cp(i,l) = cptemp[l];
        }


    }

}

//set external scalar dissipation rate
void CamResidual::setExternalScalarDissipationRate(const doublereal sr){
    throw CamError("Base class function: set scalar dissipation rate\n");
}

//
//calculate the species diffusion fluxes in kg/m2s
void CamResidual::updateDiffusionFluxes(){

    s_jk.resize(cellEnd+1,nSpc);


    std::vector<doublereal> flx;
    doublereal avgD, grad, avgRho;
    doublereal delta;


    //preperation for flux correction
    doublereal jCorr;

    for(int i=iMesh_s; i<iMesh_e; i++){

        avgRho = (m_rho[i]+m_rho[i-1])/2.0;
        delta = (dz[i]+dz[i-1])/2.0;
        flx.resize(nSpc,0.0);

        jCorr = 0.0;

        for(int l=0; l<nSpc; l++){
            grad = dydx(s_mf(i,l),s_mf(i-1,l),delta);
            avgD = (s_Diff(i-1,l)+s_Diff(i,l))/2.0;
            //if((i-1)==0)avgD = s_Diff(i-1,l);
            flx[l] = -avgD*avgRho*grad;
            //preparing for flux correction
            jCorr += flx[l];

        }

       //correction
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
    for (int i = iMesh_s; i < iMesh_e; i++) {
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

void CamResidual::extractEnergyVector(std::vector<doublereal>& vec){
    vec.resize(cellEnd,0);
    for(int i=0; i<cellEnd; i++)
        vec[i] = solvect[i*nVar+ptrT];
}

void CamResidual::extractContinuity(std::vector<doublereal>& vec){
    vec.resize(cellEnd,0);
    for(int i=0; i<cellEnd; i++)
        vec[i] = solvect[i*nVar+ptrF];
}

void CamResidual::extractSpeciesVector(std::vector<doublereal>& vec){
    vec.resize(cellEnd*nSpc,0);
    for(int i=0; i<cellEnd; i++){
        for(int l=0; l<nSpc; l++)
            vec[i*nSpc+l] = solvect[i*nVar+l];
    }
}
void CamResidual::extractSootMoments(std::vector<doublereal>& vec){
    vec.resize(cellEnd*nMoments,0);
    for(int i=0; i<cellEnd; i++){
        for(int l=0; l<nMoments; l++)
            vec[i*nMoments+l] = solvect[i*nVar+l+nSpc];
    }

}

void CamResidual::extractMomentum(std::vector<doublereal>& vec){
    vec.resize(cellEnd,0);
    for(int i=0; i<cellEnd; i++)
        vec[i] = solvect[i*nVar+ptrG];
}

void CamResidual::mergeEnergyVector(doublereal* vec){
    for(int i=0; i<cellEnd; i++)
        solvect[i*nVar+ptrT] = vec[i];
}

void CamResidual::mergeContinuity(doublereal* vec){
    for(int i=0; i<cellEnd; i++)
        solvect[i*nVar+ptrF] = vec[i];
}

void CamResidual::mergeSpeciesVector(doublereal* vec){
    for(int i=0; i<cellEnd; i++){
        for(int l=0; l<nSpc; l++)
            solvect[i*nVar+l] = vec[i*nSpc+l];
    }
}
/*
 *soot moments are arranged right below the species vector
 */
void CamResidual::mergeSootMoments(doublereal* vec){
    for(int i=0; i<cellEnd; i++){
//        for(int l=nSpc; l<(nSpc+nMomemts); l++){
//            solvect[i*nVar+l] = vec[l-nSpc];
//            //std::cout << i << "  " << solvect[i*nVar+l] << std::endl;
//        }
        //int dd; cin >> dd;
        for(int l=0; l<nMoments; l++){
            solvect[i*nVar+l+nSpc] = vec[i*nMoments+l];
            //std::cout << i << "  " << solvect[i*nVar+l+nSpc] << std::endl;
        }
        //int dd; cin >> dd;
    }

}

void CamResidual::mergeMomentum(doublereal* vec){
    for(int i=0; i<cellEnd; i++)
        solvect[i*nVar+ptrG] = vec[i];
}

void CamResidual::massMatrix(doublereal** M){
    std::cout << "Base class function: nothing implemented\n";
}

int CamResidual::eval(doublereal* y, doublereal* ydot){
    std::cout << "Base class function nothing to do\n";
    return 0;
}

//flow field evaluation
void CamResidual::calcFlowField(const doublereal& time, doublereal* y){
    throw CamError("Base class function calcFlowField nothing to do\n");
}

// external funcation call to solve the reactor mode
void CamResidual::solve(std::vector<Thermo::Mixture>& cstrs,
       const  std::vector< std::vector<doublereal> >& iniSource,
        const std::vector< std::vector<doublereal> >& fnlSource,
        Mechanism& mech, CamControl& cc,
        CamAdmin& ca, CamGeometry& cg,
        CamProfile& cp){
    throw CamError("Base class function interface call nothing implemented\n");
}

/*
 *return the species mass fraction array
 */
void CamResidual::getSpeciesMassFracs(Array2D& mf){
    mf = s_mf;
}
/*
 *return the average molar weight of the mixture to the calling program
 */
void CamResidual::getAverageMolarWeight(std::vector<doublereal>& avgMolWt_){
    avgMolWt_ = avgMolWt;
}
/*
 *return the density
 */
void CamResidual::getDensityVector(std::vector<doublereal>& density){
    density = m_rho;
}
/*
 *return the velocity
 */
void CamResidual::getVelocity(std::vector<doublereal>& vel){
    vel = m_u;
}
/*
 *return the viscosity
 */
void CamResidual::getViscosityVector(std::vector<doublereal>& viscosity){
    viscosity = m_mu;
}
/*
 *return the temperature
 */
void CamResidual::getTemperatureVector(std::vector<doublereal>& temp){
    temp = m_T;
}

/*
 *  Return the specific heat
 */
void CamResidual::getSpecificHeat(std::vector<doublereal>& spHeat){
    spHeat = m_cp;
}

/*
 *  Return the diffusion coefficient
 */
void CamResidual::getDiffusionCoefficient(Array2D& dCoeff){
    dCoeff = s_Diff;
}

/*
 *  Return the thermal conductivity
 */
void CamResidual::getThermalConductivity(std::vector<doublereal>& lambda){
    lambda = m_k;
}
/*
 *return the independant variable
 */
void CamResidual::getIndepedantVar(std::vector<doublereal>& indVar){
    indVar = reacGeom->getAxpos();
}


void CamResidual::setParticleSource(const std::vector<std::vector<doublereal> >& initial,
        const std::vector<std::vector<doublereal> >& final){
    int n = initial.size();
    int m = final.size();
    if(m!=n) throw CamError("size mismatch : initial and final particle source terms\n");

    //resize array to allocate memory
    s_ParticleBegin.resize(cellEnd,nSpc);
    s_ParticleEnd.resize(cellEnd,nSpc);

    /*
     *The source terms passed in covers only the interior mech points.
     *At the inlet and outlet the  functions contain no chemical source
     *terms. Therfore, the sources at the inlet and outllet are set to 0
     */

    for(int i=iMesh_s; i<iMesh_e; i++){
        std::vector<doublereal> s_tBegin = initial[i-1];
        std::vector<doublereal> s_tEnd = final[i-1];
        for(int l=0; l<nSpc; l++){

            s_ParticleBegin(i,l) = s_tBegin[l];
            s_ParticleEnd(i,l) = s_tEnd[l];
        }
    }


}


