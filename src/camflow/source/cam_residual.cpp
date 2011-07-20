#include "cam_residual.h"

using namespace Camflow;

CamResidual::CamResidual
(
    CamAdmin& ca,
    CamConfiguration& config,
    CamControl& cc,
    CamGeometry& cg,
    CamProfile& cp,
    CamSoot& cs,
    Mechanism& mech
)
:
  admin_(ca),
  reporter_(new CamReporter()),
  reacGeom_(cg),
  sootMom_(cs),
  control_(cc),
  camMech_(&mech),
  camMixture_(new Sprog::Thermo::Mixture(mech.Species())),
  spv_(camMixture_->Species()),
  opPre(ca.getPressure()),
  mCord(cg.getnCells()), // Need to add 2 due to the addition of 'ghost' cells in flamelet::solve
  iMesh_s(1),
  cellBegin(0),
  iMesh_e(mCord-1),
  cellEnd(mCord),
  //nMoments(6),				// Should use getNumberMoments method to set this.
  nMoments(cs.getNumMoments()),
  nSpc(mech.SpeciesCount()),
  nVar(nSpc+1+nMoments),
  ptrT(nSpc),
  nEqn(nVar*mCord),
  solverID(cc.getSolver()),
  resSp(nSpc*mCord,0.0),
  resMom(nMoments*mCord,0.0),		// ank25: moments residual
  resT(mCord,0.0),
  s_mf(mCord,nSpc),
  s_Wdot(mCord,nSpc),
  s_H(mCord,nSpc),
  s_Diff(mCord,nSpc),
  m_T(mCord,0.0),
  m_rho(mCord,0.0),
  m_cp(mCord,0.0),
  m_mu(mCord,0.0),
  m_u(mCord,0.0),
  m_k(mCord,0.0),
  dz(cg.getGeometry()),
  avgMolWt(mCord,0.0),
  moments(mCord,nMoments),			// ank25: moments analogous to s_mf
  moments_dot(mCord,nMoments)       // ank25: moments rate analogous to s_wdot
{}

CamResidual::~CamResidual()
{
    //if (camMech_ != NULL) delete camMech_;
    //if (camMixture_ != NULL) delete camMixture_;
    //if (reporter_ != NULL) delete reporter_;
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
    //moments.resize(cellEnd, nMoments);
    //moments_dot.resize(cellEnd, nMoments);

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
        camMixture_->SetMassFracs(mf);
        camMixture_->SetTemperature(temperature);
        mwt = camMixture_->getAvgMolWt();
        dens = opPre*mwt/(R*temperature);
        avgMolWt.push_back(mwt);
        camMixture_->SetMassDensity(dens);
        camMech_->Reactions().GetMolarProdRates(*camMixture_,wdot);

        m_rho.push_back(dens);
        //store the diffusion coefficient
        temp = camMixture_->getMixtureDiffusionCoeff(opPre);
        if(mom) m_mu.push_back(camMixture_->getViscosity());
        //the following properties are needed only when the
        //energy equation is solved. They may not be needed
        //for jacobian avaluation as well.
        htemp.resize(nSpc,0.0);
        cptemp.resize(nSpc,0.0);
        if(thermo){
            //store the molar enthalpy (J/mol)
            htemp = camMixture_->getMolarEnthalpy();
            //molar specific heats
            cptemp = camMixture_->getMolarSpecificHeat();
            //store the the thermal conductivity (J/m-s-K)
            m_k.push_back(camMixture_->getThermalConductivity(opPre));
            //store the specific heat capacity (J/kg K)
            m_cp.push_back(camMixture_->getSpecificHeatCapacity());
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

                slope = ( s_ParticleBegin(i,l) - s_ParticleEnd(i,l) )/control_.getMaxTime();
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
            grad = (s_mf(i,l) - s_mf(i-1,l))/delta;
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
    }

}

void CamResidual::mergeMomentum(doublereal* vec){
    for(int i=0; i<cellEnd; i++)
        solvect[i*nVar+ptrG] = vec[i];
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
    indVar = reacGeom_.getAxpos();
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



