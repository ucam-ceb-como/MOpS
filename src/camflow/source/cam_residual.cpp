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
  config_(config),
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
  //nMoments(6),    // Should use getNumberMoments method to set this.
  nMoments(cs.getNumMoments()),
  nSpc(mech.SpeciesCount()),
  nVar(nSpc+1+nMoments),
  ptrT(nSpc),
  nEqn(nVar*mCord),
  solverID(cc.getSolver()),
  resSp(nSpc*mCord,0.0),
  resMom(nMoments*mCord,0.0),   // ank25: moments residual
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
  m_G(mCord,0.0),
  dz(cg.getGeometry()),
  avgMolWt(mCord,0.0),
  moments(mCord,nMoments),      // ank25: moments analogous to s_mf
  moments_dot(mCord,nMoments),       // ank25: moments rate analogous to s_wdot
  sootComponentRatesAllCells(mCord,nMoments*4),  // ank25: used to output soot rates
  avgSootDiamMaster(mCord,0.0),    // soot properties derived from moments.
  dispersionMaster(mCord,0.0),
  sootSurfaceAreaMaster(mCord,0.0),
  sootVolumeFractionMaster(mCord,0.0),
  wdotA4Master(mCord,0.0)
{
    camMixture_->SetViscosityModel(Sprog::iChapmanEnskog);
}


CamResidual::~CamResidual()
{
    //if (camMech_ != NULL) delete camMech_;
    //if (camMixture_ != NULL) delete camMixture_;
    //if (reporter_ != NULL) delet reporter_;
}


void CamResidual::saveMixtureProp
(
    const double time,
    double* y,
    bool thermo,
    bool mom
)
{
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
    std::vector<double> mf; //mas fraction
    std::vector<double> temp, htemp, cptemp;  //diffusion coefficient
    double temperature, dens, mwt;
    /*
     *check if particle sources are present
     */
    int npSource = std::max( s_ParticleBegin.size(), s_ParticleEnd.size());
    double slope, intersect;

    for (int i=cellBegin; i< cellEnd; i++)
    {
        mf.clear();
        for (int l=0; l<nSpc; l++)
        {
            mf.push_back(y[i*nVar+l]);
        }

        temperature = y[i*nVar+ptrT];

        if (temperature < 0 || temperature > 3500)
        {
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

        if (thermo)
        {
            //store the molar enthalpy (J/mol)
            htemp = camMixture_->getMolarEnthalpy();
            //molar specific heats
            cptemp = camMixture_->getMolarSpecificHeat();
            //store the the thermal conductivity (J/m-s-K)
            m_k.push_back(camMixture_->getThermalConductivity(opPre));
            //store the specific heat capacity (J/kg K)
            m_cp.push_back(camMixture_->getSpecificHeatCapacity());
        }


        for (int l=0; l<nSpc; l++)
        {
            /*
             *if there is a particle process source present
             *add that into the species source terms to
             *account for the particle processes
             */
            double pSource = 0;
            if (npSource != 0)
            {
                 /*
                 *the particle process source terms are assumed to be
                 *piece-wise linear. At a given time, the source term
                 *is evaluated by interpolating between the initial and
                 *final value
                 */

                slope = (s_ParticleBegin(i,l) - s_ParticleEnd(i,l))
                       /control_.getMaxTime();
                intersect = s_ParticleBegin(i,l);
                pSource = slope * time + intersect;
            }

            s_Wdot(i,l) = wdot[l] + pSource;
            s_Diff(i,l) = temp[l];
            s_mf(i,l) = mf[l];

            if (mf[l] > 1.1)
            {
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
void CamResidual::updateDiffusionFluxes()
{
    s_jk.resize(cellEnd+1,nSpc);

    std::vector<double> flx;
    double avgD, grad, avgRho;
    double delta;

    //preperation for flux correction
    double jCorr;

    for (int i=iMesh_s; i<iMesh_e; i++)
    {
        avgRho = (m_rho[i]+m_rho[i-1])/2.0;
        delta = (dz[i]+dz[i-1])/2.0;
        flx.resize(nSpc,0.0);

        jCorr = 0.0;

        for (int l=0; l<nSpc; l++)
        {
            grad = (s_mf(i,l) - s_mf(i-1,l))/delta;
            avgD = (s_Diff(i-1,l)+s_Diff(i,l))/2.0;
            //if((i-1)==0)avgD = s_Diff(i-1,l);
            flx[l] = -avgD*avgRho*grad;
            //preparing for flux correction
            jCorr += flx[l];

        }

        //correction
        for (int l=0; l<nSpc; l++)
        {
            flx[l] -= s_mf(i,l)*jCorr;
            s_jk(i,l) = flx[l];
        }
    }
}

void CamResidual::updateThermo()
{
    /*
     *the conduction fluxes are stored for the
     *left face of each cell.
     */
    m_q.resize(cellEnd + 1,0.0);
    double delta;

    for (int i = iMesh_s; i < iMesh_e; i++)
    {
        delta = (dz[i]+dz[i-1])/2.0;
        double value = 0.5*(m_k[i]+m_k[i-1])*(m_T[i]-m_T[i-1])/delta;
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

void CamResidual::extractEnergyVector(std::vector<double>& vec){
    vec.resize(cellEnd,0);
    for(int i=0; i<cellEnd; i++)
        vec[i] = solvect[i*nVar+ptrT];
}

void CamResidual::extractContinuity(std::vector<double>& vec){
    vec.resize(cellEnd,0);
    for(int i=0; i<cellEnd; i++)
        vec[i] = solvect[i*nVar+ptrF];
}

void CamResidual::extractSpeciesVector(std::vector<double>& vec){
    vec.resize(cellEnd*nSpc,0);
    for(int i=0; i<cellEnd; i++){
        for(int l=0; l<nSpc; l++)
            vec[i*nSpc+l] = solvect[i*nVar+l];
    }
}
void CamResidual::extractSootMoments(std::vector<double>& vec){
    vec.resize(cellEnd*nMoments,0);
    for(int i=0; i<cellEnd; i++){
        for(int l=0; l<nMoments; l++)
            vec[i*nMoments+l] = solvect[i*nVar+l+nSpc];
    }

}

void CamResidual::extractSpeciesAndEnergyVector(std::vector<double>& vec){
    vec.resize(cellEnd*(nSpc+1),0);
    for(int i=0; i<cellEnd; i++){
        for(int l=0; l<nSpc; l++)
            vec[i*nSpc+l] = solvect[i*nVar+l];
        vec[i*nSpc+ptrT] = solvect[i*nVar+ptrT];
    }
}


void CamResidual::extractMomentum(std::vector<double>& vec){
    vec.resize(cellEnd,0);
    for(int i=0; i<cellEnd; i++)
        vec[i] = solvect[i*nVar+ptrG];
}

void CamResidual::mergeEnergyVector(double* vec){
    for(int i=0; i<cellEnd; i++)
        solvect[i*nVar+ptrT] = vec[i];
}

void CamResidual::mergeContinuity(double* vec){
    for(int i=0; i<cellEnd; i++)
        solvect[i*nVar+ptrF] = vec[i];
}

void CamResidual::mergeSpeciesVector(double* vec)
{
    for (int i=0; i<cellEnd; i++)
    {
        for (int l=0; l<nSpc; l++)
        {
            solvect[i*nVar+l] = vec[i*nSpc+l];
        }
    }
}

void CamResidual::mergeSpeciesAndEnergyVector(double* vec)
{
    for (int i=0; i<cellEnd; i++)
    {
        for (int l=0; l<nSpc; l++)
        {
            solvect[i*nVar+l] = vec[i*nSpc+l];
        }

        solvect[i*nVar+ptrT] = vec[i*nSpc+ptrT];
    }
}


/*
 *soot moments are arranged after species and temperature
 */
void CamResidual::mergeSootMoments(double* vec){
    for(int i=0; i<cellEnd; i++){
        for(int l=0; l<nMoments; l++){
            solvect[i*nVar+ptrT+1+l] = vec[i*nMoments+l];
        }
    }
}

void CamResidual::mergeMomentum(double* vec)
{
    for (int i = 0; i < cellEnd; i++)
    {
        solvect[i*nVar+ptrG] = vec[i];
    }
}

int CamResidual::eval(double* y, double* ydot){
    std::cout << "Base class function nothing to do\n";
    return 0;
}

//flow field evaluation
void CamResidual::calcFlowField(const double& time, double* y){
    throw CamError("Base class function calcFlowField nothing to do\n");
}

// external funcation call to solve the reactor mode
void CamResidual::solve
(
    std::vector<Thermo::Mixture>& cstrs,
    const std::vector< std::vector<double> >& iniSource,
    const std::vector< std::vector<double> >& fnlSource,
    Mechanism& mech, CamControl& cc,
    CamAdmin& ca, CamGeometry& cg,
    CamProfile& cp)
{
    throw CamError("Base class function interface call nothing implemented\n");
}

/*
 *return the species mass fraction array
 */
void CamResidual::getSpeciesMassFracs(Array2D& mf){
    mf = s_mf;
}

/*
 *return the soot moments array
 */
void CamResidual::getMoments(Array2D& moments_){
    moments_ = moments;
}

/*
 *return the rate of change of the soot moments array
 *ank25 added for ELFM
 */
void CamResidual::getMomentsWdot(Array2D& momentsWdot_){
    momentsWdot_ = moments_dot;
}


/*
 *return the average molar weight of the mixture to the calling program
 */
void CamResidual::getAverageMolarWeight(std::vector<double>& avgMolWt_){
    avgMolWt_ = avgMolWt;
}
/*
 *return the density
 */
void CamResidual::getDensityVector(std::vector<double>& density){
    density = m_rho;
}
/*
 *return the velocity
 */
void CamResidual::getVelocity(std::vector<double>& vel){
    vel = m_u;
}
/*
 *return the viscosity
 */
void CamResidual::getViscosityVector(std::vector<double>& viscosity){
    viscosity = m_mu;
}
/*
 *return the temperature
 */
void CamResidual::getTemperatureVector(std::vector<double>& temp){
    temp = m_T;
}
/*
 *return the average soot diameter
 */
void CamResidual::getSootAverageDiameterVector(std::vector<double>& temp){
    temp = avgSootDiamMaster;
}
/*
 *return the soot dispersion
 */
void CamResidual::getSootDispersionVector(std::vector<double>& temp){
    temp = dispersionMaster;
}
/*
 *return the soot surface area
 */
void CamResidual::getSootSurfaceAreaVector(std::vector<double>& temp){
    temp = sootSurfaceAreaMaster;
}
/*
 *return the soot volume fraction
 */
void CamResidual::getSootVolumeFractionVector(std::vector<double>& temp){
    temp = sootVolumeFractionMaster;
}
/*
 * return the rate of A4 production
 */
void CamResidual::getWdotA4interface(std::vector<double>& temp){
    temp = wdotA4Master;
}
/*
 *  Return the specific heat
 */
void CamResidual::getSpecificHeat(std::vector<double>& spHeat){
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
void CamResidual::getThermalConductivity(std::vector<double>& lambda){
    lambda = m_k;
}
/*
 *return the independant variable
 */
void CamResidual::getIndepedantVar(std::vector<double>& indVar){
    indVar = reacGeom_.getAxpos();
}


void CamResidual::setParticleSource
(
    const std::vector<std::vector<double> >& initial,
    const std::vector<std::vector<double> >& final
)
{
    if (initial.size() != final.size())
    {
        throw
            CamError
            (
                "size mismatch : initial and final particle source terms\n"
            );
    }

    //resize array to allocate memory
    s_ParticleBegin.resize(cellEnd,nSpc);
    s_ParticleEnd.resize(cellEnd,nSpc);

    /*
     *The source terms passed in covers only the interior mech points.
     *At the inlet and outlet the  functions contain no chemical source
     *terms. Therfore, the sources at the inlet and outllet are set to 0
     */
    for (int i = iMesh_s; i < iMesh_e; i++)
    {
        std::vector<double> s_tBegin = initial[i-1];
        std::vector<double> s_tEnd = final[i-1];

        for (int l = 0; l < nSpc; l++)
        {
            s_ParticleBegin(i,l) = s_tBegin[l];
            s_ParticleEnd(i,l) = s_tEnd[l];
        }
    }
}
