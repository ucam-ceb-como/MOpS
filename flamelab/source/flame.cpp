#include <string>
#include <iostream>
#include <map>
#include <vector>
#include "fl_io.h"
#include "fl_premix.h"
#include "fl_solver_control.h"
#include "fl_error_handler.h"
#include "fl_params.h"
#include "gpc.h"

using namespace std;
using namespace FlameLab;
using namespace Sprog;
int main(){
	string fChem("chem.inp");
	string fThermo("therm.dat");
	string fTrans("tran.dat");
	string fFlame("flame.xml");
	
	Reactor *flame = NULL;
	SolverManager *solver = NULL;
	SolverControl *solverControl = new SolverControl();
	Mechanism mech;
	
	
	try{
		// read chemkin input files
		IO::MechanismParser::ReadChemkin(fChem,mech,fThermo,fTrans);
		
		// create a reator object just for reading purpose		
		flame = new Reactor();
		// create flamelab object
		FlameLab::FlameLabIO *fio = new FlameLabIO();
		// read flamelab input file
		fio->readInput(fFlame,*flame,*solverControl);
		
		// create the mixture object
		Thermo::Mixture mix(mech.Species());
		//cout << boolalpha << (*flame==flame->PremixFlame) << endl;
		// if premix reactor create a premix object and assign to the base class reactor
		if( (*flame==flame->PremixFlame) || (*flame==flame->Plug) ){
			solver = new Premix(mech);
			//solver->initSolver(mech,mix,*solverControl,*flame);
			solver->solve(mech,*solverControl,*flame,*fio);
			//Premix p(mech);
			//premixFlame = &p;
			
			//premixFlame->initVariables(mech,flame->getFuelInletConditions(),mix);
		}			

		//solverControl->initSolver(*flame,mech,mix);
		
		//testing part for input parameters
		//cout << "reactor parameters" << endl;
		//cout << "aspect ratio " << flame->getAspectRatio() <<endl;
		//vector<FlameLab::real> dz = flame->getGeometry();
		//for(int i =0 ; i< flame->getnCells(); i++)
		//	cout << dz[i] << endl;						
		//cout << flame->getTemperature() << endl;
		//cout << flame->getnCells() << endl;
		//cout << flame->getLength() << endl;
		//cout << flame->getPressure() << endl;
		//cout << "fuel inlet conditions\n";
		//InitialConditions ic = flame->getFuelInletConditions();
		//cout << ic.getVelocity() << endl;
		//cout << ic.getTemperature() << endl;
		//map<std::string,real> inSpec = ic.getFraction();
		//map<std::string,real>::iterator p;
		//p = inSpec.begin();
		//while(p != inSpec.end()){
		//	cout << p->first << "  " << p->second << endl;
		//	p++;
		//}		
		//cout << "oxi inlet conditions\n";
		//ic = flame->getOxidizerInletConditions();
		//cout << ic.getVelocity() << endl;
		//cout << ic.getTemperature() << endl;
		//inSpec = ic.getFraction();		
		//p = inSpec.begin();
		//while(p != inSpec.end()){
		//	cout << p->first << "  " << p->second << endl;
		//	p++;
		//}
		//cout << "solver parameters \n";
		//cout << solverControl->getATol() << endl;
		//cout << solverControl->getRTol() << endl;
		//cout << solverControl->getIniStep() << endl;
		//cout << solverControl->getMaxStep() << endl;
		
				
	}catch(ErrorHandler rh){
		cout << rh.errorString <<endl;
	}

	delete flame;
	delete solverControl;

	cout << "Flame code finished successfully\n";	
	return 0;
}