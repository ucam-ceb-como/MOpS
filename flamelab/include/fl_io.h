#ifndef FL_IO_H
#define FL_IO_H
#include<string>
#include<vector>
#include "fl_reactor.h"
#include "fl_initial.h"
#include "fl_solver_control.h"
#include "fl_params.h"
#include "fl_premix.h"
#include "console_io.h"
#include "camxml.h"
#include "gpc_mech.h"
#include "gpc_mixture.h"
using namespace std;
namespace FlameLab{

	class Premix;//forward declaration

	class FlameLabIO{
		vector<string> monitor;
		vector<unsigned int> consoleMask;
		Console_IO flameLabConsole;
	public:
		FlameLabIO(){}
		~FlameLabIO(){}
		void readInput(const string &fileName, Reactor &reac, SolverControl &solver);
		void readGeometry(Reactor &reac, const CamXML::Element &node);
		void readOPConditions(Reactor &reac, const CamXML::Element &node);
		void readInlet(Reactor &reac, const CamXML::Element &node);
		void readNozzleConditions(Reactor &reac, InitialConditions &nozzle, const CamXML::Element &node);
		void readSolverControl(SolverControl &solver, const CamXML::Element &node);
		void readMonitor(const CamXML::Element &node);
		void prepareConsole(Sprog::Mechanism &mech, FlameLab::Premix &flame);
		void writeToConsole(Reactor &reac) const;
		
	};
}

#endif