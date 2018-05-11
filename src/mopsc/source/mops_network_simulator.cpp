 /*!
  * @file   mops_network_simulator.cpp
  * @author William Menz
  * @brief  Implementation of network simulator
  *
  *   Licence:
  *      mops is free software; you can redistribute it and/or
  *      modify it under the terms of the GNU Lesser General Public License
  *      as published by the Free Software Foundation; either version 2
  *      of the License, or (at your option) any later version.
  *
  *      This program is distributed in the hope that it will be useful,
  *      but WITHOUT ANY WARRANTY; without even the implied warranty of
  *      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  *      GNU Lesser General Public License for more details.
  *
  *      You should have received a copy of the GNU Lesser General Public
  *      License along with this program; if not, write to the Free Software
  *      Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA
  *      02111-1307, USA.
  *
  *   Contact:
  *      Prof Markus Kraft
  *      Dept of Chemical Engineering
  *      University of Cambridge
  *      New Museums Site
  *      Pembroke Street
  *      Cambridge
  *      CB2 3RA, UK
  *
  *      Email:       mk306@cam.ac.uk
  *      Website:     http://como.cheng.cam.ac.uk
  */

#include "mops_network_simulator.h"
#include "mops_solver_factory.h"
#include "loi_reduction.h"

#include <boost/functional/hash.hpp>
#include <boost/random/mersenne_twister.hpp>



//////////////////////////////////////////// aab64 ////////////////////////////////////////////
#include <iostream>
#include <string>
//////////////////////////////////////////// aab64 ////////////////////////////////////////////



namespace Mops {

/*!
 * Construct using a Mops::Simulator object
 *
 * @param sim   The Simulator to construct with
 * @param times The timesteps object to calculate over
 */
NetworkSimulator::NetworkSimulator(
        const Mops::Simulator &sim,
        const Mops::timevector &times)
: mRuns(sim.RunCount()),
  mFileBase("net"),
  mTimes(times),
  mSimInitial(NULL) {
    // Create a simulator copy
    mSimInitial = new Mops::Simulator(sim);
}

// Destructor
NetworkSimulator::~NetworkSimulator() {
    delete mSimInitial;

    // Loop over the simulator path to delete individual simulators
    // given by the second element of the pair
    for (NetworkSimulator::s_iter it=NetworkSimulator::Begin();
            it != NetworkSimulator::End(); ++it) {
        delete it->sim;
        delete it->sol;
    }
}

//! Get an iterator to the beginning of the simulator paths
NetworkSimulator::s_iter NetworkSimulator::Begin() {
    return mSimulators.begin();
}

//! Get an iterator to the end of the simulator paths
NetworkSimulator::s_iter NetworkSimulator::End() {
    return mSimulators.end();
}

/*!
 * Run the simulation.
 *
 * @param net       The network to run
 * @param solver    The Mops::Solver of choice
 * @param seed      A random seed
 */
void NetworkSimulator::Run(
        Mops::ReactorNetwork &net,
        const Mops::Solver &solver,
        std::size_t seed) {

    // Set up some storage
    NetworkSimulator::s_iter it;

    // Initialise the simulator if it is not done already
    if (mSimulators.size() < 1) Initialise(net, solver);

    // Write auxiliary file
    mSimInitial->writeAux(*(net.Mechanism()), *(&mTimes), solver);

    // Open the binary output streams for each file
    for (it=this->Begin(); it!=this->End(); ++it) {
        it->sim->openOutputFile();

        // Write the initial conditions
        it->sim->fileOutput(it->sim->m_output_step, it->sim->m_output_iter,
                *(it->reac), *(it->sol), it->sim);
    }

    // Loop over the number of runs
    for (unsigned int i(0u); i!=mRuns; ++i) {

        std::size_t iseed = seed;
        boost::hash_combine(iseed, i);
        boost::mt19937 rng(iseed);

        Mops::timevector::const_iterator iint;
        unsigned int istep(0u), global_step(0u);
        double dt, t2;

        // Fill the reactors with their initial mixtures
        net.ResetNetwork();

        // Initialise reactors and simulators
        t2 = mTimes[0].StartTime();
        for (it=this->Begin(); it!=this->End(); ++it) {
			// aab64 Initialise OMP RNGs
			it->reac->Mixture()->Setprng(iseed);

            // Set the start times
            it->reac->SetTime(t2);

            // Reset the ODE solver
            it->sol->Reset(*(it->reac));

            // Set timing parameters
            it->sim->m_cpu_start = std::clock();
            it->sim->m_runtime  = 0.0;

            // Set up some LOI data
            Mops::fvector z = Mops::fvector(net.Mechanism()->GasMech().SpeciesCount());
            if (it->sol->GetLOIStatus()) {
                it->sim->setupLOI(*(it->reac), *(it->sol));
            }


			// aab64 temporary 
			// Add headers to coagulation diagnostics file
			/*ofstream coagFile1, coagFile2;
			std::string coagfname;
			coagfname = "Coag-event-diagnostics-3micros.csv";
			coagFile1.open(coagfname.c_str());
			coagFile1 << "Time (s)" << " , " << "truek" << " , " << "majk" << " , " 
				<< "dc1_0" << " , " << "dc2_0" << " , " << "w1_0" << " , " << "w2_0" << " , "
				<< "dc1_1" << " , " << "dc2_1" << " , " << "w1_1" << " , " << "w2_1" << " , "
				<< "dca_0" << " , " << "dca_1" << "\n";
			coagFile1.close();
			coagfname = "Coag-event-diagnostics-3ms.csv";
			coagFile2.open(coagfname.c_str());
			coagFile2 << "Time (s)" << " , " << "truek" << " , " << "majk" << " , "
				<< "dc1_0" << " , " << "dc2_0" << " , " << "w1_0" << " , " << "w2_0" << " , "
				<< "dc1_1" << " , " << "dc2_1" << " , " << "w1_1" << " , " << "w2_1" << " , "
				<< "dca_0" << " , " << "dca_1" << "\n";
			coagFile2.close();*/

			// aab64 temporary 
			// Add headers to psc diagnostics file
			/*ofstream pscFile1, pscFile2;
			std::string pscfname;
			pscfname = "PSC-event-diagnostics-3micros.csv";
			pscFile1.open(pscfname.c_str());
			pscFile1 << "Time (s)" << " , " << "dc_0" << " , " << "dc_1" << " , " << "w_0" << " , " << "w_1" << "\n";
			pscFile1.close();
			pscfname = "PSC-event-diagnostics-3ms.csv";
			pscFile2.open(pscfname.c_str());
			pscFile2 << "Time (s)" << " , " << "dc_0" << " , " << "dc_1" << " , " << "w_0" << " , " << "w_1" << "\n";
			pscFile2.close();*/

			// aab64 temporary 
			// Add headers to outflow diagnostics file
			/*ofstream deathFile1;
			std::string deathfname;
			deathfname = "Outflow-event-diagnostics.csv";
			deathFile1.open(deathfname.c_str());
			deathFile1 << "Time (s)" << " , " << "Num" << " , " << "RTD" << " , " << "wt" << " , " << "age" << " , " << "mass" << "\n";
			deathFile1.close();*/

			/*ofstream pscFile1;
			std::string pscfname;
			pscfname = "sp_update_times.csv";
			pscFile1.open(pscfname.c_str());
			// t maj N Ninc Ntotal ip1flag ip1 ip2flag ip2 sp1_wt sp2_wt sp1_d sp2_d sp1_age sp2_age sp1_lut sp2_lut coag successful
			pscFile1 << "Time (s) , maj, N, Ninc, Ninc_total, ip1_flag , ip1 , ip2_flag , ip2 , sp1_wt , sp2_wt , sp1_dp , sp2_dp , sp1_age , sp2_age , sp1_lut , sp2_lut , coag successful? \n";
			pscFile1.close();*/

			//////////////////////////////////////////// aab64 ////////////////////////////////////////////
			if (it->sim->GetWriteDiagsStatus()) {
				/* Create partProc diagnostics csv file with pre/post split SV, #SPs, #events in split
				step including additions in LPDA, create gasConcFile with pre/post split concs 
			    Note that this is defintely not an elegant implentation and is only really intended 
			    to verify expected process behaviour in the TiO2 case */
				ofstream partProcFile, gasConcFile;
				unsigned int process_iter;
				std::vector<std::string> tmpPNames;
				it->reac->Mech()->ParticleMech().GetProcessNames(tmpPNames, 0);
				std::string rname (it->reac->GetName());
				std::string partfname, chemfname;
				partfname = "Part-split-diagnostics(" + rname + ").csv";
				chemfname = "Chem-split-diagnostics(" + rname + ").csv";

				// Add headers to partProc diagnostics file
				partProcFile.open(partfname.c_str());
				partProcFile << "Time (s)" << " , " << "Time out (s)" << " , " << "Step number (-)" << " , "
					<< "SV in (-)" << " , " << "SV out (-)" << " , "
					<< "SP in (-)" << " , " << "SP out (-)" << " , "
					<< "Total statistical weight pre-split (-)" << " , " << "Total statistical weight post-split (-)" << " , "
					<< "Average collision diameter pre-split (-)" << " , " << "Average collision diameter post-split (-)" << " , "
					<< "Incepting weight pre-split (-)" << " , " << "Incepting weight post-split (-)" << " , "
					<< "Incepting factor pre-split (-)" << " , " << "Incepting factor post-split (-)" << " , ";
				for (process_iter = 0; process_iter < tmpPNames.size() - 1; process_iter++) {
					partProcFile << tmpPNames[process_iter] << " , ";
				}
				partProcFile << "TransitionRegimeCoagulationTerms (kernel specific)";
				for (process_iter = tmpPNames.size(); process_iter < it->reac->Mech()->ParticleMech().GetTermCount() + 1;
					process_iter++) {
					partProcFile << " , ";
				}
				partProcFile << "FictitiousCoagulationTerms (kernel specific)";
				for (process_iter = tmpPNames.size(); process_iter < it->reac->Mech()->ParticleMech().GetTermCount() + 1;
					process_iter++) {
					partProcFile << " , ";
				}
				partProcFile << "Inflow events" << " , " << "Outflow events" << "\n";
				partProcFile.close();

				// Add headers to gasConc diagnostics file
				gasConcFile.open(chemfname.c_str());
				gasConcFile << "Time (s)" << " , " << "Time out (s)" << " , " << "Step number (-)" << " , ";
				for (process_iter = 0; process_iter < it->reac->Mech()->GasMech().Species().size() - 1; process_iter++) {
					gasConcFile << it->reac->Mech()->GasMech().Species(process_iter)->Name() << " pre-split (mol/m3)" << " , "
						<< it->reac->Mech()->GasMech().Species(process_iter)->Name() << " post-split (mol/m3)" << " , ";
				}
				gasConcFile << "TiO2 pre-split (mol/m3)" << " , " << "TiO2 post-split (mol/m3)" << " , "
					<< "Temperature pre-split (K)" << " , " << "Temperature post-split (K)" << "\n";
				gasConcFile.close();
			}
			//////////////////////////////////////////// aab64 ////////////////////////////////////////////
        }

        // Loop over time intervals
        for (iint=mTimes.begin(); iint!=mTimes.end(); ++iint) {
            // Get the step size for this interval.
            dt = (*iint).StepSize();

            // Set output parameters for this interval
            for (it=this->Begin(); it!=this->End(); ++it) {
                it->sim->m_output_step = max((int)iint->SplittingStepCount(), 0);
                it->sim->m_output_iter = max((int)it->sim->m_niter, 0);
            }

            // Loop over the steps in the interval
            for (istep=0; istep<iint->StepCount(); ++istep, ++global_step) {

                // Note incrementation of t2 here
                std::cout << "Stepping "
                        << t2 << " - " << (t2+dt) << "." << std::endl;
                t2 += dt;

                // Run the solver
                for (it=this->Begin(); it!=this->End(); ++it) {

                    it->sim->m_cpu_mark = std::clock();
                    



					//////////////////////////////////////////// aab64 ////////////////////////////////////////////
					/* If the solve function with diagnostics capacity replaces the original solve function, this
					if statement is no longer necessary */
					if (it->sim->GetWriteDiagsStatus()) {
						it->sol->Solve(*(it->reac), t2, iint->SplittingStepCount(),
							it->sim->m_niter, rng, &Mops::Simulator::fileOutput,
							(void*)(it->sim), it->sim->GetWriteDiagsStatus());
					}
					else {
						it->sol->Solve(*(it->reac), t2, iint->SplittingStepCount(),
							it->sim->m_niter, rng, &Mops::Simulator::fileOutput,
							(void*)(it->sim));
					}
					//////////////////////////////////////////// aab64 ////////////////////////////////////////////




                    std::cout << it->reac->GetName() << " done. " << std::endl;

                    // Do LOI calculation here
                    if (it->sol->GetLOIStatus()) {
                        it->sim->solveLOIJacobian(
                                *(it->reac),
                                *(it->sol),
                                istep,
                                t2);
                    }
                }

                // Set the runtime
                for (it=this->Begin(); it!=this->End(); ++it) {
                    it->sim->m_runtime +=
                            it->sim->calcDeltaCT(it->sim->m_cpu_mark);
                }

                std::cout << std::endl;

            } // (time steps)

            // Create save points and ensemble files
            for (it=this->Begin(); it!=this->End(); ++it) {
                it->sim->createSavePoint(*(it->reac), global_step, i);

                // Clear LOI data if applicable
                if (it->sol->GetLOIStatus())
                    it->reac->DestroyJac(
                            it->sim->m_loi_J,
                            it->reac->Mech()->GasMech().SpeciesCount()
                            );

                if (it->sim->m_write_ensemble_file)
                    it->sim->createEnsembleFile(*(it->reac), global_step, i);
            }
        } // (time intervals)

        std::cout << "\nFinished run " << (i+1) << " of " << mRuns << "." << std::endl;

        // Write reduced mechanism for LOI
        for (it=this->Begin(); it!=this->End(); ++it) {
            if (it->sol->GetLOIStatus()) {
                std::vector<std::string> rejects;
                rejects.clear();

                Mops::LOIReduction::RejectSpecies(
                        it->sim->m_loi_data,
                        it->sol->ReturnCompValue(),
                        it->reac->Mech(),
                        rejects,
                        it->sol->ReturnKeptSpecies());
                it->reac->Mech()->GasMech().WriteReducedMech(
                        it->sim->OutputFile() + std::string("-kept.inp"),
                        rejects);
            }
        }


    } // (n runs)


    // Close the binary output files
    for (it=this->Begin(); it!=this->End(); ++it) {
        it->sim->closeOutputFile();
    }
}

/*!
 * Create the data structure for each reactor 'node' in the network. Each
 * reactor is assigned to a Mops::Node containing the reactor, its own solver
 * and its own simulator.
 *
 * @param net       The network object to initialise a simulator for
 * @param solver    The solver to be used
 */
void NetworkSimulator::Initialise(
        Mops::ReactorNetwork &net,
        const Mops::Solver &solver) {
    // Set up some storage
    Mops::Simulator* sim;
    Mops::Solver* isol;

    // Loop over the reactor paths of the network to initialise
    // simulators and solvers for each reactor
    for (ReactorNetwork::r_iter it=net.Begin(); it!=net.End(); ++it) {
        // Create a simulator
        sim = NetworkSimulator::CreateSimulator(*it);
        sim->SetTimeVector(mTimes);

        // Create a solver and initialise it
        isol = solver.Clone();
        isol->Initialise(*(*it));

        // Turn of LOI stuff if no reactions
        if (net.Mechanism()->GasMech().ReactionCount() < 1)
            isol->SetLOIStatusFalse();

        // Now create a new node
        Mops::Node n;
        n.sim = sim;
        n.sol = isol;
        n.reac = *it;

        // Add it to our simulator paths
        mSimulators.push_back(n);
    }

    // And ensure each reactor has an outflow
    net.CheckOutflow();

    // Finally, normalise the birth/death process rates in case some streams
    // have multiple in/outflows.
    for (ReactorNetwork::r_iter it=net.Begin(); it!=net.End(); ++it) {
        (*it)->NormaliseIOProcessRates();
    }
}

/*!
 *
 * @param r     The PSR to create a simulator for
 * @return      The simulator copy
 */
Mops::Simulator* NetworkSimulator::CreateSimulator(const Mops::PSR* r) {

    // Check that the initial simulator has been specified
    if (mSimInitial == NULL)
        throw std::logic_error("Can't initialise a new simulator without a template!");

    // Now create one from the initial one, using copy constructor
    Mops::Simulator* sim = new Mops::Simulator(*mSimInitial);

    std::cout << "Creating a simulator for reactor " << r->GetName()
            << "." << std::endl;
    sim->SetOutputFile(mSimInitial->OutputFile() + std::string("(")
        + r->GetName() + std::string(")"));

    return sim;
}

/*!
 * Call to the 'sub simulators' to postprocess the network.
 */
void NetworkSimulator::PostProcess() {

    std::cout << "Postprocessing the network." << std::endl;

    // Set up some storage
    Mops::Mechanism mech;
    Mops::timevector times;
    unsigned int ncput = 0;
    vector<string> cput_head;

    // Read the auxiliary file
    mSimInitial->readAux(mech, times, ncput, cput_head);

    // Calculate number of output points.
    unsigned int npoints(1u);
    for(Mops::timevector::const_iterator i=times.begin(); i!=times.end(); ++i) {
        npoints += i->StepCount();
    }

    // Loop over the simulators
    for (NetworkSimulator::s_iter it=this->Begin(); it!=this->End(); ++it) {
        it->sim->postProcessSimulation(mech, times, npoints, ncput, cput_head);
    }
}

}
