 /*!
  * @file   mops_network_simulator.cpp
  * @author William Menz
  * @brief  Declaration of network simulator
  *
  *   About:
  *      The network simulator manages the recording of data from the network
  *      being solved. This is achieved by creating a 'Node' for each reactor,
  *      containing the reactor, its own solver and its simulator.
  *
  *      Data from the reactor is then recorded by each of these 'sub
  *      simulators' and postprocessed in the usual manner.
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
#ifndef MOPS_NETWORK_SIMULATOR_H_
#define MOPS_NETWORK_SIMULATOR_H_

#include "mops_simulator.h"
#include "mops_reactor_network.h"
#include "mops_solver_factory.h"

namespace Mops {

struct Node {

    //! A pointer to this node's reactor
    Mops::PSR* reac;

    //! A pointer to this node's simulator
    Mops::Simulator* sim;

    //! A pointer to this node's solver
    Mops::Solver* sol;

    Node(): reac(NULL), sim(NULL), sol(NULL) {}

    Node(Mops::PSR& r, Mops::Simulator& s, Mops::Solver& sl):
        reac(&r), sim(&s), sol(&sl) {}
};

class NetworkSimulator {
public:
    //! Constructor
    NetworkSimulator(
            const Mops::Simulator &sim,
            const Mops::timevector &times);

    //! Default destructor
    ~NetworkSimulator();

    //! Run the simulation!
    void Run(
            Mops::ReactorNetwork &net,
            const Mops::Solver &solver,
            std::size_t seed);

    //! Initialise the simulator for a given network
    void Initialise(
            Mops::ReactorNetwork &net,
            const Mops::Solver &sol);

    //! Postprocess the binary outputs into CSVs
    void PostProcess();

private:
    typedef std::vector<Mops::Node> SimPath;

    typedef SimPath::iterator s_iter;

    //! Get an iterator to the beginning of the simulator paths
    s_iter Begin();

    //! Get an iterator to the end of the simulator paths
    s_iter End();

    //! Creates a simulator for a PSR
    Mops::Simulator* CreateSimulator(const Mops::PSR* r);

    //! Number of runs of the network
    unsigned int mRuns;

    //! Base filename
    const std::string mFileBase;

    //! The timesteps object
    const Mops::timevector mTimes;

    //! The default Mops::Simulator object
    Mops::Simulator* mSimInitial;

    //! A list of (reactors, simulators) in the order they should be solved
    SimPath mSimulators;
};

}

#endif /* MOPS_NETWORK_SIMULATOR_H_ */
