 /*!
  * @file   mops_reactor_network.h
  * @author William Menz
  * @brief  Declaration of reactor network class
  *
  *   About:
  *      Reactors are networked in MOPS through this class. It stores a list of
  *      'static mixtures', which are used to initially fill reactors, supply
  *      inflow terms and regenerate reactors when reset.
  *
  *      The network is also responsible for setting-up flowstreams between
  *      reactors and ensuring their pointer structure is intact after resets.
  *      Of key importance in running the network is defining the order in
  *      which reactors are solved, which can be done automatically (not
  *      recommended!) or by specifying an order in the input file.
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


#ifndef MOPS_REACTOR_NETWORK_H_
#define MOPS_REACTOR_NETWORK_H_

#include "mops_psr.h"
#include "mops_reactor.h"
#include "mops_flow_stream.h"

namespace Mops {

class ReactorNetwork {
public:
    typedef std::vector<Mops::PSR*>::iterator r_iter;

    //! Constructor
    ReactorNetwork(const Mops::Mechanism &mech);

    //! Destructor
    ~ReactorNetwork();

    //! Overload of the << operator
    friend std::ostream& operator<<(
            std::ostream &os,
            const Mops::ReactorNetwork &net);

    //! Add a mixture to the list of static mixtures
    void AddMixture(const std::string &id, const Mops::Mixture &mix);

    //! Add a reactor to the reactor map
    void AddReactor(const std::string &id, Mops::PSR &reac);

    //! Fill a reactor with a mixture
    void FillReactor(const std::string &rid, const std::string &mid);

    //! Connect a flow or reactor to a reactor
    void Connect(const std::string &reac, const std::string &flow, bool in, double frac);

    //! Specify the order in which to solve the reactors
    void SpecifyPaths(const std::map<unsigned int, std::string> order);

    //! Reset the reactors to their intial state
    void ResetNetwork();

    //! Get the mechanism for the network.
    const Mops::Mechanism* Mechanism() const {return mMech;}

    //! Ensure the 'last' reactor has an outflow attached to it
    void CheckOutflow();

    //! If the reactor path is not manually specified, find it
    void FindPaths();

    //! Returns an iterator to the beginning of the reactor path
    r_iter Begin();

    //! Returns an iterator to the end of the reactor end
    r_iter End();

    //! Returns the number of reactors in the network
    unsigned int ReactorCount() const;

protected:
    template <typename T>
    static bool HasKey(const std::map<std::string, T> m, const std::string key) {
        for(typename std::map<std::string, T>::const_iterator it = m.begin(); it != m.end(); ++it)
            if (key.compare(it->first) == 0) return true;

        return false;
    }

private:
    // Hide default constructor to prevent initialisation without the mech
    ReactorNetwork();

    typedef std::map<std::string, const Mops::Mixture*> MixtureMap;

    typedef std::map<std::string, Mops::PSR*> ReactorMap;

    typedef std::map<std::string, Mops::FlowStream*> StreamMap;

    typedef std::vector<Mops::PSR*> ReactorPath;

    //! Get a pointer to a mixture given a key
    const Mops::Mixture* GetMixture(const std::string& id) {return mMixtures.at(id);}

    //! Get a pointer to a mixture given a key
    Mops::PSR* GetReactor(const std::string& id) {return mReactors.at(id);}

    //! Recursive helper function for SpecifyPaths()
    void GenerateRevPaths(Mops::PSR* r, ReactorPath &order, unsigned int &guard);

    //! Connect two reactors
    void ConnectReacAndReac(
            Mops::FlowStream &fs,
            const std::string &source,
            const std::string &sink);

    //! Connect two reactors
    void ConnectReacAndMix(
            Mops::FlowStream &fs,
            const std::string &reac,
            const std::string &mix);

    //! Connect a reactor to an empty outflow.
    void ConnectReacAndBlank(
            Mops::FlowStream &fs,
            const std::string &reac);

    //! The mechanism used for the reactor
    const Mops::Mechanism* mMech;

    //! The map of static mixtures
    MixtureMap mMixtures;

    //! The map of reactors
    ReactorMap mReactors;

    //! The map of reactor names to initial mixtures
    MixtureMap mReactorInitialMixtures;

    //! The map of flow streams
    StreamMap mStreams;

    //! A vector of pointers to reactors, defining the order of solution
    ReactorPath mReactorPath;
};

}

#endif /* MOPS_REACTOR_NETWORK_H_ */
