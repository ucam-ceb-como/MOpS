 /*!
  * @file   mops_reactor_network.cpp
  * @author William Menz
  * @brief  Implementation of reactor network class
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

#include "mops_reactor_network.h"

#include <utility>
#include <iostream>

namespace Mops {

/*!
 * Constructor, requires the mechanism object.
 *
 * @param mech  Pointer to the Mops::Mechanism
 */
ReactorNetwork::ReactorNetwork(const Mops::Mechanism &mech)
: mMech(NULL) {
    std::cout << "Creating a network." << std::endl;
    mMech = &mech;
}

/*!
 * Reactor network has ownership of the mixtures and the reactors, so
 * delete the objects associated with any of the remaining pointers.
 */
ReactorNetwork::~ReactorNetwork() {
    // Delete mixtures
    std::cout << "mops: Destructing a network." << std::endl;
    for (MixtureMap::const_iterator i = mMixtures.begin();
            i != mMixtures.end(); ++i)
        delete i->second;

    // Delete reactors
    for (ReactorMap::const_iterator i = mReactors.begin();
            i != mReactors.end(); ++i)
        delete i->second;

    // Delete flowstreams
    for (StreamMap::const_iterator i = mStreams.begin();
            i != mStreams.end(); ++i)
        delete i->second;
}

/*!
 * @param os    Output stream
 * @param net   Network object to print
 * @return      Output stream
 */
std::ostream& operator<<(
        std::ostream &os,
        const Mops::ReactorNetwork &net)
{
  os << "[ReactorNetwork]\n";
  os << "\n MIXTURES (" << net.mMixtures.size() << "):\n" ;
  for (Mops::ReactorNetwork::MixtureMap::const_iterator i = net.mMixtures.begin();
          i != net.mMixtures.end(); ++i)
      os << "   "  << i->first << ": " << *(i->second);
  os << "\n REACTORS (" << net.mReactors.size() << "):\n" ;
  for (Mops::ReactorNetwork::ReactorMap::const_iterator i = net.mReactors.begin();
          i != net.mReactors.end(); ++i)
      os << "   "  << i->first << ": " << *(i->second) << "\n";
  return os;
}

/*!
 * Adds a mixture to the network. The network takes ownership of this mixture
 * and responsibility for deletion.
 *
 * @param id    Name of mixture
 * @param mix   Mixture object to add
 */
void ReactorNetwork::AddMixture(
        const string &id,
        const Mops::Mixture &mix) {
    if (HasKey(mMixtures, id))
        throw std::runtime_error("Error, existing mixture with same ID in network.");
    mMixtures.insert(std::make_pair(id, &mix));
}

/*!
 * Adds a PSR to the network. The network takes ownership of this PSR
 * and responsibility for deletion.
 *
 * @param id    Name of PSR
 * @param mix   PSR object to add
 */
void ReactorNetwork::AddReactor(
        const string &id,
        Mops::PSR &reac) {
    if (HasKey(mReactors, id))
        throw std::runtime_error("Error, existing reactor with same ID in network.");
    mReactors.insert(std::make_pair(id, &reac));
    reac.SetName(id);
}

/*!
 * Fills a reactor using string keys
 *
 * @param rid   Key of the reactor to be filled
 * @param mid   Key of the mixture to fill the reactor with
 */
void ReactorNetwork::FillReactor(
        const string &rid,
        const string &mid) {
    if (HasKey(mReactors, rid) && HasKey(mMixtures, mid)) {
        std::cout << "Filling reactor (" << rid << ") with mixture ("
                << mid << ")." << std::endl;
        mReactors.at(rid)->Fill(*(mMixtures.at(mid)->Clone()));

        // Add an entry to mReactorInitialMixtures to record how to refill
        // the reactor
        mReactorInitialMixtures.insert(std::make_pair(rid, mMixtures.at(mid)));
    }
    else
        throw std::runtime_error("Error, unknown reactor or mixture IDs");
}


void ReactorNetwork::ConnectReacAndReac(
        Mops::FlowStream &fs,
        const std::string &source,
        const std::string &sink) {

    std::cout << "Connecting sink reactor (" << sink << ") with source reactor ("
                    << source << ")." << std::endl;

    // Connecting the inflow will also set the stream's mixture to
    // that of the inflow (sink) reactor
    fs.ConnectInflow(*(GetReactor(source)));
    fs.ConnectOutflow(*(GetReactor(sink)));

    // Create a copy of the inflow and take ownership by the reactor
    // receiving it
    // fs:out to sink:in
    GetReactor(sink)->SetInflow(fs);

    // Now join up outflow reactor connections to the cloned inflow of
    // downstream reactor
    // source:out to fs:in
    GetReactor(source)->SetOutflow(fs);

}

void ReactorNetwork::ConnectReacAndMix(
        Mops::FlowStream &fs,
        const std::string &reac,
        const std::string &mix) {

    // The inflow key is a mixture key
    std::cout << "Connecting reactor (" << reac << ") with mixture ("
                                << mix << ")." << std::endl;

    // Set only the mixture of the stream, leave the inflow alone
    fs.SetConditions(*(GetMixture(mix)));
    fs.ConnectOutflow(*(GetReactor(reac)));

    // And the connection for "this" reactor
    // fs:out to reac:in
    GetReactor(reac)->SetInflow(fs);

}

void ReactorNetwork::ConnectReacAndBlank(
        Mops::FlowStream &fs,
        const std::string &reac) {

    std::cout << "Connecting reactor (" << reac << ") to exhaust." << std::endl;

    //Set r:out to fs:in and fs:in to r:out
    GetReactor(reac)->SetOutflow(fs);
    fs.ConnectInflow(*(GetReactor(reac)));
}

void ReactorNetwork::Connect(
        const std::string &reac,
        const std::string &flow,
        bool in,
        double frac) {
    // Create a new flowstream
    Mops::FlowStream* fs = new Mops::FlowStream(*mMech);
    fs->SetFlowFraction(frac);

    // Check for blank stream
    bool blank(false);
    if (!in && flow.size() == 0) blank = true;

    if (HasKey(mReactors, reac)) {

        if (HasKey(mReactors, flow) || blank) {
            // We're connecting a reactor to a reactor
            if (in) {
                // We're connecting an inflow.
                ConnectReacAndReac(*fs, flow, reac);
            } else {
                // We're connecting an outflow
                if (blank) ConnectReacAndBlank(*fs, reac);
                else ConnectReacAndReac(*fs, reac, flow);
            }
        } else if (HasKey(mMixtures, flow)) {
            // We're connecting a reactor to a mixture
            if (in) {
                ConnectReacAndMix(*fs, reac, flow);
            } else {
                throw std::runtime_error("Error, can't connect a mixture as an outflow!");
            }
        } else {
            throw std::runtime_error(std::string("Error, couldn't find flow key ") + flow);
        }
    } else {
        throw std::runtime_error(std::string("Error, couldn't find reactor key ") + reac);
    }

    // Add it to the stream map?
    mStreams.insert(std::make_pair(flow+std::string("-")+reac, fs));
}

/*!
 * Specify the order in which to solve the network. This takes input from the
 * mops.inx file defining the network and is the preferable way of defining
 * the network.
 *
 * A map is used to provide slightly more flexibility in numbering of the
 * reactor orders than a vector. Maps always sort by their keys so simply
 * looping over the elements should do the trick.
 *
 * @param order     A map of (order, reactor id) pairs
 */
void ReactorNetwork::SpecifyPaths(
        const std::map<unsigned int, std::string> order) {
    // Loop over all specified paths
    for (std::map<unsigned int, std::string>::const_iterator it=order.begin();
            it!=order.end(); ++it) {
        mReactorPath.push_back(GetReactor(it->second));
    }
}


/*!
 * @param r         Pointer to the PSR object
 * @param order     Reverse path order vector
 * @param guard     Stops the system getting stuck in a loop
 */
void ReactorNetwork::GenerateRevPaths(
        Mops::PSR* r,
        ReactorPath &order,
        unsigned int &guard) {
    if (guard > 1000)
        throw std::runtime_error("Looped network detected, use manual reactor definition.");
    else guard++;

    order.push_back(r);
    Mops::FlowPtrVector iptrs = r->Inflows();
    for (Mops::FlowPtrVector::iterator it=iptrs.begin();
            it!=iptrs.end(); ++it) {
        if ((*it)->HasReacInflow())
            GenerateRevPaths(static_cast<PSR*>((*it)->Inflow()), order, guard);
    }
}


/*!
 * If the reactor network was correctly defined in the settings file, all
 * reactor 'nodes' should be linked to each other or mixtures via flowstream
 * objects. We should therefore be able to determine the order in which they
 * are to be solved by tracing these links to the beginning and the end.
 */
void ReactorNetwork::FindPaths() {
    // Initialise some storage
    ReactorPath rev_order;
    Mops::PSR* end(NULL);

    // Start by choosing a reactor with no out flowstream.. this is our
    // most likely point for finishing the network.
    for (ReactorMap::const_iterator it=mReactors.begin();
            it!=mReactors.end(); ++it) {
        if (!it->second->HasOutflow()) end = it->second;
    }
    // If we don't find one, choose the first
    // (some sort of recycle loop MUST be present in this case)
    if (end == NULL)
        throw std::runtime_error("Possibly a recycle loop in network, use manual definition.");

    // Now, loop over the pointers find the inflow with now
    unsigned guard(0);
    GenerateRevPaths(end, rev_order, guard);

    // If the above hasn't linked all reactors, throw error
    if(rev_order.size() != mReactors.size())
        throw std::runtime_error("Couldn't link all reactors, try manual order definition.");

    // Now create the reactor path from the order found
    for (unsigned int i(rev_order.size()); i-- > 0; ) {
        mReactorPath.push_back(rev_order[i]);
    }
}

/*!
 * Reset the reactors of the network back to their initial states. This needs
 * to be done carefully to ensure that the FlowStream pointer structure
 * remains intact.
 */
void ReactorNetwork::ResetNetwork() {

    // Loop over the reactor map
    for (MixtureMap::iterator it=mReactorInitialMixtures.begin();
            it!=mReactorInitialMixtures.end(); ++it) {

        Mops::PSR *r = mReactors.at(it->first);
        Mops::Mixture *m = it->second->Clone();

        // Fill the reactor
        r->Fill(*m, true);

        // And regenerate the pointer of the outflow to the reactor's mixture
        Mops::FlowPtrVector optrs = r->Outflows();
        for (Mops::FlowPtrVector::iterator it=optrs.begin();
                it!=optrs.end(); ++it) {
            (*it)->ConnectInflow(*r);
        }
    }

    // Once the pointer structure has been reestablished, initialise the
    // birth and death processes
    for (ReactorNetwork::r_iter it=this->Begin(); it!=this->End(); ++it) {
        (*it)->InitialiseInflows();
        (*it)->InitialiseOutflows();
    }

    // And re-normalise the flow rates in case of multiple inlets.
    for (ReactorNetwork::r_iter it=this->Begin(); it!=this->End(); ++it) {
        (*it)->NormaliseIOProcessRates();
    }
}

/*!
 * Cheks that the PSR has an outflow, if it doesn't, this method will
 * create one; connecting it to a 'blank' source.
 */
void ReactorNetwork::CheckOutflow() {
    for (ReactorNetwork::r_iter it=this->Begin(); it!=this->End(); ++it) {
        if (!(*it)->HasOutflow()) {
            Mops::FlowStream* fs = new Mops::FlowStream(*mMech);
            ConnectReacAndBlank(*fs, (*it)->GetName());
            mStreams.insert(std::make_pair((*it)->GetName()+std::string("-out"), fs));
        }
    }
}

ReactorNetwork::r_iter ReactorNetwork::Begin() {
    return mReactorPath.begin();
}

ReactorNetwork::r_iter ReactorNetwork::End() {
    return mReactorPath.end();
}

// Returns the number of reactors in the network
unsigned int ReactorNetwork::ReactorCount() const {
    return mReactors.size();
}

}
