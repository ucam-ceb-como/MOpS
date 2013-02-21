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
    std::cout << "Destructing a network." << std::endl;
    for (MixtureMap::const_iterator i = mMixtures.begin();
            i != mMixtures.end(); ++i)
        delete i->second;

    // Delete reactors
    for (ReactorMap::const_iterator i = mReactors.begin();
            i != mReactors.end(); ++i)
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


/*!
 * Connects an inflow to a reactor using string keys
 *
 * @param reac      Key of the reactor receiving the inflow
 * @param inflow    Mixture or reactor key for inflow
 */
void ReactorNetwork::Connect(
        const std::string &reac,
        const std::string &inflow) {

    if (HasKey(mReactors, reac)) {
        if (HasKey(mReactors, inflow)) {
            // The inflow key is a reactor key
            std::cout << "Connecting reactor (" << reac << ") with reactor ("
                            << inflow << ")." << std::endl;

            // Create a new flowstream
            Mops::FlowStream fs = Mops::FlowStream(*mMech);
            // Connecting the inflow will also set the stream's mixture to
            // that of the inflow reactor
            fs.ConnectInflow(*(GetReactor(inflow)));
            fs.ConnectOutflow(*(GetReactor(reac)));

            mStreams.insert(std::make_pair(inflow+std::string("-")+reac, &fs));

            // Create a copy of the inflow and take ownership by the reactor
            // receiving it
            // fs:out to reac:in
            GetReactor(reac)->SetInflow(fs);

            // Now join up outflow reactor connections to the cloned inflow of
            // downstream reactor
            // inflow:out to fs:in
            GetReactor(inflow)->SetOutflow(*(GetReactor(reac)->Inflow()));

        } else if (HasKey(mMixtures, inflow)) {
            // The inflow key is a mixture key
            std::cout << "Connecting reactor (" << reac << ") with mixture ("
                                        << inflow << ")." << std::endl;

            // Create a new flowstream
            Mops::FlowStream fs = Mops::FlowStream(*mMech);
            // Set only the mixture of the stream, leave the inflow alone
            fs.SetConditions(*(GetMixture(inflow)));
            fs.ConnectOutflow(*(GetReactor(reac)));

            mStreams.insert(std::make_pair(inflow+std::string("-")+reac, &fs));

            // And the connection for "this" reactor
            // fs:out to reac:in
            GetReactor(reac)->SetInflow(fs);

        } else
            throw std::runtime_error(std::string("Error, couldn't find inflow key ") + inflow);
    } else
        throw std::runtime_error(std::string("Error, couldn't find inflow key ") + reac);
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
 */
void ReactorNetwork::GenerateRevPaths(Mops::PSR* r, ReactorPath &order) {
    order.push_back(r);
    if (r->Inflow() != NULL) {
        if (r->Inflow()->HasReacInflow())
            GenerateRevPaths(static_cast<PSR*>(r->Inflow()->Inflow()), order);
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
    if (end == NULL) end = mReactors.begin()->second;

    // Now, loop over the pointers find the inflow with now
    GenerateRevPaths(end, rev_order);

    // If the above hasn't linked all reactors, throw error
    if(rev_order.size() != mReactors.size())
        throw std::runtime_error("Couldn't link all reactors, try manual order definition");

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
        if (r->HasOutflow()) r->Outflow()->ConnectInflow(*r);
    }

    // Once the pointer structure has been reestablished, initialise the
    // birth and death processes
    for (ReactorNetwork::r_iter it=this->Begin(); it!=this->End(); ++it) {
        (*it)->InitialiseInflow();
        (*it)->InitialiseOutflow();
    }
}


void ReactorNetwork::CheckOutflow() {
    for (ReactorNetwork::r_iter it=this->Begin(); it!=this->End(); ++it) {
        if (!(*it)->HasOutflow()) {
            std::cout << "Attaching an outflow to "
                    << (*it)->GetName() << "." << std::endl;
            Mops::FlowStream* fs = new Mops::FlowStream(*mMech);
            (*it)->SetOutflow(*fs);

            mStreams.insert(std::make_pair((*it)->GetName()+std::string("-out"), fs));
            fs->ConnectInflow(*(*it));
        }
    }
}

ReactorNetwork::r_iter ReactorNetwork::Begin() {
    return mReactorPath.begin();
}

ReactorNetwork::r_iter ReactorNetwork::End() {
    return mReactorPath.end();
}


}
