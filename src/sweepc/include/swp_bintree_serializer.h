/*!
 * @file    swp_bintree_serializer.h
 * @author  William J Menz
 * @brief   Template for serialising particle binary trees (e.g. silica)
 *
 *   Author(s):      William J Menz
 *   Project:        sweepc (population balance solver)
 *   Copyright (C) 2012 William J Menz
 *
 *   File purpose:
 *      Template for serialising particle binary trees (e.g. silica). Particle
 *      models must add the BinTreeSerializer class as a friend in order
 *      for the template class to have access to private data members.
 *
 *   Licence:
 *      This file is part of "sweepc".
 *
 *      sweepc is free software; you can redistribute it and/or
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


#ifndef SWEEP_BINTREE_SERIALIZER_H
#define SWEEP_BINTREE_SERIALIZER_H

//Includes
#include <iostream>
#include <stack>


namespace Sweep
{
namespace AggModels
{

template <class ParticleClass>
class BinTreeSerializer
{
public:
    //! Default constructor
    BinTreeSerializer() {};

    //! Default destructor
    ~BinTreeSerializer() {};

    /*!
     * @brief Serializes a particle binary tree structure
     *
     * @param[in,out]    out           Output binary stream
     * @param[in]        root          Pointer to the root node of the tree
     * @param[in,out]    duplicates    Addresses of PAHs that have already been serialised  
	 */
    void Serialize(std::ostream &out, const ParticleClass* root, void *duplicates) const
    {

        if (out.good()) {
            // Write version
            const unsigned int version = 1;
            out.write((char*)&version, sizeof(version));

            // Call the recursive Serializer!
            SerializeLoop(out, root, root, duplicates);

        }
    }

    /*!
     * @brief Helper function to write the primary state space and connectivity
     *
     * @param[in,out]    out           Output binary stream
     * @param[in]        root          Pointer to the root node of the tree
     * @param[in]        node          Pointer to the node to be investigated
     * @param[in,out]    duplicates    Addresses of PAHs that have already been serialised
     */
    void SerializeLoop(std::ostream &out, const ParticleClass* root, const ParticleClass* node, void *duplicates) const
    {
        // Serialize the state space first
        node->SerializePrimary(out, duplicates);

        // Does this primary have children?
        bool children(false);
        if (node->m_leftchild != NULL && node->m_rightchild != NULL) children = true;
        out.write((char*)&children, sizeof(children));

        // Now serialize the children
        if (children) SerializeLoop(out, root, node->m_leftchild, duplicates);
        if (children) SerializeLoop(out, root, node->m_rightchild, duplicates);

        // Now write the left/right particle connectivity
        int val(0);
        val = GetParticleIndex(root, node->m_leftparticle);
        out.write((char*)&val, sizeof(val));

        val = GetParticleIndex(root, node->m_rightparticle);
        out.write((char*)&val, sizeof(val));
    }

    /*!
     * @brief           Returns the index of a particle, given its address
     * @param node      Pointer to the current node being investigated
     * @param target    Pointer to the target node
     * @return          Value of the index
     */
    int GetParticleIndex(const ParticleClass* node, const ParticleClass* target) const
    {
        // Return zero if the target is nothing.
        if (target == NULL) return 0;

        // Otherwise, first particle has ID of 1
        int sum(1);
        bool status(false);
        GetParticleIndexLoop(node, target, &sum, &status);

        // Check the particle was found..
        if (! status) {
            throw std::runtime_error("Couldn't find particle index "
                    "(Sweep::AggModels::BinTreeSerializer::GetParticleIndex).");
        }

        return sum;
    }

    /*!
     * @brief           Recursive function to help find the node index
     *
     * This function starts at the root node and jumps to each node,
     * incrementing the counter until the desired function is found. 0 is
     * assumed to be NULL, and the root node has index 1.
     *
     * @param node      Pointer to the current node being investigated
     * @param target    Pointer to the target node
     * @param sum       Pointer to the counter
     * @param status    Flag indicating whether the node has been found
     */
    void GetParticleIndexLoop(
            const ParticleClass* node,
            const ParticleClass* target,
            int *sum,
            bool *status) const
    {
        if (target == node) {
            *status = true;
        } else {
            // Jump to next node
            (*sum)++;
            if (! *status) {
                if (node->m_leftchild != NULL)
                GetParticleIndexLoop(node->m_leftchild, target, sum, status);
            }
            if (! *status) {
                if (node->m_rightchild != NULL)
                GetParticleIndexLoop(node->m_rightchild, target, sum, status);
            }
        }
    }

    /*!
     * @brief Deserializes a binary tree structure
     *
     * @param[in,out]    in                Input binary stream
     * @param[in]        root              Pointer to the root node of the tree
     * @param[in]        model	           Particle model defining interpretation of particle data
     * @param[in,out]    pah_duplicates    Addresses of PAHs for use when reading primary particles
     *
     * @exception        runtime_error     Incorrect version number
     */
    void Deserialize(std::istream &in, ParticleClass *root, const Sweep::ParticleModel &model, void *duplicates)
    {
        if (in.good()) {
            // Read version
            unsigned int version = 0;
            in.read(reinterpret_cast<char*>(&version), sizeof(version));

            if (version == 1) {
                DeserializeLoop(in, root, root, model, duplicates);
            } else {
                throw std::runtime_error("Incorrect version number "
                        "(Sweep::AggModels::BinTreeSerializer::Deserialize).");
            }

        }
    }

    /*!
     * @brief Helper function to read the state space and connectivity
     *
     * @param[in,out]    in            Input binary stream
     * @param[in]        root          Pointer to the root node of the tree
     * @param[in]        node          Pointer to the current node being investigated
     * @param[in]        model	       Particle model defining interpretation of particle data    
     * @param[in,out]    duplicates    Addresses of PAHs for use when reading primary particles
     */
    void DeserializeLoop(std::istream &in, ParticleClass* root, ParticleClass* node, const Sweep::ParticleModel &model,
                         void *duplicates)
    {
        // Read in the state space first
        node->DeserializePrimary(in, model, duplicates);

        // Does the primary have children?
        bool children(false);
        in.read(reinterpret_cast<char*>(&children), sizeof(children));

        if (children) {
            // Create new children particles
            node->m_leftchild = new ParticleClass();
            node->m_rightchild = new ParticleClass();

            // Read in the children
            DeserializeLoop(in, root, node->m_leftchild, model, duplicates);
            node->m_leftchild->m_parent = node;
            DeserializeLoop(in, root, node->m_rightchild, model, duplicates);
            node->m_rightchild->m_parent = node;
        }
        /*
         * Now use the connectivity indices to reconnect particles
         * NOTE: we are able to do this because this section is not reached until
         * the rightmost child is created, thus *all* particles should have been
         * created before we attempt to establish pointer links between them.
         */
        int val(0);
        // Check the left particle
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        if (val != 0) node->m_leftparticle = FindParticleFromIndex(root, val);

        // Check the right particle
        in.read(reinterpret_cast<char*>(&val), sizeof(val));
        if (val != 0) node->m_rightparticle = FindParticleFromIndex(root, val);
    }

    /*!
     * @brief       Returns the address of a particle, given its index
     * @param node  Pointer to the current node being investigated
     * @param index Index of particle in tree
     * @return      Pointer to the desired node
     */
    ParticleClass* FindParticleFromIndex(ParticleClass* node, int index) {

        // Use the null pointer where 0 is given.
        if (index == 0) return NULL;

        // Otherwise, start at the top of the tree and work downwards...
        bool status(false);                 // Success flag for loop
        int sum(1);                         // Counter for tracking node number

        return FindParticleFromIndexLoop(node, &index, &sum, &status);
    }

    /*!
     * @brief           Recursive function to help find the pointer address
     * @param node      Pointer to the current node being investigated
     * @param index     Index of particle in tree
     * @param sum       Pointer to the counter
     * @param status    Flag indicating whether the node has been found
     * @return          Pointer to the desired node
     */
    ParticleClass* FindParticleFromIndexLoop(ParticleClass* node, int *index, int *sum, bool *status)
    {
        if (*index == *sum) {
            *status = true;
            return node;
        } else {
            // Jump to next node
            (*sum)++;
            ParticleClass* ans=NULL;
            if (! *status) {
                if (node->m_leftchild != NULL)
                    ans = FindParticleFromIndexLoop(node->m_leftchild, index, sum, status);
            }
            if (! *status) {
                if (node->m_rightchild != NULL)
                    ans = FindParticleFromIndexLoop(node->m_rightchild, index, sum, status);
            }
            return ans;
        }
    }

};

} // AggModels
} // Sweep

#endif
