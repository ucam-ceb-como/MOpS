/*! \file binary_tree.hpp -*- C++ -*-
 * \brief Declaration of a binary tree for particle selection and rate summation.
 * \author  Robert I A Patterson

  Licence:

    This binary tree is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Robert I A Patterson (riap@cantab.net)
*/

/***************History************************************************
 * 27 Oct 2004 begun riap
 * 03 Mar 2010 adapted for Computational Modelling Group use riap
 * *******************************************************************/


#ifndef UTILS_BINARY_TREE_HPP
#define UTILS_BINARY_TREE_HPP


#include <list>
#include <vector>
#include <cassert>
#include<iostream>

namespace Utils {

//!Structure for the efficient evaluation of sums and inversion of probability distributions.
/*!*************Notes**************************************************
 * This class consists of a tree of nodes above a list of entries.  The
 * nodes form the tree; the entries should be thought of as a list, it
 * is this list which is the contents of this container and is exposed
 * by the STL like access functions.
 *
 * Nodes and entries are numbered separately, entries are 0 based,
 * nodes are not because this would break the system for calculating
 * the parents and children of nodes.
 * The numbering system which is described in c4e 26
 * is shown schematically below for a tree with 3 levels and a
 * total capacity of 8:
 *
 *                                      1
 *  tree                 2                             3
 *                4              5              6             7
 *         --------------------------------------------------------
 *  list       0     1        2    3         4     5       6      7
 *
 * ******************************************************************
 * \todo These notes do not format properly with doxygen.
 *
 * The point of the tree is that each node contains the sum of the
 * values (not necessarily scalar) of the leaves attached to its
 * two sub-trees.  There are some terminological issues here:
 * In the definition of this class 'weight' refers to these values
 * because they are a kind of relative probability.  Value has to
 * refer to the elements of the list at the bottom of the tree
 * since this is the STL convention for referring to the items
 * one puts into a tree.
 */
template<class Weight, class Selectee> class BinaryTree
{
  public:
    //Typedefs to help convert to templates later.  These are basic
    //parts of STL containers

    //!Sum of the values in the leaves of a (sub-)tree.
    /*! Binary tree provides structured sums of this type,
     * addition must be defined for this type. */
    typedef Weight weight_type;

    //!Constant copy of sum of the values in the leaves of a (sub-)tree.
    typedef const weight_type const_weight_type; //extra to STL

    class iterator;
    class const_iterator;

    //!Difference between two iterators.
    typedef int difference_type;

    //!Type for counting numbers of entries
    typedef unsigned size_type;


    //! Point back to a particle (outside this class) associated with a leaf.
    /*!When a leaf of the binary tree is selected according to a probability
     * distribution proportional to the weights it will be necessary to give
     * the caller a means of identifying which of the particles they inserted
     * into the tree has been chosen.*/
    typedef Selectee return_pointer_type;

    //!The user should view a binary tree as a container of these.
    /*!An entry in the binary tree is associated with an externally held
     * computational particle identified by an instance of
     * return_pointer_type.  As far as the caller is concerned an
     * entry in the tree (the value_type in STL terms) is just that -
     * a weight and a way of identifying which particle in the caller
     * is associated with the entry.
     * The weight_type is constant, because there is no way for the
     * tree to stay consistent if a weight is modified - it is rather
     * like the key of a pair associative STL container.
     */
    typedef std::pair<const weight_type, return_pointer_type> value_type;

    //!The user should view a binary tree as a container of these.
    /*!
     *See \ref value_type
     */
    typedef const value_type const_value_type; //extra to STL

    //! Create an empty tree with no space for particles
    BinaryTree();

    //!Create an empty tree with space for n particles.
    BinaryTree(size_type n);

    //Container operations, as for STL

    //! Points to the first element
    iterator begin();
    //! Points to the first element
    const_iterator begin() const;
    //! Point to one past the end
    iterator end();
    //! Point to one past the end
    const_iterator end() const;

    //iterator rbegin();
    //iterator rend();

    //! Number of elements in the tree
    size_type size() const;

    //! Maximum number of elements that the tree can extended to hold
    size_type max_size() const;

    //! Maximum number of elements that the tree can hold without rebuilding
    size_type capacity() const {return mCapacity;}

    //! True if nothing is stored in a tree instance
    bool empty() const;

    //! Exchange all contents
    void swap(BinaryTree& tree);

    //! Add element to end of container
    void push_back(const_value_type& x);

    //! Remove final element
    void pop_back();

    //! Change a stored value
    void replace(iterator i, const value_type& x);

    //! Remove the element referred to by i
    iterator erase(iterator i);

    //! Remove all elements, but do not alter the capacity
    void clear();

    //! Adjust internal structure to smallest possible size to hold given number of elements
    void resize(size_type new_size);

    //! Assign from a range
    template<typename ForwardIteratorType>
      void assign(const ForwardIteratorType it_begin,
                  const ForwardIteratorType it_end);

    //! Assign the same value many times
    void assign(size_type repeat_count, const value_type& value);

    //! Sum of all weights in the tree
    const weight_type& head() const;

    //! Select an entry according to the weights
    template<typename ScalarExtractorType> const_iterator select(double d, ScalarExtractorType f) const;

    //! Select an entry according to the weights
    template<typename ScalarExtractorType> iterator select(double d, ScalarExtractorType f);

    //!Write the nodes to stdout.
    void print_tree() const;

    //! Perform internal consistency checks
    bool isValid() const;

  private:
    //!Current number of levels in the tree.
    unsigned mNumLevels;

    //!Max total capacity of tree which is \f$2^{number_of_levels}\f$.
    size_type mCapacity;

    //!Index of first space in the list.
    /*!The index of the first space will also be the number of
     * entries currently stored. */
    size_type mFirstSpace;

    //!Type of the nodes in the tree.
    class node_type;

    //!Contains the tree structure
    std::vector<node_type> mNodes;

    //!Identify the particle in the calling code associated with each leaf.
    /*!The return pointer for each leaf of the tree allows the caller
     * to work out which of its particles a leaf (value) in the tree
     * is associated with.*/
    std::vector<return_pointer_type> mIdentifiers;

    //!Largest number of entries to permit is 2^31.
    static const size_type sMaxSize = 0x80000000;

    //! Pointer to left or right part of a node
    typedef weight_type node_type::*leaf_side;

    //!Identify a leaf of the tree.
    struct leaf_pos {
      //!Index of tree node leaf is under.
      size_type node;

      //!Member pointer to show which side of node the leaf is.
      leaf_side side_ptr;
    };

    //!Find position in the tree structure of value specified by an index.
    leaf_pos place_entry(const size_type index) const;
    //!Find position in the tree structure of value referenced by an iterator.
    leaf_pos place_entry(const_iterator i) const;

    //!Find the index of the entry under a tree leaf, inverse of BinaryTree::place_entry.
    size_type entry_leaf(const leaf_pos& leaf) const;

    //Binary summation tree methods

    //!Update the sums in the tree that depend on a node.
    void ascending_recalc(const size_type start_node);

    //! Recalculate the whole tree based on the values in the leaf nodes
    void recalculate();

};

//!Holds the values at a point in the tree.
/*!
 * A node consists of two entries: the value for the left subtree
 * and the value in the for the right subtree.  The values do not
 * have to be scalar.
 */
template<class Weight, class Selectee> class BinaryTree<Weight, Selectee>::node_type {
  public:
    //!Sum of values in the sub-tree to the left of this node.
    weight_type left;
    //!Sum of values in the sub-tree to the right of this node.
    weight_type right;

    //Equality tests, use of friend is superfluous.
    /*friend bool operator==(const node_type& a, const node_type& b);
    friend bool operator!=(const node_type& a, const node_type& b);*/
    //! Compare the data in two nodes
    /*! Note that comparing floating point values for equality is dangerous
     *\param[in]    b   Value to which to compare
     */
    bool isEqual(const node_type& b);
};

//! Iterator to access BinaryTree without changing its contents
template<class Weight, class Selectee> class BinaryTree<Weight, Selectee>::const_iterator {
  public:
    //default constructor
    const_iterator();
    //copy constructor
    //const_iterator(const const_iterator& i);

    //!Construct a const version of the iterator so data cannot be changed.
    const_iterator(const iterator& i);

    //!Iterator pointing to data element n
    const_iterator(const size_type n, const BinaryTree *ptree);
    //!Prefix increment
    const_iterator& operator++();
    //!Postfix increment
    const_iterator operator++(int);
    //!Prefix decrement
    const_iterator& operator--();
    //!Postfix decrement
    const_iterator operator--(int);

    //!Dereference to get a copy of the value only.
    typename BinaryTree<Weight, Selectee>::const_value_type operator*();

    //!Multiple increment
    const_iterator& operator+=(const difference_type n);

    //! Multiple increment a copy
    const_iterator operator+(const difference_type n);

    //! Distance between iterators
    difference_type operator-(const const_iterator& it) const;
    //! Distance between iterators
    difference_type operator-(const iterator& it) const;

    //! Check if two iterators point to the same entry in the same tree
    bool operator==(const const_iterator& b) const
    {
      return (container == b.container &&
	            element_index == b.element_index);
    }

    //!Needs access to extrct the node number,
    /*! \todo try to limit this later. */
    friend class BinaryTree;
  private:
    //!Pointer to the tree this instance of const_iterator refers to.
    const BinaryTree* container;

    //!Index of the element this object points to within a binary tree.
    size_type element_index;

    //!A sort of NULL
    static const size_type singular_value = 0xFFFFFFFF;

    //An end of sequence marker
//    const size_type sentry_value;
};

//! Iterator class for BinaryTree
template<class Weight, class Selectee> class BinaryTree<Weight, Selectee>::iterator {
  public:
    //default constructor
    iterator();
    //copy constructor
    //iterator(const iterator& i);

    //Iterator pointing to data element n of BinaryTree *ptree.
    iterator(const size_type n, BinaryTree<Weight, Selectee> *ptree);

    //!Prefix increment
    iterator& operator++();
    //!Postfix increment
    iterator operator++(int);
    //!Prefix decrement
    iterator& operator--();
    //!Postfix decrement
    iterator operator--(int);

    //!Dereference to get a copy of the value only.
    typename BinaryTree<Weight, Selectee>::value_type operator*();

    //! Multiple increment
    iterator& operator+=(const difference_type n);

    //! Multiple increment a copy
    iterator operator+(const difference_type n);

    //! Distance between iterators
    difference_type operator-(const const_iterator& it) const;
    //! Distance between iterators
    difference_type operator-(const iterator& it) const;

    /*! Equality means pointing to the same element in the same tree
     * \param[in]   b   Iterator to which to compare
     */
    bool operator==(const iterator& b) const
    {
      return (container == b.container &&
	            element_index == b.element_index);
    }

    //!Needs access to extrct the node number.
    /*! \todo try to limit this later. */
    friend class BinaryTree<Weight, Selectee>;

    //!Allow access to internals of this class when constructing const iterators.
    friend BinaryTree<Weight, Selectee>::const_iterator::const_iterator(const iterator& i);
  private:
    //index of the element this object points to
    size_type element_index;
    //a pointer to the tree this is an iterator for
    BinaryTree<Weight, Selectee>* container;

    //a sort of NULL
    static const size_type singular_value = 0xFFFFFFFF;

};

/****Section 1: Constructors, destructor and assignment operator******/

/*!
 * Completely empty tree
 *
 * @deprecated  Object should be initialised to something useful, otherwise it does not need creating
 */
template<class Weight, class Selectee> BinaryTree<Weight, Selectee>::BinaryTree()
: mNumLevels(0), mCapacity(0), mFirstSpace(0)
{
    assert(isValid());
}


/*!Set up an empty tree with mCapacity for upto 2^k entries
 * where \f$ 2^{\left(k-1\right)} < n <= 2^k \f$.
 * \param[in] n Number of particles tree must hold*/
template<class Weight, class Selectee> BinaryTree<Weight, Selectee>::BinaryTree(size_type n) {

  if(n > 0) {
      //find k to meet description above
      size_type k = 1;
      while((1u << k) < n) //2^k < n
        ++k;

      //tree will need k levels to handle n entries
      mNumLevels = k;
      //maximum mCapacity of the tree is 2^k >= n
      mCapacity = 1 << k; //2^k
  }
  else {
      // Zero capacity
      mNumLevels = 0;
      mCapacity = 0;
  }

  //tree is empty
  mFirstSpace = 0;

  //Allocate space for the basic tree stucture, element 0 will in fact
  //only be used as a cache for the total sums because of the numbering
  //scheme described above which makes 1 the root node
  mNodes.resize(mCapacity, node_type());

  //Allocate space for the list beneath the tree.
  mIdentifiers.reserve(mCapacity);

  assert(isValid());
}

/********* End of section 1 ******************************************/

/****** 2 BinaryTree node implementation *****************************/
/*! See if two nodes contain the same data, but note the dangers of
 *  comparing floating point numbers for equality.
 * \param[in]   b   Value to which to compare
 */
template<class Weight, class Selectee> bool
  BinaryTree<Weight, Selectee>::node_type::isEqual(const BinaryTree<Weight, Selectee>::node_type& b) {
  return (left == b.left && right == b.right);
}

/********* End of section 2 ******************************************/

/****** 3 BinaryTree iterator implementation *************************/

/*!Default constructor; should not point to anything meaningful.*/
template<class Weight, class Selectee> BinaryTree<Weight, Selectee>::iterator::iterator() {
  element_index = singular_value;
  container = NULL; //there is no container
}

/*!Point to a specified index in a tree
 * \param[in] n index of element in tree
 * \param[in] ptree pointer to tree.*/
template<class Weight, class Selectee>
  BinaryTree<Weight, Selectee>::iterator::iterator(size_type n, BinaryTree<Weight, Selectee> *ptree) {
  element_index = n;
  container = ptree;
}

//prefix increment
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::iterator&
  BinaryTree<Weight, Selectee>::iterator::operator++()
{
   ++element_index;
   return *this;
}

//postfix increment
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::iterator
  BinaryTree<Weight, Selectee>::iterator::operator++(int)
{
  iterator i(*this); //make a copy of ourself to return
  ++element_index; //increment ourself but not the copy
  return i;
}

//prefix decrement
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::iterator&
  BinaryTree<Weight, Selectee>::iterator::operator--()
{
   --element_index;
   return *this;
}

//postfix decrement
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::iterator
  BinaryTree<Weight, Selectee>::iterator::operator--(int)
{
  iterator i(*this); //make a copy of ourself to return
  --element_index; //decrement ourself but not the copy
  return i;
}

/*!
 * @return A copy of the value to which the iterator points
 *
 * @pre    The iterator must be dereferenceable (ie point to a valid element)
 */
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::value_type
  BinaryTree<Weight, Selectee>::iterator::operator*()
{
  //put together the return value in this variable
  value_type iterator_and_weight;

  //the second part of the pair is the iterator from the ensemble
  iterator_and_weight.second = container->mIdentifiers[element_index];

  //The position of the leaf node containing the weight of the entry
  //pointed to by the iterator.
  leaf_pos weight_pos = container->place_entry(element_index);

  //Now fill up the first part of the return value
  iterator_and_weight.first =
    container->mNodes[weight_pos.node].*(weight_pos.side_ptr);

  return iterator_and_weight;
}

//!Increment an iterator by n
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::iterator&
  BinaryTree<Weight, Selectee>::iterator::operator+=(difference_type n)
{
  element_index += n;
  return *this;
}

//! Return the negation of the equality comparison
template<class Weight, class Selectee>
  bool operator!=(const typename BinaryTree<Weight, Selectee>::iterator& a,
                  const typename BinaryTree<Weight, Selectee>::iterator& b)
{
  return !(a == b);
}

/*!
 * Move iterator n places forward
 *
 * @param[in]   n       Number of places to advance
 *
 * @return      Iterator advanced by n places
 */
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::iterator
  BinaryTree<Weight, Selectee>::iterator::operator+(const difference_type n) {
    // Construct a new instance to increment and return
    return (iterator(*this) += n);
}

/*!
 * Calculate the distance between two iterators
 *
 * @pre     Both iterators refer to the same container
 *
 * @param[in]   it       Iterator from which to calculate distance
 *
 * @return      n such that *this refers to the same entry as it + n
 */
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::difference_type
  BinaryTree<Weight, Selectee>::iterator::operator-(const const_iterator& it) const {
    return element_index - it.element_index;
}

/*
 * Calculate the distance between two iterators
 *
 * @pre     Both iterators refer to the same container
 *
 * @param[in]   it       Iterator from which to calculate distance
 *
 * @return      n such that *this refers to the same entry as it + n
 */
/*template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::difference_type
  BinaryTree<Weight, Selectee>::iterator::operator-(const iterator& it) const {
    return element_index - const_iterator(it);
}*/

/***** const_iterator ************************************************/

/*!Default constructor should not point to anything meaningful.*/
template<class Weight, class Selectee> BinaryTree<Weight, Selectee>::const_iterator::const_iterator() {
  element_index = singular_value;
  container = NULL; //there is no container
}

/*! Create a version of an iterator which cannot be used to modify the
 * value it points to.
 * \param[in] i Iterator to produce a (neutered) clone of. */
template<class Weight, class Selectee> BinaryTree<Weight, Selectee>::const_iterator::const_iterator(const iterator& i) {
  element_index = i.element_index;
  container = static_cast<const BinaryTree*>(i.container);
}

/*!Iterator pointing to value at aparticular index in a tree.
 * \param[in] n Index of value to point at.
 * \param[in] ptree Pointer to tree containing value. */
template<class Weight, class Selectee> BinaryTree<Weight, Selectee>::const_iterator::const_iterator(const size_type n, \
                                          const BinaryTree<Weight, Selectee> *ptree) {
  element_index = n;
  container = ptree;
}

//prefix increment
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::const_iterator&
  BinaryTree<Weight, Selectee>::const_iterator::operator++()
{
   ++element_index;
   return *this;
}

//postfix increment
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::const_iterator
  BinaryTree<Weight, Selectee>::const_iterator::operator++(int)
{
  const_iterator i(*this); //make a copy of ourself to return
  ++element_index; //increment ourself but not the copy
  return i;
}

//prefix decrement
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::const_iterator&
  BinaryTree<Weight, Selectee>::const_iterator::operator--()
{
   --element_index;
   return *this;
}

//postfix decrement
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::const_iterator
  BinaryTree<Weight, Selectee>::const_iterator::operator--(int)
{
  const_iterator i(*this); //make a copy of ourself to return
  --element_index; //decrement ourself but not the copy
  return i;
}

/*!Constant form of the value pointed to by the iterator operated on.*/
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::const_value_type
  BinaryTree<Weight, Selectee>::const_iterator::operator*() {
  //put together the return value in this variable
  value_type iterator_and_weight;

  //the second part of the pair is the iterator from the ensemble
  iterator_and_weight.second = container->mIdentifiers[element_index];

  //The position of the leaf node containing the weight of the entry
  //pointed to by the iterator.
  leaf_pos weight_pos = container->place_entry(element_index);

  //Now fill up the first part of the return value
  iterator_and_weight.first =
    container->mNodes[weight_pos.node].*(weight_pos.side_ptr);

  return iterator_and_weight;
}

/*!
 * @param[in]   n   Amount by which to increment the iterator
 *
 * @return      The incremented iterator
 */
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::const_iterator&
  BinaryTree<Weight, Selectee>::const_iterator::operator+=(difference_type n)
{
  element_index += n;
  return *this;
}

//! Negation of an equality comparison
template<class Weight, class Selectee>
  bool operator!=(const typename BinaryTree<Weight, Selectee>::const_iterator& a,
                  const typename BinaryTree<Weight, Selectee>::const_iterator& b)
{
  return !(a == b);
}



/*!
 * Move iterator n places forward
 *
 * @param[in]   n       Number of places to advance
 *
 * @return      Iterator advanced by n places
 */
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::const_iterator
  BinaryTree<Weight, Selectee>::const_iterator::operator+(const difference_type n) {
    // Construct a new instance to increment and return
    return (const_iterator(*this) += n);
}

/*!
 * Calculate the distance between two iterators
 *
 * @pre     Both iterators refer to the same container
 *
 * @param[in]   it       Iterator from which to calculate distance
 *
 * @return      n such that *this refers to the same entry as it + n
 */
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::difference_type
  BinaryTree<Weight, Selectee>::const_iterator::operator-(const const_iterator& it) const {
    return element_index - it.element_index;
}

/*!
 * Calculate the distance between two iterators
 *
 * @pre     Both iterators refer to the same container
 *
 * @param[in]   it       Iterator from which to calculate distance
 *
 * @return      n such that *this refers to the same entry as it + n
 */
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::difference_type
  BinaryTree<Weight, Selectee>::const_iterator::operator-(const iterator& it) const {
    return element_index - it.element_index;
}


/********* End of section 3 ******************************************/

/****** 4 STL container methods **************************************/
//Class needs to have its own iterator type to prevent things like
//delete begin() which would be pretty disastrous at the moment.


//Largest possible tree
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::size_type BinaryTree<Weight, Selectee>::max_size() const
{
  return sMaxSize;
}

template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::size_type BinaryTree<Weight, Selectee>::size() const
{
  return mFirstSpace;
}

/*!
 * Remove everything stored, but keep the memory allocated and the tree structure
 * initialised.
 */
template<class Weight, class Selectee> void BinaryTree<Weight, Selectee>::clear()
{
    mFirstSpace = 0;

    // Clear all the mNodes
    mNodes.assign(mCapacity, node_type());

    mIdentifiers.clear();
    mIdentifiers.reserve(mCapacity);

    assert(isValid());
}

//Iterator for the first entry
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::iterator BinaryTree<Weight, Selectee>::begin()
{
  return iterator(0, this);
}

//const_iterator for the first entry
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::const_iterator BinaryTree<Weight, Selectee>::begin() const
{
  return const_iterator(0, this);
}

//One past the last entry
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::iterator BinaryTree<Weight, Selectee>::end()
{
  return iterator(mFirstSpace, this);
}
//const_iterator version
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::const_iterator BinaryTree<Weight, Selectee>::end() const
{
  return const_iterator(mFirstSpace, this);
}

//iterator rbegin() {
//iterator rend() {

/*!
 * Put a value into the first available space
 *
 * \param[in] x Value to put into tree.
 *
 * \pre size() < max_size()
 *
 * \todo Handle the case when size() >= capacity() (Need to build a new tree of double the size)
 */
template<class Weight, class Selectee> void BinaryTree<Weight, Selectee>::push_back(const_value_type& x) {
    assert(isValid());
    //Find the leaf node above the new entry.
    leaf_pos leaf = place_entry(mFirstSpace);

    //increment mFirstSpace now that the old space is being used
    ++mFirstSpace;

    //Put the new entry to the left or right of the node as indiated
    //by the member pointer in leaf.
    mNodes[leaf.node].*(leaf.side_ptr) = x.first;

    //store the iterator that points back into the ensemble
    mIdentifiers.push_back(x.second);

    //update the tree structre above the new entry
    ascending_recalc(leaf.node);
}

/*!
 * Remove last entry in the tree
 *
 * \pre size() > 0
 *
 */
template<class Weight, class Selectee> void BinaryTree<Weight, Selectee>::pop_back() {
    // One more space at the end now
    --mFirstSpace;
    mIdentifiers.pop_back();

    //Find the leaf node above the new space
    leaf_pos leaf = place_entry(mFirstSpace);

    // A space has 0 weight
    mNodes[leaf.node].*(leaf.side_ptr) = weight_type();
    ascending_recalc(leaf.node);

    assert(isValid());
}

/*!
 * Replace the value pointed to by i with the values in x.
 * This will be more efficient than delete(i) followed by
 * insert(x) .
 * \param[in] i Pointer to value to replace.
 * \param[in] x New value to go into tree.
 */
template<class Weight, class Selectee> void BinaryTree<Weight, Selectee>::replace(iterator i, const value_type& x) {
  //first work out the tree leaf node the iterator refers to
  leaf_pos changed_leaf = place_entry(i);

  //update the list underneath the tree
  mIdentifiers[i.element_index] = x.second;

  //now update the tree above the change
  const_cast<weight_type&>(mNodes[changed_leaf.node].*(changed_leaf.side_ptr)) = x.first;
  //finally propogate the new rates up the tree
  ascending_recalc(changed_leaf.node);
}

/*!
 * Replace the current contents of the tree with a new tree of the smallest power of 2
 * size sufficient to hold all the data in the input sequence and that is no smaller
 * than the current tree.  The value of capacity() will not be reduced, to achieve
 * this a call to \ref resize will be necessary.
 * The order of the elements in this container will be the same as their order in
 * the supplied sequence.
 *
 * @tparam  ForwardIteratorType     Multi-pass iterator with value type BinaryTree<Weight, Selectee>::value_type
 *
 * @param   it_begin        Start of sequence of values to assign to this tree
 * @param   it_end          One past end of sequence of values to assign to this tree
 *
 */
template<class Weight, class Selectee> template<class ForwardIteratorType>
  void BinaryTree<Weight, Selectee>::assign(const ForwardIteratorType it_begin,
                                            const ForwardIteratorType it_end) {

    // First count the number of elements
    const size_type inputLength = std::distance(it_begin, it_end);

    //find k to meet description above
    size_type k = 1;
    while((1u << k) < inputLength) //2^k < inputLength
        ++k;

    // Does the tree need to grow
    if(mNumLevels < k) {
      //tree will need k levels to handle n entries
      mNumLevels = k;
      //maximum capacity of the tree is 2^k >= inputLength
      mCapacity = 1 << k; //2^k

      // Make space for the additional nodes in the summation structure
      mNodes.resize(mCapacity);
    }

    // Get rid of any existing return pointers
    mIdentifiers.clear();
    mIdentifiers.reserve(mCapacity);

    // Input will go in the first available spaces so the first unused
    // space will be equal to the length of the input.
    mFirstSpace = inputLength;

    // Now read in the elements one at a time
    ForwardIteratorType it = it_begin;
    size_type i = 0;
    while(it != it_end) {
        // Store the return pointer
        mIdentifiers.push_back(it->second);

        // Now work out what to do with the weight
        leaf_pos leafPosition = place_entry(i);
        mNodes[leafPosition.node].*leafPosition.side_ptr = it->first;

        // Move on to the next particle
        ++it;
        ++i;
    }

    // Set 0 weights for the empty spaces in the tree
    while(i != mCapacity) {
      // Zero weight
      leaf_pos leafPosition = place_entry(i);
      mNodes[leafPosition.node].*leafPosition.side_ptr = weight_type();

      ++i;
    }

    // Finally initialise the binary tree summation structure
    recalculate();

    assert(isValid());
}

/*!
 * Replace the current contents of the tree repeact_coutn copies of value,
 * increasing the size of the tree to the smallest power of 2 sufficient
 * to hold this number of values if necessary.
 * Capacity will not be reduced, to achieve a reduction in capacity a call
 * to \ref resize will be necessary.
 *
 * @param   repeat_count   Number of times to repeat value in the tree
 * @param   value          Value be repeatedly inserted into the tree
 *
 */
template<class Weight, class Selectee>
  void BinaryTree<Weight, Selectee>::assign(size_type repeat_count,
                                            const value_type& value) {
    std::vector<value_type> newValues(repeat_count, value);
    assign(newValues.begin(), newValues.end());
}

/*!
 * Replace the current contents of the tree with a new tree of the smallest power of 2
 * size sufficient to hold new_size entries.  The first new_size values of the tree
 * will be preserved.  Internal re-arrangements may invalidate all iterators and change
 * capacity().
 *
 * @param   new_size    Number of elements to remain in the tree, must be >= 2
 *
 */
template<class Weight, class Selectee>
  void BinaryTree<Weight, Selectee>::resize(size_type new_size) {
    // Index of the first leaf node in mNodes
    const size_type halfCapacity = mCapacity / 2;

    // Index of last entry that should be copied into new tree.
    const size_type lastData = std::min(new_size, mFirstSpace);

    // Deliberate integer division to find the number of leaf nodes which
    // have values to both their left and their right that need copying
    // into the new tree.
    const size_type newSize2 = lastData / 2;

    // Use this vector to build up the values that will be kept in
    // the resized tree.
    typedef std::pair<weight_type, return_pointer_type> mutable_value_type;
    std::vector<mutable_value_type> newValues;
    newValues.reserve(new_size);
    // Copy the values for all the leaf nodes that have values to both left and right
    for(size_type i = 0;  i < newSize2; ++i) {
        newValues.push_back(mutable_value_type(mNodes[halfCapacity + i].left,  mIdentifiers[2 * i]));
        newValues.push_back(mutable_value_type(mNodes[halfCapacity + i].right, mIdentifiers[2 * i + 1]));
    }
    // See if there is one remaining value to copy
    if((lastData > 0) && (lastData % 2) == 0) {
        newValues.push_back(mutable_value_type(mNodes[halfCapacity + newSize2].left, mIdentifiers[lastData]));
    }

    //find smallest possible number of levels for the tree
    size_type k = 1;
    while((1u << k) < new_size) //2^k < new_size
        ++k;

    //tree will need k levels to handle new_size entries
    mNumLevels = k;

    //maximum capacity of the tree is 2^k
    mCapacity = 1u << k;

    // Maintain internal consistency with member container sizes
    mNodes.resize(mCapacity);
    mIdentifiers.reserve(mCapacity);

    // Set the new values
    assign(newValues.begin(), newValues.end());
}
/********* End of section 4 ******************************************/

/****** 5 Comparison methods *****************************************/

/*! Two binary trees are equal if they have exactly the same data in the
 * same internal order.  This says nothing about whether they have the
 * same capacities so an anomaly if possible in that if a == b then the
 * assignment a = b; may still change a.
 * \param[in] a First tree
 * \param[in] b Second tree
 */
template<class Weight, class Selectee>
  bool operator==(const BinaryTree<Weight, Selectee>& a, const BinaryTree<Weight, Selectee>& b)
{
  return a.size() == b.size() && \
                     std::equal(a.begin(), a.end(), b.begin());
}

//! Negative of equaltiy comparison
template<class Weight, class Selectee>
  bool operator!=(const BinaryTree<Weight, Selectee>& a, const BinaryTree<Weight, Selectee>& b)
{
  return !operator==(a, b);
}

/****** End of section 5 *********************************************/

/****** 6 Other methods **********************************************/

/*!
 * Return a leaf_pos structure containing the position of the leaf node
 * under which an entry at position entry in the mIdentifiers vector
 * belongs.
 *
 * @param[in]   entry   Index in mIdentifiers
 *
 * @return      Index and side of position in mNodes where the weights from entry should be stored
 */
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::leaf_pos
  BinaryTree<Weight, Selectee>::place_entry(const size_type entry) const
{
  leaf_pos lp;

  //leaf is to the right of the node if entry is an odd number
  lp.side_ptr = (entry % 2 ? &node_type::right : &node_type::left);

  //node number is 2^(mNumLevels - 1) + entry/2 (integer division)
  lp.node = (1 << (mNumLevels - 1)) + (entry >> 1);

  return lp;
}

/*!
 * Another version of place_entry to convert directly from iterators.
 *
 *@param[in]    i   Iterator from which to calculate position
 *
 *@return       Structure specifying to which side of which leaf node the iterator points
 */
template<class Weight, class Selectee>
 typename BinaryTree<Weight, Selectee>::leaf_pos
 BinaryTree<Weight, Selectee>::place_entry(const_iterator i) const
{
  return place_entry(i.element_index);
}

/*!
 * This is the inverse of \ref place_entry.
 *
 *@param[in]    leaf_pos    Structure specifying to which side of which leaf node an entry is
 *
 *@return       Index of entry in mIdentifiers
 */
template<class Weight, class Selectee>
  typename BinaryTree<Weight, Selectee>::size_type
  BinaryTree<Weight, Selectee>::entry_leaf(const leaf_pos& leaf) const
{
  //The direct inversion of place_entry would be
  //index = (leaf.node - (2^(mNumLevels - 1))) * 2
  //I have taken the *2 inside the outermost brackets.
  size_type index = (leaf.node << 1) - (1 << mNumLevels);
  if(leaf.side_ptr == &node_type::right) {
    //leaf is right hand child of node
    ++index;
  }
  return index;
}

/*!
 * Propagate the values in start_node up to the top node.
 *
 * @param[in]   start_node      Index of leaf node from which to being recalculation
 */
template<class Weight, class Selectee> void BinaryTree<Weight, Selectee>::ascending_recalc(size_type start_node) {
  size_type i, j;

  i = start_node; //start from start_node

  while(i > 1) { //1 is the top node so we can stop when we reach it.
    //in this loop i is the index of the current node, j of its parent
    j = (i >> 1); //integer division by 2 to find index of parent

    if(i & 0x1) { //i is to right of parent
      mNodes[j].right = mNodes[i].left + mNodes[i].right;
    }
    else { //i is to left of parent
      mNodes[j].left = mNodes[i].left + mNodes[i].right;
    }
    i = j; //move up to handle parent of i on next iteration
  }

  // Cache the overall sums
  mNodes[0].left = mNodes[1].left + mNodes[1].right;
}

/*!
 * @return      Sum of each of the weights in the tree
 */
template<class Weight, class Selectee>
  const typename BinaryTree<Weight, Selectee>::weight_type& BinaryTree<Weight, Selectee>::head() const
{
  return mNodes[0].left;
}

template<class Weight, class Selectee> void BinaryTree<Weight, Selectee>::print_tree() const
{
  size_type i,j;
  for(i = 0; i < mNumLevels; ++i) {
    for(j = 0; j < (1 << i); ++j) {
      std::cout << mNodes[(1 << i) + j].left << ", ";
      std::cout << mNodes[(1 << i) + j].right << "; ";
    }
    std::cout << '\n';
  }
}

/*!
 * Propagate the values in the leaf mNodes up the tree
 */
template<class Weight, class Selectee> void BinaryTree<Weight, Selectee>::recalculate() {

    // Loop over all the levels in the tree, except for the first level, which
    // is assumed to already contain the data to be summed.
    for(size_type i = mNumLevels - 1; i > 0; --i) {
        //2^(i-1) is the width of level i of the tree
        const size_type rowLength = 1u << (i - 1);

        for(size_type j = rowLength; j < 2 * rowLength; ++j) {
            mNodes[j].left  = mNodes[2 * j].left     + mNodes[2 * j].right;
            mNodes[j].right = mNodes[2 * j + 1].left + mNodes[2 * j + 1].right;
        }
    }

    // Calculate the overall sum
    mNodes[0].left = mNodes[1].left + mNodes[1].right;
}

/*!
 * Select a particle according to a caller specified weighting
 *
 * @tparam  ScalarExtractorType Univariate function that returns a real when
 *          applied to an instance of Weight.  All reasonable uses are likely
 *          to be linear functions of the properties in Weight.
 *
 * @param[in]   d   Value indicating where in the distribution to select the particle.
 *                  In general, if r is a U[0,1] random variable then d should be
 *                  generated as r * f(head())
 *
 * @param[in]   f   Univariate function object
 *
 * @return  Iterator to a selected entry, or end if nothing selected
 *
 */
template<class Weight, class Selectee> template<typename ScalarExtractorType>
  typename BinaryTree<Weight, Selectee>::const_iterator
  BinaryTree<Weight, Selectee>::select(double d, ScalarExtractorType f) const {
    // Descend the tree one level at a time
    size_type selectedIndex = 1;
    for(size_type i = 0; i < mNumLevels; ++i) {
        const double leftWeight = f(mNodes[selectedIndex].left);
        if(d <= leftWeight) {
            // Take left branch
            selectedIndex *= 2;
        }
        else {
            // This is a special case to avoid descending a zero weighted part
            // of the tree due to floating point imprecision.  This case implies
            // that the rightmost entry should be selected
            if(f(mNodes[selectedIndex].right) == 0.0) {
                // Add mCapacity because it will be subtracted below
                assert(mFirstSpace >= 1);
                selectedIndex = mFirstSpace + mCapacity - 1;
                break;
            }
            // Take right branch
            d -= leftWeight;
            selectedIndex = 2 * selectedIndex + 1;
        }
    }

    // The final left/right choice should take index beyond the leaf mNodes
    // of the tree
    assert(selectedIndex >= mCapacity);
    selectedIndex -= mCapacity;

    return const_iterator(std::min(selectedIndex, mFirstSpace) , this);
}

template<class Weight, class Selectee> template<typename ScalarExtractorType>
  typename BinaryTree<Weight, Selectee>::iterator
  BinaryTree<Weight, Selectee>::select(double d, ScalarExtractorType f) {

    // Call the const version to do the calculation
    const_iterator selectee = (static_cast<const BinaryTree<Weight, Selectee> >(this))->select(d, f);

    // Now convert to a non-const iterator
    return iterator(selectee.element_index, this);
}

/*
 * Sanity checks with error messages to describe any problems
 */
template<class Weight, class Selectee> bool BinaryTree<Weight, Selectee>::isValid() const {
    bool valid = true;

    // (1u << mNumLevels) == 2^mNumLevels
    if(!((mCapacity == (1u << mNumLevels)) || ((mCapacity == 0) && (mNumLevels == 0)))) {
        std::cerr << "Total capacity is " << mCapacity << ", but there are " << mNumLevels << "levels\n";
        valid = false;
    }

    if(!(mFirstSpace <= mCapacity)){
        std::cerr << "First space is " << mFirstSpace << ", but capacity is only " << mCapacity << '\n';
        valid = false;
    }

    if(!(mCapacity == mNodes.size())){
        std::cerr << "Number of nodes in tree (" << mNodes.size()
                  << ") inconsistent with capacity of " << mCapacity << '\n';
        valid = false;
    }

    if(mFirstSpace != mIdentifiers.size()){
        std::cerr << "Number of return pointers in tree (" << mIdentifiers.size()
                  << ") inconsistent with position of first space " << mFirstSpace << '\n';
        valid = false;
    }

    return valid;
}

} //namespace Utils
#endif //UTILS_BINARY_TREE_HPP
