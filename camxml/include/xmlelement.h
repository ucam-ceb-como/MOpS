/*
  Author(s):      Matthew Celnik (msc37)
  Project:        camxml (Cambridge c++ XML library)

  File purpose:
    Defines an XML element.
*/

#ifndef CAMXML_ELEMENT_H
#define CAMXML_ELEMENT_H

#include "xmlattribute.h"
#include <vector>
#include <string>

namespace CamXML
{
class Element
{
private:
    /* This element's tag. */
    std::string m_tag;
    /* Text stored in this element. */
    std::string m_data;
    /* Vector of element attributes. */
    std::vector<Attribute*> m_attr;
    /* Parent node of this node. */
    Element *m_parent;
    /* List of child nodes of this element. */
    std::vector<Element*> m_children;
public:
    Element(void);
    Element(const std::string &tag, const std::string &data, Element *const parent);
    ~Element(void);
public:
    /* Returns this element's tag. */
    const std::string &Tag(void) const;
    /* Returns this element's data. */
    const std::string &Data(void) const;
    /* Returns the vector of attributes. */
    const std::vector<Attribute*> &GetAttributes(void) const;
    /* Returns the attribute with the given name if exists,
       otherwise returns NULL. */
    const Attribute *const GetAttribute(const std::string &name) const;
    /* Returns the value of the attribute with the given name if exists,
       otherwise returns empty string. */
    const std::string GetAttributeValue(const std::string &name) const;
    /* Returns pointer to this elements parent if it has one,
       otherwise returns NULL. */
    Element *const Parent(void) const;
    /* Returns the vector of all child elements. */
    const std::vector<Element*> &GetChildren(void) const;
    /* Returns a vector of all child elements which have the
       given tag. */
    void GetChildren(const std::string &tag, std::vector<Element*> &children) const;
    /* Returns the first child under the current element with the given tag.
       If no child is found then returns NULL. */
    Element *const GetFirstChild(const std::string &tag) const;
    /* Returns the index of the child under this element if exists,
       otherwise return negative. */
    int GetChildIndex(const Element *const child) const;
    /* Returns true is the element contains the given child, otherwise
       false. */
    bool Contains(const Element *const child) const;
public:
    /* Sets the element's tag. */
    void SetTag(const std::string &tag);
    /* Sets the element's data string. */
    void SetData(const std::string &data);
    /* Sets the attribute with the given name to have the given
       value.  If the attribute is not in the list then it is added. */
    void SetAttribute(const std::string &name, const std::string &value);
    /* Sets the parent element of the element. Returns pointer to
       new parent on success, otherwise returns NULL. */
    const Element *const SetParent(Element *const parent);
    /* Adds a child element to the element.  Returns pointer to
       child on success, otherwise NULL. */
    const Element *const AddChild(Element *const child);
public:
    /* Deletes all attributes from the element. */
    void DeleteAttributes(void);
    /* Deletes the attribute with the give name from the element. */
    void DeleteAttribute(const std::string &name);
    /* Deletes all children of this element. */
    void DeleteChildren(void);
    /* Deletes all children with the given tag. */
    void DeleteChildren(const std::string &tag);
    /* Deletes the child at the given index from the element. */
    void DeleteChild(const unsigned int i);
public:
    /* Moves a child node to a new parent. */
    int MoveChildTo(Element *const child, Element *const newparent);
};

/* Inline get function definitions. */
inline const std::string &Element::Tag() const {return m_tag;};
inline const std::string &Element::Data() const {return m_data;};
inline const std::vector<Attribute*> &Element::GetAttributes() const {return m_attr;};
inline Element *const Element::Parent() const {return m_parent;};
inline const std::vector<Element*> &Element::GetChildren() const {return m_children;};

/* Inline set function definitions. */
inline void Element::SetTag(const std::string &tag) {m_tag = tag;};
inline void Element::SetData(const std::string &data) {m_data = data;};
};

#endif
