/*
  Author(s):      Matthew Celnik (msc37)
  Project:        camxml (Cambridge c++ XML library)
  Sourceforge:    http://sourceforge.net/projects/mopssuite
  
  Copyright (C) 2008 Matthew S Celnik.

  File purpose:
    Defines an XML element.

  Licence:
    This file is part of the "camxml" library.

    camxml is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.

  Contact:
    Dr Markus Kraft
    Dept of Chemical Engineering
    University of Cambridge
    New Museums Site
    Pembroke Street
    Cambridge
    CB2 3RA
    UK

    Email:       mk306@cam.ac.uk
    Website:     http://como.cheng.cam.ac.uk
*/

#ifndef CAMXML_ELEMENT_H
#define CAMXML_ELEMENT_H

#include "unicodef.h"
#include "stringf.h"
#include "xmlattribute.h"
#include <vector>
#include <string>

namespace CamXML
{
class Element
{
private:
    /* This element's tag. */
    std::wstring m_tag;
    /* Text stored in this element. */
    std::wstring m_data;
    /* Vector of element attributes. */
    std::vector<Attribute*> m_attr;
    /* Parent node of this node. */
    Element *m_parent;
    /* List of child nodes of this element. */
    std::vector<Element*> m_children;
public:
    Element(void);
    Element(const std::string &tag, const std::string &data, Element &parent);
    Element(const std::wstring &tag, const std::wstring &data, Element &parent);
    ~Element(void);
public:
    /* Returns this element's tag. */
    const std::string Tag(void) const;
    /* Returns this element's data. */
    const std::string Data(void) const;
    /* Returns this element's tag as wide string. */
    const std::wstring &wTag(void) const;
    /* Returns this element's data as wide string. */
    const std::wstring &wData(void) const;
    /* Returns the vector of attributes. */
    const std::vector<Attribute*> &GetAttributes(void) const;
    /* Returns the attribute with the given name if exists,
       otherwise returns NULL. */
    const Attribute *const GetAttribute(const std::string &name) const;
    /* Returns the attribute with the given name (wide string) if exists,
       otherwise returns NULL. */
    const Attribute *const GetAttribute(const std::wstring &name) const;
    /* Returns the value of the attribute with the given name (wide string) if exists,
       otherwise returns empty string. */
    const std::string GetAttributeValue(const std::string &name) const;
    /* Returns the wide string value of the attribute with the given name if exists,
       otherwise returns empty string. */
    const std::wstring wGetAttributeValue(const std::string &name) const;
    /* Returns the wide string value of the attribute with the given name if exists,
       otherwise returns empty string. */
    const std::wstring wGetAttributeValue(const std::wstring &name) const;
    /* Returns pointer to this element's parent if it has one,
       otherwise returns NULL. */
    Element *const Parent(void) const;
    /* Returns the vector of all child elements. */
    const std::vector<Element*> &GetChildren(void) const;
    /* Returns a vector of all child elements which have the
       given tag. */
    void GetChildren(const std::string &tag, std::vector<Element*> &children) const;
    /* Returns a vector of all child elements which have the
       given tag. */
    void GetChildren(const std::wstring &tag, std::vector<Element*> &children) const;
    /* Returns the child of a specified index.
       If no child is found then returns NULL. */
    Element *const GetChild(const int index) const;
    /* Returns the first child under the current element with the given tag.
       If no child is found then returns NULL. */
    Element *const GetFirstChild(const std::string &tag) const;
    /* Returns the first child under the current element with the given tag.
       If no child is found then returns NULL. */
    Element *const GetFirstChild(const std::wstring &tag) const;
    /* Returns the index of the child under this element if exists,
       otherwise return negative. */
    int GetChildIndex(const Element &child) const;
    /* Returns true is the element contains the given child, otherwise
       false. */
    bool Contains(const Element &child) const;
    /* Return XML wide string for this element */
    const std::wstring GetXMLString() const;
private:
	/* Return XML wide string for this element. Depth is needed for inserting correct
       number of spaces for each element depth. */
	const std::wstring GetXMLString(const int depth) const;
    Element *const Clone() const;
public:
    /* Sets the element's tag. */
    void SetTag(const std::string &tag);
    /* Sets the element's data string. */
    void SetData(const std::string &data);
    /* Sets the element's tag. */
    void SetTag(const std::wstring &tag);
    /* Sets the element's data wstring. */
    void SetData(const std::wstring &data);
    /* Sets the attribute with the given name to have the given
       value.  If the attribute is not in the list then it is added. */
    void SetAttribute(const std::string &name, const std::string &value);
    /* Sets the attribute with the given name to have the given
       value.  If the attribute is not in the list then it is added. */
    void SetAttribute(const std::wstring &name, const std::wstring &value);
    /* Sets the parent element of the element. Returns pointer to
       new parent on success, otherwise returns NULL. */
    const Element *const SetParent(Element &parent);
    /* Creates a new child element by specified tag and add it to the element.
       Returns pointer to child on success, otherwise NULL. */
    const Element *const AddChild(const std::string &tag);
    /* Creates a new child element by specified tag and add it to the element.
       Returns pointer to child on success, otherwise NULL. (wide string version) */
    const Element *const AddChild(const std::wstring &tag);
    /* Adds a child element to the element.  Returns pointer to
       child on success, otherwise NULL. */
    const Element *const AddChild(Element &child, bool clone=true);
public:
    /* Deletes all attributes from the element. */
    void DeleteAttributes(void);
    /* Deletes the attribute with the give name from the element. */
    void DeleteAttribute(const std::string &name);
    /* Deletes the attribute with the give name from the element. (wide string version) */
    void DeleteAttribute(const std::wstring &name);
    /* Deletes all children of this element. */
    void DeleteChildren(void);
    /* Deletes all children with the given tag. */
    void DeleteChildren(const std::string &tag);
    /* Deletes all children with the given tag.  (wide string version) */
    void DeleteChildren(const std::wstring &tag);
    /* Deletes the child at the given index from the element. */
    void DeleteChild(const unsigned int i);
public:
    /* Moves a child node to a new parent. */
    int MoveChildTo(Element &child, Element &newparent);
};

/* Inline get function definitions. */
inline const std::string Element::Tag() const {return ComoUnicode::WStringToString(m_tag);};
inline const std::string Element::Data() const {return ComoUnicode::WStringToString(m_data);};
inline const std::wstring &Element::wTag() const {return m_tag;};
inline const std::wstring &Element::wData() const {return m_data;};
inline const std::vector<Attribute*> &Element::GetAttributes() const {return m_attr;};
inline Element *const Element::Parent() const {return m_parent;};
inline const std::vector<Element*> &Element::GetChildren() const {return m_children;};

/* Inline set function definitions. */
inline void Element::SetTag(const std::string &tag) {m_tag = ComoUnicode::StringToWString(tag);};
inline void Element::SetData(const std::string &data) {m_data = ComoUnicode::StringToWString(data);};
inline void Element::SetTag(const std::wstring &tag) {m_tag = tag;};
inline void Element::SetData(const std::wstring &data) {m_data = data;};
};

#endif
