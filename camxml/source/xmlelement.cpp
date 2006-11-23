#include "xmlelement.h"

using namespace CamXML;
using namespace std;

Element::Element(void)
{
    m_tag = "";
    m_data = "";
    m_attr.clear();
    m_parent = NULL;
    m_children.clear();
}

Element::Element(const std::string &tag, const std::string &data, Element *const parent)
{
    m_tag = tag;
    m_data = data;
    m_parent = parent;
}

Element::~Element(void)
{
    DeleteAttributes();
    DeleteChildren();
}

/* Definitions of property get functions. */

const Attribute *const Element::GetAttribute(const std::string &name) const
{
    vector<Attribute*>::const_iterator i;
    for (i=m_attr.begin(); i!=m_attr.end(); i++) {
        if ((*i)->GetName().compare(name) == 0) {
            return (*i);
        }
    }
    return NULL;
}

const std::string Element::GetAttributeValue(const std::string &name) const
{
    const Attribute *const pattr = GetAttribute(name);
    if (pattr) {
        return pattr->GetValue();
    } else {
        return "";
    }
}

void Element::GetChildren(const std::string &tag, vector<Element *> &children) const
{
    // Clear any elements already in the child list.
    children.clear();

    // Loop through child elements, and save pointers to those
    // with matching tags.
    vector<Element*>::const_iterator i;
    for (i=m_children.begin(); i!=m_children.end(); i++) {
        if ((*i)->Tag().compare(tag) == 0) {
            children.push_back(*i);
        }
    }
}

Element *const Element::GetFirstChild(const std::string &tag) const
{
    // Loop through child elements, and find that with the given tag.
    vector<Element*>::const_iterator i;
    for (i=m_children.begin(); i!=m_children.end(); i++) {
        if ((*i)->Tag().compare(tag) == 0) {
            return *i;
        }
    }
    return NULL;
}

int Element::GetChildIndex(const CamXML::Element *const child) const
{
    int k;
    vector<Element*>::const_iterator i;
    for (i=m_children.begin(),k=0; i!=m_children.end(); i++,k++) {
        if ((*i) == child) {
            return k;
        }
    }
    return -1;
}

bool Element::Contains(const CamXML::Element *const child) const
{
    vector<Element*>::const_iterator i;
    for (i=m_children.begin(); i!=m_children.end(); i++) {
        if ((*i) == child) {
            return true;
        }
    }
    return false;
}

/* Definitions of property set functions. */

void Element::SetAttribute(const std::string &name, const std::string &value)
{
    vector<Attribute*>::iterator i;
    for (i=m_attr.begin(); i!=m_attr.end(); i++) {
        if ((*i)->GetName().compare(name) == 0) {
            (*i)->SetValue(value);
            return;
        }
    }

    // Got here because the element does not yet contain this attribute.
    m_attr.push_back(new Attribute(name, value));
}

const Element *const Element::SetParent(CamXML::Element *const parent)
{
    // Check that the new parent does not already contain this
    // element.
    if (!parent->Contains(this)) {
        if (m_parent!=NULL) {
            // Element already has a parent, so move it from
            // one to the other.
            m_parent->MoveChildTo(this, parent);
        } else {
            // This element does not have a parent, so just add
            // it to the new one.
            parent->AddChild(this);
        }
    }
    m_parent = parent;
    return m_parent;
}

const Element *const Element::AddChild(CamXML::Element *const child)
{
    // Check that the child is not already under this element.
    if (!Contains(child)) {
        m_children.push_back(child);
        child->SetParent(this);
    }
    return child;
}

/* Definitions of functions which delete elements and parts thereof. */

void Element::DeleteAttributes()
{
    vector<Attribute*>::iterator i;
    for (i=m_attr.begin(); i!=m_attr.end(); i++) {
        delete *i;
    }
    m_attr.clear();
}

void Element::DeleteAttribute(const std::string &name)
{
    vector<Attribute*>::iterator i;
    int k;
    for (i=m_attr.begin(),k=0; i!=m_attr.end(); i++,k++) {
        if ((*i)->GetName().compare(name) == 0) {
            delete *i;
            m_attr.erase(i);
            return;
        }
    }
}

void Element::DeleteChildren()
{
    vector<Element*>::iterator i;
    for(i=m_children.begin(); i!=m_children.end(); i++) {
        delete *i;
    }
    m_children.clear();
}

void Element::DeleteChildren(const std::string &tag)
{
    vector<Element*>::iterator i = m_children.begin();
    while(i!=m_children.end()) {
        if ((*i)->Tag().compare(tag) == 0) {
            delete *i;
            i = m_children.erase(i);
        } else {
            i++;
        }
    }
}

void Element::DeleteChild(const unsigned int i)
{
    if (i<(unsigned int)m_children.size()) {
        delete m_children[i];
        m_children.erase(m_children.begin()+i);
    }
}

/* Definitions of other Element functions. */

int Element::MoveChildTo(CamXML::Element *const child, CamXML::Element *const newparent)
{
    if (this!=newparent) {
        int i = this->GetChildIndex(child);
        if ((i>=0) && (newparent->GetChildIndex(child) < 0)) {
            // We have checked that we are not moving the child to the same
            // parent element, that the child belongs to this element and that
            // the new parent does not contain the child already.
            newparent->AddChild(child);
            m_children.erase(m_children.begin()+i);
            return 0;
        } else {
            return -1;
        }
    } else {
        return -2;
    }
}