#include "xmlelement.h"

using namespace CamXML;
using namespace std;

Element::Element(void)
{
    m_tag = L"";
    m_data = L"";
    m_attr.clear();
    m_parent = NULL;
    m_children.clear();
}

Element::Element(const std::string &tag, const std::string &data, Element &parent)
{
    m_tag = ComoUnicode::StringToWString(tag);
    m_data = ComoUnicode::StringToWString(data);
    m_parent = &parent;
}

Element::Element(const std::wstring &tag, const std::wstring &data, Element &parent)
{
    m_tag = tag;
    m_data = data;
    m_parent = &parent;
}

Element::~Element(void)
{
    DeleteAttributes();
    DeleteChildren();
}

/* Definitions of property get functions. */

const Attribute *const Element::GetAttribute(const std::string &name) const
{
    return GetAttribute(ComoUnicode::StringToWString(name));
}

const Attribute *const Element::GetAttribute(const std::wstring &name) const
{
    vector<Attribute*>::const_iterator i;
    for (i=m_attr.begin(); i!=m_attr.end(); i++) {
        if ((*i)->wGetName().compare(name) == 0) {
            return (*i);
        }
    }
    return NULL;
}

const std::string Element::GetAttributeValue(const std::string &name) const
{
    return ComoUnicode::WStringToString(wGetAttributeValue(name));
}

const std::wstring Element::wGetAttributeValue(const std::string &name) const
{
    return wGetAttributeValue(ComoUnicode::StringToWString(name));
}

const std::wstring Element::wGetAttributeValue(const std::wstring &name) const
{
    const Attribute *const pattr = GetAttribute(name);
    if (pattr) {
        return pattr->wGetValue();
    } else {
        return L"";
    }
}

void Element::GetChildren(const std::string &tag, vector<Element *> &children) const
{
    GetChildren(ComoUnicode::StringToWString(tag), children);
}

void Element::GetChildren(const std::wstring &tag, vector<Element *> &children) const
{
    // Clear any elements already in the child list.
    children.clear();

    // Loop through child elements, and save pointers to those
    // with matching tags.
    vector<Element*>::const_iterator i;
    for (i=m_children.begin(); i!=m_children.end(); i++) {
        if ((*i)->wTag().compare(tag) == 0) {
            children.push_back(*i);
        }
    }
}

Element *const Element::GetFirstChild(const std::string &tag) const
{
    return GetFirstChild(ComoUnicode::StringToWString(tag));
}

Element *const Element::GetFirstChild(const std::wstring &tag) const
{
    // Loop through child elements, and find that with the given tag.
    vector<Element*>::const_iterator i;
    for (i=m_children.begin(); i!=m_children.end(); i++) {
        if ((*i)->wTag().compare(tag) == 0) {
            return *i;
        }
    }
    return NULL;
}

int Element::GetChildIndex(const CamXML::Element &child) const
{
    int k;
    vector<Element*>::const_iterator i;
    for (i=m_children.begin(),k=0; i!=m_children.end(); i++,k++) {
        if ((*i) == &child) {
            return k;
        }
    }
    return -1;
}

bool Element::Contains(const CamXML::Element &child) const
{
    vector<Element*>::const_iterator i;
    for (i=m_children.begin(); i!=m_children.end(); i++) {
        if ((*i) == &child) {
            return true;
        }
    }
    return false;
}
 
const std::wstring Element::GetXMLString(const int depth) const {
    std::wstring m_data_temp = ComoString::trim(m_data);
    
    std::wstring xml_temp = ComoString::insertSpace(depth) + L"<" + m_tag;
    // Add Attributes to this element XML Tag
    vector<Attribute*>::const_iterator i;
    for (i=m_attr.begin(); i!=m_attr.end(); i++) {
        xml_temp.append(L" " + (*i)->wGetName() + L"=\"" + (*i)->wGetValue() + L"\"");
    }
    if ((m_children.size() == 0) && (m_data_temp.compare(L"") == 0)) {
        xml_temp.append(L"/>");
    } else {
        xml_temp.append(L">");
        vector<Element*>::const_iterator j;
        for (j=m_children.begin(); j!=m_children.end(); j++) {
            xml_temp.append(L"\n" + (*j)->GetXMLString(depth+1));
        }
        if (m_data_temp.compare(L"") != 0) {
            if (m_children.size() == 0) {
                xml_temp.append(m_data_temp);
                xml_temp.append(L"</" + m_tag + L">");
            } else {
                xml_temp.append(L"\n" + ComoString::insertSpace(depth+1) + m_data_temp);
                xml_temp.append(L"\n" + ComoString::insertSpace(depth) + L"</" + m_tag + L">");
            }
        } else {
            xml_temp.append(L"\n" + ComoString::insertSpace(depth) + L"</" + m_tag + L">");
        }
        
    }
    return xml_temp;
}

const std::wstring Element::GetXMLString() const {
	return GetXMLString(0);
}


/* Definitions of property set functions. */

void Element::SetAttribute(const std::string &name, const std::string &value)
{
    SetAttribute(ComoUnicode::StringToWString(name),ComoUnicode::StringToWString(value));
}

void Element::SetAttribute(const std::wstring &name, const std::wstring &value)
{
    vector<Attribute*>::iterator i;
    for (i=m_attr.begin(); i!=m_attr.end(); i++) {
        if ((*i)->wGetName().compare(name) == 0) {
            (*i)->SetValue(value);
            return;
        }
    }

    // Got here because the element does not yet contain this attribute.
    m_attr.push_back(new Attribute(name, value));
}

const Element *const Element::SetParent(CamXML::Element &parent)
{
    // Check that the new parent does not already contain this
    // element.
    if (!parent.Contains(*this)) {
        if (m_parent!=NULL) {
            // Element already has a parent, so move it from
            // one to the other.
            m_parent->MoveChildTo(*this, parent);
        } else {
            // This element does not have a parent, so just add
            // it to the new one.
            parent.AddChild(*this, false);
        }
    }
    m_parent = &parent;
    return m_parent;
}

const Element *const Element::AddChild(const std::string &tag)
{
    return AddChild(ComoUnicode::StringToWString(tag));
}

const Element *const Element::AddChild(const std::wstring &tag)
{
    // Check that the child is not already under this element.
    Element * el = new Element(tag, L"", *this->m_parent);
    return AddChild(*el, false);
}

const Element *const Element::AddChild(CamXML::Element &child, bool clone)
{
    Element * a_child;
    if (clone) {
        a_child = child.Clone();
    } else {
        a_child = &child;
    }
        // Check that the child is not already under this element.
        if (!Contains(*a_child)) {
            m_children.push_back(a_child);
            a_child->SetParent(*this);
        }
    return a_child;
}

Element *const Element::Clone() const
{
    Element * el = new Element();
    el->m_tag = this->wTag();
    el->m_data = this->wData();
    el->m_parent = NULL;
    // Allocate new memories for Attributes, copy values and push them into m_attr
    size_t n_attr = this->m_attr.size();
    el->m_attr.clear();
    for(size_t i = 0; i < n_attr; i++) {
        Attribute * attr = new Attribute(this->m_attr.at(i)->wGetName(), this->m_attr.at(i)->wGetValue());
        el->m_attr.push_back(attr);
    }
    // Allocate new memories for Elements, deep copy elements and push them into m_children
    size_t n_ch = this->m_children.size();
    el->m_children.clear();
    for(size_t i = 0; i < n_ch; i++) {
        Element * ch_el = this->m_children.at(i)->Clone();
        el->m_children.push_back(ch_el);
    }

    return el;
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
    DeleteAttribute(ComoUnicode::StringToWString(name));
}

void Element::DeleteAttribute(const std::wstring &name)
{
    vector<Attribute*>::iterator i;
    int k;
    for (i=m_attr.begin(),k=0; i!=m_attr.end(); i++,k++) {
        if ((*i)->wGetName().compare(name) == 0) {
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
    DeleteChildren(ComoUnicode::StringToWString(tag));
}

void Element::DeleteChildren(const std::wstring &tag)
{
    vector<Element*>::iterator i = m_children.begin();
    while(i!=m_children.end()) {
        if ((*i)->wTag().compare(tag) == 0) {
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

int Element::MoveChildTo(CamXML::Element &child, CamXML::Element &newparent)
{
    if (this!=&newparent) {
        int i = this->GetChildIndex(child);
        if ((i>=0) && (newparent.GetChildIndex(child) < 0)) {
            // We have checked that we are not moving the child to the same
            // parent element, that the child belongs to this element and that
            // the new parent does not contain the child already.
            newparent.AddChild(child, false);
            m_children.erase(m_children.begin()+i);
            return 0;
        } else {
            return -1;
        }
    } else {
        return -2;
    }
}