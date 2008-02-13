#include "xmlattribute.h"

using namespace CamXML;

Attribute::Attribute(void)
{
}

Attribute::Attribute(const std::string &name, const std::string &value)
{
    m_name = ComoUnicode::StringToWString(name);
    m_value = ComoUnicode::StringToWString(value);
}
Attribute::Attribute(const std::wstring &name, const std::wstring &value)
{
    m_name = name;
    m_value = value;
}

Attribute::~Attribute(void)
{
}
