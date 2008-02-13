/*
  Author(s):      Matthew Celnik (msc37)
  Project:        camxml (Cambridge c++ XML library)

  File purpose:
    Defines an XML attribute.
*/

#ifndef CAMXML_ATTRIBUTE_H
#define CAMXML_ATTRIBUTE_H

#include "unicodef.h"
#include <vector>
#include <string>

namespace CamXML
{
class Attribute
{
private:
    /* Attribute name (bit before '=' sign). */
    std::wstring m_name;
    /* Attribute value (bit after '=' sign). */
    std::wstring m_value;
public:
    Attribute(void);
    Attribute(const std::string &name, const std::string &value);
    Attribute(const std::wstring &name, const std::wstring &value);
    ~Attribute(void);
public:
    /* Returns this attribute's name. */
    const std::string GetName(void) const;
    /* Returns this attribute's value. */
    const std::string GetValue(void) const;
    /* Returns this attribute's name as wide string. */
    const std::wstring &wGetName(void) const;
    /* Returns this attribute's value as wide string */
    const std::wstring &wGetValue(void) const;
public:
    /* Sets the attribute name. */
    void SetName(const std::string &name);
    /* Sets the attribute value. */
    void SetValue(const std::string &value);
    /* Sets the attribute name. (wide string version) */
    void SetName(const std::wstring &name);
    /* Sets the attribute value. (wide string version) */
    void SetValue(const std::wstring &value);
};

/* Inline function definitions. */
inline const std::string Attribute::GetName() const {return ComoUnicode::WStringToString(m_name);};
inline const std::string Attribute::GetValue() const {return ComoUnicode::WStringToString(m_value);};
inline const std::wstring &Attribute::wGetName() const {return m_name;};
inline const std::wstring &Attribute::wGetValue() const {return m_value;};
inline void Attribute::SetName(const std::string &name) {m_name = ComoUnicode::StringToWString(name);};
inline void Attribute::SetValue(const std::string &value) {m_value = ComoUnicode::StringToWString(value);};
inline void Attribute::SetName(const std::wstring &name) {m_name = name;};
inline void Attribute::SetValue(const std::wstring &value) {m_value = value;};
};

#endif
