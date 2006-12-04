/*
  Author(s):      Matthew Celnik (msc37)
  Project:        camxml (Cambridge c++ XML library)

  File purpose:
    Defines an XML attribute.
*/

#ifndef CAMXML_ATTRIBUTE_H
#define CAMXML_ATTRIBUTE_H

#include <vector>
#include <string>

namespace CamXML
{
class Attribute
{
private:
    /* Attribute name (bit before '=' sign). */
    std::string m_name;
    /* Attribute value (bit after '=' sign). */
    std::string m_value;
public:
    Attribute(void);
    Attribute(const std::string &name, const std::string &value);
    ~Attribute(void);
public:
    /* Returns this attribute's name. */
    const std::string &GetName(void) const;
    /* Returns this attribute's value. */
    const std::string &GetValue(void) const;
public:
    /* Sets the attribute name. */
    void SetName(const std::string &name);
    /* Sets the attribute value. */
    void SetValue(const std::string &value);
};

/* Inline function definitions. */
inline const std::string &Attribute::GetName() const {return m_name;};
inline const std::string &Attribute::GetValue() const {return m_value;};
inline void Attribute::SetName(const std::string &name) {m_name = name;};
inline void Attribute::SetValue(const std::string &value) {m_value = value;};
};

#endif
