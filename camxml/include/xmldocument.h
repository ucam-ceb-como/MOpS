/*
  Author(s):      Matthew Celnik (msc37)
  Project:        camxml (Cambridge c++ XML library)

  File purpose:
    Defines an XML document.
*/

#ifndef CAMXML_DOCUMENT_H
#define CAMXML_DOCUMENT_H

#include "unicodef.h"
#include "xmlattribute.h"
#include "xmlelement.h"
#include <string>
#include <fstream>
#include <wchar.h>
#include <sstream>
#include <tchar.h>

using namespace CamXML;
using namespace std;

namespace CamXML
{
class Document
{
private:
    enum STATUS {Outside,Begin,End,ReadTag,ReadEndTag1,ReadEndTag,
                 ReadComment,ReadComment1,ReadComment2,
                 ReadAttr,ReadDocInfo,ReadData,CloseEmptyElement,
                 CloseElement,CloseDocInfo,EndComment,
                 Fail=-1};
private:
    /* Root element in document. */
    Element *m_root;
public:
    Document(void);
    ~Document(void);
public:
    /* Loads a document from file, building the data structure
       in the process. Returns 0 on success, otherwise negative. */
    int Load(const std::string &filename);
    /* Loads a document from file, building the data structure
       in the process. Returns 0 on success, otherwise negative. */
    int Load(const std::wstring &filename);
    /* Saves the XML data to file.  Overwrites whatever was previously
       in the file.  Returns 0 on success, otherwise negative.  (wide string version) */
    int Save(const std::string &filename);
    /* Returns the root element of this document. */
    int Save(const std::wstring &filename);
    /* Returns the root element of this document. (wide string version) */
    inline Element *const Root(void) {return m_root;};
    inline Element *const SetRoot(Element &root) {m_root->~Element(); m_root = &root; return m_root;};
private:
    Element *const parseElement(std::wistringstream &fin, STATUS &st);
    Attribute *const parseAttr(std::wistringstream &fin, STATUS &st);
    bool isWhiteSpace(const wchar_t c);
    bool isLetter(const wchar_t c);
    bool isLetterOrNum(const wchar_t c);
};

inline bool Document::isLetterOrNum(const wchar_t c)
{
    return isLetter(c) ||
           (c==L'0') || (c==L'1') || (c==L'2') || (c==L'3') || (c==L'4') ||
           (c==L'5') || (c==L'6') || (c==L'7') || (c==L'8') || (c==L'9');
}
inline bool Document::isWhiteSpace(const wchar_t c)
{
    return (c==L' ') || (c==L'\n') || (c==L'\r') || (c==L'\t');
}

inline bool Document::isLetter(const wchar_t c)
{
    return (c==L'a') || (c==L'b') || (c==L'c') || (c==L'd') || (c==L'e') ||
           (c==L'f') || (c==L'g') || (c==L'h') || (c==L'i') || (c==L'j') ||
           (c==L'k') || (c==L'l') || (c==L'm') || (c==L'n') || (c==L'o') ||
           (c==L'p') || (c==L'q') || (c==L'r') || (c==L's') || (c==L't') ||
           (c==L'u') || (c==L'v') || (c==L'w') || (c==L'x') || (c==L'y') ||
           (c==L'z') ||
           (c==L'A') || (c==L'B') || (c==L'C') || (c==L'D') || (c==L'E') ||
           (c==L'F') || (c==L'G') || (c==L'H') || (c==L'I') || (c==L'J') ||
           (c==L'K') || (c==L'L') || (c==L'M') || (c==L'N') || (c==L'O') ||
           (c==L'P') || (c==L'Q') || (c==L'R') || (c==L'S') || (c==L'T') ||
           (c==L'U') || (c==L'V') || (c==L'W') || (c==L'X') || (c==L'Y') ||
           (c==L'Z');
}


};

#endif
