/*
  Author(s):      Matthew Celnik (msc37)
  Project:        camxml (Cambridge c++ XML library)

  File purpose:
    Defines an XML document.
*/

#ifndef CAMXML_DOCUMENT_H
#define CAMXML_DOCUMENT_H

#include "xmlelement.h"
#include "xmlattribute.h"
#include <string>
#include <fstream>

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
    /* Saves the XML data to file.  Overwrites whatever was previously
       in the file.  Returns 0 on success, otherwise negative. */
    int Save(const std::string &filename);
    /* Returns the root element of this document. */
    inline Element *const Root(void) {return m_root;};
private:
    Element *const parseElement(std::ifstream &fin, STATUS &stat);
    Attribute *const parseAttr(std::ifstream &fin, STATUS &stat);
    bool isWhiteSpace(const char c);
    bool isLetter(const char c);
    bool isLetterOrNum(const char c);
};

inline bool Document::isWhiteSpace(const char c)
{
    return (c==' ') || (c=='\n') || (c=='\t') || (c=='\r');
}

inline bool Document::isLetter(const char c)
{
    return (c=='a') || (c=='b') || (c=='c') || (c=='d') || (c=='e') ||
           (c=='f') || (c=='g') || (c=='h') || (c=='i') || (c=='j') ||
           (c=='k') || (c=='l') || (c=='m') || (c=='n') || (c=='o') ||
           (c=='p') || (c=='q') || (c=='r') || (c=='s') || (c=='t') ||
           (c=='u') || (c=='v') || (c=='w') || (c=='x') || (c=='y') ||
           (c=='z') ||
           (c=='A') || (c=='B') || (c=='C') || (c=='D') || (c=='E') ||
           (c=='F') || (c=='G') || (c=='H') || (c=='I') || (c=='J') ||
           (c=='K') || (c=='L') || (c=='M') || (c=='N') || (c=='O') ||
           (c=='P') || (c=='Q') || (c=='R') || (c=='S') || (c=='T') ||
           (c=='U') || (c=='V') || (c=='W') || (c=='X') || (c=='Y') ||
           (c=='Z');
}

inline bool Document::isLetterOrNum(const char c)
{
    return isLetter(c) ||
           (c=='0') || (c=='1') || (c=='2') || (c=='3') || (c=='4') ||
           (c=='5') || (c=='6') || (c=='7') || (c=='8') || (c=='9');
}

};

#endif