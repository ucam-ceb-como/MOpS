#include "xmldocument.h"
#include "xmlelement.h"
#include <fstream>
#include <string>

using namespace CamXML;
using namespace std;

Document::Document(void)
{
    m_root = NULL;
}

Document::~Document(void)
{
    if (m_root!=NULL) delete m_root;
}

int Document::Load(const std::string &filename)
{
    // Open file for reading.
    ifstream fin(filename.c_str(), ios::in);

    if (fin.good()) {
        STATUS st = Outside;

        try {
            m_root = parseElement(fin, st);
        } catch (exception &e) {
            fin.close();
            throw e;
        }
        
        fin.close();
        return 0;
    } else {
        // Failed to open file.
        throw invalid_argument(string("CamXML could not open file: ").append(filename));
    }
}

Element *const Document::parseElement(std::ifstream &fin, STATUS &st)
{
    string::size_type pos=0;
    char c;
    string tag, data, comment;
    Element *el = NULL;
    bool iscomment=false;

    while((st!=End) && (st!=Fail) && (fin.good())) {
        fin.get(c);

        switch (st) {
            case Outside:
                if (c=='<') {
                    // We are now beginning to read a tag.
                    st = Begin;
                }
                break;
            case Begin:
                if (c=='!') {
                    // This begins a comment (maybe).
                    st = ReadComment1;
                } else if (c=='/') {
                    // This is a closing tag.
                    st = ReadEndTag1;
                } else if (c=='?') {
                    st = ReadDocInfo;
                } else if (isLetter(c)) {
                    // Reading element tag name.
                    el = new Element();
                    st = ReadTag;
                    tag = c;
                } else {
                    // Something has gone wrong.
                    st = Fail;
                    if (el!=NULL) delete el;
                    throw range_error("Invalid first character in opening tag.");
                }
                break;
            case ReadDocInfo:
                if (c=='?') {
                    st = CloseDocInfo;
                }
                break;


            /* Reading tag name. */

            case ReadTag:
                if (isLetterOrNum(c) || (c=='-') || (c=='_') || (c==':')) {
                    // Valid tag name character.
                    tag.append(&c,1);
                } else if (c=='/') {
                    // Empty element, now closing.
                    st = CloseEmptyElement;
                } else if (c=='>') {
                    // Now read element data, and perhaps sub-elements.
                    st = ReadData;
                    data.clear();
                } else if (isWhiteSpace(c)) {
                    // End of tag, now read attributes.
                    st = ReadAttr;
                } else {
                    // Invalid tag name character.
                    st = Fail;
                    if (el!=NULL) delete el;
                    throw range_error("Invalid character in opening tag name.");
                }
                break;


            /* Reading closing tags. */

            case ReadEndTag1:
                if (isLetter(c)) {
                    // Valid tag name first character.
                    st = ReadEndTag;
                } else {
                    // Invalid tag name character.
                    st = Fail;
                    if (el!=NULL) delete el;
                    throw range_error("Invalid first character in closing tag.");
                }
                break;

            case ReadEndTag:
                if (isLetterOrNum(c) || (c=='-') || (c=='_') || (c==':')) {
                } else if (c=='>') {
                    // That's the end!
                    st = End;
                } else if (isWhiteSpace(c)) {
                    // Now waiting for end.
                    st = CloseElement;
                } else {
                    // Invalid tag name character.
                    st = Fail;
                    if (el!=NULL) delete el;
                    throw range_error("Invalid character in closing tag name.");
                }
                break;


            /* Reading element data and sub-elements. */

            case ReadData:
                if (c=='<') {
                    // Either opening a sub-element or closing the
                    // current one.  Better check.
                    fin.get(c);
                    if (c=='/') {
                        // Closing current element.
                        st = ReadEndTag1;
                    } else {
                        // Parsing a child element.
                        fin.unget();
                        STATUS childstat = Begin;

                        try {
                            Element *child = parseElement(fin, childstat);
                            if ((child!=NULL) && (childstat != Fail)) {
                                el->AddChild(child);
                            } else if (childstat!=EndComment) {
                                st = Fail;
                                if (el!=NULL) delete el;
                                if (child!=NULL) delete child;
                                throw exception(string("Failed to read child of element: ").append(tag).c_str());
                            }
                        } catch (exception &e) {
                            throw e;
                        }
                    }
                } else {
                    // Add character to data.
                    data.append(&c,1);
                }
                break;


            /* Reading element attributes. */

            case ReadAttr:
                if (c=='>') {
                    // Finished reading opening tag, now read data.
                    st = ReadData;
                } else if (c=='/') {
                    // Closing an empty element.
                    st = CloseEmptyElement;
                } else if (isLetter(c)) {
                    // Valid start character of attribute name.
                    fin.unget();
                    STATUS attrstat = Begin;

                    try {
                        Attribute *attr = parseAttr(fin, attrstat);
                        if ((attr!=NULL) && (attrstat!=Fail)) {
                            el->SetAttribute(attr->GetName(), attr->GetValue());
                            delete attr;
                        } else {
                            st = Fail;
                            if (el!=NULL) delete el;
                            throw exception(string("Failed to read attribute of element: ").append(tag).c_str());
                        }
                    } catch (exception &e) {
                        throw e;
                    }

                } else if (isWhiteSpace(c)) {
                    // White space is ok between attributes.
                } else {
                    // Invalid attribute character!
                    st = Fail;
                    if (el!=NULL) delete el;
                    throw range_error(string("Invalid first character of attribute name of element: ").append(tag));
                }
                break;


            /* Closing the element. */

            case CloseElement:
                if (c=='>') {
                    // This is the end.
                    st = End;
                }
                break;
            case CloseEmptyElement:
                // Expecting ">" immediately.
                if (c=='>') {
                    st = End;
                } else {
                    st = Fail;
                    if (el!=NULL) delete el;
                    throw range_error(string("Expected '>' not found when closing empty element: ").append(tag));
                }
                break;
            case CloseDocInfo:
                // Expecting ">" immediately.
                if (c=='>') {
                    // That was just the document info, now try
                    // to read an element.
                    st = Outside;
                } else {
                    st = Fail;
                    if (el!=NULL) delete el;
                    throw range_error(string("Expected '>' not found when closing doc info."));
                }
                break;

            /* Reading XML comments. */

            case ReadComment1:
                // Expect a "-" next.
                if (c=='-') {
                    st = ReadComment2;
                } else {
                    st = Fail;
                    if (el!=NULL) delete el;
                    throw range_error("Expected first '-' not found when opening a comment.");
                }
                break;
            case ReadComment2:
                // Expect second "-" next.
                if (c=='-') {
                    st = ReadComment;
                    data.clear();
                    iscomment = true;
                    delete el;
                } else {
                    st = Fail;
                    if (el!=NULL) delete el;
                    throw range_error("Expected second '-' not found when opening a comment.");
                }
                break;
            case ReadComment:
                if (c=='-') {
                    if (data[data.length()-1]=='-') {
                        // Check that next character is ">", otherwise
                        // this is an error.
                        if (fin.peek()!='>') {
                            st = Fail;
                            if (el!=NULL) delete el;
                            throw range_error("Expected '>' not found when closing a comment.");
                        }
                    }
                    data.append(&c,1);
                } else if (c=='>') {
                    if ((data.length()>1) && (data.substr(data.length()-2,2)=="--")) {
                        data.erase(data.length()-2,2);
                        st = End;
                    } else {
                        data.append(&c);
                    }
                } else {
                    data.append(&c,1);
                }
                break;
        }
    }

    if (st==Fail) {
        if (el!=NULL) delete el;
        throw exception("Unhandled error occured.");
    } else if ((st==End) && !iscomment) {
        el->SetTag(tag);
        el->SetData(data);
    } else if ((st==End) && iscomment) {
        if (el!=NULL) delete el;
        el = NULL;
        st = EndComment;
    }
    return el;
}

Attribute *const Document::parseAttr(std::ifstream &fin, CamXML::Document::STATUS &st)
{
    Attribute *a=NULL;
    string name="", data="";
    bool waitQuote=true;
    char c;

    while((st!=End) && (st!=Fail) && (fin.good())) {
        fin.get(c);
        switch (st) {
            case Begin:
                if (isLetter(c)) {
                    st = ReadAttr;
                    a = new Attribute();
                    name.append(&c,1);
                }
                break;
            case ReadAttr:
                if (c=='=') {
                    st = ReadData;
                    waitQuote = true;
                } else {
                    name.append(&c,1);
                }
                break;
            case ReadData:
                if (waitQuote) {
                    if (c=='"') {
                        waitQuote = false;
                    } else {
                        // Quote mark missing after equals sign.
                        st = Fail;
                    }
                } else {
                    if (c=='"') {
                        // We have read all the data for this attribute.
                        st = End;
                    } else {
                        data.append(&c,1);
                    }
                }
                break;
        }
    }

    if ((st!=Fail)) {
        a->SetName(name);
        a->SetValue(data);
    } else {
        if (a!=NULL) delete a;
        a = NULL;
    }        
    return a;
}