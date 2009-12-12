/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        camxml (Cambridge c++ XML library)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    Implementation of the functions declared in the unicodef.h
    header file.

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

#include "unicodef.h"
#include <stdexcept>

namespace ComoUnicode {

    bool is_valid_header(std::istream &is, const unsigned char BOM[], int n) {
        // Check if n < 1, return true as there is no checking
        if (n < 1) { return true; }
        // Reading the first bufferSize characters
        const int bufferSize = n;
        // new char Buffer
        char * buffer = new char [bufferSize];
        // Bom status
        bool is_bom = false;
        // Current position of get pointer in istream
        int current_pos = is.tellg();
        // Reset position of get pointer
        is.seekg (0, std::ios::beg);
        // Read the first bufferSize bytes from file
        is.read((char*) buffer, bufferSize*sizeof(char));
        if (is.gcount() == bufferSize ) {
            is_bom = (static_cast<unsigned char>(0xFF & buffer[0]) == BOM[0]);
            for (int i = 1; i < bufferSize; i++) {
                is_bom = is_bom && (static_cast<unsigned char>(0xFF & buffer[i]) == BOM[i]);
            }
        }
        // Reset position of get pointer back to normal
        is.seekg (current_pos);
        // free buffer
        delete [] buffer;
        return is_bom;
    }

    bool is_utf8wBOM(std::istream &is) {
        // Define BOM charaters
        const unsigned char BOM[] = {0xEF, 0xBB, 0xBF};
        return is_valid_header(is, BOM, 3);
    }

    bool is_utf8wBOM(const std::string &filename) {
        // Define BOM charaters
        std::ifstream is(filename.c_str());
        bool is_bom = false;
        // Check if file is openned
        if (is.good()) {
            is_bom = is_utf8wBOM(is);
        }
        else {
            // throw exception when file could not be openned
            throw std::invalid_argument(std::string("CamXML could not open file: ").append(filename));
        }
        is.close();
        return is_bom;
    }

    bool is_utf16BE(std::istream &is) {
        // Define BOM charaters
        const unsigned char BOM[] = {0xFE, 0xFF};
        return is_valid_header(is, BOM, 2);
    }

    bool is_utf16BE(const std::string &filename) {
        // Define BOM charaters
        std::ifstream is(filename.c_str());
        bool is_bom = false;
        // Check if file is openned
        if (is.good()) {
            is_bom = is_utf16BE(is);
        }
        else {
            // throw exception when file could not be openned
            throw std::invalid_argument(std::string("CamXML could not open file: ").append(filename));
        }
        is.close();
        return is_bom;
    }

    bool is_utf16LE(std::istream &is) {
        // Define BOM charaters
        const unsigned char BOM[] = {0xFF, 0xFE};
        return is_valid_header(is, BOM, 2);
    }

    bool is_utf16LE(const std::string &filename) {
        // Define BOM charaters
        std::ifstream is(filename.c_str());
        bool is_bom = false;
        // Check if file is openned
        if (is.good()) {
            is_bom = is_utf16LE(is);
        }
        else {
            // throw exception when file could not be openned
            throw std::invalid_argument(std::string("CamXML could not open file: ").append(filename));
        }
        is.close();
        return is_bom;
    }

    void utf8decode(std::wstring &dest, const std::string &src) {
	    size_t i = 0;
	    unsigned char *s = (unsigned char *) src.c_str();

	    while (i < src.size())
	    {
		    const wchar_t c = s[i++];
		    // Decode U-0 to U-7F 
		    if ((c&0x80) == 0x00)
		    {
			    dest += c;
			    continue;
		    }
		    // Decode U-80 to U-7FF
		    if ((c&0xE0) == 0xC0)
		    {
			    if (i<src.size())
			    {
				    const wchar_t d = s[i++];
				    dest += (c&0x1f)<<6 | (d&0x3f);
				    continue;
			    }
		    }
		    // Decode U-800 to U-FFFF
		    if ((c&0xF0) == 0xE0)
		    {
			    if (i+1<src.size())
			    {
				    const wchar_t d = s[i++];
				    const wchar_t e = s[i++];
				    dest += (c&0x0f)<<12 | (d&0x3f)<<6 | (e&0x3f);
				    continue;
			    }
		    }
	    }
    }

    void utf8encode(std::string &dest, const std::wstring &src) {
	    wchar_t *s = (wchar_t *) src.c_str();

	    for (size_t i = 0; i < src.length(); i++)
	    {
		    const wchar_t c = s[i];
		    // Encode U-0 to U-7F 
		    if (c < 0x80)
		    {
			    dest += static_cast<unsigned char> (c);
			    continue;
		    }
		    // Decode U-80 to U-7FF
		    if ((c >= 0x80)&&(c < 0x800))
		    {
			    //if (i<src.size())
			    {
				    dest += static_cast<unsigned char>((c >> 6)|0xc0);
				    dest += static_cast<unsigned char>((c&0x3f)|0x80);
				    continue;
			    }
		    }
		    // Decode U-800 to U-FFFF
		    if (((c >= 0x800)&&(c < 0xD800))||
			    ((c >= 0xE000)&&(c < 0x10000)))
		    {
			    //if (i+1<src.size())
			    {
				    dest += static_cast<unsigned char>((c >> 12)|0xe0);
				    dest += static_cast<unsigned char>(((c&0xfff) >> 6)|0x80);
				    dest += static_cast<unsigned char>((c&0x3f)|0x80);
				    continue;
			    }
		    }
	    }
    }
    
    std::wstring StringToWString(const std::string &src) {
        std::wstring dest;
	    dest.resize(src.size());
	    for (size_t i=0; i<src.size(); i++)
		    dest[i] = static_cast<unsigned char>(src[i]);
        return dest;
    }

    void StringToWString(std::wstring &dest,const std::string &src) {
	    dest.resize(src.size());
	    for (size_t i=0; i<src.size(); i++)
		    dest[i] = static_cast<unsigned char>(src[i]);
    }
    
    std::string WStringToString(const std::wstring &src){
        std::string dest;
	    dest.resize(src.size());
	    for (size_t i=0; i<src.size(); i++)
		    // Ternary Conditional
		    // if src[i] < 256 then src[i] else '?'
		    dest[i] = src[i] < 256 ? src[i] : '?';
        return dest;
    }

    void WStringToString(std::string &dest,const std::wstring &src){
	    dest.resize(src.size());
	    for (size_t i=0; i<src.size(); i++)
		    // Ternary Conditional
		    // if src[i] < 256 then src[i] else '?'
		    dest[i] = src[i] < 256 ? src[i] : '?';
    }

    void WriteUTF8(const std::string &filename, const std::wstring &src) {
        std::ofstream fout(filename.c_str());
        if (fout.good()) {
            try {
                std::string s2write;
                utf8encode(s2write, src);
                fout << s2write;
            }
            catch (std::exception &e) {
                fout.close();
                throw;
            }
            fout.close();
        }
        else {
            // throw exception when file could not be openned to write
            throw std::invalid_argument(std::string("CamXML could not open file: ").append(filename));
        }
    }

    void WriteUTF8wBOM(const std::string &filename, const std::wstring &src) {
        std::ofstream fout(filename.c_str());
        if (fout.good()) {
            try {
                std::string s2write;
                const unsigned char BOM[] = {0xEF, 0xBB, 0xBF};
                utf8encode(s2write, src);
                fout << BOM;
                fout << s2write;
            }
            catch (std::exception &e) {
                fout.close();
                throw;
            }
            fout.close();
        }
        else {
            // throw exception when file could not be openned to write
            throw std::invalid_argument(std::string("CamXML could not open file: ").append(filename));
        }
    }

    void is2str(std::string &dest, std::istream &is) {
        // Current position of get pointer in istream
        int current_pos = is.tellg();
        // Reset position of get pointer to beginning of stream
        is.seekg (0, std::ios::beg);
        const int bufferSize = 1024;
        char buffer[bufferSize];
        while (is.good()) {
            is.read(buffer,bufferSize*sizeof(char));
            int count = is.gcount();
            dest.append(buffer,count);
        }
        // Reset position of get pointer back to normal
        is.seekg (current_pos);
    }

    void StreamToWString(std::istream &is, std::wstring &wstr) {
        // Current position of get pointer in istream
        int current_pos = is.tellg();
        // Reset position of get pointer to beginning of stream
        is.seekg (0, std::ios::beg);
            if (is_utf16BE(is)) {
                // Not yet support
                throw std::invalid_argument(std::string("UTF-16BE input detected and not yet support by CamXML."));
            }
            else if (is_utf16LE(is)) {
                // Not yet support
                throw std::invalid_argument(std::string("UTF-16LE input detected and not yet support by CamXML."));
            }
            else {
                // If cannot find particular format, use utf8
                // Read istream to std::string
                std::string str;
                str.clear();
                is2str(str, is);
                // remove BOM
                if (is_utf8wBOM(is))
                    str.erase(0,3);
                // Convert str ascii to unicode
                utf8decode(wstr, str);
            }
        // Reset position of get pointer back to normal
        is.seekg (current_pos);
    }

    void ReadFile(const std::string &filename, std::wstring &wstr) {
        std::ifstream is(filename.c_str());
        // If file is opened without fails
        if (is.good()) {
            StreamToWString(is, wstr);
            is.close();
        }
        else {
            // Failed to open file.
            throw std::invalid_argument(std::string("CamXML could not open file: ").append(filename));
        }
    }

}
