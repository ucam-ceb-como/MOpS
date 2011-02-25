/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        camxml (Cambridge c++ XML library)
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    Contains functions to read, write, encode, decode
    and check unicode string and file. This is created to work
    for both Windows and Linux operating system. Due to the fact that
    wchar_t is 2 byte long on Windows and 4 byte long on Linux, all
    functions in this code are implemented for 2 bytes system. So it can
    encode/decode upto about first 65000 unicode coding point.

  Future Revisions:
    Adds utf16decode, utf16encode, uft32decode, uft32encode
    is_utf32BE, is_utf32LE.
    Modifies StreamToWString to automatically read above additional
    file type.
    Check the size of wchar_t and implement the encoder/decoder to
    work with any system automatically. Current work for both system but
    does not make use of 4 byte long wchar_t

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


#ifndef COMO_UNICODE_FUNCTIONS_H
#define COMO_UNICODE_FUNCTIONS_H

#include <string>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <stdio.h>

namespace ComoUnicode {
    // Check the first n bytes of istream for BOM header.
    bool is_valid_header(std::istream &is, const unsigned char BOM[], int n);
    // Check the first 3 bytes of istream for utf8 with BOM header.
    bool is_utf8wBOM(std::istream &is);
    // Check the first 3 bytes of file for utf8 with BOM header.
    bool is_utf8wBOM(const std::string &filename);
    // Check the first 3 bytes of file for utf16 with BE BOM header.
    bool is_utf16BE(std::istream &is);
    // Check the first 3 bytes of file for utf16 with BE BOM header.
    bool is_utf16BE(const std::string &filename);
    // Check the first 3 bytes of file for utf16 with LE BOM header.
    bool is_utf16LE(std::istream &is);
    // Check the first 3 bytes of file for utf16 with LE BOM header.
    bool is_utf16LE(const std::string &filename);
    // Decode mutibyte utf8 upto 3 bytes and stores coding point
    // of unicode in wstring. This works on both size of wchar_t (2 and 4 bytes).
    void utf8decode(std::wstring &dest, const std::string &src);
    // Encode Unicode coding point into multibyte uft8 and store in a
    // string ready for write to file.
    void utf8encode(std::string &dest, const std::wstring &src);
    // Convert std::string to std::wstring, no lost
    std::wstring StringToWString(const std::string &src);
    // Convert std::string to std::wstring, no lost
    void StringToWString(std::wstring &dest,const std::string &src);
    // Convert std::wstring to std::string, Unicode above 0xFF are discarded
    std::string WStringToString(const std::wstring &src);
    // Convert std::wstring to std::string, Unicode above 0xFF are discarded
    void WStringToString(std::string &dest,const std::wstring &src);
    // Write a wstring to UTF-8 file
    void WriteUTF8(const std::string &filename, const std::wstring &src);
    // Write a wstring to UTF-8 with BOM
    void WriteUTF8wBOM(const std::string &filename, const std::wstring &src);
    // Convert istream to string
    void is2str(std::string &dest, std::istream &is);
    // convert istream of any unicode to wstring
    void StreamToWString(std::istream &is, std::wstring &wstr);
    // Read any file format to wstring
    void ReadFile(const std::string &filename, std::wstring &wstr);
};

#endif
