/*
  Author(s):      Weerapong Phadungsukanan (wp214)
  Project:        sprog (gas-phase chemical kinetics).
  Sourceforge:    http://sourceforge.net/projects/mopssuite

  Copyright (C) 2008 Weerapong Phadungsukanan.

  File purpose:
    This is the main include file for other projects which use
    mops.  All components of sprog are included here.

  Licence:
    This file is part of "sprog".

    sprog is free software; you can redistribute it and/or
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
#include "gpc_string.h"

namespace Sprog {
    namespace IO {
        namespace StringFunc {
            ck_pos getKeyPos(const std::string &key, std::string &ckstr) {
                ck_pos pos;
                // convert to capital
                std::string chemkinCapstr = Strings::convertToCaps(ckstr);
                pos.begin = chemkinCapstr.find(key, 0);
                // try smaller keyword if not found then throw an exception
                if (pos.begin == std::string::npos) {
                    pos.begin = chemkinCapstr.find(key.substr(0,4), 0);
                }
                pos.end = chemkinCapstr.find(END_KEYWORD, pos.begin);
                return pos;
            }

            std::string extract_CK_elements_str(std::string &ckstr) {
                // get position of all keywords
                ck_pos el_pos, sp_pos, rt_pos, tm_pos;
                el_pos = getKeyPos(EL_KEYWORD, ckstr);
                if (el_pos.begin == std::string::npos) throw std::invalid_argument("ELEM or ELEMENTS keyword not found in chemkin input file.");
                if (el_pos.end   == std::string::npos) throw std::invalid_argument("END keyword missing after ELEM or ELEMENTS keywords in chemkin input file.");
                // copy only part of element string from chemkin string
                std::string ck_el_str = ckstr.substr(el_pos.begin, el_pos.end - el_pos.begin + END_KEYWORD.length());
                // check if other keywords are found in the chemkin element string
                el_pos = getKeyPos(EL_KEYWORD, ck_el_str);
                sp_pos = getKeyPos(SP_KEYWORD, ck_el_str);
                tm_pos = getKeyPos(TM_KEYWORD, ck_el_str);
                rt_pos = getKeyPos(RT_KEYWORD, ck_el_str);
                if ((el_pos.begin != 0) ||
                    (sp_pos.begin != std::string::npos) ||
                    (tm_pos.begin != std::string::npos) ||
                    (rt_pos.begin != std::string::npos)) {
                        // found other keywords
                        throw std::invalid_argument("Keywords found within ELEMENTS definition.");
                }
                return ck_el_str;
            }

            std::string extract_CK_species_str(std::string &ckstr) {
                // get position of all keywords
                ck_pos el_pos, sp_pos, rt_pos, tm_pos;
                sp_pos = getKeyPos(SP_KEYWORD, ckstr);
                if (sp_pos.begin == std::string::npos) throw std::invalid_argument("SPEC or SPECIES keyword not found in chemkin input file.");
                if (sp_pos.end   == std::string::npos) throw std::invalid_argument("END keyword missing after SPEC or SPECIES keywords in chemkin input file.");
                // copy only part of species string from chemkin string
                std::string ck_sp_str = ckstr.substr(sp_pos.begin, sp_pos.end - sp_pos.begin + END_KEYWORD.length());
                // check if other keywords are found in the chemkin species string
                el_pos = getKeyPos(EL_KEYWORD, ck_sp_str);
                sp_pos = getKeyPos(SP_KEYWORD, ck_sp_str);
                tm_pos = getKeyPos(TM_KEYWORD, ck_sp_str);
                rt_pos = getKeyPos(RT_KEYWORD, ck_sp_str);
                if ((el_pos.begin != std::string::npos) ||
                    (sp_pos.begin != 0) ||
                    (tm_pos.begin != std::string::npos) ||
                    (rt_pos.begin != std::string::npos)) {
                        // found other keywords
                        throw std::invalid_argument("Keywords found within SPECIES definition.");
                }
                return ck_sp_str;
            }

            std::string extract_CK_reactions_str(std::string &ckstr) {
                // get position of all keywords
                ck_pos el_pos, sp_pos, rt_pos, tm_pos;
                rt_pos = getKeyPos(RT_KEYWORD, ckstr);
                if (rt_pos.begin == std::string::npos) throw std::invalid_argument("REAC or REACTIONS keyword not found in chemkin input file.");
                if (rt_pos.end   == std::string::npos) throw std::invalid_argument("END keyword missing after REAC or REACTIONS keywords in chemkin input file.");
                // copy only part of reactions string from chemkin string
                std::string ck_rt_str = ckstr.substr(rt_pos.begin, rt_pos.end - rt_pos.begin + END_KEYWORD.length());
                // check if other keywords are found in the chemkin reactions string
                el_pos = getKeyPos(EL_KEYWORD, ck_rt_str);
                sp_pos = getKeyPos(SP_KEYWORD, ck_rt_str);
                tm_pos = getKeyPos(TM_KEYWORD, ck_rt_str);
                rt_pos = getKeyPos(RT_KEYWORD, ck_rt_str);
                if ((el_pos.begin != std::string::npos) ||
                    (sp_pos.begin != std::string::npos) ||
                    (tm_pos.begin != std::string::npos) ||
                    (rt_pos.begin != 0)) {
                        // found other keywords
                        throw std::invalid_argument("Keywords found within REACTIONS definition.");
                }
                return ck_rt_str;
            }

            std::string extract_CK_thermo_str(std::string &ckstr) {
                // get position of all keywords
                ck_pos el_pos, sp_pos, rt_pos, tm_pos;
                tm_pos = getKeyPos(RT_KEYWORD, ckstr);
                if (tm_pos.begin == std::string::npos) throw std::invalid_argument("THER or THERMO keyword not found in chemkin input file.");
                if (tm_pos.end   == std::string::npos) throw std::invalid_argument("END keyword missing after THER or THERMO keywords in chemkin input file.");
                // copy only part of thermo string from chemkin string
                std::string ck_tm_str = ckstr.substr(tm_pos.begin, tm_pos.end - tm_pos.begin + END_KEYWORD.length());
                // check if other keywords are found in the chemkin reactions string
                el_pos = getKeyPos(EL_KEYWORD, ck_tm_str);
                sp_pos = getKeyPos(SP_KEYWORD, ck_tm_str);
                tm_pos = getKeyPos(TM_KEYWORD, ck_tm_str);
                rt_pos = getKeyPos(RT_KEYWORD, ck_tm_str);
                if ((el_pos.begin != std::string::npos) ||
                    (sp_pos.begin != std::string::npos) ||
                    (tm_pos.begin != 0) ||
                    (rt_pos.begin != std::string::npos)) {
                        // found other keywords
                        throw std::invalid_argument("Keywords found within THERMO definition.");
                }
                return ck_tm_str;
            }

            std::string CK_is2str(std::istream &isin) {
                // Current position of get pointer in istream
                int current_pos = isin.tellg();
                // Reset position of get pointer to beginning of stream
                isin.clear();
                isin.seekg (0, std::ios::beg);

                // Convert chemkin file stream to a string and store in chemkinstr
                char c;
                std::string chemkinstr;
                // must try to get the frist character in the stream first because the stream might be empty
                isin.get(c);
                while (isin.good()) {
                    // skiping comment until new line or return character is found
                    if  (c=='!') {
                        do {
                            isin.get(c);
                        } while (isin.good() && (c!='\r') && (c!='\n'));
                        if (!isin.good()) break;
                    }
                    if (c=='\r') {
                        // replace return character by new line character, if some system doesn't read \r\n as \n
                        c = '\n';
                    } else if (c=='\t') {
                        // replace tab character by a space bar character
                        c = ' ';
                    }
                    chemkinstr.append(1,c);
                    // get next character
                    isin.get(c);
                };

                // Reset position of get pointer back to normal
                isin.clear();
                isin.seekg (current_pos);
                return chemkinstr;
            }

            void remove_CK_keyword(std::string &ckstr, const std::string &key) {
                size_t pos = ckstr.find(key);
                if (pos != std::string::npos) {
                    ckstr.erase(pos,key.length());
                } else {
                    pos = ckstr.find(key.substr(0,4));
                    if (pos != std::string::npos) {
                        ckstr.erase(pos,key.substr(0,4).length());
                    }
                }
            }
        }
    }
}
