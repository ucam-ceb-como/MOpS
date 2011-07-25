/*
 * stringFunctions.h
 *
 *  Created on: Jun 24, 2011
 *      Author: Laurence R. McGlashan
 */

#ifndef STRINGFUNCTIONS_H_
#define STRINGFUNCTIONS_H_

#include "boost/regex.hpp"
#include <string>
#include <iostream>
#include <vector>
#include <sstream>
#include <fstream>

namespace IO
{

    //! Output a vector.
    template <class T>
    std::ostream& operator<< (std::ostream& output, const std::vector<T>& v)
    {
        output << "[ Size = " << v.size() << "\n";
        for (typename std::vector<T>::const_iterator ii = v.begin(); ii != v.end(); ++ii)
        {
            output << *ii << std::endl;
        }
        output << "]";
        return output;
    }

    void checkNumberFormat(std::string& t);

    template <class T>
    T from_string
    (
        const std::string& s,
        std::ios_base& (*f)(std::ios_base&) = std::dec
    )
    {
      std::string check = s;
      checkNumberFormat(check);
      std::istringstream iss(check);
      T t;
      iss >> f >> t;
      return t;
    }

    //! Read a file into a std::string.
    std::string fileToString(const std::string& fileName);

    /**
     * Private self-contained function. It does exactly as it is named.
     *
     * @param file
     * @return
     */
    std::vector<std::string> fileToStrings(const std::string file);

    std::string regex_escape
    (
        const std::string& string_to_escape
    );

    std::string replaceComments(const std::string& stringToReplace);

    std::string convertToCaps(const std::string &str);
    std::vector<std::string> convertToCaps(const std::vector<std::string>& str);

} // namespace IO

#endif /* STRINGFUNCTIONS_H_ */
