/*
 * stringFunctions.h
 *
 *  Created on: Jun 24, 2011
 *      Author: Laurence R. McGlashan
 */

#ifndef STRINGFUNCTIONS_H_
#define STRINGFUNCTIONS_H_

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


    template <class T>
    T from_string
    (
        const std::string& s,
        std::ios_base& (*f)(std::ios_base&) = std::dec
    )
    {
      std::istringstream iss(s);
      T t;
      iss >> f >> t;
      return t;
    }

    //! Read a file into a std::string.
    const std::string fileToString(const std::string& fileName);

    /**
     * Private self-contained function. It does exactly as it is named.
     *
     * @param file
     * @return
     */
    const std::vector<std::string> fileToStrings(const std::string file);

    std::string regex_escape
    (
        const std::string& string_to_escape
    );

    std::string replaceComments(std::string stringToReplace);

    std::string convertToCaps(const std::string &str);

    std::string trim(const std::string &str);

} // namespace IO

#endif /* STRINGFUNCTIONS_H_ */
