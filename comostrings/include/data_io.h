/* 
 * File:   data_io.h
 * Author: vj231
 *
 * Created on 12 February 2009, 16:04
 */

#ifndef _DATA_IO_H
#define	_DATA_IO_H
#include <string>
#include <vector>
#include <fstream>
#include <ios>

class DataIO{
    
    std::fstream dataStream;

public:
    DataIO(){}
    ~DataIO(){}
    DataIO(const std::string &name, bool old=false);


    void open(const std::string &name, bool old);
    void close();

    void write(const std::vector<std::string>& data);
    void write(const std::vector<float>& data);
    void write(const std::vector<double>& data);
    void write(const std::vector<long double> &data);
    
    template<class T>
    void writeData(const std::vector<T> &data);

};

#endif	/* _DATA_IO_H */

