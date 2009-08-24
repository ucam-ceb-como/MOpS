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
using namespace std;
class DataIO{
    
    fstream dataStream;

public:
    DataIO(){}
    ~DataIO(){}
    DataIO(const string &name, bool old=false);


    void open(const string &name, bool old);
    void close();

    void write(const vector<string>& data);
    void write(const vector<float>& data);
    void write(const vector<double>& data);
    void write(const vector<long double> &data);
    
    template<class T>
    void writeData(const vector<T> &data);

};

#endif	/* _DATA_IO_H */

