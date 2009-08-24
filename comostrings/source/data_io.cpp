
#include <ios>


#include "data_io.h"
#include "string_functions.h"


//parameterised constructor
DataIO::DataIO(const string& name, bool old){
    
    open(name, old);
}

void DataIO::close(){
    
    dataStream.close();
}
//open the file

void DataIO::open(const string& name, bool old){

    if (old){
        dataStream.open(name.c_str(),ios_base::out | ios_base::in | ios_base::app);
    }else{
        dataStream.open(name.c_str(),ios_base::out | ios_base::in | ios_base::trunc);
    }

    
}

void DataIO::write(const vector<string>& data){

    writeData<string>(data);
    //writeHeader(data);
}

void DataIO::write(const vector<float>& data){
    writeData<float>(data);
}

void DataIO::write(const vector<double>& data){
    writeData<double>(data);
}

void DataIO::write(const vector<long double>& data){
    writeData<long double>(data);
}


template<class T>
void DataIO::writeData(const vector<T>& data)
{


    //dataStream.setf(ios::left);
    dataStream.setf(ios::scientific);
    if(dataStream.is_open() && dataStream.good() && (data.size()>0) ){
        typename vector<T>::const_iterator p = data.begin();
        //dataStream << (*p++);
        while(p!=data.end()){
            
            dataStream.width(12);
            dataStream.setf(ios::left);;
            dataStream <<  (*p++) << "  ";
        }
        dataStream << "\n";
        dataStream.flush();
    }

}

