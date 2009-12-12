
#include <ios>


#include "data_io.h"
#include "string_functions.h"


//parameterised constructor
DataIO::DataIO(const std::string& name, bool old){
    
    open(name, old);
}

void DataIO::close(){
    
    dataStream.close();
}
//open the file

void DataIO::open(const std::string& name, bool old){

    if (old){
        dataStream.open(name.c_str(),std::ios_base::out | std::ios_base::in | std::ios_base::app);
    }else{
        dataStream.open(name.c_str(),std::ios_base::out | std::ios_base::in | std::ios_base::trunc);
    }

    
}

void DataIO::write(const std::vector<std::string>& data){

    writeData<std::string>(data);
    //writeHeader(data);
}

void DataIO::write(const std::vector<float>& data){
    writeData<float>(data);
}

void DataIO::write(const std::vector<double>& data){
    writeData<double>(data);
}

void DataIO::write(const std::vector<long double>& data){
    writeData<long double>(data);
}


template<class T>
void DataIO::writeData(const std::vector<T>& data)
{


    //dataStream.setf(ios::left);
    dataStream.setf(std::ios::scientific);
    if(dataStream.is_open() && dataStream.good() && (data.size()>0) ){
        typename std::vector<T>::const_iterator p = data.begin();
        //dataStream << (*p++);
        while(p!=data.end()){
            
            dataStream.width(12);
            dataStream.setf(std::ios::left);;
            dataStream <<  (*p++) << "  ";
        }
        dataStream << "\n";
        dataStream.flush();
    }

}

