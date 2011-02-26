/*
  Author(s):      Matthew Celnik (msc37)
  Project:        comostrings_test

  File purpose:
    Main test routine for comostrings library.
*/

#include "../../../src/io/comostrings/include/csv_io.h"
#include "../../../src/io/comostrings/include/string_functions.h"
#include <string>
#include <vector>
#include <iostream>

using namespace std;

int main(void)
{
    // Open a csv file for writing.
    CSV_IO io;
    io.Open("test.csv", true);

    // Write a header row.
    vector<string> head;
    head.push_back("Time");
    head.push_back("V1");
    head.push_back("V2");
    io.Write(head);

    // Write a data row.
    vector<double> data;
    data.push_back(0.0);
    data.push_back(1.1110383453454543523453245e-12);
    data.push_back(2.2222);
    io.Write(data);
    io.Write(data);
    io.Write(data);
    io.Write(data);

    // Close the file.
    io.Close();

    // Now attempt to read the file.
    vector<string> headin;
    vector<long double> datain;
    io.Open("test.csv");
    io.Read(headin);
    io.Read(datain);
    io.Close();

    // Print the data.
    cout << headin.size() << ' ' << headin[0] << ' ' << 
            headin[1] << ' ' << headin[2] << '\n';
    cout << datain.size() << ' ' << datain[0] << ' ' << 
            datain[1] << ' ' << datain[2] << '\n';

    return 0;

}
