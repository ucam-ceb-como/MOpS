/*! Little utility to generate a camflow grid.
 *  Laurence McGlashan (lrm29@cam.ac.uk)
 */

#include "boost/program_options.hpp"
#include <boost/math/distributions/normal.hpp> // for normal_distribution

#include <fstream>
#include <iomanip>

using namespace boost::program_options;
using boost::math::normal;
using namespace std;

int main(int argc, char *argv[])
{

    cout << argv[0] << " by Laurence McGlashan\n" << endl;

    string outputFile = "grid.inp";
    double stMixFrac = -1;
    int numberOfCells = -1;

    {
        // Parse arguments.
        options_description desc("Allowed options for program");
        desc.add_options()
            ("help", "Show this help message.")
            ("stoich", value<double>(), "Stoichiometric mixture fraction.")
            ("outputFile", value<string>(), "Output file for grid.")
            ("cells", value<int>(), "Number of cells.");
            
        variables_map vm;
        store(parse_command_line(argc, argv, desc), vm);
        notify(vm);    

        if (vm.count("help")) {
            cout << desc << "\n";
            return 1;
        }

        if (vm.count("stoich")) {
            stMixFrac = vm["stoich"].as<double>();
            if (stMixFrac <= 0) 
                throw std::logic_error("Stochiometric mixture fraction must be greater than 0.\n");
            cout<< "Stochiometric mixture fraction is " << stMixFrac << ".\n";
        } else {
            throw std::logic_error("Stochiometric mixture fraction was not provided.\n");
        }
        
        if (vm.count("cells")) {
            numberOfCells = vm["cells"].as<int>();           
            if (numberOfCells <= 0) 
                throw std::logic_error("Number of cells must be greater than 2.\n");
            if (numberOfCells%2 != 0) 
                throw std::logic_error("Number of cells must be even.\n");
            cout<< "Number of cells is " << numberOfCells << ".\n";
        } else {
            throw std::logic_error("Number of cells was not provided.\n");
        }
        
        if (vm.count("outputFile")) {
            cout<< "Grid will be output to " 
                << vm["outputFile"].as<string>() << ".\n";
        } else {
            cout<< "outputFile was not provided. Use default of grid.inp.\n";
        }   
    }
    
    // Calculate the grid here. The values are the cell edges.
    vector<double> grid;
    grid.push_back(0.0);
    grid.push_back(stMixFrac);
    grid.push_back(1.0);
    
    // We will add this many points to the grid.
    int numberOfPoints = numberOfCells - 2;
    

    //////////// Generate grid here ///////////

    normal s(stMixFrac,2.0/3.0);
    for (size_t i=0; i<=numberOfCells; ++i)
    {
        cout << cdf(s, i*2.0/numberOfCells) << endl;
    }
    
    //////////// End Generate grid here ///////
    
    // Sort the values.
    sort(grid.begin(), grid.end());
    
    // Output the grid.
    ofstream out;
    out.open(outputFile.c_str(), ios::trunc);
    if (out.good())
    {
        for (size_t i=0; i<grid.size(); ++i)
        {
            out << setprecision(10) << grid[i] << endl;
            //cout << grid[i] << endl;
        }
    }
    
    
    
    cout << "\nProgram End." << endl;
    
    return 0;

}
