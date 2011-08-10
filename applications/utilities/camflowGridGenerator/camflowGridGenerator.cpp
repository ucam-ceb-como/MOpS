/*! Little utility to generate a camflow grid.
 *  Laurence McGlashan (lrm29@cam.ac.uk)
 */

#include "boost/program_options.hpp"
#include <boost/math/distributions/triangular.hpp>
#include <boost/math/distributions/normal.hpp>
#include <fstream>
#include <iomanip>
#include <numeric>

using namespace boost::program_options;
using namespace boost::math;
using namespace std;

int main(int argc, char *argv[])
{

    cout << argv[0] << " by Laurence McGlashan\n" << endl;

    string outputFile = "grid.inp";
    string distType;
    double stMixFrac = -1;
    int numberOfCells, numberOfCellsLower, numberOfCellsUpper = -1;

    {
        // Parse arguments.
        options_description desc("Allowed options for program");
        desc.add_options()
            ("help", "Show this help message.")
            ("stoich", value<double>(), "Stoichiometric mixture fraction.")
            ("outputFile", value<string>(), "Output file for grid.")
            ("cells", value<int>(), "Number of cells.")
            ("distribution", value<string>(), "Type of distribution to weight points against (triangular|normal).");

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

        if (vm.count("distribution")) {
            distType = vm["distribution"].as<string>();
            if (distType != "triangular" && distType != "normal")
                    throw std::logic_error("Distribution " + distType + " not available.\n");
            cout<< "Distribution used will be " << distType << ".\n";
        } else {
            throw std::logic_error("Distribution was not provided.\n");
        }

        if (vm.count("outputFile")) {
            outputFile = vm["outputFile"].as<string>();
            cout<< "Grid will be output to "
                << outputFile << ".\n";
        } else {
            cout<< "outputFile was not provided. Use default of grid.inp.\n";
        }
    }

    // Calculate the grid here. The values are the cell edges.
    vector<double> grid;
    grid.push_back(0.0);
    grid.push_back(stMixFrac);
    grid.push_back(1.0);

    numberOfCellsLower = numberOfCells/2.0;
    numberOfCellsUpper = numberOfCells/2.0;

    //////////// Generate grid here ///////////

    // Construct uniform grid.
    for (size_t i=1; i<numberOfCellsLower; ++i)
    {
            grid.push_back(2.0*i*stMixFrac/numberOfCells);
    }
    for (size_t i=1; i<numberOfCellsUpper; ++i)
    {
            grid.push_back(stMixFrac + 2.0*i*(1.0-stMixFrac)/numberOfCells);
    }

    // Sort the values.
    sort(grid.begin(), grid.end());

    if (distType == "triangular")
    {
        triangular s(0.0,stMixFrac,1.0);
        double max = pdf(s,stMixFrac);
        for (size_t i=3; i<=numberOfCells; ++i)
        {
            grid[i] = stMixFrac + (grid[i]-stMixFrac)*(1.0-pdf(s,grid[i])/max);
        }
    }

    if (distType == "normal")
    {
        normal s(stMixFrac,0.5);
        vector<double> spacingFactor, pdfSaved;
        
        for (size_t i=0; i<=numberOfCells; ++i) 
            pdfSaved.push_back(pdf(s,grid[i]));
            
        for (size_t i=0; i<numberOfCells; ++i)
            spacingFactor.push_back(abs(pdfSaved[i+1] - pdfSaved[i]));

        double normalise = accumulate(spacingFactor.begin(),spacingFactor.end(),0.0);
        for (size_t i=1; i<numberOfCells; ++i)
            grid[i] = grid[i-1] + abs(pdfSaved[i]-pdfSaved[i-1])/normalise;
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
        }
    }

    cout << "\nProgram End." << endl;

    return 0;

}

