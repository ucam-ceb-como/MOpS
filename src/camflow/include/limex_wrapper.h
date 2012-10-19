#ifndef LIMEX_WRAPPER_H
#define	LIMEX_WRAPPER_H

#include <vector>
#include "cam_residual.h"




//using namespace Cratl;



namespace Gadgets{
    class LimexWrapper{
    
    std::vector<double> fcn;          //Residual function
    std::vector<double> jacobian;     //Jacobian of residual function

    Camflow::CamResidual *reactorPtr;
    int eqnSize;
    double y;
    double resNorm;                    //residual norm
    double resTol;                     //residual tolerence
    double tInitial;                  //initial integration time
  

    public:
        //Default constructor
        LimexWrapper();
        //Destructor
        ~LimexWrapper(){}

        //Initialize limex
        void initLimex( const int nEq, std::vector<double> y, const double tStart=0);

        //set the tolerance parameters
        
        void setTolerances(const std::vector<double> rTols, const std::vector<double> aTols);

        //set the initial step size
        void setIntialialStep(const double step =0 );

        //set the upper and lower diagonal
        void setBandWidth(const int upper, const int lower);

        //start the integration
        void integrate(const double tFinal);
       
        //start the DAE integration
        void integrateDAE(const double tFinal);

        void calcResidualNorm();        

    private:

        int n;                            //number of equations
        std::vector<double> v_y;      //vector of dependent variables
        std::vector<double> v_yPrime; //derivatives of dependent variables
        
        std::vector<double> v_aTol;//   //Vector of absolute tolerance
        std::vector<double> v_rTol;//   //Vector of relative tolerance
                       
        double iniStep;               //initial step
        double hMax;                  //max allowed step size
        double tBegin;                //initial integration time
        double tEnd;                  //final integration time
        double tEnding;               //final integration time

        std::vector<int> iOpt;            //Integration control
        std::vector<int> iPos;            //constraint for +veness
        std::vector<double> rOpt;     //Integration control
        std::vector<int> iFail;           //Error indication

        int bandUpper;                    //upper diagonal
        int bandLower;                    //lower diagonal

    };
}


#endif	/* LIMEX_WRAPPER_H */
