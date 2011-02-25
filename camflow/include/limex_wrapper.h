#ifndef LIMEX_WRAPPER_H
#define	LIMEX_WRAPPER_H

#include <vector>
#include "cam_residual.h"

#define doublereal double



//using namespace Cratl;



namespace Gadgets{
    class LimexWrapper{
    
    std::vector<doublereal> fcn;          //Residual function
    std::vector<doublereal> jacobian;     //Jacobian of residual function

    Camflow::CamResidual *reactorPtr;
    int eqnSize;
    doublereal y;
    doublereal resNorm;                    //residual norm
    doublereal resTol;                     //residual tolerence
    doublereal tInitial;                  //initial integration time
  

    public:
        //Default constructor
        LimexWrapper();
        //Destructor
        ~LimexWrapper(){}

        //Initialize limex
        void initLimex( const int nEq, std::vector<doublereal> y, const doublereal tStart=0);

        //set the tolerance parameters
        
        void setTolerances(const std::vector<doublereal> rTols, const std::vector<doublereal> aTols);

        //set the initial step size
        void setIntialialStep(const doublereal step =0 );

        //set the upper and lower diagonal
        void setBandWidth(const int upper, const int lower);

        //start the integration
        void integrate(const doublereal tFinal);
       
        //start the DAE integration
        void integrateDAE(const doublereal tFinal);

        void calcResidualNorm();        

    private:

        int n;                            //number of equations
        std::vector<doublereal> v_y;      //vector of dependent variables
        std::vector<doublereal> v_yPrime; //derivatives of dependent variables
        
        std::vector<doublereal> v_aTol;//   //Vector of absolute tolerance
        std::vector<doublereal> v_rTol;//   //Vector of relative tolerance
                       
        doublereal iniStep;               //initial step
        doublereal hMax;                  //max allowed step size
        doublereal tBegin;                //initial integration time
        doublereal tEnd;                  //final integration time
        doublereal tEnding;               //final integration time

        std::vector<int> iOpt;            //Integration control
        std::vector<int> iPos;            //constraint for +veness
        std::vector<doublereal> rOpt;     //Integration control
        std::vector<int> iFail;           //Error indication

        int bandUpper;                    //upper diagonal
        int bandLower;                    //lower diagonal

    };
}


#endif	/* LIMEX_WRAPPER_H */
