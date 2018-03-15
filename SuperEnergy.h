//The Initial wave must be coherent with the QWave order!!!!.

#include <Eigen/Core>
#include <SymEigsSolver.h>
#include <GenEigsSolver.h>
#include <MatOp/SparseGenMatProd.h>
#include <iostream>
#include "SingleSub.h"

using namespace Spectra;

#ifndef SUPERENERGY_H
#define SUPERENERGY_H
class SuperEnergy
{
public:

        //SuperEnergy(){};
        SuperEnergy(Parameter&para,SpMat& sup)
        {
                
                int a(6);

		SparseGenMatProd<double> op(sup);
                SymEigsSolver<double, SMALLEST_ALGE, SparseGenMatProd<double>> eigs(&op, 2, a);
                eigs.init();
                eigs.compute(1000000);
                if (eigs.info() == SUCCESSFUL)
                {
                        para.Energy = eigs.eigenvalues()(0);cout.precision(15);
			cout<<eigs.eigenvalues()(0)<<endl<<eigs.eigenvalues()(1)<<endl;
//				<<eigs.eigenvalues()(2)<<endl;
                        //std::cout << eigs.num_iterations() << std::endl;
                }

                
        };
        

        
        
};



#endif
