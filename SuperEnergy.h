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
private:
	double _energy1, _energy2;
public:
	const double& energy1(){return _energy1;};
	const double& energy2(){return _energy2;};
        //SuperEnergy(){};
        SuperEnergy(SpMat& sup)
        {
                
                int a(6);

		SparseGenMatProd<double> op(sup);
                SymEigsSolver<double, SMALLEST_ALGE, SparseGenMatProd<double>> eigs(&op, 2, a);
                eigs.init();
                eigs.compute(1000000);
                if (eigs.info() == SUCCESSFUL)
                {
                        //para.Energy = eigs.eigenvalues()(0);cout.precision(15);
			
			_energy1=eigs.eigenvalues()(0)<eigs.eigenvalues()(1)?eigs.eigenvalues()(0):eigs.eigenvalues()(1);

			_energy2=eigs.eigenvalues()(0)>eigs.eigenvalues()(1)?eigs.eigenvalues()(0):eigs.eigenvalues()(1);
                        
                }

                
        };
        

        
        
};



#endif
