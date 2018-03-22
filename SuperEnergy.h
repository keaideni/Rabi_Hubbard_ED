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

	VectorXd _state1, _state2;
public:
	const double& energy1()const{return _energy1;};
	const double& energy2()const{return _energy2;};
	const VectorXd& state1()const{return _state1;};
	const VectorXd& state2()const{return _state2;};
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

			_state1=eigs.eigenvalues()(0)<eigs.eigenvalues()(1)?eigs.eigenvectors(2).col(0):eigs.eigenvectors(2).col(1);
			_energy2=eigs.eigenvalues()(0)>eigs.eigenvalues()(1)?eigs.eigenvalues()(0):eigs.eigenvalues()(1);
			//vectorXd haha(eigs.eigenvectors(1))
			_state2=eigs.eigenvalues()(0)>eigs.eigenvalues()(1)?eigs.eigenvectors(2).col(0):eigs.eigenvectors(2).col(1);
//			cout<<eigs.eigenvectors().cols()<<endl;
                }

                
        };
        

        
        
};



#endif
