//for the M and N single site. Because the matrix is sparse, so construct a class especially.
#ifndef SINGLE_SUB_H
#define SINGLE_SUB_H
#include "Parameter.h"
#include <Eigen/Sparse>
using namespace std;
using namespace Eigen;

typedef Eigen::SparseMatrix<double> SpMat;

class SingleSub
 {

 private:
        SpMat _System;
        SpMat _SysA;
        SpMat _SysAdag;
        SpMat _SysA1;
        SpMat _SysAdag1;
	SpMat _SysEye;
        void Kron(SpMat& ab, const SpMat& a, const SpMat& b);

 public:
        const SpMat& System()const{return _System;};
        const SpMat& SysA()const{return _SysA;};
        const SpMat& SysAdag()const{return _SysAdag;};
        const SpMat& SysA1()const{return _SysA1;};
        const SpMat& SysAdag1()const{return _SysAdag1;};
	const SpMat& SysEye()const{return _SysEye;};
         SingleSub(){};
         ~SingleSub(){};
         SingleSub(const Parameter& para);
	 SingleSub(const Parameter& para, const SingleSub& a, const SingleSub& b);
         
 }; 


#endif // SINGLE_SUB_H
