#include "SingleSub.h"
#include <cmath>
#define PI 3.14159265358979323846
void SingleSub::Kron(SpMat& ab, const SpMat& a, const SpMat& b)
{
        ab.setZero();
        ab.resize(a.rows()*b.rows(), a.cols()*b.cols());

        for (int k=0; k<a.outerSize(); ++k)
                for(SpMat::InnerIterator it(a,k); it; ++it)
                {
                        for(int l=0; l<b.outerSize(); ++l)
                        {
                                for(SpMat::InnerIterator itt(b,l); itt; ++itt)
                                {
                                        ab.insert(it.row()*b.rows()+itt.row(), 
                                                it.col()*b.cols()+itt.col())
                                        =it.value()*itt.value();
                                }
                        }
                }

}

SingleSub::SingleSub(const Parameter& para):
_System(para.nmax()*2, para.nmax()*2),
_SysA(para.nmax()*2, para.nmax()*2),
_SysAdag(para.nmax()*2, para.nmax()*2),
_SysA1(para.nmax()*2, para.nmax()*2),
_SysAdag1(para.nmax()*2, para.nmax()*2),
_SysEye(para.nmax()*2, para.nmax()*2),
_ParticleNo(para.nmax()*2, para.nmax()*2),
_Parity(para.nmax()*2, para.nmax()*2)
{
        SpMat tempA(para.nmax()+1, para.nmax()+1), 
        tempAdag(para.nmax()+1, para.nmax()+1),
        tempEye(para.nmax()+1, para.nmax()+1),
	tempParity1(para.nmax()+1, para.nmax()+1);
        for(int i=0; i<para.nmax(); ++i)
        {
                tempA.insert(i, i+1)=sqrt(i+1);
                tempAdag.insert(i+1, i)=sqrt(i+1);
                tempEye.insert(i, i)=1;
		tempParity1.insert(i, i)=cos(PI*i);
        }//cout<<tempA<<endl;
	tempParity1.insert(para.nmax(), para.nmax())=cos(PI*(para.nmax()));
        tempEye.insert(para.nmax(), para.nmax())=1;
        SpMat Sigmamin(2,2), Sigmaplu(2,2), Sigmaeye(2,2), tempParity2(2,2);
        Sigmaeye.insert(0,0)=1; Sigmaeye.insert(1,1)=1;
        Sigmamin.insert(0,1)=1; Sigmaplu.insert(1,0)=1;
	
        tempParity2.insert(0,0)=cos(PI*0); tempParity2.insert(1,1)=cos(PI*1);

        Kron(_SysA, tempA, Sigmaeye);Kron(_SysAdag, tempAdag, Sigmaeye);
	Kron(_SysEye, tempEye, Sigmaeye);
	_SysA1=_SysA; _SysAdag1=_SysAdag;

        Kron(_System, tempAdag*tempA, Sigmaeye);
        //test=================================
        //_System*=-1;
        //=====================================

        SpMat temp;
        Kron(temp, tempEye, Sigmaplu*Sigmamin); _System+=temp; _ParticleNo=_System;
        Kron(temp, tempAdag, Sigmamin); temp*=para.gr(); _System+=temp;
        Kron(temp, tempA, Sigmaplu); temp*=para.gr(); _System+=temp;

        Kron(temp, tempAdag, Sigmaplu); temp*=para.gcr(); _System+=temp;
        Kron(temp, tempA, Sigmamin); temp*=para.gcr(); _System+=temp;

	Kron(_Parity, tempParity1, tempParity2);
	MatrixXd haha(_Parity);
	//cout<<haha<<endl;

}

SingleSub::SingleSub(const Parameter& para, const SingleSub& SubL, const SingleSub& SubR)
{
        Kron(_System, SubL.System(), SubR.SysEye());SpMat temp;
        Kron(temp, SubL.SysEye(), SubR.System());_System+=temp;
        Kron(temp, SubL.SysAdag(), SubR.SysA1());temp*=para.Jr();temp*=-1;_System+=temp;
        Kron(temp, SubL.SysA(), SubR.SysAdag1());temp*=para.Jr();temp*=-1;_System+=temp;
        Kron(temp, SubL.SysAdag(), SubR.SysAdag1());temp*=para.Jcr();temp*=-1;_System+=temp;
        Kron(temp, SubL.SysA(), SubR.SysA1());temp*=para.Jcr();temp*=-1;_System+=temp;

        Kron(_SysA, SubL.SysEye(), SubR.SysA());
        Kron(_SysAdag, SubL.SysEye(), SubR.SysAdag());
        Kron(_SysA1, SubL.SysA1(), SubR.SysEye());
        Kron(_SysAdag1, SubL.SysAdag1(), SubR.SysEye());

        Kron(_SysEye, SubL.SysEye(), SubR.SysEye());

	Kron(_ParticleNo, SubL._ParticleNo, SubR._SysEye);
	Kron(_Parity, SubL._Parity, SubR._Parity);
}


