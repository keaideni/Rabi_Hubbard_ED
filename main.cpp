//#include "test.h"
#include "SuperEnergy.h"



int main(void)
{
        Parameter para;
	
	SingleSub a0(para);
	SingleSub a1(para, a0, a0);

	SingleSub a2(para, a1, a1);//for the N=4 case.

	//SingleSub a3(para, a2, a1);

	SpMat temp(a2.SysAdag1()*a2.SysA()+a2.SysAdag()*a2.SysA1());

	temp*=para.Jr();temp*=-1;

	SpMat System(a2.System());
	System+=temp;

	SuperEnergy sup(para, System);
	cout<<para.Energy<<endl;
}
