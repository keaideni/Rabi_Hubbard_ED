//#include "test.h"
#include "SuperEnergy.h"
#include <omp.h>
#include <fstream>


template<class T>
struct Energy{
	T coupling;
	T energy1;T energy2;
};



int main(int argc, char* argv[])
{
        
	
	vector<Parameter> vecpara;
	int nworks(10);
	int num_thread(2);

	omp_set_num_threads(num_thread);
	
	for(int i=0; i<nworks; ++i)
	{
		Parameter para;
		para.ChangeJr(0.01*i);
		vecpara.push_back(para);
	}
	vector<vector<Energy<double>>> a(num_thread);
	#pragma omp parallel for
	for(int i=0; i<nworks; ++i)
	{
		SingleSub a0(vecpara.at(i));
		SingleSub a1(vecpara.at(i), a0, a0);

		SingleSub a2(vecpara.at(i), a1, a1);//for the N=4 case.

		//SingleSub a3(vecpara.at(i), a2, a1);

		SpMat temp(a2.SysAdag1()*a2.SysA()+a2.SysAdag()*a2.SysA1());

		temp*=vecpara.at(i).Jr();temp*=-1;

		SpMat System(a2.System());
		System+=temp;

		SuperEnergy sup(System);

		Energy<double> tempdata={vecpara.at(i).Jr(),sup.energy1(), sup.energy2()};
		int j=omp_get_thread_num();
		a.at(j).push_back(tempdata);

	}

	ofstream outfile("result");
	outfile.precision(15);

	for(auto it=a.begin(); it!=a.end();++it)
	{
		for(auto itt=it->begin(); itt!=it->end(); ++itt)
		{
			outfile<<(*itt).coupling<<"\t"<<(*itt).energy1<<"\t"<<(*itt).energy2<<endl;
		}
	}

	outfile.close();

	return 0;

	//for(auto it=vecpara.begin(); it!=vecpara.end(); ++it)
	{
		//cout<<it->gcr()<<endl;
	}

	//#pragma_omp_parallel_for

	//MPI_Finalize();
	//cout<<para.Energy<<endl;
}
