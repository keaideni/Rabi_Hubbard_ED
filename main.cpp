//#include "test.h"
#include "SuperEnergy.h"
#include "Wave.h"
#include "mpi.h"
//#include <omp.h>





int main(int argc, char* argv[])
{
        
	
	vector<Parameter> vecpara;
	//int nworks(2);
	//int num_thread(2);

	//omp_set_num_threads(num_thread);
	
	int myid, numprocess;
	int groupn(2);

	int nproper(9);//to be change as the number of properties change.
	for(int i=0; i<groupn; ++i)
	{
		Parameter para;
		para.ChangeJr(0.075*i);
		vecpara.push_back(para);
	}
	//#pragma omp parallel for
	

	MPI_Status status;
	MPI_Init(&argc, &argv);
	

	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocess);
	//myid=0;numprocess=1;

	int everygroup(groupn/numprocess);

	if(myid==0)
	{
				
		vector<vector<Energy<double>>> a(1);
		
		for(int i=0; i<everygroup; ++i)
		{
			a.at(myid).push_back(calcuwave(vecpara, myid, i, everygroup));
		}

		for(int id=1; id<numprocess; ++id)
		{
			MatrixXRd b(MatrixXd::Zero(everygroup, nproper));

			MPI_Recv(&b(0, 0), everygroup*nproper, MPI_DOUBLE, id, id, MPI_COMM_WORLD, &status);
	
			vector<Energy<double>> temp(mat2vec(b));
			a.push_back(temp);
			
		}
		ofstream outfile("result");
		outfile.precision(15);

		outfile<<"coupling"<<"\t"<<"ground"<<"\t"<<"excited"<<"\t"<<"orderparameter"<<"\t"<<"ParticleNo"<<"\t"<<"AParticleNo"<<"\t"<<"Parity"<<"\t"<<"Asqure"<<"\t"<<"SecondOrder"<<endl;

		for(auto it=a.begin(); it!=a.end();++it)
		{
			for(auto itt=it->begin(); itt!=it->end(); ++itt)
			{
				outfile<<(*itt).coupling<<"\t"<<(*itt).energy1<<"\t"<<(*itt).energy2<<"\t"<<(*itt).orderparameter<<"\t"<<(*itt).ParticleNo<<"\t"<<(*itt).AParticleNo<<"\t"<<(*itt).Parity<<"\t"<<(*itt).Asqure<<"\t"<<(*itt).SecondOrder<<endl;
			}//this should be changed if add new properties.
		}

		outfile.close();

	}

	else
	{
				
		vector<Energy<double>> a;
		for(int i=0; i<everygroup; ++i)
		{
			a.push_back(calcuwave(vecpara, myid, i, everygroup));
			//cout<<"haha"<<endl;
		}
		
		//cout<<"haha"<<endl;
		MatrixXRd temp(vec2mat(a));
		//cout<<temp<<endl;
		MPI_Send(&temp(0, 0), everygroup*nproper, MPI_DOUBLE, 0, myid, MPI_COMM_WORLD);


	}


	MPI_Finalize();

	return 0;

	//for(auto it=vecpara.begin(); it!=vecpara.end(); ++it)
	{
		//cout<<it->gcr()<<endl;
	}

	//#pragma_omp_parallel_for

	//MPI_Finalize();
	//cout<<para.Energy<<endl;
}
