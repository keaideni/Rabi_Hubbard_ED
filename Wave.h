/*************************************************************************
    > File Name: Wave.h
    > Author: ma6174
    > Mail: ma6174@163.com 
    > Created Time: 2018年03月22日 星期四 14时12分29秒
 ************************************************************************/

#include<iostream>
#include<algorithm>
#include<cmath>
#include<sstream>
#include <fstream>
#include <iomanip>
#include"SuperEnergy.h"
using namespace std;

typedef Matrix<double, Dynamic, 9> MatrixXRd;//the second number should be changed by the properties calculated.
struct amplitude
{
	double amp;
	vector<int> Q;
	int N;
};
bool comp(const amplitude& a, const amplitude& b);
bool comp(const amplitude& a, const amplitude& b)
{
	return abs(a.amp)>abs(b.amp);
}

template<class T>
struct Energy{
	T coupling;
	T energy1;T energy2;
	T orderparameter;
	T ParticleNo;
	T AParticleNo;
	T Parity;
	T Asqure;
	T SecondOrder;
};//save all the properties calculated. And should be change if there is more properties to calculate.

vector<amplitude> wave(const VectorXd& wave, const double& err, const int& nmax);
vector<amplitude> wave(const VectorXd& wave, const double& err, const int& nmax)
{
	vector<amplitude> tempwave;
	int it(0);
	for(int a1=0; a1<=nmax;++a1)
	{
		for(int q1=0;q1<=1;++q1)
		{
			
			for(int a2=0; a2<=nmax;++a2)
			{
				for(int q2=0;q2<=1;++q2)
				{

	for(int a3=0; a3<=nmax;++a3)
	{
		for(int q3=0;q3<=1;++q3)
		{

	for(int a4=0; a4<=nmax;++a4)
	{
		for(int q4=0; q4<=1; ++q4)
		{
			int temp[]={a1, q1, a2, q2, a3, q3, a4, q4};
			int N=a1+a2+a3+a4+q1+q2+q3+q4;
			vector<int> tempvec(temp, temp+sizeof(temp)/sizeof(int));
			amplitude tempamp={wave(it++), tempvec, N};
			tempwave.push_back(tempamp);
		}
	}
		}
	}
				}
			}
		}
	}

	sort(tempwave.begin(), tempwave.end(), comp);
	vector<amplitude> reamp;
	double total(0);it=0;
	while(1-total>err)
	{
		total+=pow(tempwave.at(it).amp, 2);
		reamp.push_back(tempwave.at(it++));
	}

	return reamp;
}

string itos(const int& i)
{

	stringstream s;
	s<<i;
	return s.str();
}

Energy<double> calcuwave(const vector<Parameter>& vecpara, const int& myid, const int i, const int& everygroup);
Energy<double> calcuwave(const vector<Parameter>& vecpara, const int& myid, const int i, const int& everygroup)
{

	SingleSub a0(vecpara.at(myid*everygroup+i));
	SingleSub a1(vecpara.at(myid*everygroup+i), a0, a0);
	SingleSub a2(vecpara.at(myid*everygroup+i), a1, a1);//for the N=4 case.
	//SingleSub a3(vecpara.at(myid*evergroup+i), a2, a0);

	SpMat temp(a2.SysAdag1()*a2.SysA()+a2.SysAdag()*a2.SysA1());

	temp*=vecpara.at(myid*everygroup+i).Jr();temp*=-1;

	SpMat System(a2.System());
	System+=temp;

	SuperEnergy sup(System);

	double orderpara(abs(sup.state1().transpose()*a2.SysA()*sup.state1()));
	//cout<<"haha"<<endl;
	double ParticleNo(sup.state1().transpose()*(a2.ParticleNo()*sup.state1()));
	double AParticleNo(sup.state1().transpose()*a2.SysAdag()*a2.SysA()*sup.state1());
	double Parity(sup.state1().transpose()*a2.Parity()*sup.state1());
	double Asqure(sup.state1().transpose()*a2.SysA()*a2.SysA()*sup.state1());
	double SecondOrder((sup.state1().transpose()*a2.SysAdag()*a2.SysAdag()*a2.SysA()*a2.SysA()*sup.state1()));
	
	Energy<double> tempdata={vecpara.at(myid*everygroup+i).Jr(),sup.energy1(), sup.energy2(), orderpara, ParticleNo, AParticleNo, Parity, Asqure, SecondOrder/pow(AParticleNo, 2)};//this should be changed if there is something more to calculate.
		
	string filename1("./wave/ground");
	string filename2("./wave/excited");
	filename1+=itos(myid*everygroup+i);
	filename2+=itos(myid*everygroup+i);
	ofstream outfile1(filename1);
	ofstream outfile2(filename2);

	outfile1<<vecpara.at(myid*everygroup+i).gr()<<"\t"<<vecpara.at(myid*everygroup+i).gcr()<<"\t"<<vecpara.at(myid*everygroup+i).Jr()<<"\t";
	outfile2<<vecpara.at(myid*everygroup+i).gr()<<"\t"<<vecpara.at(myid*everygroup+i).gcr()<<"\t"<<vecpara.at(myid*everygroup+i).Jr()<<"\t";
		
	outfile1.precision(15);
	outfile2.precision(15);
	vector<amplitude> ground(wave(sup.state1(), 0.0001, vecpara.at(myid*everygroup+i).nmax()));
	vector<amplitude> excited(wave(sup.state2(), 0.0001, vecpara.at(myid*everygroup+i).nmax()));

	outfile1<<ground.size()<<endl;
	outfile2<<excited.size()<<endl;

	for(auto it=ground.begin(); it!=ground.end(); ++it)
	{
		outfile1<<setw(20)<<(*it).amp<<"\t";
		for(auto itt=(*it).Q.begin(); itt!=(*it).Q.end();++itt)
		{
			outfile1<<(*itt)<<"\t";
		}
		outfile1<<(*it).N<<endl;
	}
	outfile1.close();
	for(auto it=excited.begin(); it!=excited.end(); ++it)
	{
		outfile2<<setw(20)<<(*it).amp<<"\t";
		for(auto itt=(*it).Q.begin(); itt!=(*it).Q.end();++itt)
		{
			outfile2<<(*itt)<<"\t";
		}
		outfile2<<(*it).N<<endl;
	}
	outfile2.close();
	return tempdata;
}

MatrixXRd vec2mat(const vector<Energy<double>>& vec);
MatrixXRd vec2mat(const vector<Energy<double>>& vec)
{

	MatrixXRd a(MatrixXd::Zero(vec.size(), 9));//cout<<vec.size()<<endl;
	int nrow(0);
	for(auto it=vec.begin(); it!=vec.end(); ++it)
	{
		VectorXd temp(9);
		temp<<it->coupling, it->energy1, it->energy2, it->orderparameter, it->ParticleNo, it->AParticleNo, it->Parity, it->Asqure, it->SecondOrder;//here should be changed if the properties calculated added.
		//cout<<temp<<endl;
		a.row(nrow++)=temp.transpose();
	}

	return a;
}

vector<Energy<double>> mat2vec(const MatrixXRd& a);
vector<Energy<double>> mat2vec(const MatrixXRd& a)
{
	vector<Energy<double>> temp;
	for(int i=0; i<a.rows(); ++i)
	{
		temp.push_back({a(i, 0), a(i, 1), a(i, 2), a(i, 3), a(i, 4), a(i, 5), a(i, 6), a(i, 7), a(i, 8)});
	}

	return temp;
}
