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
#include"SuperEnergy.h"
using namespace std;

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


vector<amplitude> wave(const VectorXd& wave, const double& err);
vector<amplitude> wave(const VectorXd& wave, const double& err)
{
	vector<amplitude> tempwave;
	int it(0);
	for(int a1=0; a1<=6;++a1)
	{
		for(int q1=0;q1<=1;++q1)
		{
			
			for(int a2=0; a2<=6;++a2)
			{
				for(int q2=0;q2<=1;++q2)
				{

	for(int a3=0; a3<=6;++a3)
	{
		for(int q3=0;q3<=1;++q3)
		{

	for(int a4=0; a4<=6;++a4)
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
