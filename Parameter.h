#ifndef PARAMETER_H
#define PARAMETER_H
#include<fstream>
#include<string>
#include<iostream>

class Parameter
{
private:
        double _gr, _gcr, _Jr, _Jcr;
        int _D;
        int _LatticeSize;
        int _nmax;
public:
        const double& gr()const{return _gr;};
        const double& gcr()const{return _gcr;};
        const double& Jr()const{return _Jr;};
        const double& Jcr()const{return _Jcr;};
        const int& LatticeSize()const{return _LatticeSize;};
        const int& D()const{return _D;};
        const int& nmax()const{return _nmax;};
	
	template<class T>
        void ChangeJcr(const T& dd){_Jcr=dd;};
	template<class T>
        void Changegcr(const T& dd){_gcr=dd;};
	template<class T>
        void Changegr(const T& dd){_gr=dd;};
	template<class T>
        void ChangeJr(const T& dd){_Jr=dd;};
        //Parameter();
        ~Parameter(){};
        double Energy;

        Parameter()
        {
                std::ifstream infile("Parameter");
                if(!infile)
                {
                        std::cerr<<"The file Parameter doesn't exit!"<<std::endl;
                        exit(true);
                }

                std::string temp;
                infile>>temp>>_nmax>>temp>>_D>>temp>>_LatticeSize>>temp>>_gr>>temp>>_gcr
                >>temp>>_Jr>>temp>>_Jcr;
        };
        
};


#endif // PARAMETER_H
