#ifndef CALCU_H
#define CALCU_H
#define pi 3.1415926

#include <cmath>
#include "DMRG.h"
MatrixXd Kron(const MatrixXd& a, const MatrixXd& b);
void ReadTruncM(MatrixXd& A, const int& logo);

double ParticleNo(const Parameter& para);
double ParticleNo(const Parameter& para)
{
        Sub Sys, m(para,1);
        vector<MatrixXd> ParticleMat;
        ParticleMat.push_back(m.SysAdag()*m.SysA());

        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        ParticleMat.at(j)=(Kron(ParticleMat.at(j), m.SysEye()));
                        if(i==para.LatticeSize()/2)continue;
                        MatrixXd tempTrunc;
                        ReadTruncM(tempTrunc, i);
                        ParticleMat.at(j)=tempTrunc.adjoint()*ParticleMat.at(j)*tempTrunc;

                }
                Sys.Read(i-1);
                ParticleMat.push_back(Kron(Sys.SysEye(), m.SysAdag()*m.SysA()));
                if(i==para.LatticeSize()/2)break;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, i);
                ParticleMat.at(i-1)=tempTrunc.adjoint()*ParticleMat.at(i-1)*tempTrunc;
        }
        ParticleMat.push_back(m.SysAdag()*m.SysA());
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        ParticleMat.at(para.LatticeSize()/2+j)=
                        (Kron(ParticleMat.at(para.LatticeSize()/2+j), m.SysEye()));
                        MatrixXd tempTrunc;
                        if(i==para.LatticeSize()/2)continue;
                        ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
                        ParticleMat.at(para.LatticeSize()/2+j)
                        =tempTrunc.adjoint()*ParticleMat.at(para.LatticeSize()/2+j)*tempTrunc;

                }
                Sys.Read(para.LatticeSize()+2-i);
                ParticleMat.push_back(Kron(Sys.SysEye(), m.SysAdag()*m.SysA()));
                if(i==para.LatticeSize()/2)break;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
                ParticleMat.at(para.LatticeSize()/2+i-1)
                =tempTrunc.adjoint()*ParticleMat.at(para.LatticeSize()/2+i-1)*tempTrunc;
        }


        

        MatrixXd wave;
        ReadTruncM(wave, 10000);

        vector<double> Particle;
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                Particle.push_back((wave.adjoint()*ParticleMat.at(i)*wave).trace());
        }for(int i=para.LatticeSize()/2; i<para.LatticeSize(); ++i)
        {
                Particle.push_back((wave*(ParticleMat.at(i)).transpose()*wave.adjoint()).trace());
        }
        ofstream outfile("./result/ParticleNo");outfile.precision(20);
        for(int i=0; i<para.LatticeSize(); ++i)
        {
                outfile<<"i= "<<(i+1)<<" ,ParticleNo= "<<Particle.at(i)<<endl;
        }
        outfile.close();

        return Particle.at(para.LatticeSize()/2);


}





double SigmaParticleNo(const Parameter& para);
double SigmaParticleNo(const Parameter& para)
{
			
	Sub Sys, m(para, 1);
	Matrix2d tempsigma;
	tempsigma<<0,0,0,1;
	MatrixXd tempeye(MatrixXd::Identity(m.nmax+1, m.nmax+1));
	MatrixXd Sigma(Kron(tempeye, tempsigma));


	vector<MatrixXd> SigmaParticleMat;
	SigmaParticleMat.push_back(Sigma);
	for(int i=2; i<=para.LatticeSize()/2; ++i)
	{
		for(int j=0; j<i-1; ++j)
		{
			SigmaParticleMat.at(j)=Kron(SigmaParticleMat.at(j), m.SysEye());
			if(i==para.LatticeSize()/2)continue;
			MatrixXd tempTrunc;
			ReadTruncM(tempTrunc, i);
			SigmaParticleMat.at(j)=tempTrunc.adjoint()*SigmaParticleMat.at(j)*tempTrunc;
		}
		Sys.Read(i-1);
		SigmaParticleMat.push_back(Kron(Sys.SysEye(), Sigma));
		if(i==para.LatticeSize()/2)break;
		MatrixXd tempTrunc;
		ReadTruncM(tempTrunc, i);
		SigmaParticleMat.at(i-1)=tempTrunc.adjoint()*SigmaParticleMat.at(i-1)*tempTrunc;

	}
	SigmaParticleMat.push_back(Sigma);
	for(int i=2; i<=para.LatticeSize()/2;++i)
	{

		for(int j=0; j<i-1; ++j)
		{
							
			SigmaParticleMat.at(para.LatticeSize()/2+j)=
				Kron(SigmaParticleMat.at(para.LatticeSize()/2+j), m.SysEye());
			if(i==para.LatticeSize()/2)continue;
			MatrixXd tempTrunc;
			ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
			SigmaParticleMat.at(para.LatticeSize()/2+j)=
			tempTrunc.adjoint()*SigmaParticleMat.at(para.LatticeSize()/2+j)*tempTrunc;

		}
		Sys.Read(para.LatticeSize()+2-i);
		SigmaParticleMat.push_back(Kron(Sys.SysEye(), Sigma));
		if(i==para.LatticeSize()/2)break;
		MatrixXd tempTrunc;
		ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
		SigmaParticleMat.at(para.LatticeSize()/2+i-1)
		=tempTrunc.adjoint()*SigmaParticleMat.at(para.LatticeSize()/2+i-1)*tempTrunc;
	}

        MatrixXd wave;
        ReadTruncM(wave, 10000);

        vector<double> SigmaParticle;
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                SigmaParticle.push_back((wave.adjoint()*SigmaParticleMat.at(i)*wave).trace());
        }for(int i=para.LatticeSize()/2; i<para.LatticeSize(); ++i)
        {
                SigmaParticle.push_back
			((wave*(SigmaParticleMat.at(i)).transpose()*wave.adjoint()).trace());
        }
        ofstream outfile("./result/SigmaParticleNo");outfile.precision(20);
        for(int i=0; i<para.LatticeSize(); ++i)
        {
                outfile<<"i= "<<(i+1)<<" ,ParticleNo= "<<SigmaParticle.at(i)<<endl;
        }
        outfile.close();

        return SigmaParticle.at(para.LatticeSize()/2);
}





double OrderParameter(const Parameter& para);
double OrderParameter(const Parameter& para)
{
        Sub Sys, m(para,1);
        vector<MatrixXd> OrderMat;
        OrderMat.push_back(m.SysA());

        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        OrderMat.at(j)=(Kron(OrderMat.at(j), m.SysEye()));
                        if(i==para.LatticeSize()/2)continue;
                        MatrixXd tempTrunc;
                        ReadTruncM(tempTrunc, i);
                        OrderMat.at(j)=tempTrunc.adjoint()*OrderMat.at(j)*tempTrunc;

                }
                Sys.Read(i-1);
                OrderMat.push_back(Kron(Sys.SysEye(), m.SysA()));
                if(i==para.LatticeSize()/2)break;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, i);
                OrderMat.at(i-1)=tempTrunc.adjoint()*OrderMat.at(i-1)*tempTrunc;
        }
        OrderMat.push_back(m.SysA());
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        OrderMat.at(para.LatticeSize()/2+j)=
                        (Kron(OrderMat.at(para.LatticeSize()/2+j), m.SysEye()));
                        MatrixXd tempTrunc;
                        if(i==para.LatticeSize()/2)continue;
                        ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
                        OrderMat.at(para.LatticeSize()/2+j)
                        =tempTrunc.adjoint()*OrderMat.at(para.LatticeSize()/2+j)*tempTrunc;

                }
                Sys.Read(para.LatticeSize()+2-i);
                OrderMat.push_back(Kron(Sys.SysEye(), m.SysA()));
                if(i==para.LatticeSize()/2)break;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
                OrderMat.at(para.LatticeSize()/2+i-1)
                =tempTrunc.adjoint()*OrderMat.at(para.LatticeSize()/2+i-1)*tempTrunc;
        }


        

        MatrixXd wave;
        ReadTruncM(wave, 10000);

        vector<double> OrderNo;
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                OrderNo.push_back((wave.adjoint()*OrderMat.at(i)*wave).trace());
        }for(int i=para.LatticeSize()/2; i<para.LatticeSize(); ++i)
        {
                OrderNo.push_back((wave*(OrderMat.at(i)).transpose()*wave.adjoint()).trace());
        }
        ofstream outfile("./result/Order");outfile.precision(20);
        for(int i=0; i<para.LatticeSize(); ++i)
        {
                outfile<<"i= "<<(i+1)<<" ,Order= "<<OrderNo.at(i)<<endl;
        }
        outfile.close();

        return OrderNo.at(para.LatticeSize()/2);
}

void Correlation(const Parameter& para);
void Correlation(const Parameter& para)
{
        Sub m(para, 0);
        vector<MatrixXd> CorrMat;
        MatrixXd Adag1(m.SysAdag());

        CorrMat.push_back(m.SysAdag()*m.SysA());
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        CorrMat.at(j)=Kron(CorrMat.at(j), m.SysEye());
                        
                        if(i==para.LatticeSize()/2)continue;
                        
                        MatrixXd tempTrunc;
                        ReadTruncM(tempTrunc, i);
                        CorrMat.at(j)=tempTrunc.adjoint()*CorrMat.at(j)*tempTrunc;
                        
                }
                CorrMat.push_back(Kron(Adag1, m.SysA()));
                Adag1=Kron(Adag1, m.SysEye());
                if(i==para.LatticeSize()/2)continue;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, i);
                CorrMat.at(i-1)=tempTrunc.adjoint()*CorrMat.at(i-1)*tempTrunc;
                Adag1=tempTrunc.adjoint()*Adag1*tempTrunc;
        }

        vector<double> Corr;
        MatrixXd wave;
        ReadTruncM(wave, 10000);
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                Corr.push_back((wave.transpose()*CorrMat.at(i)*wave).trace());
        }
        Sub Sys;
        Sys.Read(para.LatticeSize()/2+2);
        MatrixXd A(Kron(Sys.SysA1(), m.SysEye()));
        Corr.push_back((wave.adjoint()*Adag1*wave*A.transpose()).trace());

        ofstream outfile("./result/Correlation");
        for(int i=0; i<Corr.size(); ++i)
        {
                outfile<<"R= "<<i<<" ,Corr(R)= "<<Corr.at(i)<<endl;
        }
        outfile.close();
}



void SigmaCorrelation(const Parameter& para);
void SigmaCorrelation(const Parameter& para)
{
        Sub m(para, 0);
        vector<MatrixXd> CorrMat;
	Matrix2d tempsigmaplu, tempsigmamin;
	tempsigmaplu<<0,0,1,0;        tempsigmamin<<0,1,0,0;
	MatrixXd tempeye(MatrixXd::Identity(m.nmax+1, m.nmax+1));
        MatrixXd Adag1(Kron(tempeye, tempsigmaplu)), A(Kron(tempeye, tempsigmamin));

        CorrMat.push_back(Adag1*A);
        for(int i=2; i<=para.LatticeSize()/2; ++i)
        {
                for(int j=0; j<i-1; ++j)
                {
                        CorrMat.at(j)=Kron(CorrMat.at(j), m.SysEye());
                        
                        if(i==para.LatticeSize()/2)continue;
                        
                        MatrixXd tempTrunc;
                        ReadTruncM(tempTrunc, i);
                        CorrMat.at(j)=tempTrunc.adjoint()*CorrMat.at(j)*tempTrunc;
                        
                }
                CorrMat.push_back(Kron(Adag1, A));
                Adag1=Kron(Adag1, m.SysEye());
                if(i==para.LatticeSize()/2)continue;
                MatrixXd tempTrunc;
                ReadTruncM(tempTrunc, i);
                CorrMat.at(i-1)=tempTrunc.adjoint()*CorrMat.at(i-1)*tempTrunc;
                Adag1=tempTrunc.adjoint()*Adag1*tempTrunc;
        }

        vector<double> Corr;
        MatrixXd wave;
        ReadTruncM(wave, 10000);
        for(int i=0; i<para.LatticeSize()/2; ++i)
        {
                Corr.push_back((wave.transpose()*CorrMat.at(i)*wave).trace());
        }
        Sub Sys;
        Sys.Read(para.LatticeSize()/2+2);
        MatrixXd A1(Kron(Sys.SysA1(), m.SysEye()));
        Corr.push_back((wave.adjoint()*Adag1*wave*A1.transpose()).trace());

        ofstream outfile("./result/SigmaCorrelation");
        for(int i=0; i<Corr.size(); ++i)
        {
                outfile<<"R= "<<i<<" ,Corr(R)= "<<Corr.at(i)<<endl;
        }
        outfile.close();
}



MatrixXd Kron(const MatrixXd& a, const MatrixXd& b)
{
        MatrixXd ab(MatrixXd::Zero(a.rows()*b.rows(),a.cols()*b.cols()));

        int sizer(b.rows()),sizec(b.cols());
        for(int i=0; i<a.rows(); ++i)
        {
                for(int j=0; j<a.cols(); ++j)
                {
                        int startr(i*b.rows()),startc(j*b.cols());

                        ab.block(startr,startc,sizer, sizec)=a(i,j)*b;
                }
        }

        return ab;
}

double secondcorrelation(const Parameter& para);
double secondcorrelation(const Parameter& para)
{
	Sub a,m(para,0);
	a.Read(para.LatticeSize()/2-1);
	MatrixXd wave;
	ReadTruncM(wave, 10000);
	MatrixXd Adag(Kron(a.SysEye(),m.SysAdag())), A(Kron(a.SysEye(),m.SysA()));
	double fenzi((wave.transpose()*Adag*Adag*A*A*wave).trace());
	double fenmu(pow((wave.transpose()*Adag*A*wave).trace(),2));

	return fenzi/fenmu;

}

double Parity(const Parameter& para);
double Parity(const Parameter& para)
{

        Sub m(para, 0);
	Matrix2d tempsigma;
	tempsigma<<0,0,0,1;
	MatrixXd tempeye(MatrixXd::Identity(m.nmax+1, m.nmax+1));
        MatrixXd particleno(m.SysAdag()*m.SysA()+Kron(tempeye, tempsigma));
	
	//std::cout<<particleno<<endl<<"haha"<<endl;
	for(int i=0; i<particleno.rows(); ++i)
	{
		particleno(i,i)=cos(particleno(i,i)*pi);
	}
	//std::cout<<particleno<<endl;
	
	MatrixXd paritymatrix1(particleno);
	for(int i=2; i<para.LatticeSize()/2; ++i)
	{
		//std::cout<<paritymatrix1<<endl<<"haha"<<endl;
		paritymatrix1=Kron(paritymatrix1, particleno);

		//std::cout<<paritymatrix1<<endl<<"haha"<<endl;
		MatrixXd tempTrunc;
		ReadTruncM(tempTrunc, i);
		//MatrixXd temp,temp1,temp2;
		paritymatrix1=tempTrunc.adjoint()*paritymatrix1*tempTrunc;
		/*temp1=tempTrunc.adjoint()*tempTrunc;
		temp2=MatrixXd::Identity(tempTrunc.cols(), tempTrunc.cols());
		double sum(0);
		for(int irow=0; irow<tempTrunc.cols(); ++irow)
		{

			for(int icol=0; icol<tempTrunc.cols(); ++icol)
			{
				sum+=abs(temp2(irow, icol)-temp1(irow, icol));
			}
		}
		std::cout<<sum<<endl<<"haha"<<endl;
		paritymatrix1=temp;
		
		sum=0;
		for(int irow=0; irow<tempTrunc.cols(); ++irow)
		{

			for(int icol=0; icol<tempTrunc.cols(); ++icol)
			{
				if(irow==icol)continue;
				sum+=abs(paritymatrix1(irow, icol));
			}
		}
		std::cout<<sum<<endl<<"haha"<<endl;*/
	}
	paritymatrix1=Kron(paritymatrix1, particleno);

	
	MatrixXd paritymatrix2(particleno);
	for(int i=2; i<para.LatticeSize()/2; ++i)
	{
		paritymatrix2=Kron(paritymatrix2, particleno);
		MatrixXd tempTrunc;
		ReadTruncM(tempTrunc, para.LatticeSize()-i+1);
		paritymatrix2=tempTrunc.adjoint()*paritymatrix2*tempTrunc;
	}
	paritymatrix2=Kron(paritymatrix2, particleno);

	
        MatrixXd wave;
        ReadTruncM(wave, 10000);

	return (wave.adjoint()*paritymatrix1*wave*paritymatrix2.transpose()).trace();


}

#endif // CALCU_H
