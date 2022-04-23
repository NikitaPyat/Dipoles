#include "pch.h"
#include <iostream>
#include <fstream>
#include <complex>
#include <math.h>
#include "DipoleFields.h"
using namespace std;

    const double PI = 3.14159265359;
	const double mu = 4 * PI * 1.0E-7;
	const double eps0 = 8.8541878e-12;
	double lam0om1 = 2 * PI / sqrt(eps0 * mu); // для частоты omega = 1

	double eps1;
	double eps2;
	double h;
	double ro;
	double z0;
	double z;
	double w;
	double omega;

	double k0k0;
	double k1k1;
	double k2k2;
	double k10;
	double k21;
	ofstream ofs;
	
	// передаются нормированные параметры из ModelParams
	void Init_ModelParams(double _eps1, double _eps2, double _h,
		double _ro, double _z0, double _z, double _omega)
	{
		eps1 = _eps1;
		eps2 = _eps2;
		omega = _omega;
		w = lam0om1 / omega;

		h = _h * w;
		ro = _ro * w;
		z0 = _z0 * w;
		z = _z * w;
			
		k0k0 = 4 * PI * PI / w / w;
		k1k1 = k0k0 * eps1;
		k2k2 = k0k0 * eps2;

		k10 = eps1;
		k21 = eps2 / eps1;
	};

	extern "C"  _declspec(dllexport)
	void Test(double _eps1, double _eps2, double _h,
		      double _ro, double _z0, double _z, double _wave,
		      int nlam, double* lam, double* res)
	{
		
		Init_ModelParams(_eps1, _eps2, _h, _ro, _z0, _z, _wave);
		
		double fValues[2];

		for (int j = 0; j < nlam; j++)
		{
			double lambda = lam[j];
			W(lambda, 2, fValues);
			res[j * 2] = fValues[0];
			res[j * 2 + 1] = fValues[1];
		}
	}

	void W(double lambda, int nf, double* fValues)
	{
		/*ofs.open("testR1.txt", ios::app);
	
		ofs << "PI = " << PI << endl;
		ofs << "mu = " << mu << endl;
		ofs << "eps0 = " << eps0 << endl;
		ofs << "lam0om1 = " << lam0om1 << endl;

					ofs << "eps1 = " << eps1 << endl;
					ofs << "eps2 = " << eps2 << endl;
					ofs << "omega = " << omega << endl;
					ofs << "w = " << w << endl;
					ofs << "h = " << h << endl;
					ofs << "z = " << z << endl;
					ofs << "z0 = " << z0 << endl;
					ofs << "ro = " << ro << endl;
					ofs << "k0k0 = " << k0k0 << endl;
					ofs << "k1k1 = " << k1k1 << endl;
					ofs << "k2k2 = " << k2k2 << endl;
					ofs << "k10 = " << k10 << endl;
					ofs << "k21 = " << k21 << endl;

					ofs << "lambda = " << lambda << endl;*/

		complex<double> tmp = complex<double> (lambda * lambda - k0k0, 0);
		complex<double> nu0 = sqrt(tmp);

		tmp = complex<double>(lambda * lambda - k1k1, 0);
		complex<double> nu1 = sqrt(tmp);
		complex<double> nu1k21 = k21 * nu1;

		tmp = complex<double>(lambda * lambda - k2k2, 0);
		complex<double> nu2 = sqrt(tmp);

		complex<double> e1 = exp(-2.0 * nu1 * h);
		complex<double> F1 = (nu1k21 - nu2) / (nu1k21 + nu2) * e1;

		complex<double> G0;

		if (abs(lambda * lambda - k1k1) > 1.0E-7)
		{
		 complex<double> G0N = k10 * (F1 + 1.0);
		 complex<double> G0D = nu0 * G0N - nu1 * (F1 - 1.0);
		 G0 = 2.0 * G0N / G0D;
		}
    	else
		{
		  complex<double> G0N = k10 * (k21 + h * nu2);
		  complex<double> G0D = 2.0 * G0N * nu0 + nu2 * (1.0 - F1);
		  G0 = 4.0 * G0N / G0D;
		}
		complex<double> e0 = exp(-nu0 * (z + z0));
		complex<double>	WF = G0 * e0;
	    fValues[0] = WF.real();
		fValues[1] = WF.imag();
		//ofs.close();
	}



//	complex<double> getU1(double lambda) {
//
//		ofs.open("testR1.txt", ios::app);
//
//		complex<double> nu0 = sqrt(complex<double>(lambda * lambda - k0k0, 0));
//		complex<double> nu1 = sqrt(complex<double>(lambda * lambda - k1k1, 0));
//		complex<double> nu2 = sqrt(complex<double>(lambda * lambda - k2k2, 0));
//
//		complex<double> nu12 = (nu1 - nu2) / (nu1 + nu2);
//		complex<double> nu01 = (nu0 - nu1) / (nu0 + nu1);
//
//		if (lambda > 1.0E-7)
//		{
//			complex<double> R1 = nu12 * exp(-2.0 * nu1 * h);
//			complex<double> R1p = 1.0 + R1;
//			complex<double> R1m = 1.0 - R1;
//			complex<double> tmp1 = 2.0 * R1p / (nu1 * R1m + nu0 * R1p);
//			complex<double> tmp2 = mu * exp(-nu0 * (z + z0));
//			ofs << "lambda = " << lambda << endl;
//			ofs << "R1 = " << R1 << endl;
//			ofs << "tmp1 = " << tmp1 << endl;
//			ofs.close();
//			return tmp1 * tmp2;
//		}
//		else
//		{
//			complex<double> tmp1 = 1.0 + nu2 * h;
//			complex<double> tmp2 = nu0 * tmp1 + nu2 - nu1 - nu1 * nu2 * h;
//			complex<double> tmp3 = tmp1 / tmp2;
//			complex<double> tmp4 = mu * exp(-nu0 * (z + z0));
//
//			ofs << "lambda = " << lambda << endl;
//			ofs << "lambda * lambda = " << lambda * lambda << endl;
//			ofs << "k2k2 = " << k2k2 << endl;
//			ofs << "lambda * lambda - k2k2 = " << lambda * lambda - k2k2 << endl;
//
//			ofs << "nu0 = " << nu2 << endl;
//			ofs << "nu1 = " << nu2 << endl;
//			ofs << "nu2 = " << nu2 << endl;
//
//			ofs << "tmp1 = " << tmp1 << endl;
//			ofs << "tmp2 = " << tmp2 << endl;
//			ofs << "tmp3 = " << tmp3 << endl;
//			ofs << "tmp4 = " << tmp4 << endl;
//			ofs.close();
//			return tmp3 * tmp4;
//		}
//	}
//};