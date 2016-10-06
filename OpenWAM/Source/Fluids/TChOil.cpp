//#include "stdafx.h"
#include "TChOil.h"
#include <cmath>
//#include <cmath>
//#include "CoolPropLib.h"

TChOil::TChOil(std::string nombre, double RU, double PMA) {
	setName(nombre);
	_RU = RU;
	_PM = PMA;
	_R = _RU / _PM;
	_HV = 0.;
	_OC = 0.;
	_HC = 0.;
	_FE = 0.;
	_CAHFL = 0.;
	_CBHFL = 0.;
	_DENSF = 1.2;
}

TChOil::~TChOil() {
	// TODO Auto-generated destructor stub
}

double TChOil::FunCv(double T) {

	double RaizdeT = sqrt(T);

	return FunCp(T) - _R;

}

double TChOil::FunCp(double T) {

	VectorXd A(5);
	A(0) = -34040.12013;
	A(1) = -71.24978467;
	A(2) = 0.107361063;
	A(3) = -6.63966E-05;
	A(4) = 8670.542844;
	VectorXd vecT(5);
	vecT(0) = 1;
	vecT(1) = T;
	vecT(2) = pow2(T);
	vecT(3) = pow3(T);
	vecT(4) = log(T);

	return A.transpose() * vecT;

}

double TChOil::FunU(double T) {

	return FunCv(T) * T;
}

double TChOil::FunH(double T) {

	return FunU(T) + FunR() * T;
}

double TChOil::Funk(double T) {

	VectorXd A(5);
	A(0) = -13.78897843;
	A(1) = -0.028315003;
	A(2) = 3.80794E-05;
	A(3) = -2.25321E-08;
	A(4) = 3.437962885;
	VectorXd vecT(5);
	vecT(0) = 1;
	vecT(1) = T;
	vecT(2) = pow2(T);
	vecT(3) = pow3(T);
	vecT(4) = log(T);

	return A.transpose() * vecT;
}

double TChOil::FunVisc(double T) {

	return ViscCoef[0] * exp(ViscCoef[1] / (T - ViscCoef[2]));

}

double TChOil::FunPr(double T) {

	return FunVisc(T) * FunCp(T) / Funk(T);
}

double TChOil::FunRHO(double T) {
	VectorXd A(5);
	A(0) = -16572.74497;
	A(1) = -38.71871842;
	A(2) = 0.054641551;
	A(3) = -3.45327E-05;
	A(4) = 4398.316892;
	VectorXd vecT(5);
	vecT(0) = 1;
	vecT(1) = T;
	vecT(2) = pow2(T);
	vecT(3) = pow3(T);
	vecT(4) = log(T);

	return A.transpose() * vecT;

}
