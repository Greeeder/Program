//#include "stdafx.h"
#include "TChH2Ol.h"
#include <cmath>
//#include <cmath>
//#include "CoolPropLib.h"

TChH2Ol::TChH2Ol(std::string nombre, double RU, double PMA) {
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
	_DENSF = 1000;
}

TChH2Ol::~TChH2Ol() {
	// TODO Auto-generated destructor stub
}

double TChH2Ol::FunCv(double T) {

	double RaizdeT = sqrt(T);

	return FunCp(T) - _R;

}

double TChH2Ol::FunCp(double T) {

	double RaizdeT = sqrt(T);

	return (22.605 - 0.09067 * RaizdeT + (-826.53 * RaizdeT + 13970.1 - 82114 / RaizdeT) / T) * _R;

}

double TChH2Ol::FunU(double T) {

	return FunCv(T) * T;
}

double TChH2Ol::FunH(double T) {

	return FunU(T) + FunR() * T;
}

double TChH2Ol::Funk(double T) {

	return ((9.496332E-09 * T - 1.707697E-05) * T + 9.183462E-03) * T - 8.626578E-01;
}

double TChH2Ol::FunVisc(double T) {

	return ((-2.632351E-09 * T + 2.737629E-06) * T - 9.530709E-04) * T + 1.114642E-01;

}

double TChH2Ol::FunPr(double T) {

	return ((-2.022269E-05 * T + 2.106518E-02) * T - 7.340298E+00) * T + 8.581110E+02;
}

double TChH2Ol::FunRHO(double T) {

	return (-3.6002E-03*T + 1.9045E+00) * T + 7.4865E+02;
}
