//#include "stdafx.h"
#include "TChH2Ov.h"
#include <cmath>
//#include <cmath>
//#include "CoolPropLib.h"

TChH2Ov::TChH2Ov(std::string nombre, double RU, double PMA) {
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

TChH2Ov::~TChH2Ov() {
	// TODO Auto-generated destructor stub
}

double TChH2Ov::FunCv(double T) {

	double RaizdeT = sqrt(T);

	return FunCp(T) - _R;

}

double TChH2Ov::FunCp(double T) {

	double RaizdeT = sqrt(T);

	return (22.605 - 0.09067 * RaizdeT + (-826.53 * RaizdeT + 13970.1 - 82114 / RaizdeT) / T) * _R;

}

double TChH2Ov::FunU(double T) {

	return FunCv(T) * T;
}

double TChH2Ov::FunH(double T) {

	return FunU(T) + FunR() * T;
}

