//#include "stdafx.h"
#include "TChCO2.h"
#include <cmath>
//#include <cmath>
//#include "CoolPropLib.h"

TChCO2::TChCO2(std::string nombre, double RU, double PMA) {
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

TChCO2::~TChCO2() {
	// TODO Auto-generated destructor stub
}

double TChCO2::FunCv(double T) {

	double RaizdeT = sqrt(T);

	return FunCp(T) - _R;

}

double TChCO2::FunCp(double T) {

	double RaizdeT = sqrt(T);

	return (12.019 - 0.03566 * RaizdeT + (-142.34 * RaizdeT - 163.7 + 9470 / RaizdeT) / T) * _R;

}

double TChCO2::FunU(double T) {

	return FunCv(T) * T;
}

double TChCO2::FunH(double T) {

	return FunU(T) + FunR() * T;
}

