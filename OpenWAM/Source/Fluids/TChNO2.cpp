//#include "stdafx.h"
#include "TChNO2.h"
#include <cmath>
//#include <cmath>

TChNO2::TChNO2(std::string nombre, double RU, double PMA) {
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

TChNO2::~TChNO2() {
	// TODO Auto-generated destructor stub
}

double TChNO2::FunCv(double T) {

	return FunCp(T) - _R;

}

double TChNO2::FunCp(double T) {

	// TODO Añadir las propiedades del NO2

	return 1000;

}

double TChNO2::FunU(double T) {

	return FunCv(T) * T;
}

double TChNO2::FunH(double T) {

	return FunU(T) + FunR() * T;
}

