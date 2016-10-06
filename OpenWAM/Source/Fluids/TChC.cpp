//#include "stdafx.h"
#include "TChC.h"
#include <cmath>
//#include <cmath>
//#include "CoolPropLib.h"

TChC::TChC(std::string nombre, double RU, double PMA) {
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

TChC::~TChC() {
	// TODO Auto-generated destructor stub
}

double TChC::FunCv(double T) {

	return FunCp(T) - _R;

}

double TChC::FunCp(double T) {

	// TODO Añadir las propiedades del CO

	return 1000;

}

double TChC::FunU(double T) {

	return FunCv(T) * T;
}

double TChC::FunH(double T) {

	return FunU(T) + FunR() * T;
}

