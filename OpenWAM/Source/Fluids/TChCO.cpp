//#include "stdafx.h"
#include "TChCO.h"
#include <cmath>
//#include <cmath>
//#include "CoolPropLib.h"

TChCO::TChCO(std::string nombre, double RU, double PMA) {
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

TChCO::~TChCO() {
	// TODO Auto-generated destructor stub
}

double TChCO::FunCv(double T) {

	return FunCp(T) - _R;

}

double TChCO::FunCp(double T) {

	// TODO Añadir las propiedades del CO

	return 1000;

}

double TChCO::FunU(double T) {

	return FunCv(T) * T;
}

double TChCO::FunH(double T) {

	return FunU(T) + FunR() * T;
}

