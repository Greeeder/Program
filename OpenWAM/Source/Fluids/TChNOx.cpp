//#include "stdafx.h"
#include "TChNO.h"
#include <cmath>
//#include <cmath>

TChNO::TChNO(std::string nombre, double RU, double PMA) {
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

TChNO::~TChNO() {
	// TODO Auto-generated destructor stub
}

double TChNO::FunCv(double T) {

	return FunCp(T) - _R;

}

double TChNO::FunCp(double T) {

	// TODO Añadir las propiedades del NO

	return 1000;

}

double TChNO::FunU(double T) {

	return FunCv(T) * T;
}

double TChNO::FunH(double T) {

	return FunU(T) + FunR() * T;
}

