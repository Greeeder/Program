/**
 * @file TChN2.h
 * @author F.J. Arnau <farnau@mot.upv.es>
 * @version 0.1.0
 *
 * @section LICENSE
 *
 * This file is part of OpenWAM.
 *
 * OpenWAM is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * OpenWAM is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with OpenWAM.  If not, see <http://www.gnu.org/licenses/>.
 *
 * @section DESCRIPTION
 * The TChN2 class represents the chemical component N2.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of the N2.
 *
 */

//#include "stdafx.h"
#include "TChN2.h"
#include <cmath>
//#include <cmath>
//#include "CoolPropLib.h"

TChN2::TChN2(std::string nombre, double RU, double PMA) {
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

TChN2::~TChN2() {
	// TODO Auto-generated destructor stub
}

double TChN2::FunCv(double T) {

	double RaizdeT = sqrt(T);

	return FunCp(T) - _R;

}

double TChN2::FunCp(double T) {

	double RaizdeT = sqrt(T);

	return (12.531 - 0.05932 * RaizdeT + (-352.3 * RaizdeT + 5279.1 - 27358 / RaizdeT) / T) * _R;

}

double TChN2::FunU(double T) {

	return FunCv(T) * T;
}

double TChN2::FunH(double T) {

	return FunU(T) + FunR() * T;
}

double TChN2::FunVisc(double T){
	// From CoolProp
	return ((-1.204057883911676e-11 * T) + 4.885447828536354e-08) * T + 4.620770796830561e-06;
}