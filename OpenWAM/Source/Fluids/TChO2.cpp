/**
 * @file TChO2.h
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
 * The TChO2 class represents the chemical component O2.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of the O2.
 *
 */

//#include "stdafx.h"
#include "TChO2.h"
#include <cmath>
//#include <cmath>
//#include "CoolPropLib.h"

TChO2::TChO2(std::string nombre, double RU, double PMA) {
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

TChO2::~TChO2() {
	// TODO Auto-generated destructor stub
}

double TChO2::FunCv(double T) {

	double RaizdeT = sqrt(T);

	return FunCp(T) - _R;

}

double TChO2::FunCp(double T) {

	double RaizdeT = sqrt(T);

	return (-0.112 + 0.0479 * RaizdeT + (195.42 * RaizdeT - 4426.1 + 32538 / RaizdeT) / T) * _R;

}

double TChO2::FunU(double T) {

	return FunCv(T) * T;
}

double TChO2::FunH(double T) {

	return FunU(T) + FunR() * T;
}

double TChO2::FunVisc(double T){
	// From CoolProp
	return ((-1.47482922732296e-11 * T) + 5.91566665017342e-08) * T + 4.63698691859697e-06;
}

