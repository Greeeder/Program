/**
 * @file TChGasoil.cpp
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
 * The TChGasoil class represents the chemical component Gasoil.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of the Gasoil.
 *
 */

//#include "stdafx.h"
#include "TChGasoil.h"
#include "Math_wam.h"

TChGasoil::TChGasoil(std::string nombre, double CAHFL, double CBHFL, double RU, double PM, double HC, double OC,
					 double HV, double DENSF, double PMA, double XO2) {
	setName(nombre);
	_CAHFL = CAHFL;
	_CBHFL = CBHFL;
	_R = RU / PM;
	_PM = PM;

	_HC = HC;
	_OC = OC;
	_HV = HV;
	_DENSF = DENSF;

	_FE = 1 / ((PMA / XO2) * ((1 + HC / 4 - OC / 2) / (12 + HC + 16 * OC)));

}

TChGasoil::~TChGasoil() {
	// TODO Auto-generated destructor stub
}

double TChGasoil::FunCv(double T) {

	return -256.4 + 6.95372 * T - 0.00404715 * pow2<double>(T) + 9.10259E-07 * pow3<double>(T) + 1458487 *
		   (1 / pow2<double>(T));
}

double TChGasoil::FunCp(double T) {

	//double RaizdeT = sqrt(T);

	//return  641.154 + T * (0.43045 + T * (-0.0001125 + T * 8.979e-9));

	return 0.;
}

double TChGasoil::FunU(double T) {
	return -1234157.8 - 256.4 * T + 3.47686 * pow(T, 2) - 0.00134905 * pow(T, 3) + 2.27565E-07 * pow(T,
			4) - 1458487 * pow(T, (-1));
}

double TChGasoil::FunH(double T) {
	//return _CAHFL + _CBHFL * (_TF + _C_TIY);
	return _CAHFL + _CBHFL * T;
}

