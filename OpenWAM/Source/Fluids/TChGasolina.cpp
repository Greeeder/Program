/**
 * @file TChGasolina.cpp
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
 * The TChGasolina class represents the chemical component gasoline.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of the Gasoline.
 *
 */

//#include "stdafx.h"
#include "TChGasolina.h"
#include <cmath>
#include "Math_wam.h"

TChGasolina::TChGasolina(std::string nombre, double CAHFL, double CBHFL, double RU, double PM, double HC, double OC,
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

TChGasolina::~TChGasolina() {
	// TODO Auto-generated destructor stub
}

double TChGasolina::FunCv(double T) {

	if(T < 1000) {
		return -12291002.2 * (1 / pow2<double>(T)) + 227580.041 * (1 / T) - 1545.51269 + 10.8382363 * T - 0.00837844 *
			   pow2<double>(T) + 3.2557E-06 * pow3<double>(T) - 4.0429E-10 * pow4<double>(T)
			   - 72.78128689;
	} else {
		return 984559799 * (1 / pow2<double>(T)) - 3394060.95 * (1 / T) + 5673.52925 + 1.036209 * T - 0.00036926 * pow2<double>
			   (T) + 5.2754E-08 * pow3<double>(T) - 2.7797E-12 * pow4<double>(T)
			   - 72.78128689 + 0.0312498;
	}
}

double TChGasolina::FunCp(double T) {

	//double RaizdeT = sqrt(T);

	//return  641.154 + T * (0.43045 + T * (-0.0001125 + T * 8.979e-9));

	return 0.;
}

double TChGasolina::FunU(double T) {
	if(T < 1000) {
		return (-1 * 3230332.559) + 12291002.16 * pow(T, (-1)) + 227580.0409 * log(T) - 1618.293972 * T + 5.41911816 * pow(T,
				2) - 0.002792812 * pow(T, 3) + 8.13916E-07 * pow(T, 4)
			   - 8.08583E-11 * pow(T, 5);
	} else {
		return (18516600.1 + 72.34589655) - 984559798.9 * pow(T,
				(-1)) - 3394060.946 * log(T) + 5600.74796 * T + 0.5181045 * pow(T, 2) - 0.00012309 * pow(T, 3) + 1.31188E-08 * pow(T, 4)
			   - 5.5593E-13 * pow(T, 5);
	}
}

double TChGasolina::FunH(double T) {
	//return _CAHFL + _CBHFL * (_TF + _C_TIY);

	return _CAHFL + _CBHFL * T;
}

