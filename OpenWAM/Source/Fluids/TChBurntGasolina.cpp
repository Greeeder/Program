/**
 * @file TChBurntGasolina.cpp
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
 * The TChBurntGasolina class represents the chemical component burnt gas produced by a
 * combustion with gasoline.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of a burnt gas.
 *
 */

//#include "stdafx.h"
#include "TChBurntGasolina.h"
#include <cmath>
#include "Math_wam.h"

TChBurntGasolina::TChBurntGasolina(std::string nombre, double RU, double PM) {
	setName(nombre);
	_R = RU / PM;
	_PM = PM;
	_HV = 0.;
	_OC = 0.;
	_HC = 0.;
	_FE = 0.;
	_CAHFL = 0.;
	_CBHFL = 0.;
	_DENSF = 1.2;
}

TChBurntGasolina::~TChBurntGasolina() {
	// TODO Auto-generated destructor stub
}

double TChBurntGasolina::FunCv(double T) {

	return 655.1103865 + 0.401687993 * T - 6.98091E-05 * pow2<double>(T) - 9.31055E-09 * pow3<double>
		   (T) + 2.53355E-12 * pow4<double>(T);
}

double TChBurntGasolina::FunCp(double T) {

	//double RaizdeT = sqrt(T);

	//return  641.154 + T * (0.43045 + T * (-0.0001125 + T * 8.979e-9));

	return FunCv(T) + FunR();
}

double TChBurntGasolina::FunU(double T) {
	return -3205539.609 + 660.7221242 * T + 0.290764161 * pow(T, 2) - 0.000133742 * pow(T, 3) + 3.25919E-08 * pow(T, 4);
}

double TChBurntGasolina::FunH(double T) {
	return FunU(T) + FunR() * T;
}

