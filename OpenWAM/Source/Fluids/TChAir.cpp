/**
 * @file TChAir.h
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
 * The TChAir class represents the chemical component Air.
 *
 * This class include different method to determine the chemical and
 * thermodynamics properties of the Air.
 *
 */

//#include "stdafx.h"
#include "TChAir.h"
#include <cmath>
//#include <cmath>
//#include "CoolPropLib.h"

TChAir::TChAir(std::string nombre, double RU, double PMA) {
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

TChAir::~TChAir() {
	// TODO Auto-generated destructor stub
}

double TChAir::FunCv(double T) {

	double RaizdeT = sqrt(T);

	return -10.4199 * RaizdeT + 2522.88 + (-67227.1 * RaizdeT + 917124.4 - 4174853.6 / RaizdeT) / T;

	//return PropsSI("Cvmass","T",T+273,"P",1e5,"Air");

	//return -10.4199 * pow(T, 0.5) + 2522.88 - 67227.1 * pow(T, (-0.5)) + 917124.4 * pow(T, (-1)) - 4174853.6 * pow(T, (-1.5));
}

double TChAir::FunCp(double T) {

	//double RaizdeT = sqrt(T);

	//return -10.4199 * RaizdeT + 2522.88	+ (-67227.1 * RaizdeT + 917124.4 - 4174853.6 / RaizdeT) / T;

	//return PropsSI("Cvmass","T",T+273,"P",1e5,"Air");

	return FunCv(T) + FunR();
}

double TChAir::FunU(double T) {
	double RaizdeT = sqrt(T);

	return -4193697.9 - 6.94661 * T * RaizdeT + 2522.881 * T - 134454.16 * RaizdeT + 917124.39 * log(T) + 8349707.14 *
		   (1 / RaizdeT);
}

double TChAir::FunH(double T) {

	//hainy = Ua + RU / v.PMA * (v.TF + v.C_TIY)
	//RU / v.PMA es la R del aire
	return FunU(T) + FunR() * T;
}

double TChAir::Funk(double T) {

	return (-8.39061e-09 * T + 7.05256e-05) * T + 6.51528e-03;

}

double TChAir::FunVisc(double T) {

	return 1.4615e-6 * pow150(T) / (T + 110.4);

}

