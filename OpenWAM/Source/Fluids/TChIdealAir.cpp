/**
 * @file TChIdealAir.cpp
 * @author L.M. Garcia-Cuevas <luiga12@mot.upv.es>
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
 * The TChIdealAir class represents the chemical component Air.
 *
 * This class defines different methods to determine the chemical and
 * thermodynamics properties of ideal air.
 *
 */

#include "TChIdealAir.h"
#include <cmath>

TChIdealAir::TChIdealAir(std::string nombre, double RU, double PMA) {
	setName(nombre);
	_RU = __R::Universal;
	_PM = __PM::IdealAir;
	_R = _RU / _PM;
	_HV = 0.;
	_OC = 0.;
	_HC = 0.;
	_FE = 0.;
	_CAHFL = 0.;
	_CBHFL = 0.;
	_DENSF = 1.2;
}

TChIdealAir::~TChIdealAir() {
	// TODO Auto-generated destructor stub
}

double TChIdealAir::FunCv(double T) {
	return __Gamma::Cv;
}

double TChIdealAir::FunCp(double T) {
	return __Gamma::Cp;
}

double TChIdealAir::FunU(double T) {
	//! \todo Internal energy of ideal air.
	double RaizdeT = sqrt(T);

	return -4193697.9 - 6.94661 * T * RaizdeT + 2522.881 * T - 134454.16 * RaizdeT + 917124.39 * log(T) + 8349707.14 *
		   (1 / RaizdeT);
}

double TChIdealAir::FunH(double T) {

	//hainy = Ua + RU / v.PMA * (v.TF + v.C_TIY)
	//RU / v.PMA es la R del aire
	return FunU(T) + FunR() * T;
}

double TChIdealAir::Funk(double T) {
	return 2.5317018834285707E-2;
}

double TChIdealAir::FunVisc(double T) {
	return 1.789223325167428E-5;
}

